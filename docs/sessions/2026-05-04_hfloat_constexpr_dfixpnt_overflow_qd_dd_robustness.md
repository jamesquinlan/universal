# Session: hfloat constexpr, dfixpnt wide-instantiation overflow, qd/dd to_digits robustness

**Date:** 2026-05-04
**Branches:** `main`, `feat/issue-732-hfloat-constexpr` (merged), `fix/issue-804-dfixpnt-wide-overflow` (merged), `fix/issue-801-to-digits-renormalize` (PR #808 open)
**Build directories:** `build_ci/` (gcc), `build_ci_clang/` (clang), `build_ubsan/` (clang + UBSan, new)

## Objectives

1. Verify PR #805 (dfloat constexpr) post-merge state vs CodeRabbit review findings
2. Resolve issue #732 -- bring `hfloat` to full constexpr
3. Resolve issue #804 -- dfixpnt overflow in `to_int64()` and `operator=(double)` for `ndigits >= 20`
4. Resolve issue #801 -- `to_digits()` robustness against non-canonical limb inputs across `qd`, `dd`, and `floatcascade<N>`
5. Document a local UBSan workflow after CI caught a regression that the gcc/clang Release builds missed

## Summary of Changes

### PR #805 verification (no-op)

Started the session investigating CodeRabbit's review on the merged PR #805 (dfloat
constexpr). Discovered all six findings (2 Critical, 4 Major) had already been
addressed in a fixup commit `ca043f2a` that was squash-merged into `0d4375b2`.
Verified each fix in the current `main`:

| # | Finding | Location |
|---|---------|----------|
| 1 | `pow10_64_safe` chunked rescale | `dfloat_impl.hpp:1099-1127` (chunked 10^19 multiply/divide) |
| 2 | `convert_signed/unsigned` 64-bit boundary | `dfloat_impl.hpp:1234-1271` (`-(v+1)+1ull`, `setbits(abs_v)`) |
| 3 | Plain `char` via `convert_signed` | `dfloat_impl.hpp:248` |
| 4 | `bit_cast` guarded for -0.0 detection | `dfloat_impl.hpp:1076-1083` (`if constexpr (sw::is_bit_cast_constexpr_v)`) |
| 5 | `BLOCKBINARY_THROW_ARITHMETIC_EXCEPTION` wiring | `dfloat.hpp:44-47` (forwarding removed, comment) |
| 6 | `constexpr.cpp` value-initialize | `static/float/dfloat/api/constexpr.cpp:137-145` |

CodeRabbit's threads still appeared unresolved on the GitHub UI; the actual code
was correct. Created and immediately deleted `fix/dfloat-postmerge-cr-805`.

### hfloat full constexpr (#732, PR #806)

Moderate-scope promotion paralleling dfloat #805 / dfixpnt #803 / qd #800. hfloat
is self-contained (no `blockbinary` / `blocktriple` dependencies), so the
promotion was localized to 3 files.

**Files modified:**
- `include/sw/universal/number/hfloat/hfloat_impl.hpp` -- ~355 line diff
- `include/sw/universal/number/hfloat/hfloat_fwd.hpp` -- 4-line diff
- `static/float/hfloat/api/constexpr.cpp` -- new (240 lines)

**Key transformations:**
- All trivial helpers, selectors, ctors, assignments, arithmetic (+= -= *= /= ++ --),
  comparisons (== != < <= > >=), unary minus, conversion-out (operator float/double/
  long double), abs, fabs marked `constexpr`
- `convert_ieee754(double)` rewritten without `std::isnan/isinf/fabs/frexp/ldexp`:
  - NaN via `x != x`
  - Infinity via `numeric_limits<double>::max()` bracket
  - `if constexpr (sw::is_bit_cast_constexpr_v)` extracts IEEE 754 fields directly;
    runtime fallback retains `std::frexp/ldexp`
  - Mantissa-with-hidden-bit shifting replaces `ldexp(frac, shift)` with integer math
- `convert_to_double()` rewritten with constexpr power-of-2 scaling loop;
  `is_constant_evaluated()` dispatch keeps runtime on `std::ldexp`
- Plain `char` routed via `convert_signed(static_cast<int>(rhs))` -- applied the
  CodeRabbit lesson from PR #805 preemptively
- `operator/=` divide-by-zero throw fenced under `!std::is_constant_evaluated()`

**Tests added in `constexpr.cpp`:**
- Issue #732 acceptance form
- All four binary operators + compound forms (equality via `(a-b).iszero()` to
  handle HFP wobbling precision)
- All six comparison operators
- `operator double()` / `operator float()` constexpr verification
- Increment/decrement on value-initialized `HShort{}`
- SpecificValue construction with HFP saturation invariants
  (`infpos == maxpos`, `qnan == zero`)
- `isinf()` / `isnan()` always return false
- `1/0` and `0/0` saturate to zero in constant-eval

**CR follow-up commit `ca043f2a`-style:**
- CR Critical: silent-zero in `convert_ieee754` runtime fallback when
  `sw::is_bit_cast_constexpr_v` was false AND constant-evaluated. My `else { return
  *this; }` path produced a wrong result instead of a compile-time error.
  Fix: removed the `if constexpr` guard entirely; bit_cast is called
  unconditionally. On platforms with constexpr bit_cast it works at constant-eval;
  on platforms without, the constexpr caller naturally fails to compile pointing
  at the bit_cast call -- the correct behavior
- CR Major: `operator==` / `operator<` were lossy via `double()` round-trip
  (`hfloat_long` has 56-bit fraction, double has 53). Replaced with
  sign/exponent/fraction tuple compare via constexpr `unpack`. After
  `normalize_and_pack`, equal values have identical tuples; ordering well-defined
  per HFP semantics
- CR Major: `convert_signed/unsigned` round-tripped via double. Added
  `pack_uint64()` direct integer packing for `fbits <= 56` (highest bit position
  -> hex_exp -> fraction shift); for `fbits >= 64` (hfloat_extended) keeps double
  round-trip since uint64 can't hold a 112-bit fraction. INT64_MIN-safe via
  `-(v+1) + 1ull` identity
- Tests added: 2^53 vs 2^53+1 distinction in hfloat_long, INT64_MIN sign
  preservation
- Squash-merged at `a0443659`

### dfixpnt wide-instantiation overflow (#804, PR #807)

Two pre-existing UB bugs surfaced by CodeRabbit on PR #803, both reachable now
that conversion paths are marked constexpr. UB silently produces garbage at
runtime; in a constant expression the compiler hard-fails.

**Fix 1: `to_int64()` (`dfixpnt_impl.hpp:559`)**

LSD-first accumulator with `scale *= 10` overflowed `long long` for `idigits >= 19`
since 10^19 > LLONG_MAX (~9.22e18). `dfixpnt<25, 5>` (idigits=20) hit this -- scale
reached 10^20 -> UB.

Replaced with MSD-first Horner over `unsigned long long` with per-step overflow
detection; clamps to `[LLONG_MIN, LLONG_MAX]`. Matches `blockdecimal::to_long_long`
pattern.

**Fix 2: `operator=(double)` (`dfixpnt_impl.hpp:146`)**

Materialized `scaled` (potentially > UINT64_MAX for `ndigits >= 20`) into uint64_t
-- UB per C++20 [conv.fpint]. `max_magnitude = 10^ndigits - 1` exceeds UINT64_MAX
for `ndigits >= 20`, so the existing clamp didn't protect the cast.

Replaced uint64 materialization with FP-domain digit extraction:
- `q = v / 10.0`
- `q_floor = (q < 2^53) ? cast/back : q` (10.0 is exactly representable; for q
  in dense-integer range cast truncates correctly; above that q is its own floor)
- `digit = v - 10*q_floor` clamped to `[0, 9]` against FP noise

Constexpr-safe for arbitrary `ndigits` within double's exponent range (~308 digits).

**Tests added** to `static/fixpnt/decimal/api/constexpr.cpp`:
- `dfixpnt<25, 5>` integer round-trip + LLONG_MAX/MIN clamp at maxpos/maxneg
- `dfixpnt<20, 5>` 15-digit boundary preservation + maxpos saturation
- `dfixpnt<25, 5>(1e15)` exercises wide `to_int64`

**CR follow-up:**
- CR Major: `radix > 308` or `ndigits > 308` would overflow the
  `for (i=0; i<radix; ++i) scaled *= 10.0` loop to +inf, and `inf >= max_magnitude
  (inf)` would silently coerce ANY input to all-nines. Added `static_assert(ndigits
  <= max_exponent10)` at the head of `operator=(double)` -- type remains usable
  for wider configurations via non-FP construction paths

Squash-merged at `6fe64ba3`.

### qd/dd/floatcascade to_digits robustness (#801, PR #808 open)

`to_digits()` in `qd`, `dd`, and the shared `floatcascade<N>` (used by
`dd_cascade` / `td_cascade` / `qd_cascade`) assumed canonical limb form
(high-magnitude-first, `|x[i+1]| <= ulp(x[i])/2`). Non-canonical inputs constructed
via the public raw-limb constructors drove the iterative digit-extraction loop to
NaN, surfacing as either spurious stderr warnings or, before PR #800's defensive
guard, as C++20 [conv.fpint] UB at `static_cast<int>(r[0])`.

**Issue #801 reproduction:**
```cpp
sw::universal::qd a(1.0, 0.0, 0.0, std::pow(2.0, -196.0));
std::cout << a;   // pre-fix: stderr warning + garbage output
```

**Initial fix:**
- `qd_impl.hpp::to_digits()` -- `r.renorm()` after the abs; reverted PR #800's
  defensive NaN guard (premature -- see below)
- `dd_impl.hpp` -- added `dd::renorm()` member (mirrors `qd::renorm`); call
  `r.renorm()` in `to_digits`
- `floatcascade.hpp::to_digits()` -- `r = expansion_ops::renormalize(r)` after the
  abs

**dd::renorm() inf-guard discovery:**

Initial dd::renorm using `quick_two_sum(hi, lo, lo)` corrupted `dd::maxpos`:
- `hi = DBL_MAX`, `lo` positive
- `s = hi + lo` overflows to +inf
- `new_lo = lo - (s - hi) = lo - inf = -inf`
- Result: `(inf, -inf)` -- broken

Found by diffing pre/post-fix outputs of `dd_api_api`. Added inf-guard:
```cpp
if (sw::universal::is_inf_cx(hi + lo)) return;  // boundary preserves canonical
```

**CodeRabbit review pass on PR #808 (3 Major findings):**

| # | Finding | Fix |
|---|---------|-----|
| 1 | `dd::renorm()` `quick_two_sum` requires `\|hi\| >= \|lo\|`; fails for `dd(0.0, 1e100)` | Switched to `two_sum` (order-independent + finite-result safety) |
| 2 | `dd::to_digits` computed `e` from pre-renorm `hi` | Moved frexp + iszero check AFTER renorm |
| 3 | `qd::to_digits` same issue | Same pattern |

After the CR fixes, `dd(1.0, 1e100)` correctly prints as `1.000...e+100`
(promote-lo case handled). `dd(0.0, 1e100)` directly via `dd::renorm()` now
produces canonical `(1e100, 0)`. Note: `dd(0.0, 1e100)` printed via `operator<<`
still returns `0` because `dd::to_string()` short-circuits via `dd::iszero()`
(which inspects only `this->hi`) BEFORE `to_digits` is reached. That is the
broader "Item 2" documentation/contract issue identified in #801; not in scope
for this PR's focused `to_digits` robustness fix.

**UBSan regression and second CR pass:**

CI UBSan caught: `qd_impl.hpp:1229: runtime error: -nan is outside the range of
representable values of type 'int'`. The renorm at to_digits entry IS sufficient
for the issue's reproduction case, but the iterative subtraction/multiplication in
the digit-extraction loop can still drift `r[0]` to NaN for extreme input
magnitudes (subnormal-dominant qd values exercised by `api.cpp`'s `<6,7>`
precision-progression cases at line 266).

The acceptance criterion in #801 ("the defensive NaN guard in PR #800 can be
reverted once the algorithm is robust at the entry point") was over-optimistic:
the loop iterations are also a source of drift, not just the entry. Re-added the
`(v != v) ? 0 : static_cast<int>(v)` guard at the cast site, mirrored to dd and
floatcascade for defense-in-depth.

**Local UBSan build directory:**

After CI caught the regression, set up `build_ubsan/` locally with `clang + Debug
+ -DUNIVERSAL_ENABLE_UBSAN=ON` mirroring `.github/workflows/sanitizers.yml`. Local
UBSan run is now part of the verification loop for changes touching numeric
conversion paths. Documented the workflow in `docs/build/local-sanitizer.md`
(committed by the user).

**PR title fix:**

Original title `fix(qd,dd,floatcascade): ...` failed `lint-pr-title` (the
amannn/action-semantic-pull-request action requires a single scope from a fixed
list, and `floatcascade` was not in the allowlist). Changed to `fix(qd):
to_digits robust to non-canonical limbs (qd/dd/floatcascade)` -- single scope
in the parenthetical, full coverage in the description.

PR #808 still open at session close, awaiting full CI cycle on doc-fix commit
`d8ae6bd6`.

## Key Decisions

### hfloat: skip raw-limb constructor change (option b)

Issue #732 mentioned the open question of whether to renormalize in the raw-limb
constructor. Chose to fix only at `to_digits` entry (effectively option a applied
locally) rather than touching all constructors. Rationale: the constructor is on
the hot path for EFT-output values (which are already canonical); fixing
locally at I/O is precision-equivalent for canonical inputs and corrects
non-canonical ones at the moment they matter.

### qd::to_digits: keep PR #800 NaN guard

The acceptance criterion in #801 read "the defensive NaN guard in PR #800 can be
reverted once the algorithm is robust at the entry point." Initial fix removed
the guard; UBSan caught a residual NaN drift INSIDE the loop. Restored the guard
and mirrored it to dd and floatcascade. The guard is `O(1)` per iteration and
cheap; the entry-point renorm covers the canonicalization failure mode and the
guard covers the iterative drift mode. Both are needed.

### Local UBSan workflow

After CI caught the qd loop drift that gcc/clang Release builds passed silently,
created `build_ubsan/` configured exactly as the CI sanitizers job. Now part of
the pre-push verification for changes touching numeric conversion. Documented in
`docs/build/local-sanitizer.md` (committed by the user as `8a500b34`).

## Challenges & Solutions

### dd_api_api warning regression bisection

After the initial dd::renorm fix, `dd_api_api` warning count went from 17 (pre-fix
baseline) to 25 (post-fix). Diffed pre/post-fix outputs to identify the regression:
`dd::maxpos` printed as `. e+00` (garbage). Traced to `quick_two_sum(DBL_MAX, lo,
lo)` overflowing to `(inf, -inf)`. Added `is_inf_cx(hi + lo)` guard. Back to 17
baseline.

### CR finding triage on PR #808

CR's Critical finding was about `dd::renorm` accepting `|hi| < |lo|` inputs (not
addressed by the maxpos guard). Switched primitive from `quick_two_sum` (requires
ordering) to `two_sum` (order-independent). Verified `dd(0.0, 1e100).renorm()` now
correctly canonicalizes to `(1e100, 0)`. Kept the inf-guard for the maxpos
boundary -- two_sum still overflows when `hi + lo` exceeds DBL_MAX.

### lint-pr-title constraint

`amannn/action-semantic-pull-request@v5` is configured with a fixed allowlist of
scopes. `floatcascade` was not in the list, and multi-scope `(qd,dd,floatcascade)`
was rejected. Looking at PRs #800/#803/#805/#806/#807: all use single scope.
Adopted the same pattern: `fix(qd): ... (qd/dd/floatcascade)` with the
parenthetical clarifying the wider impact.

## Performance Characteristics

- hfloat constexpr promotion: zero runtime cost (`std::is_constant_evaluated()`
  dispatch keeps runtime path on `std::ldexp`/`std::frexp`)
- dfixpnt FP-domain digit extraction: comparable to the pre-fix uint64
  materialization for `ndigits <= 19`; slightly more iterations for wider
  configurations but bounded by `ndigits` (constant per call)
- qd/dd/floatcascade `to_digits` renormalize: one renormalization per
  `to_string` call. `to_string` already does many high-precision operations,
  so the perf impact is negligible

## Test Results

| Target | gcc | clang | UBSan | notes |
|--------|-----|-------|-------|-------|
| hfloat full suite (9 targets) | PASS | PASS | n/a | new constexpr.cpp + 8 existing |
| dfixpnt full suite (8 targets) | PASS | PASS | n/a | new wide-instantiation tests |
| qd_api_api | PASS | PASS | PASS | 6 pre-existing subnormal warnings unchanged |
| dd_api_api | PASS | PASS | PASS | 17 pre-existing subnormal warnings unchanged |
| ddc_api_api / tdc_api_api / qdc_api_api | PASS | PASS | PASS | 0 warnings |

## Next Steps

- Wait for PR #808 CI to complete (~20 min from `d8ae6bd6`); merge if green
- Item 2 from issue #801 (raw-limb constructor canonicalization OR contract
  documentation for `iszero`/`sign` on non-canonical inputs) deferred -- separate
  scope, requires either a wider semantic change or a focused docs commit
- Subnormal-magnitude `to_digits` warnings (qd=6, dd=17) remain as a separate
  pre-existing edge case (algorithm struggles near 2^-1022); not addressed by
  this PR

## References

- Issue #732 (hfloat constexpr) -- closed by PR #806
- Issue #804 (dfixpnt wide overflow) -- closed by PR #807
- Issue #801 (qd/dd to_digits robustness) -- PR #808 open
- Epic #723 (constexpr support across Universal number systems) -- parent for
  hfloat work
- PR #800 (qd constexpr) -- introduced the NaN band-aid that #801 partially
  addresses, then reinstates
- PR #803 (dfixpnt constexpr) -- the constexpr promotion that surfaced #804
- `.github/workflows/sanitizers.yml` -- CI UBSan / ASan job that caught the
  qd loop drift regression

## Appendix: Local UBSan Quick Reference

```bash
# One-time configure
mkdir build_ubsan && cd build_ubsan
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=clang \
  -DCMAKE_CXX_COMPILER=clang++ -DUNIVERSAL_BUILD_CI=ON \
  -DUNIVERSAL_ENABLE_UBSAN=ON

# Build a target (with safety check)
pgrep -a make
cmake --build . --target qd_api_api -j4

# Run with CI-matching options
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ./static/highprecision/qd/qd_api_api

# Full ctest under UBSan
UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ctest --output-on-failure
```

See `docs/build/local-sanitizer.md` for the full workflow.
