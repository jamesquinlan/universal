// logic.cpp: functional tests for logic tests on classic cfloats
//
// Copyright (C) 2017-2021 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <universal/utility/directives.hpp>

#include <universal/number/cfloat/cfloat.hpp>
#include <universal/verification/test_status.hpp>

namespace sw::universal {

	template<typename TestType>
	int VerifyCfloatLogicEqual() {
		constexpr size_t max = TestType::nbits > 16 ? 16 : TestType::nbits;
		size_t NR_TEST_CASES = (size_t(1) << max);
		int nrOfFailedTestCases = 0;
		for (unsigned i = 0; i < NR_TEST_CASES; i++) {
			TestType a;
			a.setbits(i);
			for (unsigned j = 0; j < NR_TEST_CASES; j++) {
				TestType b;
				b.setbits(j);

				// set the golden reference

					// initially, we thought this would be the same behavior as IEEE floats
					// ref = double(a) == double(b);
					// but we have found that some compilers (MSVC) take liberty with NaN
					// \fp:fast		floating point model set to fast
					//	NaN == NaN  : IEEE = true    cfloat = true  because we have unique encodings for +-NaN
					//	NaN == real : IEEE = true    cfloat = false
					// \fp:strict	floating point model set to strict
					//	NaN == NaN  : IEEE = false    cfloat = true
					//	NaN == real : IEEE = false    cfloat = false
					// and thus we can't rely on IEEE float as reference

				// instead, use the bit pattern as reference
				bool ref = (i == j);

				bool result = (a == b);
				if (ref != result) {
					nrOfFailedTestCases++;
					std::cout << a << " != " << b << " fails: reference is " << ref << " actual is " << result << std::endl;
				}
			}
		}
		return nrOfFailedTestCases;
	}

	template<typename TestType>
	int VerifyCfloatLogicNotEqual() {
		constexpr size_t max = TestType::nbits > 16 ? 16 : TestType::nbits;
		size_t NR_TEST_CASES = (size_t(1) << max);
		int nrOfFailedTestCases = 0;
		for (unsigned i = 0; i < NR_TEST_CASES; i++) {
			TestType a;
			a.setbits(i);
			for (unsigned j = 0; j < NR_TEST_CASES; j++) {
				TestType b;
				b.setbits(j);

				// set the golden reference

					// initially, we thought this would be the same behavior as IEEE floats
					// ref = double(a) == double(b);
					// but we have found that some compilers (MSVC) take liberty with NaN
					// \fp:fast		floating point model set to fast
					//	NaN == NaN  : IEEE = true    cfloat = true  because we have unique encodings for +-NaN
					//	NaN == real : IEEE = true    cfloat = false
					// \fp:strict	floating point model set to strict
					//	NaN == NaN  : IEEE = false    cfloat = true
					//	NaN == real : IEEE = false    cfloat = false
					// and thus we can't rely on IEEE float as reference

				// instead, use the bit pattern as reference
				bool ref = (i != j);

				bool result = (a != b);
				if (ref != result) {
					nrOfFailedTestCases++;
					std::cout << a << " != " << b << " fails: reference is " << ref << " actual is " << result << std::endl;
				}
			}
		}
		return nrOfFailedTestCases;
	}

	template<typename TestType>
	int VerifyCfloatLogicLessThan() {
		constexpr size_t max = TestType::nbits > 16 ? 16 : TestType::nbits;
		size_t NR_TEST_CASES = (size_t(1) << max);
		int nrOfFailedTestCases = 0;
		for (unsigned i = 0; i < NR_TEST_CASES; i++) {
			TestType a;
			a.setbits(i);
			for (unsigned j = 0; j < NR_TEST_CASES; j++) {
				TestType b;
				b.setbits(j);

				// set the golden reference

					// initially, we thought this would be the same behavior as IEEE floats
					// ref = double(a) == double(b);
					// but we have found that some compilers (MSVC) take liberty with NaN
					// \fp:fast		floating point model set to fast
					//	NaN == NaN  : IEEE = true    cfloat = true  because we have unique encodings for +-NaN
					//	NaN == real : IEEE = true    cfloat = false
					// \fp:strict	floating point model set to strict
					//	NaN == NaN  : IEEE = false    cfloat = true
					//	NaN == real : IEEE = false    cfloat = false
					// and thus we can't rely on IEEE float as reference

				// since this function is only useful for small cfloat<>s, we can depend on the double conversion
				bool ref = (double(a) < double(b));

				bool result = (a < b);
				if (ref != result) {
					nrOfFailedTestCases++;
					std::cout << a << " < " << b << " fails: reference is " << ref << " actual is " << result << std::endl;
				}
			}
		}
		return nrOfFailedTestCases;
	}

	template<typename TestType>
	int VerifyCfloatLogicLessOrEqualThan() {
		constexpr size_t max = TestType::nbits > 16 ? 16 : TestType::nbits;
		size_t NR_TEST_CASES = (size_t(1) << max);
		int nrOfFailedTestCases = 0;
		for (unsigned i = 0; i < NR_TEST_CASES; i++) {
			TestType a;
			a.setbits(i);
			for (unsigned j = 0; j < NR_TEST_CASES; j++) {
				TestType b;
				b.setbits(j);

				// set the golden reference

					// initially, we thought this would be the same behavior as IEEE floats
					// ref = double(a) == double(b);
					// but we have found that some compilers (MSVC) take liberty with NaN
					// \fp:fast		floating point model set to fast
					//	NaN == NaN  : IEEE = true    cfloat = true  because we have unique encodings for +-NaN
					//	NaN == real : IEEE = true    cfloat = false
					// \fp:strict	floating point model set to strict
					//	NaN == NaN  : IEEE = false    cfloat = true
					//	NaN == real : IEEE = false    cfloat = false
					// and thus we can't rely on IEEE float as reference

				// since this function is only useful for small cfloat<>s, we can depend on the double conversion
				bool ref = (double(a) <= double(b));

				bool result = (a <= b);
				if (ref != result) {
					nrOfFailedTestCases++;
					std::cout << a << " < " << b << " fails: reference is " << ref << " actual is " << result << std::endl;
				}
			}
		}
		return nrOfFailedTestCases;
	}

	template<typename TestType>
	int VerifyCfloatLogicGreaterThan() {
		constexpr size_t max = TestType::nbits > 16 ? 16 : TestType::nbits;
		size_t NR_TEST_CASES = (size_t(1) << max);
		int nrOfFailedTestCases = 0;
		for (unsigned i = 0; i < NR_TEST_CASES; i++) {
			TestType a;
			a.setbits(i);
			for (unsigned j = 0; j < NR_TEST_CASES; j++) {
				TestType b;
				b.setbits(j);

				// set the golden reference

					// initially, we thought this would be the same behavior as IEEE floats
					// ref = double(a) == double(b);
					// but we have found that some compilers (MSVC) take liberty with NaN
					// \fp:fast		floating point model set to fast
					//	NaN == NaN  : IEEE = true    cfloat = true  because we have unique encodings for +-NaN
					//	NaN == real : IEEE = true    cfloat = false
					// \fp:strict	floating point model set to strict
					//	NaN == NaN  : IEEE = false    cfloat = true
					//	NaN == real : IEEE = false    cfloat = false
					// and thus we can't rely on IEEE float as reference

				// since this function is only useful for small cfloat<>s, we can depend on the double conversion
				bool ref = (double(a) > double(b));

				bool result = (a > b);
				if (ref != result) {
					nrOfFailedTestCases++;
					std::cout << a << " < " << b << " fails: reference is " << ref << " actual is " << result << std::endl;
				}
			}
		}
		return nrOfFailedTestCases;
	}

	template<typename TestType>
	int VerifyCfloatLogicGreaterOrEqualThan() {
		constexpr size_t max = TestType::nbits > 16 ? 16 : TestType::nbits;
		size_t NR_TEST_CASES = (size_t(1) << max);
		int nrOfFailedTestCases = 0;
		for (unsigned i = 0; i < NR_TEST_CASES; i++) {
			TestType a;
			a.setbits(i);
			for (unsigned j = 0; j < NR_TEST_CASES; j++) {
				TestType b;
				b.setbits(j);

				// set the golden reference

					// initially, we thought this would be the same behavior as IEEE floats
					// ref = double(a) == double(b);
					// but we have found that some compilers (MSVC) take liberty with NaN
					// \fp:fast		floating point model set to fast
					//	NaN == NaN  : IEEE = true    cfloat = true  because we have unique encodings for +-NaN
					//	NaN == real : IEEE = true    cfloat = false
					// \fp:strict	floating point model set to strict
					//	NaN == NaN  : IEEE = false    cfloat = true
					//	NaN == real : IEEE = false    cfloat = false
					// and thus we can't rely on IEEE float as reference

				// since this function is only useful for small cfloat<>s, we can depend on the double conversion
				bool ref = (double(a) >= double(b));

				bool result = (a >= b);
				if (ref != result) {
					nrOfFailedTestCases++;
					std::cout << a << " < " << b << " fails: reference is " << ref << " actual is " << result << std::endl;
				}
			}
		}
		return nrOfFailedTestCases;
	}

} // namespace sw::universal

// Regression testing guards: typically set by the cmake configuration, but MANUAL_TESTING is an override
#define MANUAL_TESTING 0
// REGRESSION_LEVEL_OVERRIDE is set by the cmake file to drive a specific regression intensity
// It is the responsibility of the regression test to organize the tests in a quartile progression.
//#undef REGRESSION_LEVEL_OVERRIDE
#ifndef REGRESSION_LEVEL_OVERRIDE
#undef REGRESSION_LEVEL_1
#undef REGRESSION_LEVEL_2
#undef REGRESSION_LEVEL_3
#undef REGRESSION_LEVEL_4
#define REGRESSION_LEVEL_1 1
#define REGRESSION_LEVEL_2 1
#define REGRESSION_LEVEL_3 0
#define REGRESSION_LEVEL_4 0
#endif

int main()
try {
	using namespace sw::universal;

	std::string test_suite  = "cfloat<> logic operator validation";
	std::string test_tag    = "logic";
//	bool reportTestCases    = false;
	int nrOfFailedTestCases = 0;

#if MANUAL_TESTING

	float b = 1.0f;

	std::cout << "correct pattern as defined by IEEE-754 is:\nF F F F F T   <--- correct pattern\n";
	{
		float test = std::numeric_limits<float>::quiet_NaN();

		printf("%c %c %c %c %c %c\n",
			(test < test) ? 'T' : 'F',
			(test <= test) ? 'T' : 'F',
			(test == test) ? 'T' : 'F',
			(test > test) ? 'T' : 'F',
			(test >= test) ? 'T' : 'F',
			(test != test) ? 'T' : 'F'
		);

		printf("%c %c %c %c %c %c\n",
			(test < b) ? 'T' : 'F',
			(test <= b) ? 'T' : 'F',
			(test == b) ? 'T' : 'F',
			(test > b) ? 'T' : 'F',
			(test >= b) ? 'T' : 'F',
			(test != b) ? 'T' : 'F'
		);
	}

	{
		float test = std::numeric_limits<float>::signaling_NaN();

		printf("%c %c %c %c %c %c\n",
			(test < test) ? 'T' : 'F',
			(test <= test) ? 'T' : 'F',
			(test == test) ? 'T' : 'F',
			(test > test) ? 'T' : 'F',
			(test >= test) ? 'T' : 'F',
			(test != test) ? 'T' : 'F'
		);

		printf("%c %c %c %c %c %c\n",
			(test < b) ? 'T' : 'F',
			(test <= b) ? 'T' : 'F',
			(test == b) ? 'T' : 'F',
			(test > b) ? 'T' : 'F',
			(test >= b) ? 'T' : 'F',
			(test != b) ? 'T' : 'F'
		);
	}

	std::cout << "cfloat\n";
	{
		cfloat<16, 5> test = std::numeric_limits<cfloat<16, 5>>::quiet_NaN();

		printf("%c %c %c %c %c %c\n",
			(test < test) ? 'T' : 'F',
			(test <= test) ? 'T' : 'F',
			(test == test) ? 'T' : 'F',
			(test > test) ? 'T' : 'F',
			(test >= test) ? 'T' : 'F',
			(test != test) ? 'T' : 'F'
		);

		printf("%c %c %c %c %c %c\n",
			(test < b) ? 'T' : 'F',
			(test <= b) ? 'T' : 'F',
			(test == b) ? 'T' : 'F',
			(test > b) ? 'T' : 'F',
			(test >= b) ? 'T' : 'F',
			(test != b) ? 'T' : 'F'
		);
	}

	{
		cfloat<16, 5> test = std::numeric_limits<cfloat<16, 5>>::signaling_NaN();

		printf("%c %c %c %c %c %c\n",
			(test < test) ? 'T' : 'F',
			(test <= test) ? 'T' : 'F',
			(test == test) ? 'T' : 'F',
			(test > test) ? 'T' : 'F',
			(test >= test) ? 'T' : 'F',
			(test != test) ? 'T' : 'F'
		);

		printf("%c %c %c %c %c %c\n",
			(test < b) ? 'T' : 'F',
			(test <= b) ? 'T' : 'F',
			(test == b) ? 'T' : 'F',
			(test > b) ? 'T' : 'F',
			(test >= b) ? 'T' : 'F',
			(test != b) ? 'T' : 'F'
		);
	}

#else

	cfloat<16, 5> a;

	std::cout << "Logic: operator==()\n";

	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat< 4, 2> >(), "cfloat< 4,2>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat< 5, 2> >(), "cfloat< 5,2>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat< 6, 2> >(), "cfloat< 6,2>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat< 7, 2> >(), "cfloat< 7,2>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat< 8, 2> >(), "cfloat< 8,2>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat< 9, 2> >(), "cfloat< 9,2>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat<10, 2> >(), "cfloat<10,2>", "==");

	if (!(a == 0)) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> == 0", "== int literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> == 0", "== int literal");
	}
	if (!(a == 0.0f)) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> == 0.0f", "== float literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> == 0.0f", "== float literal");
	}
	if (!(a == 0.0)) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> == 0.0", "== double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> == 0.0", "== double literal");
	}
#if LONG_DOUBLE_SUPPORT
	if (!(a == 0.0l)) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> == 0.0l", "== long double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> == 0.0l", "== long double literal");
	}
#endif

	std::cout << "Logic: operator!=()\n";
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat< 4, 1, uint8_t, true, true, false> >(), "cfloat< 4,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat< 5, 1, uint8_t, true, true, false> >(), "cfloat< 5,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat< 6, 1, uint8_t, true, true, false> >(), "cfloat< 6,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat< 7, 1, uint8_t, true, true, false> >(), "cfloat< 7,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat< 8, 1, uint8_t, true, true, false> >(), "cfloat< 8,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat< 9, 1, uint8_t, true, true, false> >(), "cfloat< 9,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat<10, 1, uint8_t, true, true, false> >(), "cfloat<10,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat<12, 1, uint8_t, true, true, false> >(), "cfloat<12,1>", "!=");

	if (a != 0) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> != 0", "!= int literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> != 0", "!= int literal");
	}
	if (a != 0.0f) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> != 0.0f", "!= float literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> != 0.0f", "!= float literal");
	}
	if (a != 0.0) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> != 0.0", "!= double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> != 0.0", "!= double literal");
	}
#if LONG_DOUBLE_SUPPORT
	if (a != 0.0l) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> != 0.0l", "!= long double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> != 0.0l", "!= long double literal");
	}
#endif

	std::cout << "Logic: operator<()\n";
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 4, 1, uint8_t, true, true, false> >(), "cfloat< 4,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 5, 1, uint8_t, true, true, false> >(), "cfloat< 5,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 6, 1, uint8_t, true, true, false> >(), "cfloat< 6,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 7, 1, uint8_t, true, true, false> >(), "cfloat< 7,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 8, 1, uint8_t, true, true, false> >(), "cfloat< 8,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 9, 1, uint8_t, true, true, false> >(), "cfloat< 9,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat<10, 1, uint8_t, true, true, false> >(), "cfloat<10,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat<12, 1, uint8_t, true, true, false> >(), "cfloat<12,1>", "<");

	if (a < 0) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0", "< int literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0", "< int literal");
	}
	if (a < 0.0f) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0.0f", "< float literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0.0f", "== float literal");
	}
	if (a < 0.0) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0.0", "< double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0.0", "< double literal");
	}
#if LONG_DOUBLE_SUPPORT
	if (a < 0.0l) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0.0l", "< long double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0.0l", "== long double literal");
	}
#endif

	std::cout << "Logic: operator<=()\n";
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 4, 1, uint8_t, true, true, false> >(), "cfloat< 4,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 5, 1, uint8_t, true, true, false> >(), "cfloat< 5,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 6, 1, uint8_t, true, true, false> >(), "cfloat< 6,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 7, 1, uint8_t, true, true, false> >(), "cfloat< 7,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 8, 1, uint8_t, true, true, false> >(), "cfloat< 8,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat< 9, 1, uint8_t, true, true, false> >(), "cfloat< 9,1>", "<");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicLessThan< cfloat<10, 1, uint8_t, true, true, false> >(), "cfloat<10,1>", "<");

	if (a < 0) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0", "< int literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0", "< int literal");
	}
	if (a < 0.0f) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0.0f", "< float literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0.0f", "== float literal");
	}
	if (a < 0.0) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0.0", "< double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0.0", "< double literal");
	}
#if LONG_DOUBLE_SUPPORT
	if (a < 0.0l) {
		nrOfFailedTestCases += ReportTestResult(1, "cfloat<16,1> < 0.0l", "< long double literal");
	}
	else {
		ReportTestResult(0, "cfloat<16,1> < 0.0l", "== long double literal");
	}
#endif


#if STRESS_TESTING
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat<12, 1> >(), "cfloat<12,1>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat<14, 1> >(), "cfloat<14,1>", "==");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicEqual< cfloat<16, 1> >(), "cfloat<16,1>", "==");

	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat<12, 1> >(), "cfloat<12,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat<14, 1> >(), "cfloat<14,1>", "!=");
	nrOfFailedTestCases += ReportTestResult(VerifyCfloatLogicNotEqual< cfloat<16, 1> >(), "cfloat<16,1>", "!=");

#endif  // STRESS_TESTING

#endif  // MANUAL_TESTING

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << "Caught ad-hoc exception: " << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::universal_arithmetic_exception& err) {
	std::cerr << "Caught unexpected universal arithmetic exception : " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::universal_internal_exception& err) {
	std::cerr << "Caught unexpected universal internal exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const std::runtime_error& err) {
	std::cerr << "Caught runtime exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
