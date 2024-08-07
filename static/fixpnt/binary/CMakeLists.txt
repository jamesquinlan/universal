file(GLOB API_SRC        "api/*.cpp")
file(GLOB CONVERSION_SRC "conversion/*.cpp")
file(GLOB LOGIC_SRC      "logic/*.cpp")
file(GLOB ARITHMETIC_SRC "arithmetic/*.cpp")
file(GLOB COMPLEX_SRC    "complex/*.cpp")
file(GLOB MATH_SRC       "math/*.cpp")

compile_all("true" "fixpnt" "Number Systems/static/fixed-point/binary/fixpnt/api" "${API_SRC}")
compile_all("true" "fixpnt" "Number Systems/static/fixed-point/binary/fixpnt/conversion" "${CONVERSION_SRC}")
compile_all("true" "fixpnt" "Number Systems/static/fixed-point/binary/fixpnt/logic" "${LOGIC_SRC}")
compile_all("true" "fixpnt" "Number Systems/static/fixed-point/binary/fixpnt/arithmetic" "${ARITHMETIC_SRC}")
compile_all("true" "fixpnt" "Number Systems/static/fixed-point/binary/fixpnt/math" "${MATH_SRC}")

# Compiler specific environments
message(STATUS "CMAKE_CXX_COMPILER ID is -${CMAKE_CXX_COMPILER_ID}-")

# exclude XCode15
if(CMAKE_CXX_COMPILER_ID MATCHES  "AppleClang")
    message(STATUS "Ignoring complex tests for AppleClang")

else()
    message(STATUS "Adding complex test for all other non-AppleClang environments")
    compile_all("true" "fixpnt" "Number Systems/static/fixed-point/binary/fixpnt/complex" "${COMPLEX_SRC}")

endif()

