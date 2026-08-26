if(NOT DEFINED BENCHMARK OR NOT DEFINED INSTANCE)
  message(FATAL_ERROR "BENCHMARK and INSTANCE are required")
endif()

set(repetitions 10)
set(seed 123456789)

execute_process(
  COMMAND "${BENCHMARK}" "${INSTANCE}" "${repetitions}" "${seed}" "1"
  RESULT_VARIABLE first_result
  OUTPUT_VARIABLE first_output
  ERROR_VARIABLE first_error
)
if(NOT first_result EQUAL 0)
  message(FATAL_ERROR "First seeded run failed: ${first_error}")
endif()

execute_process(
  COMMAND "${BENCHMARK}" "${INSTANCE}" "${repetitions}" "${seed}" "4"
  RESULT_VARIABLE second_result
  OUTPUT_VARIABLE second_output
  ERROR_VARIABLE second_error
)
if(NOT second_result EQUAL 0)
  message(FATAL_ERROR "Second seeded run failed: ${second_error}")
endif()

# Runtime naturally varies; every solution-related field must be byte-identical.
string(REGEX REPLACE "elapsed_ms=[^\r\n]*" "" first_solution "${first_output}")
string(REGEX REPLACE "elapsed_ms=[^\r\n]*" "" second_solution "${second_output}")
string(REGEX REPLACE "max_threads=[^\r\n]*" "" first_solution "${first_solution}")
string(REGEX REPLACE "max_threads=[^\r\n]*" "" second_solution "${second_solution}")

if(NOT first_solution STREQUAL second_solution)
  message(FATAL_ERROR
    "Fixed seed changed between one and four threads.\n"
    "First:\n${first_output}\nSecond:\n${second_output}")
endif()
