# Find the numdiff program. Sets:
#  NUMDIFF_FOUND - true if numdiff was found
#  NUMDIFF_EXE - the location of the numdiff program
#
# Part of this code adapted from the deal.II library:
#   https://github.com/dealii/dealii
include(FindPackageHandleStandardArgs)

if (DEFINED ENV{NUMDIFF_DIR})
  set(NUMDIFF_DIR "$ENV{NUMDIFF_DIR}")
endif()

find_program(NUMDIFF_EXE
  NAMES numdiff
  PATH_SUFFIXES bin
  HINTS ${NUMDIFF_DIR}
  )

find_package_handle_standard_args(NUMDIFF
  REQUIRED_VARS NUMDIFF_EXE
  )
mark_as_advanced(NUMDIFF_EXE)

## check that we have numdiff
if (NOT NUMDIFF_FOUND)
  set(NUMDIFF_DIR "${NUMDIFF_DIR}" CACHE STRING "Directory containing the numdiff program (for tests).")
  return()
endif()

## check that it executes correctly
execute_process(COMMAND ${NUMDIFF_EXE} "-v" TIMEOUT 4 OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE excode)
if (NOT "${excode}" STREQUAL "0")
  message(STATUS "numdiff did not run or timed out")
  unset(NUMDIFF_FOUND)
  return()
endif()

## check that it gives reasonable results
string(RANDOM _suffix)
set(_first_test_file_name "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/numdiff-test-${_suffix}-1.txt")
set(_second_test_file_name "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/numdiff-test-${_suffix}-2.txt")
file(WRITE "${_first_test_file_name}" "0.99999999998\n2.0\n1.0\n")
file(WRITE "${_second_test_file_name}" "1.00000000001\n2.0\n1.0\n")

execute_process(COMMAND ${NUMDIFF_EXE}
  "-r" "1.0e-8" "--" "${_first_test_file_name}" "${_second_test_file_name}"
  TIMEOUT 4 # seconds
  OUTPUT_QUIET
  ERROR_QUIET
  RESULT_VARIABLE _numdiff_tolerance_test_status
  )
file(REMOVE ${_first_test_file_name})
file(REMOVE ${_second_test_file_name})

if(NOT "${_numdiff_tolerance_test_status}" STREQUAL "0")
  message(STATUS "numdiff did not work correctly")
  unset(NUMDIFF_FOUND)
  return()
endif()


