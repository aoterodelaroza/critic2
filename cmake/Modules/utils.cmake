## List all variables
## sakra at https://stackoverflow.com/questions/9298278/cmake-print-out-all-accessible-variables-in-a-script
macro(list_variables)
  get_cmake_property(_variableNames VARIABLES)
  list (SORT _variableNames)
  foreach (_variableName ${_variableNames})
      message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endmacro()

## Keep a bundled/vendored numerical library optimized even in Debug
## builds by compiling it with -O2 in every configuration. lang is the
## language of the target's sources (C, Fortran, ...).
function(c2_optimize_vendored tgt lang)
  if (CMAKE_${lang}_COMPILER_ID MATCHES "GNU|Clang|Intel")
    target_compile_options(${tgt} PRIVATE -O2)
  endif()
endfunction()
