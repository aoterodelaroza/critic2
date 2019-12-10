## List all variables
## sakra at https://stackoverflow.com/questions/9298278/cmake-print-out-all-accessible-variables-in-a-script
macro(list_variables)
  get_cmake_property(_variableNames VARIABLES)
  list (SORT _variableNames)
  foreach (_variableName ${_variableNames})
      message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endmacro()
