set(SOURCES json_file_module.F90 json_kinds.F90 json_module.F90 json_parameters.F90
    json_string_utilities.F90 json_value_module.F90)

add_library(json-fortran OBJECT ${SOURCES})
