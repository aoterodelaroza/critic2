cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(critic2 LANGUAGES C CXX Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(critic2_URL "https://aoterodelaroza.github.io/critic2/")
set(critic2_EMAIL "aoterodelaroza@gmail.com")
set(critic2_LICENSE "GNU/GPL version 3")
set(critic2_DESCRIPTION "Analysis of quantum mechanical calculations in molecules and solids")

## list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
option(ENABLE_GUI "Enables the critic2 GUI" OFF)
include(GNUInstallDirs)

find_package(OpenMP MODULE COMPONENTS Fortran)

add_subdirectory(src)
add_subdirectory(dat)
if (ENABLE_GUI)
  add_subdirectory(src-gui)
endif()

## get_cmake_property(_variableNames VARIABLES)
## list (SORT _variableNames)
## foreach (_variableName ${_variableNames})
##     message(STATUS "${_variableName}=${${_variableName}}")
## endforeach()
