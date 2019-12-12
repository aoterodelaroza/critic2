## Find the readline library. Sets:
## READLINE_FOUND - whether the following items were found
## READLINE_LIBRARY - the readline library
## READLINE_INCLUDE_DIRS - the directory containing the header files

include(FindPackageHandleStandardArgs)

if (DEFINED ENV{READLINE_DIR})
  set(READLINE_DIR "$ENV{READLINE_DIR}")
endif()

find_library(READLINE_LIBRARY
    NAMES readline
    PATH_SUFFIXES lib readline
    HINTS ${READLINE_DIR}
)

find_path(READLINE_INCLUDE_DIRS
  NAMES readline/readline.h
  PATH_SUFFIXES include
  HINTS ${READLINE_DIR}
  )

find_package_handle_standard_args(READLINE
  REQUIRED_VARS READLINE_LIBRARY READLINE_INCLUDE_DIRS
  )

mark_as_advanced(READLINE_LIBRARY READLINE_INCLUDE_DIRS)
if (NOT READLINE_FOUND)
  set(READLINE_DIR "" CACHE STRING "Directory containing the readline library.")
endif()
