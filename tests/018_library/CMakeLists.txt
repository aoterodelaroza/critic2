## C API test for the critic2 shared library (libcritic2). Unlike the
## other test bundles, this is not driven by a .cri input: it builds a
## small C program against the library and runs it. Only added when
## BUILD_LIBRARY is ON (see ../CMakeLists.txt).

add_executable(libcritic2_ctest libcritic2_ctest.c)
target_link_libraries(libcritic2_ctest PRIVATE libcritic2)
## for libcritic2.h, which lives next to the library sources
target_include_directories(libcritic2_ctest PRIVATE "${CMAKE_SOURCE_DIR}/src")
if (NOT WIN32)
  target_link_libraries(libcritic2_ctest PRIVATE m)
endif()

add_test(
  NAME "-LIB-:018_library/c_api"
  COMMAND libcritic2_ctest
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
set_property(TEST "-LIB-:018_library/c_api" PROPERTY LABELS "regression quick nodata")
set_property(TEST "-LIB-:018_library/c_api" PROPERTY ENVIRONMENT CRITIC_HOME=${CMAKE_SOURCE_DIR})
