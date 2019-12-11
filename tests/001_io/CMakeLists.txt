set(TESTS
  001_simple_arithmetics
  )

get_filename_component(curdir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

foreach(test ${TESTS})
  ## run test
  add_test(NAME "${curdir}/${test}-run" COMMAND critic2 -q "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_BINARY_DIR}/${test}.cro")

  ## result test
  if (NUMDIFF_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
    include("${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
    list(LENGTH _numdiff_opts _nfile)
    math(EXPR _nfile "${_nfile} - 1")
    foreach(i RANGE "${_nfile}")
      list(GET _numdiff_opts ${i} _opts)
      list(GET _list_files ${i} _file)
      add_test(NAME "${curdir}/${test}-result" COMMAND "${NUMDIFF_EXE}" "-s='\\n\ ,()\\t'" "-q" "${_opts}"
      	"${CMAKE_CURRENT_SOURCE_DIR}/${_file}" "${CMAKE_CURRENT_BINARY_DIR}/${_file}")
      set_tests_properties("${curdir}/${test}-result" PROPERTIES DEPENDS "${curdir}/${test}-run")
    endforeach()
  endif()
endforeach()
