## enable testing

## list of test bundles
set(TESTDIRS 
  001_io      ## input and output; arithmetics
  002_crystal ## reading a crystal structure (CRYSTAL)
  )

## define the macro that will run the tests
macro (runtests)
  get_filename_component(curdir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  foreach(test ${ARGN})
    ## gen test
    add_test(
      NAME "${curdir}/${test}-gen" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      CONFIGURATIONS gen
      )
    set_tests_properties("${curdir}/${test}-gen" PROPERTIES LABELS "gen")

    ## run test
    add_test(
      NAME "${curdir}/${test}-run" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_BINARY_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
      )

    ## result tests
    if (NUMDIFF_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
      include("${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
      list(LENGTH _numdiff_opts _nfile)
      math(EXPR _nfile "${_nfile} - 1")
      if (_list_files)
	foreach(i RANGE "${_nfile}")
	  list(GET _numdiff_opts ${i} _opts)
	  list(GET _list_files ${i} _file)
	  add_test(
	    NAME "${curdir}/${test}-result" 
	    COMMAND sh -c "${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' -q ${_opts} ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}"
	    WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
	    )
	  set_tests_properties("${curdir}/${test}-result" PROPERTIES DEPENDS "${curdir}/${test}-run")
	endforeach()
      endif()
    endif()
  endforeach()
endmacro()

## process all subdirectories
foreach(dir ${TESTDIRS})
  add_subdirectory(${dir})
endforeach()
