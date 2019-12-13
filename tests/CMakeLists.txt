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
    ## gather info about the test
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" check_list REGEX "^ *#.*check *: ")
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" delete_list REGEX "^ *#.*delete *:")

    ## gen test
    add_test(
      NAME "${curdir}/${test}-gen" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      CONFIGURATIONS gen
      )
    set_tests_properties("${curdir}/${test}-gen" PROPERTIES LABELS "gen")

    ## gen-clean test
    if (delete_list)
      set(_delete_list "")
      foreach (i IN LISTS delete_list)
        string(REGEX REPLACE "^ *#.*delete *:" "" i "${i}")
        string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(APPEND _delete_list "${llist}")
      endforeach()
      list(TRANSFORM _delete_list PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
      string(REPLACE ";" " " delete_string "${_delete_list}")
      message("delete: ${delete_string}")

      add_test(
        NAME "${curdir}/${test}-gen_clean" 
        COMMAND sh -c "${CMAKE_COMMAND} -E remove ${delete_string}"
        WORKING_DIRECTORY ""
        CONFIGURATIONS gen
        )
      set_tests_properties("${curdir}/${test}-gen_clean" PROPERTIES 
        LABELS "gen" 
        DEPENDS "${curdir}/${test}-gen"
        )
    endif()

    ## run test
    add_test(
      NAME "${curdir}/${test}-run" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_BINARY_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
      )

    ## check tests
    if (NUMDIFF_FOUND AND check_list)
      foreach (i IN LISTS check_list)
        string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
        string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
        list(GET llist 0 file)
        list(GET llist 1 numdiff_opts)
        add_test(
          NAME "${curdir}/${test}-check_${file}" 
          COMMAND sh -c "${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' -q ${numdiff_opts} ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BUILD_DIR}"
          )
        set_tests_properties("${curdir}/${test}-check_${file}" PROPERTIES DEPENDS "${curdir}/${test}-run")
      endforeach()
    endif()
  endforeach()
endmacro()

## process all subdirectories
foreach(dir ${TESTDIRS})
  add_subdirectory(${dir})
endforeach()
