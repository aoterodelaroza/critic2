## enable testing

## list of test bundles
set(TESTDIRS 
  001_io      ## input and output; arithmetics
  002_crystal ## reading a crystal structure (CRYSTAL)
  )

## unpack data for the tests (source)
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/zz_source/")
  message(STATUS "Unpacking the data for tests in ${CMAKE_CURRENT_SOURCE_DIR}/")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar "xJf" "zz_source.tar.xz" WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

## create symlink for the test data (build)
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/zz_source/")
  message(STATUS "Creating symlink for tests in ${CMAKE_CURRENT_BINARY_DIR}/")
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/zz_source" "${CMAKE_CURRENT_BINARY_DIR}/zz_source")
endif()

## define the macro that will run the tests
macro (runtests)
  get_filename_component(curdir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/ref")

  foreach(test ${ARGN})
    ## gather info about the test
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" check_list REGEX "^ *#.*check *: ")
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" delete_list REGEX "^ *#.*delete *:")

    ## gen test
    add_test(
      NAME "GEN:${curdir}/${test}" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      CONFIGURATIONS gen
      )
    set_tests_properties("GEN:${curdir}/${test}" PROPERTIES LABELS "gen")

    ## gen_copy test
    set(_check_list "")
    if (check_list)
      foreach (i IN LISTS check_list)
        string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
        string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
        list(GET llist 0 file)
	list(APPEND _check_list "${file}")
      endforeach()
      list(TRANSFORM _check_list PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
      string(REPLACE ";" " " check_string "${_check_list}")

      add_test(
        NAME "GEN-COPY-:${curdir}/${test}" 
        COMMAND sh -c "${CMAKE_COMMAND} -E copy ${check_string} ${CMAKE_CURRENT_SOURCE_DIR}/ref/"
        WORKING_DIRECTORY ""
        CONFIGURATIONS gen
        )
      set_tests_properties("GEN-COPY-:${curdir}/${test}" PROPERTIES 
        LABELS "gen" 
        DEPENDS "GEN:${curdir}/${test}"
        )
    endif()

    ## gen_clean test
    if (_check_list OR delete_list)
      set(_delete_list "")
      foreach (i IN LISTS delete_list)
        string(REGEX REPLACE "^ *#.*delete *:" "" i "${i}")
        string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(APPEND _delete_list "${llist}")
      endforeach()
      list(TRANSFORM _delete_list PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
      list(APPEND _delete_list "${_check_list}")
      string(REPLACE ";" " " delete_string "${_delete_list}")

      add_test(
        NAME "GEN-CLEAN-:${curdir}/${test}" 
        COMMAND sh -c "${CMAKE_COMMAND} -E remove ${delete_string}"
        WORKING_DIRECTORY ""
        CONFIGURATIONS gen
        )
      set_tests_properties("GEN-CLEAN-:${curdir}/${test}" PROPERTIES 
        LABELS "gen" 
        DEPENDS "GEN-COPY:${curdir}/${test}"
        )
    endif()

    ## run test
    add_test(
      NAME "-RUN-:${curdir}/${test}" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_BINARY_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      )

    ## check tests
    if (NUMDIFF_FOUND AND check_list)
      foreach (i IN LISTS check_list)
        string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
        string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
        list(GET llist 0 file)
        list(GET llist 1 numdiff_opts)
        add_test(
          NAME "CHECK:${curdir}/ref/${file}" 
          COMMAND sh -c "${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' ${numdiff_opts} ${CMAKE_CURRENT_SOURCE_DIR}/ref/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
          )
        set_tests_properties("CHECK:${curdir}/ref/${file}" PROPERTIES DEPENDS "-RUN-:${curdir}/${test}")
      endforeach()
    endif()
  endforeach()
endmacro()

## process all subdirectories
foreach(dir ${TESTDIRS})
  add_subdirectory(${dir})
endforeach()
