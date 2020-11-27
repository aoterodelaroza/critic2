## enable testing

## test data package info
set(TEST_PACKAGE_NAME "zz_source.tar.xz")
set(TEST_PACKAGE_TARGET "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PACKAGE_NAME}")
#set(TEST_URL 
#  http://schooner.chem.dal.ca/downloads/critic2/critic2_test_package.tar.xz
#  )
set(TEST_SHA1_HASH e86d633baa2624f39f4e881b089fa19fcd057f17)

## list of test bundles
set(TESTDIRS 
  001_io        ## input and output; arithmetics
  002_crystal   ## reading a crystal structure (CRYSTAL)
  003_molecule  ## reading a molecular structure (MOLECULE)
  004_load      ## loading a scalar field (LOAD)
  005_plot      ## plotting keywords (POINT,LINE,PLANE,CUBE)
  006_write     ## structure writers (WRITE)
  007_xtaltools ## crystallography tools (NEWCELL,COMPARE,...)
  008_libxc     ## tests for the libxc library
  009_intgrid   ## grid integration
  010_auto      ## critical point search (AUTO,...)
  011_mesh      ## mesh integrations (MOLCALC,...)
  012_libcint   ## tests for the libcint library
  )

## define the macro that will run the tests
macro (runtests)
  get_filename_component(curdir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  foreach(test ${ARGN})
    ## gather info about the test
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
      include("${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
    else()
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" check_list REGEX "^ *#.*check *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" delete_list REGEX "^ *#.*delete *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" label_list REGEX "^ *#.*labels *: ")
    endif()
    string(REGEX REPLACE "^ *#.*labels *:" "" label_list "${label_list}")

    ## gen test
    add_test(
      NAME "GEN:${curdir}/${test}" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      CONFIGURATIONS gen
      )
    set_property(TEST "GEN:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")

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
      set_property(TEST "GEN-COPY-:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")
      set_property(TEST "GEN-COPY-:${curdir}/${test}" PROPERTY DEPENDS "GEN:${curdir}/${test}")
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
      set_property(TEST "GEN-CLEAN-:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")
      set_property(TEST "GEN-CLEAN-:${curdir}/${test}" PROPERTY DEPENDS "GEN-COPY:${curdir}/${test}")
    endif()

    ## run test
    add_test(
      NAME "-RUN-:${curdir}/${test}" 
      COMMAND critic2 -q -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_BINARY_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      )
    set_property(TEST "-RUN-:${curdir}/${test}" PROPERTY LABELS "${label_list}")

    ## check tests
    if (NUMDIFF_FOUND AND check_list)
      foreach (i IN LISTS check_list)
        string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
        string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(LENGTH llist nn)
        list(GET llist 0 file)
	if ("${nn}" GREATER 1)
          list(GET llist 1 numdiff_opts)
	endif()
        add_test(
          NAME "CHECK:${curdir}/ref/${file}" 
          COMMAND sh -c "${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' ${numdiff_opts} ${CMAKE_CURRENT_SOURCE_DIR}/ref/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}"
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
          )
	set_property(TEST "CHECK:${curdir}/ref/${file}" PROPERTY LABELS "${label_list}")
	set_property(TEST "CHECK:${curdir}/ref/${file}" PROPERTY DEPENDS "-RUN-:${curdir}/${test}")
      endforeach()
    endif()
  endforeach()
endmacro()

## if the file exists, calculate the hash and verify
if (EXISTS "${TEST_PACKAGE_TARGET}")
  file(SHA1 "${TEST_PACKAGE_TARGET}" _hash)
  if ("${_hash}" STREQUAL "${TEST_SHA1_HASH}")
    set(TEST_PACKAGE_FOUND "YES")
  endif()
endif()

## fetch the test data package, if necessary - verify the hash
if(NOT "${TEST_PACKAGE_FOUND}")
  message(FATAL_ERROR "A package file with the correct sha1sum could not be found.")
  # foreach(url ${TEST_URL})
  #   if (EXISTS "${TEST_PACKAGE_TARGET}")
  #     file(REMOVE "${TEST_PACKAGE_TARGET}")
  #   endif()
  #   message(STATUS "package: ${TEST_PACKAGE_TARGET}")
  #   file(DOWNLOAD "${url}" "${TEST_PACKAGE_TARGET}"
  #     TIMEOUT 60 
  #     TLS_VERIFY ON
  #     SHOW_PROGRESS
  #     STATUS status
  #     )
  #   list(GET status 0 error_code)
  #   if(NOT "${error_code}")
  #     file(SHA1 "${TEST_PACKAGE_TARGET}" _hash)
  #     if ("${_hash}" STREQUAL "${TEST_SHA1_HASH}")
  # 	set(TEST_PACKAGE_FOUND "YES")
  # 	break()
  #     endif()
  #   endif()
  # endforeach()
endif()

## unpack data for the tests (source)
if ("${TEST_PACKAGE_FOUND}")
  if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/zz_source/")
    message(STATUS "Unpacking the data for tests in ${CMAKE_CURRENT_SOURCE_DIR}/")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar "xJf" "${TEST_PACKAGE_TARGET}" WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  ## create symlink for the test data (build)
  if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/zz_source/")
    message(STATUS "Creating symlink for tests in ${CMAKE_CURRENT_BINARY_DIR}/")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/zz_source" "${CMAKE_CURRENT_BINARY_DIR}/zz_source")
  endif()

  ## process all subdirectories
  foreach(dir ${TESTDIRS})
    add_subdirectory(${dir})
  endforeach()
else()
  message(STATUS "Skipping tests because test data package could not be found or downloaded.")
endif()
