## enable testing

## Directives recognized inside a test .cri file (lines beginning with "#"):
##   ## check:  FILE [numdiff_opts]   output file to compare against ref/FILE
##   ## delete: FILE ...              extra files to remove after generation
##   ## labels: LABEL ...             ctest labels (e.g. regression quick nodata)
##   ## grep:   FILE 'REGEX'          only keep lines of FILE matching REGEX
##                                    (extended regex) before comparing
##   ## sed:    FILE 'PROGRAM'        run "sed -e PROGRAM" on FILE before
##                                    comparing (line ranges, deletions, ...)
## grep:/sed: select the portion of FILE that must match for the test to pass;
## the same filter is applied when generating ref/FILE and when checking. A
## file may carry both (grep runs first, then sed). FILE must also appear in a
## check: directive. Single quotes inside REGEX/PROGRAM are not supported.

## test data package info
set(TEST_PACKAGE_NAME "zz_source.tar.xz")
set(TEST_PACKAGE_TARGET "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_PACKAGE_NAME}")

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
  013_misctools ## miscellaneous tools
  014_bisect    ## bisection keywords
  015_grdplot   ## gradient path tracing plots
  016_qtree     ## qtree keyword
  )

## Build the filter pipeline (grep and/or sed) for a checked file, given the
## per-file key. Returns the empty string in ${outvar} if the file has no
## grep:/sed: directive. The pipeline is meant to be fed via stdin, e.g.
##   cat FILE | <filter> > OUT
function (make_filter outvar key)
  set(_cmd "")
  if (DEFINED grepfor_${key})
    set(_cmd "grep -E -- '${grepfor_${key}}'")
  endif()
  if (DEFINED sedfor_${key})
    if (_cmd)
      set(_cmd "${_cmd} | sed -e '${sedfor_${key}}'")
    else()
      set(_cmd "sed -e '${sedfor_${key}}'")
    endif()
  endif()
  set(${outvar} "${_cmd}" PARENT_SCOPE)
endfunction()

## define the macro that will run the tests
macro (runtests)
  get_filename_component(curdir "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  foreach(test ${ARGN})
    ## gather info about the test (reset first so directives never leak from a
    ## previous test, e.g. via the .cmake-override path below)
    set(sed_list "")
    set(grep_list "")
    set(check_list "")
    set(delete_list "")
    set(label_list "")
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
      include("${CMAKE_CURRENT_SOURCE_DIR}/${test}.cmake")
    else()
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" sed_list REGEX "^ *#.*sed *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" grep_list REGEX "^ *#.*grep *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" check_list REGEX "^ *#.*check *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" delete_list REGEX "^ *#.*delete *: ")
      file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" label_list REGEX "^ *#.*labels *: ")
    endif()
    string(REGEX REPLACE "^ *#.*labels *:" "" label_list "${label_list}")

    ## build the per-file grep:/sed: filter map (keyed by test+file so filters
    ## never leak across tests)
    ## each grep:/sed: line is "FILE PROGRAM"; split on the first blank (the
    ## program may itself contain spaces, so a negated-space regex is unsafe in
    ## cmake) and strip one optional pair of surrounding single quotes
    foreach (i IN LISTS grep_list)
      string(REGEX REPLACE "^ *#.*grep *:" "" i "${i}")
      string(STRIP "${i}" i)
      string(FIND "${i}" " " sp)
      if (sp GREATER 0)
	string(SUBSTRING "${i}" 0 ${sp} gfile)
	string(SUBSTRING "${i}" ${sp} -1 gpat)
	string(STRIP "${gpat}" gpat)
	string(REGEX REPLACE "^'(.*)'$" "\\1" gpat "${gpat}")
	string(MAKE_C_IDENTIFIER "${test}_${gfile}" gkey)
	set(grepfor_${gkey} "${gpat}")
      endif()
    endforeach()
    foreach (i IN LISTS sed_list)
      string(REGEX REPLACE "^ *#.*sed *:" "" i "${i}")
      string(STRIP "${i}" i)
      string(FIND "${i}" " " sp)
      if (sp GREATER 0)
	string(SUBSTRING "${i}" 0 ${sp} sfile)
	string(SUBSTRING "${i}" ${sp} -1 sprog)
	string(STRIP "${sprog}" sprog)
	string(REGEX REPLACE "^'(.*)'$" "\\1" sprog "${sprog}")
	string(MAKE_C_IDENTIFIER "${test}_${sfile}" skey)
	set(sedfor_${skey} "${sprog}")
      endif()
    endforeach()

    ## skip tests if we have no package unless they are self-contained
    if (NOT "${TEST_DIR_FOUND}" AND NOT ${label_list} MATCHES nodata)
      continue()
    endif()

    ## gen test
    add_test(
      NAME "GEN:${curdir}/${test}"
      COMMAND critic2 -q -t -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      CONFIGURATIONS gen
      )
    set_property(TEST "GEN:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")
    set_property(TEST "GEN:${curdir}/${test}" PROPERTY ENVIRONMENT CRITIC_HOME=${CMAKE_SOURCE_DIR})

    ## gen_copy test: copy each checked file into ref/, or (when it has a
    ## grep:/sed: filter) write only the selected portion into ref/
    set(_check_list "")
    if (check_list)
      set(copy_cmds "")
      foreach (i IN LISTS check_list)
	string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
	string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(GET llist 0 file)
	string(STRIP "${file}" file)
	list(APPEND _check_list "${file}")
	string(MAKE_C_IDENTIFIER "${test}_${file}" fkey)
	make_filter(flt "${fkey}")
	if (flt)
	  set(one "cat ${CMAKE_CURRENT_SOURCE_DIR}/${file} | ${flt} > ${CMAKE_CURRENT_SOURCE_DIR}/ref/${file}")
	else()
	  set(one "${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_CURRENT_SOURCE_DIR}/ref/")
	endif()
	if (copy_cmds)
	  set(copy_cmds "${copy_cmds} && ${one}")
	else()
	  set(copy_cmds "${one}")
	endif()
      endforeach()
      list(TRANSFORM _check_list PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

      add_test(
	NAME "GEN-COPY-:${curdir}/${test}"
	COMMAND sh -c "${copy_cmds}"
	WORKING_DIRECTORY ""
	CONFIGURATIONS gen
	)
      set_property(TEST "GEN-COPY-:${curdir}/${test}" PROPERTY LABELS "gen ${label_list}")
      set_property(TEST "GEN-COPY-:${curdir}/${test}" PROPERTY DEPENDS "GEN:${curdir}/${test}")
      set_property(TEST "GEN-COPY-:${curdir}/${test}" PROPERTY ENVIRONMENT CRITIC_HOME=${CMAKE_SOURCE_DIR})
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
      set_property(TEST "GEN-CLEAN-:${curdir}/${test}" PROPERTY ENVIRONMENT CRITIC_HOME=${CMAKE_SOURCE_DIR})
    endif()

    ## run test
    add_test(
      NAME "-RUN-:${curdir}/${test}"
      COMMAND critic2 -q -t -l "${CMAKE_CURRENT_SOURCE_DIR}/${test}.cri" "${CMAKE_CURRENT_BINARY_DIR}/${test}.cro"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      )
    set_property(TEST "-RUN-:${curdir}/${test}" PROPERTY LABELS "${label_list}")
    set_property(TEST "-RUN-:${curdir}/${test}" PROPERTY ENVIRONMENT CRITIC_HOME=${CMAKE_SOURCE_DIR})

    ## check tests
    if (NUMDIFF_FOUND AND check_list)
      foreach (i IN LISTS check_list)
	string(REGEX REPLACE "^ *#.*check *:" "" i "${i}")
	string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)"  llist "${i}")
	list(LENGTH llist nn)
	list(GET llist 0 file)
	string(STRIP "${file}" file)
	set(numdiff_opts "")
	if ("${nn}" GREATER 1)
	  list(GET llist 1 numdiff_opts)
	endif()
	## if the file has a grep:/sed: filter, compare ref/FILE (already
	## filtered at generation) against the filtered build output
	string(MAKE_C_IDENTIFIER "${test}_${file}" fkey)
	make_filter(flt "${fkey}")
	if (flt)
	  set(check_cmd "cat ${CMAKE_CURRENT_BINARY_DIR}/${file} | ${flt} > ${CMAKE_CURRENT_BINARY_DIR}/${file}.cmp && ${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' -r 1e-7 ${numdiff_opts} ${CMAKE_CURRENT_SOURCE_DIR}/ref/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}.cmp")
	else()
	  set(check_cmd "${NUMDIFF_EXE} -s ' \t\n=,:;<>[](){}^' -r 1e-7 ${numdiff_opts} ${CMAKE_CURRENT_SOURCE_DIR}/ref/${file} ${CMAKE_CURRENT_BINARY_DIR}/${file}")
	endif()
	add_test(
	  NAME "CHECK:${curdir}/ref/${file}"
	  COMMAND sh -c "${check_cmd}"
	  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
	  )
	set_property(TEST "CHECK:${curdir}/ref/${file}" PROPERTY LABELS "${label_list}")
	set_property(TEST "CHECK:${curdir}/ref/${file}" PROPERTY DEPENDS "-RUN-:${curdir}/${test}")
	set_property(TEST "CHECK:${curdir}/ref/${file}" PROPERTY ENVIRONMENT CRITIC_HOME=${CMAKE_SOURCE_DIR})
      endforeach()
    endif()
  endforeach()
endmacro()

## unpack the package if necessary
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/zz_source/" AND EXISTS "${TEST_PACKAGE_TARGET}")
  message(STATUS "Unpacking the data for tests in ${CMAKE_CURRENT_SOURCE_DIR}/")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar "xJf" "${TEST_PACKAGE_TARGET}" WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

## create symlink for the test data (build)
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/zz_source/")
  message(STATUS "Creating symlink for tests in ${CMAKE_CURRENT_BINARY_DIR}/")
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/zz_source" "${CMAKE_CURRENT_BINARY_DIR}/zz_source")
endif()

## set the flag
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/zz_source")
  set(TEST_DIR_FOUND "YES")
endif()

## process all subdirectories
foreach(dir ${TESTDIRS})
  add_subdirectory(${dir})
endforeach()
