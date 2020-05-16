if (LIBXC_FOUND)

  set(TESTS
    001_utils         ## libxc keyword, utils
    002_molcalc_libxc ## molcalc using libxc and a wfn file
    )

  runtests(${TESTS})
endif()
