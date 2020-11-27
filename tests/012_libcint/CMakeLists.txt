if (LIBCINT_FOUND)

  set(TESTS
    001_molcalc_hf  ## MOLCALC keyword with HF option
    )

  runtests(${TESTS})
endif()
