set(TESTS )
if (LIBCINT_FOUND)
set(TESTS ${TESTS}
  001_sigmahole    ## sigmahole keyword
  )
endif()

set(TESTS
  ${TESTS}
  002_stm_current  ## stm, constant current  
  003_stm_height   ## stm, constant height  
  004_xdm_qe       ## xdm, qe keyword
  005_xdm_postg    ## xdm, postg keyword
  006_xdm_wfn      ## xdm, wavefunction calculation
  007_xdm_grid     ## xdm, grid keyword
  )

runtests(${TESTS})
