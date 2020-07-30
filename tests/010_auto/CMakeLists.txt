set(TESTS
  001_auto_simple  ## a simple search for critical points (AUTO)
  002_auto_opts    ## AUTO options
  003_auto_chk     ## AUTO: chk generation
  004_auto_seedobj ## AUTO: seedobj option
  )

runtests(${TESTS})
