set(TESTS
  001_environ_simple  ## environ keyword, simple
  002_environ_opts    ## environ keyword, options
  003_environ_atomby  ## environ keyword, atom and by keywords
  )

runtests(${TESTS})
