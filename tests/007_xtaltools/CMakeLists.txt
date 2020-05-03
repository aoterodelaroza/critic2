set(TESTS
  001_environ_simple  ## environ keyword, simple
  002_environ_opts    ## environ keyword, options
  003_environ_atomby  ## environ keyword, atom and by keywords
  004_ewald           ## ewald keyword
  005_ewald_pot       ## ewald keyword, potential
  006_packing         ## packing keyword
  007_packing_vdw     ## packing keyword, vdw option
  )

runtests(${TESTS})
