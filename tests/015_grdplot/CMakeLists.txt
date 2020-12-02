set(TESTS
  001_fluxprint_basic  ## fluxprint keyword, basic
  002_fluxprint_opts   ## fluxprint keyword, options
  003_grdvec_crys      ## grdvec, in a crystal
  004_grdvec_mols      ## grdvec, in molecules
  )

runtests(${TESTS})
