set(TESTS
  001_thermo_gamma           ## thermo keyword, acoustic modes at gamma left out (commensurate, unshifted mesh)
  002_thermo_cutoff          ## thermo keyword, default cutoff (only zero or negative modes left out)
  003_born_realspace         ## born keyword, real-space dipole-dipole sum (Lambda independence, NOREAL)
  )

runtests(${TESTS})
