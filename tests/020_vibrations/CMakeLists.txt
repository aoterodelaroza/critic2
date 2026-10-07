set(TESTS
  001_thermo_gamma           ## thermo keyword, acoustic modes at gamma left out (commensurate, unshifted mesh)
  002_thermo_cutoff          ## thermo keyword, default cutoff (only zero or negative modes left out)
  003_born_realspace         ## born keyword, real-space dipole-dipole sum (Lambda independence, NOREAL)
  004_load_dynmat            ## load_dynmat keyword, force constants from dynamical matrices at commensurate q
  005_fc2_supercell          ## load_fc2 keyword, the supercell as NEWCELL integers, PHONOPY order, a file; default file name
  006_fc2_write              ## write_fc2 keyword, compact and full formats, units of the calculator, round trips
  007_calcq_forms            ## calcq keyword, one-line form, lists, LINE, MESH, APPEND; print summary
  008_load_qe                ## load keyword, QE matdyn.modes, matdyn.eig, dyn; thermo on stored frequencies
  009_load_phonopy           ## load keyword, phonopy anime.ascii and mesh.yaml
  011_load_crystal           ## load keyword, CRYSTAL output; thermo on a gamma-only crystal
  012_load_gaussian          ## load keyword, Gaussian log and fchk; thermo on molecules (CO2 linear)
  013_clear                  ## clear keyword, born, freq, fc2 and everything
  014_create_displacements   ## create_displacements keyword, integers, NICE, SUPERCELL file; QE, aims and VASP inputs
  015_create_displacements_random ## create_displacements keyword, RANDOM, SEED, PLUSMINUS
  016_read_forces_fd         ## read_forces keyword, finite differences, template and list (exact round trip)
  017_read_forces_random     ## read_forces keyword, random displacements, NORIDGE, PAIRS, RIDGE, cross-validation; memory
  018_load_dynmat_skewed     ## load_dynmat keyword, a skewed supercell matrix gives the same frequencies (shortest images in a reduced cell)
  050_err_nocalc             ## errors: load_fc2 without a calculator
  051_err_fc2units           ## errors: load_fc2 with a unit word
  052_err_calcname           ## errors: unknown calculator
  053_err_thermo_meshrk      ## errors: thermo MESH and RKLENGTH
  054_err_thermo_shift       ## errors: thermo SHIFT without a mesh
  055_err_thermo_dos         ## errors: thermo DOS without a mesh
  056_err_thermo_commmesh    ## errors: thermo COMMENSURATE and MESH
  057_err_thermo_nofc2       ## errors: thermo MESH without force constants
  058_err_thermo_tmax        ## errors: thermo TMAX below TMIN
  059_err_thermo_nofreq      ## errors: thermo with no frequencies
  060_err_cd_seed            ## errors: create_displacements SEED without RANDOM
  061_err_cd_niceints        ## errors: create_displacements NICE and a supercell matrix
  062_err_cd_nosc            ## errors: create_displacements without a supercell
  063_err_rf_ridge           ## errors: read_forces RIDGE and NORIDGE
  064_err_rf_nodataset       ## errors: read_forces with a missing dataset
  065_err_load_unknown       ## errors: load of an unknown format
  066_err_load_fc2           ## errors: load of a FORCE_CONSTANTS file
  067_err_dynmat_nosc        ## errors: load_dynmat without the supercell
  068_err_unknown            ## errors: unknown vibrations keyword
  069_err_write_units        ## errors: write_fc2 with a unit word
  070_err_thermo_multiq      ## errors: thermo on stored frequencies at several q-points
  )

## the phonopy hdf5 reader needs the HDF5 library
if (USE_HDF5 AND HDF5_FOUND)
  set(TESTS ${TESTS}
    010_load_phonopy_hdf5    ## load keyword, phonopy mesh.hdf5 (needs hdf5)
    )
endif()

runtests(${TESTS})
