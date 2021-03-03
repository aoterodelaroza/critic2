set(TESTS
  001_write_xyz        ## write xyz file, with options
  002_write_xyz_border ## write xyz: border option
  003_write_xyz_nmer   ## write xyz: nmer option
  004_write_gjf        ## write gjf files
  005_write_cml        ## write cml files
  006_write_obj        ## write obj files
  007_write_ply        ## write ply files
  008_write_off        ## write off files
  009_write_qein       ## write QE inputs
  010_write_tess       ## write tessel input files
  011_write_cri        ## write critic2 input files
  012_write_vasp       ## write VASP input files (CONTCAR,POSCAR)
  013_write_abin       ## write abinit input files
  014_write_elk        ## write elk input files
  015_write_gau        ## write Gaussian input files (periodic)
  016_write_cif        ## write cif files, nosym option
  017_write_cif_long   ## write cif files, whole space group
  018_write_m          ## write octave/escher files
  019_write_db         ## write dcp/db files
  020_write_gulp       ## write gulp files, dreiding option
  021_write_lammps     ## write lammps files
  022_write_siesta     ## write siesta files (fdf,STRUCT_IN)
  023_write_dftb       ## write dftb+ files (gen,hsd)
  024_write_pyscf      ## write pyscf files, crystals and molecules
  025_write_d12        ## write d12 files
  026_write_res        ## write res files
  027_write_d12_withsym ## write d12 file, with symmetry
  028_write_fhi        ## write FHIaims input files
  )

runtests(${TESTS})
