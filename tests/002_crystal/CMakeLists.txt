set(TESTS
  001_crystalenv-simple    ## simple use of the crystal environment
  002_crystalenv-arithm    ## crystal environment with arithmetics
  003_crystalenv-opts      ## options in the crystal environment
  004_crystalenv-tricky1   ## tricky cryst. env.: almost hexagonal
  005_crystalenv-tricky2   ## tricky cryst. env.: atom in a very large box
  006_crystalenv-tricky3   ## tricky cryst. env.: almost hexagonal (2)
  007_crystal_abinit       ## crystal structure from abinit files
  008_crystal_crystal      ## from crystal files
  009_crystal_elk          ## from elk files
  010_crystal_cifallspg    ## from cif files: whole spacegroup collection
  011_crystal_cifone       ## from cif files: only one crystal
  012_crystal_cifallblock  ## from cif files: all blocks from a cif file
  013_crystal_cifoneblock  ## from cif files: one block from a cif file
  014_crystal_xsf          ## from xsf files
  015_crystal_axsf         ## from axsf files; include nudging
  016_crystal_library      ## from the crystal structure library
  017_crystal_shelx        ## from shelx files
  018_crystal_cube         ## from a cube file
  019_crystal_bincube      ## from a binary cube file
  020_crystal_wien         ## from wien2k files
  021_crystal_vasp         ## from VASP files
  022_crystal_dftb+        ## from DFTB+ files
  023_crystal_qein_simple  ## from QE inputs
  024_crystal_qeout_simple ## from QE outputs; include step number in opts
  025_crystal_pwc          ## from pwc files
  026_crystal_siesta       ## from siesta files
  027_crystal_16           ## from a dmacrys fort.16 file
  028_crystal_21           ## from a dmacrys/neighcrys fort.21 file
  029_crystal_cifnoxyz     ## cif without xyz symops
  030_crystal_dmain        ## from a dmacrys input file
  031_crystal_aimsin       ## from an FHIaims input file
  032_crystal_aimsout      ## from an FHIaims output file
  033_crystal_aimsmol      ## from an FHIaims molecule
  034_crystal_tinkerfrac   ## from a TINKER frac file
  035_crystal_forceformat  ## force a crystal file format
  036_crystal_cell         ## from a CASTEP cell file
  037_crystal_geom         ## from a CASTEP geom file
  038_crystal_fploout      ## from an FPLO output file
  039_crystal_magres       ## from a magres file
  040_crystal_akaikkr      ## from an akaikkr input file
  041_crystal_xband        ## from an xband sys file
  )

runtests(${TESTS})
