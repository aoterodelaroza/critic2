set(TESTS
  001_load-wfn             ## load scalar file from wfn file
  002_load-wfx             ## from wfx file
  003_load-cube            ## from cube file
  004_load-cubeopts        ## from cube file; check interpolation options
  005_load-cubezpsp        ## from cube file; check core augmentation
  006_load-bincube         ## from binary cube file
  007_load-abinit          ## from abinit files
  008_load-vasp            ## from VASP files
  009_load-xsf             ## from xsf files
  010_load-wien2k          ## from WIEN2k files
  011_load-elk             ## from elk files (STATE.OUT)
  012_load-pi              ## from aiPI files
  013_load-dftb+           ## from DFTB+ files
  014_load-siesta          ## from siesta files
  015_load-elkgrid         ## from elk grid files (.OUT)
  016_load-fchk            ## from Gaussian fchk files
  017_load-molden          ## from psi4 molden files
  018_load-gaussian_hf     ## from Gaussian files; restricted HF case
  019_load-gaussian_uhf    ## from Gaussian files; unrestricted HF case
  020_load-gaussian_mp2    ## from Gaussian files; MP2 case
  021_load-psi4_hf         ## from psi4 files; RHF case
  022_load-psi4_uhf        ## from psi4 files; UHF case
  023_load-fchk_virtuals   ## from fchk files; including virtual orbitals
  024_load-molden_virtuals ## from psi4 files; including virtual orbitals
  025_load-pwc             ## from pwc files
  026_load-vaspspin        ## from VASP files; choose block from CHGCAR
  027_load-promolecular    ## load promolecular density, including fragments
  028_load-copy            ## copy fields
  029_load-force           ## force-load with a specific reader
  030_load-newgrid         ## load new grids from other fields
  031_load-promolecular2   ## load a promolecular grid from fragments
  032_load-core            ## load core density
  033_load-gridcalc        ## load with grid operations (lap, grad, etc.)
  034_load-wiencalc        ## load using WIEN2k operations (add/sub)
  035_load-normalize       ## load grids with normalization
  036_load-typnuc          ## typnuc option to load
  037_load-wiennorm        ## WIEN2k normalization option to load
  038_load-numericalderivs ## ghost loading and numerical derivatives
  039_benchmark            ## benchmark keyword
  040_grid-keywords        ## grid keywords (count, etc.)
  041_unload               ## unload keyword
  042_setfield             ## setfield keyword
  043_reference            ## reference keyword
  044_load-gaussian-rohf   ## from Gaussian files; restricted-open HF case
  045_load-pwc-contrib     ## read pwc file; only some contributions
  046_load-orca_rhf        ## from orca files; RHF case
  047_load-orca_uhf        ## from orca files; UHF case
  048_check_molwfn_lmax    ## molecular wavefunctions; check max. angular momentum
  049_check_molwfn_wfntyp  ## molecular wavefunctions; check wavefunction type
  050_load-corekeyw        ## core, nocore, zpsp keywords
  051_load-fmt             ## load grid from CASTEP fmt file
  052_load-fplogrid        ## load grid from FPLO grid file (001)
  053_load-moldenxtb       ## load molden file from xtb
  )

runtests(${TESTS})
