set(TESTS
  001_molecule-xyz       ## read molecule structure from xyz files; including options
  002_molecule-gaussian  ## from Gaussian files (log,wfn,wfx,fchk)
  003_molecule-cube      ## from Gaussian cube files
  004_molecule-bincube   ## from binary cube files
  005_molecule-psi4      ## from psi4 files (output & molden)
  006_molecule-adf       ## from ADF files (molden)
  007_molecule-dftb+     ## from DFTB+ files
  008_moleculeenv-simple ## read molecule structure using a molecule environment; includes options
  009_molecule-library   ## from the molecule structure library
  010_molcell            ## molcell keyword
  )

runtests(${TESTS})
