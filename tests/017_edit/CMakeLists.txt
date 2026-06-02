set(TESTS
  001_edit_delete    ## EDIT: delete atoms/molecules/hydrogens
  002_edit_move      ## EDIT: move atoms (relative/absolute/bohr/nneq)
  003_edit_cellmove  ## EDIT: cellmove cell lengths/angles/volume
  )

runtests(${TESTS})
