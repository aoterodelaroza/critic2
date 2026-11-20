set(TESTS
  001_edit_delete    ## EDIT: delete atoms/molecules/hydrogens
  002_edit_move      ## EDIT: move atoms (relative/absolute/bohr/nneq)
  003_edit_cellmove  ## EDIT: cellmove cell lengths/angles/volume
  004_edit_vacuum    ## EDIT: vacuum along an axis (center/bottom/cut/orthogonal/keepmol)
  )

runtests(${TESTS})
