set(TESTS
  001_environ_simple     ## environ keyword, simple
  002_environ_opts       ## environ keyword, options
  003_environ_atomby     ## environ keyword, atom and by keywords
  004_ewald              ## ewald keyword
  005_ewald_pot          ## ewald keyword, potential
  006_packing            ## packing keyword
  007_packing_vdw        ## packing keyword, vdw option
  008_newcell_supercell  ## newcell keyword, create a supercell
  009_newcell_manual     ## newcell keyword, manual specification of the new cell
  010_newcell_spglib     ## newcell keyword, spglib transformations
  011_newcell_reduced    ## newcell keyword, reduced cells
  012_newcell_refine     ## newcell keyword, refine keyword
  013_powder_basic       ## powder keyword
  014_powder_options     ## powder keyword, options
  )

runtests(${TESTS})
