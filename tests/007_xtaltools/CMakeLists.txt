set(TESTS
  001_environ_simple         ## environ keyword, simple
  002_environ_opts           ## environ keyword, options
  003_environ_atomby         ## environ keyword, atom and by keywords
  004_ewald                  ## ewald keyword
  005_ewald_pot              ## ewald keyword, potential
  006_packing                ## packing keyword
  007_packing_vdw            ## packing keyword, vdw option
  008_newcell_supercell      ## newcell keyword, create a supercell
  009_newcell_manual         ## newcell keyword, manual specification of the new cell
  010_newcell_spglib         ## newcell keyword, spglib transformations
  011_newcell_reduced        ## newcell keyword, reduced cells
  012_newcell_refine         ## newcell keyword, refine keyword
  013_powder_basic           ## powder keyword
  014_powder_options         ## powder keyword, options
  015_rdf_basic              ## rdf keyword
  016_rdf_options            ## rdf keyword, options
  017_compare_crystals       ## compare keyword, crystals
  018_compare_crystals_opts  ## compare keyword, crystals with options
  019_compare_molecules      ## compare keyword, molecules
  020_compare_molecules_opts ## compare keyword, molecules with options
  021_atomlabel              ## atomlabel keyword
  022_bondfactor             ## bondfactor keyword
  023_coord                  ## coord keyword
  024_coord_opts             ## coord keyword, with options
  025_econ                   ## econ keyword
  026_identify               ## identify keyword
  027_identify_file          ## identify keyword, from file
  028_spg                    ## spg keyword
  029_sym_basic              ## sym keyword
  030_sym_large              ## sym keyword, with large structure
  031_sym_opts               ## sym keyword, options
  032_sym_refine             ## sym keyword, refine keyword
  033_vdw                    ## vdw keyword
  034_sym_wholemols          ## sym keyword, wholemols keyword
  035_compare_reduce         ## compare keyword, reduce option
  036_radii                  ## radii keyword
  037_polyhedra              ## polyhedra keyword
  038_molreorder_mol         ## molreorder keyword, molecule target
  039_molreorder_cry         ## molreorder keyword, crystal target
  040_molreorder_moveatoms   ## molreorder keyword, moveatoms
  041_kpoints                ## kpoints keyword
  042_bz                     ## bz keyword
  043_molmove                ## molmove keyword
  044_amd_basic              ## amd keyword
  )

runtests(${TESTS})
