set(TESTS
  001_point_crys        ## point keyword
  002_pointprop1        ## pointprop: naming, clear, listing
  003_pointprop2        ## pointprop: arithmetic expressions, field selectors
  004_pointprop3        ## pointprop: shortcuts, MOs
  005_point_mol         ## point keyword in molecules; read points from a file
  006_line_crys         ## line keyword in crystals
  007_line_mol          ## line: in molecules
  008_plane_crys        ## plane keyword in crystals
  009_plane_contour     ## plane: contour plots
  010_plane_colormap    ## plane: colormap plots
  011_plane_relief      ## plane: relief plots
  012_plane_mols        ## plane: in molecules
  013_cube_simple       ## cube keyword in crystals, field selectors
  014_cube_header       ## cube: write the header only
  015_cube_cell         ## cube: cell keyword
  016_cube_grid         ## cube: grid keyword
  017_cube_files        ## cube: different file formats
  018_cube_mlwf         ## cube: mlwf keyword
  019_cube_mlwf_spin    ## cube: mlwf keyword, with spinpolarization
  020_cube_wannier      ## cube: wannier keyword
  021_cube_wannier_spin ## cube: wannier keyword, with spinpolarization
  022_cube_unk          ## cube: unk keyword
  023_cube_unk_sym      ## cube: unk keyword, with symmetry
  024_cube_unk_spin     ## cube: unk keyword, with spinpolarization
  )

runtests(${TESTS})
