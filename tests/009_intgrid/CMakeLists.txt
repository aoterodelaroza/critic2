set(TESTS
  001_yt_basic             ## yt keyword, basic
  002_yt_options           ## yt keyword, options
  003_yt_discard           ## yt keyword, discard
  004_yt_plot              ## yt keyword, plotting and output options
  005_bader_basic          ## bader keyword, basic
  006_bader_options        ## bader keyword, options
  007_bader_discard        ## bader keyword, discard
  008_bader_plot           ## bader keyword, plotting and output options
  009_integrable_basic     ## integrable keyword, basic
  010_integrable_multipole ## integrable keyword, multipole
  011_dis_wnr              ## integrable keyword, dis, with wnr
  012_dis_psink            ## integrable keyword, dis, with psink
  013_dis_wnr_feo          ## integrable keyword, dis, with wnr (spinpol)
  014_dis_psink_feo        ## integrable keyword, dis, with psink (spinpol)
  015_dis_psink_na         ## integrable keyword, dis, with psink (metal)
  016_hirshfeld_grid       ## hirshfeld keyword, using grids
  017_hirshfeld_nogrid_crystal ## hirshfeld keyword, using meshes in a crystal
  018_hirshfeld_nogrid_mol ## hirshfeld keyword, using meshes in a molecule
  019_isosurface           ## isosurface keyword
  020_grid_only            ## only option in YT, BADER
  )

runtests(${TESTS})
