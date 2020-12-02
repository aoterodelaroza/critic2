set(TESTS
  001_basinplot_basic  ## basinplot keyword, basic
  002_basinplot_opts   ## basinplot keyword, options
  003_bundleplot_basic ## bundleplot keyword, basic
  004_bundleplot_opts  ## bundleplot keyword, options
  005_integrals_basic  ## integrals keyword, basic
  006_integrals_radial ## integrals keyword, with int_radial
  007_sphereintegrals  ## sphereintegrals keyword
  )

runtests(${TESTS})
