set(TESTS
  001_yt_basic      ## yt keyword, basic
  002_yt_options    ## yt keyword, options
  003_yt_discard    ## yt keyword, discard
  004_yt_plot       ## yt keyword, plotting and output options
  005_bader_basic   ## bader keyword, basic
  006_bader_options ## bader keyword, options
  007_bader_discard ## bader keyword, discard
  008_bader_plot    ## bader keyword, plotting and output options
  )

runtests(${TESTS})
