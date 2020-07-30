set(TESTS
  001_auto_simple      ## a simple search for critical points (AUTO)
  002_auto_opts        ## AUTO options
  003_auto_chk         ## AUTO: chk generation
  004_auto_seedobj     ## AUTO: seedobj option
  005_auto_clip_cube   ## AUTO: clip cube option
  006_auto_clip_sphere ## AUTO: clip cube option
  007_auto_seed        ## AUTO: seed options
  )

runtests(${TESTS})
