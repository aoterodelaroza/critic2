set(TESTS
  001_auto_simple      ## a simple search for critical points (AUTO)
  002_auto_opts        ## AUTO options
  003_auto_chk         ## AUTO: chk generation
  004_auto_seedobj     ## AUTO: seedobj option
  005_auto_clip_cube   ## AUTO: clip cube option
  006_auto_clip_sphere ## AUTO: clip cube option
  007_auto_seed        ## AUTO: seed options
  008_cpreport         ## cpreport keyword
  009_cpreport_xyz     ## cpreport keyword, xyz
  010_cpreport_obj     ## cpreport keyword, obj,ply,off
  011_cpreport_graph   ## cpreport keyword, graph option
  )

runtests(${TESTS})
