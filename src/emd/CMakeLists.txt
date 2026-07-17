set(SOURCES emd.c)
add_library(emd OBJECT ${SOURCES})

## keep this vendored code optimized even in Debug builds
c2_optimize_vendored(emd C)
