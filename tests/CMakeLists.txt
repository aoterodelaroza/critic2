## enable testing

## list of test bundles
set(TESTDIRS 
  001_io ## input and output; arithmetics
  )

foreach(dir ${TESTDIRS})
  add_subdirectory(${dir})
endforeach()
