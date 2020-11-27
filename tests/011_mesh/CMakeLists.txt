set(TESTS
  001_meshtype      ## basic use of the MESHTYPE keyword
  002_molcalc_basic ## basic use of the MOLCALC keyword
  003_molcalc_peach ## MOLCALC/PEACH keyword 
  )

runtests(${TESTS})
