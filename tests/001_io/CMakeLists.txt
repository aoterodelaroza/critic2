set(TESTS
  001_simple_arithmetics ## a simple sum
  002_arithmetic1        ## more complex arithmetic operations, with no structure
  003_dosinput           ## same as 002 with CRLF line terminators
  004_macinput           ## same as 002 with CR line terminators (except first line)
  005_root_caps          ## root keyword and keyword capitalization
  006_echo               ## echo keyword
  007_reset_end          ## reset and end keywords
  008_clear              ## clear keyword
  009_library            ## library keyword
  010_list               ## list keyword
  011_run_system         ## run and system keywords
  )

runtests(${TESTS})
