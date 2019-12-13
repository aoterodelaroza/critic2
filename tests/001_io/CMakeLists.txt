set(TESTS
  001_simple_arithmetics ## a simple sum
  002_arithmetic1        ## more complex arithmetic operations, with no structure
  003_dosinput           ## same as 002 with CRLF line terminators
  004_macinput           ## same as 002 with CR line terminators (except first line)
  005_root_caps          ## root keyword and keyword capitalization
  )

runtests(${TESTS})
