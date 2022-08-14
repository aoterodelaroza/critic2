  ! -*-F90-*-

  !xx! tinycthread constants
  integer(c_int), bind(C,name="const_thrd_error") :: thrd_error
  integer(c_int), bind(C,name="const_thrd_success") :: thrd_success
  integer(c_int), bind(C,name="const_thrd_timedout") :: thrd_timedout
  integer(c_int), bind(C,name="const_thrd_busy") :: thrd_busy
  integer(c_int), bind(C,name="const_thrd_nomem") :: thrd_nomem
  ! Mutex types
  integer(c_int), bind(C,name="const_mtx_plain") :: mtx_plain
  integer(c_int), bind(C,name="const_mtx_timed") :: mtx_timed
  integer(c_int), bind(C,name="const_mtx_recursive") :: mtx_recursive
