program nlopt_test
  implicit none

  include 'nlopt.f'

  integer*8 :: opt

  opt = 0
  call nlo_create(opt, NLOPT_LD_SLSQP, 2)
  call nlo_destroy(opt)

end program nlopt_test
