program libxc_test
  use xc_f90_lib_m
  implicit none

  integer :: iver(3)

  call xc_f90_version(iver(1),iver(2),iver(3))

end program libxc_test
