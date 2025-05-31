program libxc_test
  use xc_f03_lib_m
  implicit none

  integer :: iver(3)

  call xc_f03_version(iver(1),iver(2),iver(3))

end program libxc_test
