!! Link test for the xtb library: critic2 binds the xtb C API directly from
!! Fortran (see src/energy@proc.F90), so all that is needed is that these
!! symbols resolve in the library found by FindXTB.cmake. The program is
!! linked but never run. The entry points named are the ones an xtb too old
!! for critic2 does not have: xtb_loadGFNFF (GFN-FF) and xtb_getVirial
!! (periodic stress).
program xtb_test
  use iso_c_binding, only: c_ptr, c_double
  implicit none

  interface
     function c_xtb_new_environment() bind(c,name="xtb_newEnvironment") result(h)
       import c_ptr
       type(c_ptr) :: h
     end function c_xtb_new_environment
     subroutine c_xtb_load_gfnff(env,mol,calc,filename) bind(c,name="xtb_loadGFNFF")
       import c_ptr
       type(c_ptr), value :: env, mol, calc, filename
     end subroutine c_xtb_load_gfnff
     subroutine c_xtb_get_virial(env,res,virial) bind(c,name="xtb_getVirial")
       import c_ptr, c_double
       type(c_ptr), value :: env, res
       real(c_double), intent(out) :: virial(*)
     end subroutine c_xtb_get_virial
  end interface

  type(c_ptr) :: env
  real(c_double) :: virial(9)

  env = c_xtb_new_environment()
  call c_xtb_load_gfnff(env,env,env,env)
  call c_xtb_get_virial(env,env,virial)

end program xtb_test
