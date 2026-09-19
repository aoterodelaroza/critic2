!! Link test for the tblite library: critic2 binds the tblite C API directly
!! from Fortran (see src/energy@proc.F90), so all that is needed is that these
!! symbols resolve in the library found by FindTBLITE.cmake. The program is
!! linked but never run. tblite_update_structure_geometry and
!! tblite_get_result_virial are the entry points a tblite too old for critic2
!! does not have.
program tblite_test
  use iso_c_binding, only: c_ptr, c_double
  implicit none

  interface
     function c_tblite_new_context() bind(c,name="tblite_new_context") result(h)
       import c_ptr
       type(c_ptr) :: h
     end function c_tblite_new_context
     subroutine c_tblite_get_singlepoint(ctx,mol,calc,res) bind(c,name="tblite_get_singlepoint")
       import c_ptr
       type(c_ptr), value :: ctx, mol, calc, res
     end subroutine c_tblite_get_singlepoint
     subroutine c_tblite_update_structure_geometry(err,mol,positions,lattice) &
        bind(c,name="tblite_update_structure_geometry")
       import c_ptr
       type(c_ptr), value :: err, mol, positions, lattice
     end subroutine c_tblite_update_structure_geometry
     subroutine c_tblite_get_result_virial(err,res,virial) bind(c,name="tblite_get_result_virial")
       import c_ptr, c_double
       type(c_ptr), value :: err, res
       real(c_double), intent(out) :: virial(*)
     end subroutine c_tblite_get_result_virial
  end interface

  type(c_ptr) :: ctx
  real(c_double) :: virial(9)

  ctx = c_tblite_new_context()
  call c_tblite_update_structure_geometry(ctx,ctx,ctx,ctx)
  call c_tblite_get_singlepoint(ctx,ctx,ctx,ctx)
  call c_tblite_get_result_virial(ctx,ctx,virial)

end program tblite_test
