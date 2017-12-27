  subroutine assert_non_null(ptr,routine,ptrname)
    use iso_c_binding, only: c_associated
    use tools_io, only: ferror, faterr
    character*(*), intent(in) :: routine, ptrname
    type(c_ptr), intent(in) :: ptr

    if (.not.c_associated(ptr)) &
       call ferror(trim(routine),"Unexpected NULL pointer - " // trim(ptrname),faterr)

  end subroutine assert_non_null

