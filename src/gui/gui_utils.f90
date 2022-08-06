! Some utilities for building the GUI (e.g. wrappers around ImGui routines).

module gui_utils
  use iso_c_binding
  implicit none

  private

  public :: igIsItemHovered_delayed

  ! module procedure interfaces
  interface
     module function igIsItemHovered_delayed(flags,thr)
       integer(c_int), value :: flags
       real(c_float), intent(in) :: thr
       logical(c_bool) :: igIsItemHovered_delayed
     end function igIsItemHovered_delayed
  end interface

end module gui_utils
