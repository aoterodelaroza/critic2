! Some utilities for building the GUI (e.g. wrappers around ImGui routines).

submodule (gui_utils) proc
  use iso_c_binding
  implicit none

contains

  ! Returns true if the last item has been hovered for at least thr
  ! seconds.
  module function igIsItemHovered_delayed(flags,thr)
    use gui_main, only: g
    use gui_interfaces_cimgui, only: igIsItemHovered
    integer(c_int), value :: flags
    real(c_float), intent(in) :: thr
    logical(c_bool) :: igIsItemHovered_delayed

    igIsItemHovered_delayed = igIsItemHovered(flags) .and. (g%HoveredIdTimer >= thr)

  end function igIsItemHovered_delayed

end submodule proc
