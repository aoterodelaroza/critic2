module gui_keybindings
  use iso_c_binding, only: c_int
  implicit none

  private

  public :: erase_bind
  public :: set_bind
  public :: set_default_keybindings
  public :: is_bind_event

  ! Public list of binds
  integer, parameter, public :: BIND_QUIT = 1 ! quit the program
  integer, parameter, public :: BIND_NUM = 1 ! total number of binds
  ! #define BIND_CLOSE_LAST_DIALOG 1  // Closes the last window
  ! #define BIND_CLOSE_ALL_DIALOGS 2  // Closes all windows
  ! #define BIND_VIEW_ALIGN_A_AXIS 3  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_B_AXIS 4  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_C_AXIS 5  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_X_AXIS 6  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_Y_AXIS 7  // Align view with a axis
  ! #define BIND_VIEW_ALIGN_Z_AXIS 8  // Align view with a axis
  ! #define BIND_NAV_ROTATE        9  // Rotate the camera (navigation)
  ! #define BIND_NAV_TRANSLATE     10 // Camera pan (navigation)
  ! #define BIND_NAV_ZOOM          11 // Camera zoom (navigation)
  ! #define BIND_NAV_RESET         12 // Reset the view (navigation)

  ! module procedure interfaces
  interface
     module subroutine erase_bind(key, mod, group)
       integer(c_int), intent(in) :: key, mod
       integer, intent(in) :: group
     end subroutine erase_bind
     module subroutine set_bind(bind, key, mod)
       use tools_io, only: ferror, faterr
       integer, intent(in) :: bind
       integer(c_int), intent(in) :: key, mod
     end subroutine set_bind
     module subroutine set_default_keybindings()
     end subroutine set_default_keybindings
     module function is_bind_event(bind,held)
       integer, intent(in) :: bind
       logical, intent(in), optional :: held
       logical :: is_bind_event
     end function is_bind_event
  end interface

end module gui_keybindings
