! Copyright (c) 2019 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! The class to handle ImGui windows.
module gui_window
  use iso_c_binding
  use param, only: isformat_unknown
  implicit none

  private

  ! user data for the file open dialog
  type, bind(c) :: opendialog_userdata
     integer(c_int) :: mol = -1 ! -1 = auto, 0 = crystal, 1 = molecule
     logical(c_bool) :: showhidden = .false._c_bool ! show hidden files
     integer(c_int) :: isformat = isformat_unknown ! force structure format
  end type opendialog_userdata

  ! Wrapper class to handle ImGui windows
  type window
     ! global window parameters
     logical :: isinit = .false. ! whether this window has been initialized
     logical(c_bool) :: isopen ! whether the window is open
     integer :: type ! the window type
     integer(c_int) :: id ! internal ID for this window
     integer(c_int) :: flags ! window flags
     character(kind=c_char,len=:), allocatable :: name ! name of the window
     type(c_ptr) :: ptr ! ImGuiWindow*/ImGuiFileDialog* pointer (use only after Begin())
     ! tree table parameters
     integer :: table_selected = 1 ! the system selected in a table (input to iord)
     integer, allocatable :: iord(:) ! table order
     integer(c_int) :: table_sortcid = 0 ! sort table by this column id
     integer(c_int) :: table_sortdir = 1 ! sort table with this direction
     logical :: forceresize = .false. ! make true to force resize of columns
     logical :: forcesort = .false. ! make true to force a sort of the tree
     logical :: forceupdate = .false. ! make true to force an update of the tree
     logical :: forceinit = .false. ! make true to force an initialization of the systems
     integer :: forceremove = 0 ! make an integer to remove one of the systems
     ! opendialog parameters
     type(opendialog_userdata) :: od_data ! for the side pane callback
   contains
     procedure :: init => window_init
     procedure :: end => window_end
     procedure :: draw => window_draw
     ! tree procedures
     procedure :: draw_tree
     procedure :: update_tree
     procedure :: sort_tree
     ! opendialog procedures
     procedure :: draw_opendialog
  end type window
  public :: window

  ! the window stack and named windows
  integer, public :: nwin
  type(window), allocatable, target, public :: win(:)
  integer, public :: iwin_tree
  integer, public :: iwin_console
  integer, public :: iwin_view

  ! window types
  integer, parameter, public :: wintype_tree = 1
  integer, parameter, public :: wintype_view = 2
  integer, parameter, public :: wintype_console = 3
  integer, parameter, public :: wintype_opendialog = 4

  ! routines to manipulate the window stack
  public :: stack_create_window

  interface
     module function stack_create_window(type,isopen)
       integer, intent(in) :: type
       logical, intent(in) :: isopen
       integer :: stack_create_window
     end function stack_create_window
     module subroutine window_init(w,type,isopen)
       class(window), intent(inout) :: w
       integer, intent(in) :: type
       logical, intent(in) :: isopen
     end subroutine window_init
     module subroutine window_end(w)
       class(window), intent(inout) :: w
     end subroutine window_end
     module subroutine window_draw(w)
       class(window), intent(inout), target :: w
     end subroutine window_draw
     module subroutine draw_tree(w)
       class(window), intent(inout), target :: w
     end subroutine draw_tree
     module subroutine update_tree(w)
       class(window), intent(inout) :: w
     end subroutine update_tree
     module subroutine sort_tree(w,cid,dir)
       class(window), intent(inout) :: w
       integer(c_int), intent(in) :: cid, dir
     end subroutine sort_tree
     module subroutine draw_opendialog(w)
       class(window), intent(inout), target :: w
     end subroutine draw_opendialog
  end interface

end module gui_window
