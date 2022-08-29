! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
  use global, only: rborder_def
  use param, only: isformat_unknown
  implicit none

  private

  ! user data for the file open dialog
  type, bind(c) :: dialog_userdata
     type(c_ptr) :: ptr = c_null_ptr ! the pointer for the file dialog
     integer(c_int) :: mol = -1 ! -1 = auto, 0 = crystal, 1 = molecule
     logical(c_bool) :: showhidden = .false._c_bool ! show hidden files
     integer(c_int) :: isformat = isformat_unknown ! force structure format
     logical(c_bool) :: readlastonly = .false._c_bool ! read only the last structure
     integer(c_int) :: purpose ! the purpose of the dialog
     logical(c_bool) :: molcubic = .false. ! whether to read the cell as cubic in a molecule
     real(c_float) :: rborder = rborder_def ! border for the cell in a molecule
  end type dialog_userdata

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
     integer, allocatable :: forceremove(:) ! make an integer to remove one of the systems
     ! dialog parameters
     integer :: dialog_purpose ! purpose of the dialog (open, save,...)
     type(dialog_userdata) :: dialog_data ! for the side pane callback
     ! input console parameters
     integer :: inpcon_selected = 1 ! the system selected in the input console
   contains
     procedure :: init => window_init ! initialize the window
     procedure :: end => window_end ! finalize the window
     procedure :: draw => window_draw ! draw the window, calls one of the draw commands below
     ! tree procedures
     procedure :: draw_tree ! draw a tree
     procedure :: update_tree ! update the system information shown by the tree
     procedure :: sort_tree ! sort the systems in the tree
     ! dialog procedures
     procedure :: draw_dialog ! draw an open/save dialog
     ! input console procedures
     procedure :: draw_ci ! draw the input console
     procedure :: run_commands_ci ! run the commands in the input buffer
     procedure :: block_gui_ci ! block the GUI while the input buffer commands are run
     procedure :: read_output_ci ! read the generated output and maybe create a command i/o
     procedure :: fill_input_ci ! fill the input buffer with the given string
     procedure :: get_input_details_ci ! get the system & field strings for current input
     ! output console procedures
     procedure :: draw_co ! draw the output console
  end type window
  public :: window

  ! the window stack and named windows
  integer, public :: nwin
  type(window), allocatable, target, public :: win(:)
  integer, public :: iwin_tree
  integer, public :: iwin_console_input
  integer, public :: iwin_console_output
  integer, public :: iwin_view

  ! window types
  integer, parameter, public :: wintype_tree = 1
  integer, parameter, public :: wintype_view = 2
  integer, parameter, public :: wintype_console_input = 3
  integer, parameter, public :: wintype_console_output = 4
  integer, parameter, public :: wintype_dialog = 5

  ! window purposes
  integer, parameter, public :: wpurp_unknown = 0
  integer, parameter, public :: wpurp_dialog_openfiles = 1
  integer, parameter, public :: wpurp_dialog_savelogfile = 2

  ! routines to manipulate the window stack
  public :: stack_create_window

  interface
     module function stack_create_window(type,isopen,purpose)
       integer, intent(in) :: type
       logical, intent(in) :: isopen
       integer, intent(in), optional :: purpose
       integer :: stack_create_window
     end function stack_create_window
     module subroutine window_init(w,type,isopen,purpose)
       class(window), intent(inout) :: w
       integer, intent(in) :: type
       logical, intent(in) :: isopen
       integer, intent(in), optional :: purpose
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
     module subroutine draw_dialog(w)
       class(window), intent(inout), target :: w
     end subroutine draw_dialog
     module subroutine draw_ci(w)
       class(window), intent(inout), target :: w
     end subroutine draw_ci
     module subroutine run_commands_ci(w)
       class(window), intent(inout), target :: w
     end subroutine run_commands_ci
     module subroutine block_gui_ci(w)
       class(window), intent(inout), target :: w
     end subroutine block_gui_ci
     module function read_output_ci(w,iscom)
       class(window), intent(inout), target :: w
       logical, intent(in) :: iscom
       logical :: read_output_ci
     end function read_output_ci
     module subroutine fill_input_ci(w,str)
       class(window), intent(inout), target :: w
       character(len=*), intent(in) :: str
     end subroutine fill_input_ci
     module subroutine get_input_details_ci(w,csystem,cfield)
       class(window), intent(inout), target :: w
       character(len=:), allocatable, intent(inout) :: csystem, cfield
     end subroutine get_input_details_ci
     module subroutine draw_co(w)
       class(window), intent(inout), target :: w
     end subroutine draw_co
  end interface

end module gui_window
