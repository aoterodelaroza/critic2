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
module windows
  use iso_c_binding
  use scenes, only: representation, scene
  use interfaces_cimgui, only: ImVec2
  use global, only: rborder_def
  use param, only: isformat_unknown
  implicit none

  private

  ! the buffer for the output console
  character(kind=c_char,len=:), allocatable, target :: outputb
  integer(c_size_t) :: lob = 0
  integer(c_size_t), parameter :: maxlob = 2000000

  ! the buffer for the input console
  character(kind=c_char,len=:), allocatable, target :: inputb
  integer(c_size_t), parameter :: maxlib = 40000

  ! command input and output
  integer, parameter :: command_inout_empty = 0
  integer, parameter :: command_inout_used = 1
  type command_inout
     integer :: id ! unique command ID
     integer :: status = command_inout_empty ! status of this command
     integer(c_size_t) :: size = 0 ! size of the output
     real(c_float) :: scrolly = 0._c_float ! scrolling position in console output view
     character(len=:,kind=c_char), allocatable :: tooltipinfo ! tooltip info for command (c_null_char term)
     character(len=:,kind=c_char), allocatable :: input ! command input (c_null_char term)
     character(len=:,kind=c_char), allocatable :: output ! command output (c_null_char term)
   contains
     procedure :: end => command_end
  end type command_inout

  ! command inout stack
  integer :: ncomid = 0 ! unique command ID generator (always incremented)
  integer :: ncom = 0 ! number of commands (com(:))
  integer :: nicom = 0 ! number of ordered commands (icom(:))
  integer :: idcom = 0 ! current output shown (0 = all)
  integer, allocatable :: icom(:) ! order in which to show the commands
  type(command_inout), allocatable, target :: com(:) ! the actual commands
  integer(c_size_t), parameter :: maxcomout = 20000000 ! maximum command size

  ! view mouse behavior parameters
  integer, parameter :: MB_navigation = 1

  ! user data for the file open dialog
  type, bind(c) :: dialog_userdata
     type(c_ptr) :: dptr = c_null_ptr ! the pointer for the file dialog
     integer(c_int) :: mol = -1 ! -1 = auto, 0 = crystal, 1 = molecule
     logical(c_bool) :: showhidden = .false._c_bool ! show hidden files
     integer(c_int) :: isformat = isformat_unknown ! force structure format
     logical(c_bool) :: readlastonly = .false._c_bool ! read only the last structure
     integer(c_int) :: purpose ! the purpose of the dialog
     logical(c_bool) :: molcubic = .false. ! whether to read the cell as cubic in a molecule
     real(c_float) :: rborder = real(rborder_def,c_float) ! border for the cell in a molecule
  end type dialog_userdata

  ! Wrapper class to handle ImGui windows
  type window
     ! global window parameters
     logical :: isinit = .false. ! whether this window has been initialized
     logical :: firstpass = .true. ! true if this is the first draw() call
     logical(c_bool) :: isopen ! whether the window is open
     integer :: type ! the window type
     integer(c_int) :: id ! internal ID for this window (index from win(:))
     integer :: idparent = 0 ! internal ID (from win(:)) for the caller window
     integer(c_int) :: flags ! window flags
     character(kind=c_char,len=:), allocatable :: name ! name of the window
     character(kind=c_char,len=:), allocatable :: errmsg ! error message in the window
     type(c_ptr) :: ptr ! ImGuiWindow* pointer (use only after Begin())
     type(c_ptr) :: dptr ! ImGuiFileDialog* pointer for dialogs
     integer :: isys = 1 ! the system on which the window operates
     real*8 :: timelastupdate ! time the window data was last updated
     ! tree table parameters
     integer :: table_selected = 1 ! the system selected in a table (input to iord)
     integer, allocatable :: iord(:) ! table order
     integer(c_int) :: table_sortcid = 0 ! sort table by this column id
     integer(c_int) :: table_sortdir = 1 ! sort table with this direction
     logical :: forceresize = .false. ! make true to force resize of columns
     logical :: forcesort = .false. ! make true to force a sort of the tree
     logical :: forceupdate = .false. ! make true to force an update of the tree
     logical :: forceinit = .false. ! make true to force an initialization of the systems
     integer, allocatable :: forceremove(:) ! enter integers to remove one or more systems
     integer :: forceselect = 0 ! make the tree select this system in the next pass
     ! view parameters
     type(scene), pointer :: sc ! pointer to the view scene
     integer(c_int) :: FBO ! draw framebuffer
     integer(c_int) :: FBOtex ! draw framebuffer -> texture
     integer(c_int) :: FBOdepth ! draw framebuffer -> depth buffer
     integer(c_int) :: FBOpick ! object pick framebuffer
     integer(c_int) :: FBOdepthp ! object pick framebuffer -> depth buffer
     integer(c_int) :: FBOrgba ! object pick framebuffer -> rgba buffer
     integer(c_int) :: FBOside ! side of the render texture (pixels)
     type(ImVec2) :: v_rmin, v_rmax ! view image rectangle
     integer :: view_selected = 1 ! the system selected in the view window
     logical :: forcerender = .true. ! force render of the scene
     integer :: view_mousebehavior = MB_navigation ! mouse behavior in the view
     integer :: idexportwin = 0 ! the ID for the export window
     ! dialog parameters
     integer :: dialog_purpose ! purpose of the dialog (open, save,...)
     type(dialog_userdata) :: dialog_data ! for the side pane callback
     logical :: forcequitdialog = .false. ! for the dialog to quit
     ! input console parameters
     integer :: inpcon_selected = 1 ! the system selected in the input console
     ! new structure from library & save image
     character(kind=c_char,len=:), allocatable :: okfile ! ok file
     character(kind=c_char,len=:), allocatable :: okfilter ! ok filter
     logical :: okfile_set = .false. ! whether the library file has been set by the user
     logical :: okfile_read = .false. ! whether the structure list should be re-read from the lib
     ! load field parameters
     ! scf plot and tree plot parameters
     real(c_double) :: ymin, ymax ! y-end of the plot
     integer :: plotn ! number of plot data
     real(c_double), allocatable :: plotx(:), ploty(:) ! plot data
     ! edit representation parameters
     type(representation), pointer :: rep => NULL() ! the representation on which the e.r. window operates
     ! export image parameters
     integer :: idsave = 0 ! window ID for the save file dialog (export image)
     integer(c_int) :: nsample ! number of samples for anti-aliasing
     integer(c_int) :: jpgquality ! jpg quality
     logical(c_bool) :: exportview ! export viewport or whole texture
     integer(c_int) :: npixel ! number of pixels in the export buffer
     logical(c_bool) :: transparentbg ! transparent background
     ! rebond parameters
   contains
     procedure :: init => window_init ! initialize the window
     procedure :: end => window_end ! finalize the window
     procedure :: focused => window_focused ! whether the root window is focused (even if not current)
     procedure :: draw => window_draw ! draw the window, calls one of the draw commands below
     ! tree procedures
     procedure :: draw_tree ! draw a tree
     procedure :: update_tree ! update the system information shown by the tree
     procedure :: sort_tree ! sort the systems in the tree
     procedure :: draw_treeplot
     ! view procedures
     procedure :: draw_view ! draw a view
     procedure :: create_texture_view ! create the texture for the view
     procedure :: delete_texture_view ! delete the texture for the view
     procedure :: process_events_view ! process the mouse events in the view
     procedure :: select_view
     procedure :: draw_selection_tooltip ! draw the measure selection tooltip
     procedure :: mousepos_to_texpos ! mouse position to texture position
     procedure :: texpos_to_mousepos ! texture position to mouse position
     procedure :: getpixel ! get depth or color from texture position
     procedure :: view_to_texpos ! view coordinates to texture position
     procedure :: texpos_to_view ! texture position to view coordinates
     procedure :: world_to_texpos ! world coordinates to texture position
     procedure :: texpos_to_world ! texture position to world coordinates
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
     ! about window
     procedure :: draw_about ! draw the about window
     ! new structure procedures
     procedure :: draw_new_struct
     ! new structure from library procedures
     procedure :: draw_new_struct_from_library
     ! load field
     procedure :: draw_load_field
     ! scf plot
     procedure :: update_scfplot
     procedure :: draw_scfplot
     ! edit representation
     procedure :: update_editrep
     procedure :: draw_editrep
     procedure :: draw_editrep_atoms
     procedure :: draw_editrep_unitcell
     ! export image
     procedure :: draw_exportimage
     ! rebond
     procedure :: draw_rebond
     ! preferences
     procedure :: draw_preferences
  end type window
  public :: window

  ! the window stack and named windows
  integer, public :: nwin
  type(window), allocatable, target, public :: win(:)
  integer, public :: iwin_tree
  integer, public :: iwin_console_input
  integer, public :: iwin_console_output
  integer, public :: iwin_view
  integer, public :: iwin_about

  ! window types
  integer, parameter, public :: wintype_tree = 1
  integer, parameter, public :: wintype_view = 2
  integer, parameter, public :: wintype_console_input = 3
  integer, parameter, public :: wintype_console_output = 4
  integer, parameter, public :: wintype_dialog = 5
  integer, parameter, public :: wintype_new_struct = 6
  integer, parameter, public :: wintype_new_struct_library = 7
  integer, parameter, public :: wintype_load_field = 8
  integer, parameter, public :: wintype_scfplot = 9
  integer, parameter, public :: wintype_editrep = 10
  integer, parameter, public :: wintype_exportimage = 11
  integer, parameter, public :: wintype_about = 12
  integer, parameter, public :: wintype_rebond = 13
  integer, parameter, public :: wintype_preferences = 14
  integer, parameter, public :: wintype_treeplot = 15

  ! window purposes
  integer, parameter, public :: wpurp_unknown = 0
  integer, parameter, public :: wpurp_dialog_openfiles = 1
  integer, parameter, public :: wpurp_dialog_savelogfile = 2
  integer, parameter, public :: wpurp_dialog_openlibraryfile = 3
  integer, parameter, public :: wpurp_dialog_openfieldfile = 4
  integer, parameter, public :: wpurp_dialog_openonefilemodal = 5
  integer, parameter, public :: wpurp_dialog_saveimagefile = 6

  ! routines to manipulate the window stack
  public :: stack_realloc_maybe
  public :: stack_create_window
  public :: update_window_id

  !xx! Interfaces
  interface
     module subroutine command_end(c)
       class(command_inout), intent(inout) :: c
     end subroutine command_end
     module subroutine stack_realloc_maybe()
     end subroutine stack_realloc_maybe
     module function stack_create_window(type,isopen,purpose,isys,irep,idcaller)
       integer, intent(in) :: type
       logical, intent(in) :: isopen
       integer, intent(in), optional :: purpose
       integer, intent(in), optional :: isys
       integer, intent(in), optional :: irep
       integer, intent(in), optional :: idcaller
       integer :: stack_create_window
     end function stack_create_window
     module subroutine update_window_id(id,changed)
       integer, intent(inout) :: id
       integer, intent(out), optional :: changed
     end subroutine update_window_id
     module subroutine window_init(w,type,isopen,id,purpose,isys,irep,idcaller)
       class(window), intent(inout), target :: w
       integer, intent(in) :: type
       logical, intent(in) :: isopen
       integer, intent(in) :: id
       integer, intent(in), optional :: purpose
       integer, intent(in), optional :: isys
       integer, intent(in), optional :: irep
       integer, intent(in), optional :: idcaller
     end subroutine window_init
     module subroutine window_end(w)
       class(window), intent(inout), target :: w
     end subroutine window_end
     module function window_focused(w)
       class(window), intent(inout) :: w
       logical :: window_focused
     end function window_focused
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
     module subroutine draw_treeplot(w)
       class(window), intent(inout), target :: w
     end subroutine draw_treeplot
     module subroutine draw_view(w)
       class(window), intent(inout), target :: w
     end subroutine draw_view
     module subroutine create_texture_view(w,atex)
       class(window), intent(inout), target :: w
       integer, intent(in) :: atex
     end subroutine create_texture_view
     module subroutine delete_texture_view(w)
       class(window), intent(inout), target :: w
     end subroutine delete_texture_view
     module subroutine select_view(w,isys)
       class(window), intent(inout), target :: w
       integer, intent(in) :: isys
     end subroutine select_view
     module subroutine process_events_view(w,hover,idx)
       class(window), intent(inout), target :: w
       logical, intent(in) :: hover
       integer(c_int), intent(in) :: idx(4)
     end subroutine process_events_view
     module subroutine draw_selection_tooltip(w,idx)
       class(window), intent(inout), target :: w
       integer(c_int), intent(in) :: idx(4)
     end subroutine draw_selection_tooltip
     module subroutine mousepos_to_texpos(w,pos)
       class(window), intent(inout), target :: w
       type(ImVec2), intent(inout) :: pos
     end subroutine mousepos_to_texpos
     module subroutine texpos_to_mousepos(w,pos)
       class(window), intent(inout), target :: w
       type(ImVec2), intent(inout) :: pos
     end subroutine texpos_to_mousepos
     module subroutine getpixel(w,fb,pos,depth,rgba)
       class(window), intent(inout), target :: w
       integer(c_int), intent(in) :: fb
       type(ImVec2), intent(inout) :: pos
       real(c_float), intent(out), optional :: depth
       real(c_float), intent(out), optional :: rgba(4)
     end subroutine getpixel
     module subroutine view_to_texpos(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine view_to_texpos
     module subroutine texpos_to_view(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine texpos_to_view
     module subroutine world_to_texpos(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine world_to_texpos
     module subroutine texpos_to_world(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine texpos_to_world
     module subroutine draw_dialog(w)
       class(window), intent(inout), target :: w
     end subroutine draw_dialog
     module subroutine draw_ci(w)
       class(window), intent(inout), target :: w
     end subroutine draw_ci
     module subroutine run_commands_ci(w)
       class(window), intent(inout), target :: w
     end subroutine run_commands_ci
     module subroutine block_gui_ci(w,allsys)
       class(window), intent(inout), target :: w
       logical, intent(in) :: allsys
     end subroutine block_gui_ci
     module function read_output_ci(w,iscom,cominfo)
       class(window), intent(inout), target :: w
       logical, intent(in) :: iscom
       character*(*), intent(in), optional :: cominfo
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
     module subroutine draw_about(w)
       class(window), intent(inout), target :: w
     end subroutine draw_about
     module subroutine draw_preferences(w)
       class(window), intent(inout), target :: w
     end subroutine draw_preferences
     module subroutine draw_new_struct(w)
       class(window), intent(inout), target :: w
     end subroutine draw_new_struct
     module subroutine draw_new_struct_from_library(w)
       class(window), intent(inout), target :: w
     end subroutine draw_new_struct_from_library
     module subroutine draw_load_field(w)
       class(window), intent(inout), target :: w
     end subroutine draw_load_field
     module subroutine update_scfplot(w)
       class(window), intent(inout), target :: w
     end subroutine update_scfplot
     module subroutine draw_scfplot(w)
       class(window), intent(inout), target :: w
     end subroutine draw_scfplot
     module subroutine update_editrep(w)
       class(window), intent(inout), target :: w
     end subroutine update_editrep
     module subroutine draw_editrep(w)
       class(window), intent(inout), target :: w
     end subroutine draw_editrep
     module function draw_editrep_atoms(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical(c_bool) :: changed
     end function draw_editrep_atoms
     module function draw_editrep_unitcell(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical(c_bool) :: changed
     end function draw_editrep_unitcell
     module subroutine draw_exportimage(w)
       class(window), intent(inout), target :: w
     end subroutine draw_exportimage
     module subroutine draw_rebond(w)
       class(window), intent(inout), target :: w
     end subroutine draw_rebond
  end interface

end module windows
