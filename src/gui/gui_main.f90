! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! Structure class and routines for basic crystallography computations
module gui_main
  use iso_c_binding, only: c_ptr, c_float, c_int, c_null_ptr, c_bool
  use interfaces_cimgui, only: ImGuiIO, ImGuiContext, ImVec4, ImVec2, ImGuiViewport,&
     ImFontAtlas
  use windows, only: window
  use systemmod, only: system
  use crystalseedmod, only: crystalseed
  use scenes, only: scene
  use global, only: bondfactor_def
  use param, only: maxzat0, atmcov0
  implicit none

  private

  ! variables to GUI's structures & data
  type(ImGuiIO), pointer, public :: io ! pointer to ImGui's IO object
  type(ImFontAtlas), pointer, public :: fonts ! pointer to IO%Fonts
  type(ImGuiContext), pointer, public :: g ! pointer to ImGui's context
  type(ImGuiViewport), pointer, public :: mainvwp ! pointer to main viewport
  type(c_ptr), public :: rootwin ! the root window pointer (GLFWwindow*)
  type(ImVec2), public :: fontsize ! font size (sensitive to scaling)
  real(c_float), parameter, public :: fontbakesize = 16._c_float ! normal bake size (for GUI)
  real(c_float), parameter, public :: fontbakesize_large = 128._c_float ! large bake size (for rendering)
  type(c_ptr), public :: font_normal, font_large ! GUI and rendering font pointers

  ! GUI control parameters
  ! integer(c_int), parameter, public :: ms_samples = 1 ! number of samples in multisamples
  logical, public :: tooltip_enabled = .true. ! whether tooltips are enabled
  real(c_float), public :: tooltip_delay = 0.5_c_float ! tooltip delay, in seconds
  real(c_float), public :: tooltip_wrap_factor = 40._c_float ! tooltip wrap factor (fontsize)
  logical, parameter, public :: reuse_mid_empty_systems = .false. ! whether to reuse the empty systems in the middle
  logical, public :: tree_select_updates_inpcon = .true. ! selecting in tree chooses system in input console
  logical, public :: tree_select_updates_view = .true. ! selecting in tree chooses system in view
  integer, public :: lockbehavior = 1 ! 0=no-lock, 1=only SCF, 2=all systems

  ! GUI colors, defaults
  real(c_float), parameter, public :: ColorTableCellBg_def(4,0:8) = reshape((/& ! tree table name cell
     0.80_c_float,0.00_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic
     1.00_c_float,0.43_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic, layered
     1.00_c_float,0.79_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic, chain
     0.58_c_float,1.00_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic, molecular
     0.00_c_float,1.00_c_float,0.96_c_float,0.4_c_float,&  ! slab
     0.00_c_float,0.28_c_float,1.00_c_float,0.4_c_float,&  ! chain
     0.51_c_float,0.51_c_float,0.51_c_float,0.4_c_float,&  ! molecule in a box
     0.56_c_float,0.00_c_float,1.00_c_float,0.4_c_float,&  ! single molecule
     1.00_c_float,0.00_c_float,0.66_c_float,0.4_c_float/),&! molecular cluster
     shape(ColorTableCellBg_def))
  real(c_float), parameter, public :: ColorHighlightScene_def(4) = (/1.0_c_float,1.0_c_float,0.5_c_float,0.7_c_float/) ! hover highlight of atoms in view
  real(c_float), parameter, public :: ColorHighlightSelectScene_def(4) = (/0.2_c_float,0.64_c_float,0.9_c_float,0.7_c_float/) ! edit geometry selection color
  real(c_float), parameter, public :: ColorMeasureSelect_def(4,4) = reshape((/&
     1._c_float,  0.4_c_float, 0.4_c_float, 0.5_c_float,&
     0.4_c_float, 1._c_float,  0.4_c_float, 0.5_c_float,&
     0.4_c_float, 0.4_c_float, 1._c_float, 0.5_c_float,&
     0.9_c_float, 0.7_c_float, 0.4_c_float, 0.5_c_float/),(/4,4/))

  ! GUI colors, actual values
  real(c_float), public :: ColorTableCellBg(4,0:8)
  type(ImVec4), parameter, public :: ColorDialogDir = ImVec4(0.9, 0.9, 0.5, 1.0) ! directories in the dialog
  type(ImVec4), parameter, public :: ColorDialogFile = ImVec4(1.0, 1.0, 1.0, 1.0) ! files in the dialog
  type(ImVec4), parameter, public :: ColorHighlightText = ImVec4(0.2, 0.64, 0.9, 1.0) ! highlighted text
  type(ImVec4), parameter, public :: ColorDangerButton = ImVec4(0.50, 0.08, 0.08, 1.0) ! important button
  type(ImVec4), parameter, public :: ColorDangerText = ImVec4(0.80, 0.08, 0.08, 1.0) ! important text
  type(ImVec4), parameter, public :: ColorWaitBg = ImVec4(0.80, 0.80, 0.80, 0.6) ! dim the background while waiting
  type(ImVec4), parameter, public :: ColorFrameBgAlt = ImVec4(0.29,0.16,0.48,0.54) ! alternate framebg
  type(ImVec4), parameter, public :: ColorFieldSelected = ImVec4(0.91,1.00,0.00,0.31) ! selected field
  type(ImVec4), parameter, public :: ColorTableHighlightRow = ImVec4(1._c_float,0.8_c_float,0.1_c_float,0.5_c_float) ! selectable highlight color
  real(c_float), public :: ColorHighlightScene(4)
  real(c_float), public :: ColorHighlightSelectScene(4)
  real(c_float), public :: ColorMeasureSelect(4,4)

  ! element colors
  real(c_float), public :: ColorElement(3,0:maxzat0)

  ! system status (from lower to higher initialization level)
  integer, parameter, public :: sys_empty = 0 ! not in use
  integer, parameter, public :: sys_loaded_not_init = 1 ! the seed is available, waiting for initialization
  integer, parameter, public :: sys_initializing = 2 ! the system is initializing
  integer, parameter, public :: sys_ready = 3 ! the data is ready but thread is still working, so not initialized yet
  integer, parameter, public :: sys_init = 4 ! the system is initialized

  ! system configuration type
  type :: sysconf
     ! system ID and properties
     integer :: id ! ID for this system
     integer :: status = sys_empty ! current status
     logical :: hidden = .false. ! whether it is hidden in the tree view (filter)
     logical :: showfields = .false. ! whether to show the fields in the tree view
     type(crystalseed) :: seed ! generating seed
     logical :: has_field = .false. ! true if the seed has a field
     logical :: has_vib = .false. ! true if the seed has vibrational data
     integer :: iperiod = 0 ! periodicity (see iperiod_*)
     integer :: collapse ! 0 if independent, -1 if master-collapsed, -2 if master-extended, <n> if dependent on n
     type(c_ptr) :: thread_lock = c_null_ptr ! the lock for initialization of this system
     ! system name
     character(len=:), allocatable :: fullname ! full-path name
     logical :: renamed = .false. ! true if the system has been renamed
     ! scene
     type(scene) :: sc ! scene for the system in the main view
     ! bonding
     real*8 :: atmcov(0:maxzat0) = atmcov0 ! covalent radii for bonding
     real*8 :: bondfactor = bondfactor_def ! bond factor for bonding calculation
     ! highlights
     real(c_float), allocatable :: highlight_rgba(:,:) ! highlight colors
     real(c_float), allocatable :: highlight_rgba_transient(:,:) ! transient highlight colors
     logical :: highlight_transient_set = .false. ! set to false at beginning of main loop; clears transient highlight at end of loop
     ! time
     real*8 :: timelastchange_geometry = 0d0   ! time system last changed geometry
     real*8 :: timelastchange_rebond = 0d0     ! time system last was rebonded
     real*8 :: timelastchange_buildlists = 0d0 ! time system last required a list rebuild
     real*8 :: timelastchange_render = 0d0     ! time system last required a render
   contains
     procedure :: highlight_atoms
     procedure :: highlight_clear
     procedure :: remove_highlighted_atoms
     procedure :: set_timelastchange
  end type sysconf

  ! list of changes to the system, in order of severity
  integer, parameter, public :: lastchange_render = 0     ! system needs a new render
  integer, parameter, public :: lastchange_buildlists = 1 ! system needs building new lists
  integer, parameter, public :: lastchange_rebond = 2     ! system has been rebonded
  integer, parameter, public :: lastchange_geometry = 3   ! system geometry has changed

  ! system arrays
  integer, public :: nsys = 0
  type(system), allocatable, target, public :: sys(:)
  type(sysconf), allocatable, target, public :: sysc(:)

  ! flags to control main's behavior
  integer, public :: force_run_commands = 0 ! execute commands from the input console (0=no,1=only selected,2=all)
  logical, public :: force_quit_threads = .false. ! set to true to force all threads to quit as soon as possible

  ! public procedures
  public :: gui_start
  public :: launch_initialization_thread
  public :: kill_initialization_thread
  public :: are_threads_running
  public :: system_shorten_names
  public :: add_systems_from_seeds
  public :: add_systems_from_name
  public :: remove_system
  public :: remove_systems
  public :: reread_system_from_file
  public :: duplicate_system
  public :: set_default_ui_settings
  public :: regenerate_system_pointers
  public :: ok_system

  interface
     module subroutine gui_start()
     end subroutine gui_start
     module subroutine launch_initialization_thread()
     end subroutine launch_initialization_thread
     module subroutine kill_initialization_thread()
     end subroutine kill_initialization_thread
     module function are_threads_running()
       logical :: are_threads_running
     end function are_threads_running
     module subroutine system_shorten_names()
     end subroutine system_shorten_names
     module subroutine add_systems_from_seeds(nseed,seed,collapse,iafield,iavib,forceidx)
       integer, intent(in) :: nseed
       type(crystalseed), allocatable, intent(in) :: seed(:)
       logical, intent(in), optional :: collapse
       integer, intent(in), optional :: iafield, iavib
       integer, intent(in), optional :: forceidx
     end subroutine add_systems_from_seeds
     module subroutine add_systems_from_name(name,mol,isformat,readlastonly,rborder,molcubic,&
        forceidx)
       character(len=*), intent(in) :: name
       integer, intent(in) :: mol
       integer, intent(in) :: isformat
       logical, intent(in) :: readlastonly
       real*8, intent(in) :: rborder
       logical, intent(in) :: molcubic
       integer, intent(in), optional :: forceidx
     end subroutine add_systems_from_name
     recursive module subroutine remove_system(idx,kill_dependents_if_extended)
       integer, intent(in) :: idx
       logical, intent(in), optional :: kill_dependents_if_extended
     end subroutine remove_system
     module subroutine remove_systems(idx)
       integer, intent(in) :: idx(:)
     end subroutine remove_systems
     recursive module subroutine reread_system_from_file(idx)
       integer, intent(in) :: idx
     end subroutine reread_system_from_file
     recursive module subroutine duplicate_system(idx)
       integer, intent(in) :: idx
     end subroutine duplicate_system
     module subroutine set_default_ui_settings()
     end subroutine set_default_ui_settings
     module subroutine regenerate_system_pointers()
     end subroutine regenerate_system_pointers
     module function ok_system(isys,level)
       integer, intent(in) :: isys, level
       logical :: ok_system
     end function ok_system
     ! sysconf
     module subroutine highlight_atoms(sysc,transient,idx,type,rgba)
       class(sysconf), intent(inout) :: sysc
       logical, intent(in) :: transient
       integer, intent(in) :: idx(:)
       integer, intent(in) :: type
       real(c_float), intent(in) :: rgba(:,:)
     end subroutine highlight_atoms
     module subroutine highlight_clear(sysc,transient,idx,type)
       class(sysconf), intent(inout) :: sysc
       logical, intent(in) :: transient
       integer, intent(in), optional :: idx(:)
       integer, intent(in), optional :: type
     end subroutine highlight_clear
     module subroutine remove_highlighted_atoms(sysc)
       class(sysconf), intent(inout) :: sysc
     end subroutine remove_highlighted_atoms
     module subroutine set_timelastchange(sysc,level)
       class(sysconf), intent(inout) :: sysc
       integer, intent(in) :: level
     end subroutine set_timelastchange
  end interface

end module gui_main
