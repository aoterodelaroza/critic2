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
  implicit none

  private

  ! variables to GUI's structures & data
  real*8, public :: time ! the time
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
  logical(c_bool), public :: tooltip_enabled = .true. ! whether tooltips are enabled
  real(c_float), public :: tooltip_delay = 0.5_c_float ! tooltip delay, in seconds
  real(c_float), public :: tooltip_wrap_factor = 40._c_float ! tooltip wrap factor (fontsize)
  logical, parameter, public :: reuse_mid_empty_systems = .false. ! whether to reuse the empty systems in the middle
  logical(c_bool), public :: tree_select_updates_inpcon = .true. ! selecting in tree chooses system in input console
  logical(c_bool), public :: tree_select_updates_view = .true. ! selecting in tree chooses system in view
  integer, public :: lockbehavior = 1 ! 0=no-lock, 1=only SCF, 2=all systems

  ! GUI colors
  type(ImVec4), public :: ColorTableCellBg_Mol     = ImVec4(0.43,0.8 ,0.  ,0.2)  ! tree table name cell, molecule
  type(ImVec4), public :: ColorTableCellBg_MolClus = ImVec4(0.0 ,0.8 ,0.43,0.2)  ! tree table name cell, molecular cluster
  type(ImVec4), public :: ColorTableCellBg_MolCrys = ImVec4(0.8 ,0.43,0.0 ,0.2)  ! tree table name cell, molecular crystal
  type(ImVec4), public :: ColorTableCellBg_Crys3d  = ImVec4(0.8 ,0.  ,0.0 ,0.2)  ! tree table name cell, 3d crystal
  type(ImVec4), public :: ColorTableCellBg_Crys2d  = ImVec4(0.8 ,0.  ,0.43,0.2)  ! tree table name cell, 2d crystal
  type(ImVec4), public :: ColorTableCellBg_Crys1d  = ImVec4(0.8 ,0.43,0.43,0.2)  ! tree table name cell, 1d crystal
  type(ImVec4), parameter, public :: ColorDialogDir = ImVec4(0.9, 0.9, 0.5, 1.0) ! directories in the dialog
  type(ImVec4), parameter, public :: ColorDialogFile = ImVec4(1.0, 1.0, 1.0, 1.0) ! files in the dialog
  type(ImVec4), parameter, public :: ColorHighlightText = ImVec4(0.2, 0.64, 0.9, 1.0) ! highlighted text
  type(ImVec4), parameter, public :: ColorDangerButton = ImVec4(0.50, 0.08, 0.08, 1.0) ! important button
  type(ImVec4), parameter, public :: ColorDangerText = ImVec4(0.80, 0.08, 0.08, 1.0) ! important text
  type(ImVec4), parameter, public :: ColorWaitBg = ImVec4(0.80, 0.80, 0.80, 0.6) ! dim the background while waiting
  type(ImVec4), parameter, public :: ColorFrameBgAlt = ImVec4(0.29,0.16,0.48,0.54) ! alternate framebg
  type(ImVec4), parameter, public :: ColorFieldSelected = ImVec4(0.91,1.00,0.00,0.31) ! alternate framebg

  ! systems arrays
  integer, parameter, public :: sys_empty = 0
  integer, parameter, public :: sys_loaded_not_init = 1
  integer, parameter, public :: sys_initializing = 2
  integer, parameter, public :: sys_init = 3
  type :: sysconf
     integer :: id ! ID for this system
     integer :: status = sys_empty ! current status
     logical :: hidden = .false. ! whether it is hidden in the tree view (filter)
     logical :: showfields = .false. ! whether to show the fields in the tree view
     type(crystalseed) :: seed ! generating seed
     logical :: has_field = .false. ! true if the seed has a field
     integer :: collapse ! 0 if independent, -1 if master-collapsed, -2 if master-extended, <n> if dependent on n
     type(c_ptr) :: thread_lock = c_null_ptr ! the lock for initialization of this system
     character(len=:), allocatable :: fullname ! full-path name
     logical :: renamed = .false. ! true if the system has been renamed
     integer :: idwin_plotscf = 0 ! window ID for the scf plot
     type(scene) :: sc ! scene for the system in the main view
  end type sysconf
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
  public :: set_default_ui_settings

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
     module subroutine add_systems_from_seeds(nseed,seed,collapse,iafield)
       integer, intent(in) :: nseed
       type(crystalseed), allocatable, intent(in) :: seed(:)
       logical, intent(in), optional :: collapse
       integer, intent(in), optional :: iafield
     end subroutine add_systems_from_seeds
     module subroutine add_systems_from_name(name,mol,isformat,readlastonly,rborder,molcubic)
       character(len=*), intent(in) :: name
       integer, intent(in) :: mol
       integer, intent(in) :: isformat
       logical, intent(in) :: readlastonly
       real*8, intent(in) :: rborder
       logical, intent(in) :: molcubic
     end subroutine add_systems_from_name
     recursive module subroutine remove_system(idx)
       integer, intent(in) :: idx
     end subroutine remove_system
     module subroutine set_default_ui_settings()
     end subroutine set_default_ui_settings
  end interface

end module gui_main
