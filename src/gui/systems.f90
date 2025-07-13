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

! This module handles the systems and system configurations.
module systems
  use iso_c_binding
  use scenes, only: scene
  use systemmod, only: system
  use crystalseedmod, only: crystalseed
  use global, only: bondfactor_def
  use types, only: thread_info
  use param, only: maxzat0, atmcov0
  implicit none

  private

  ! threads in execution
  integer, parameter, public :: nthread = 1
  type(c_ptr), target, allocatable, public :: thread(:)
  type(thread_info), target, allocatable, public :: thread_ti(:)

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
  public :: regenerate_system_pointers
  public :: ok_system

  !xx! Interfaces
  interface
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
     module subroutine reread_system_from_file(idx)
       integer, intent(in) :: idx
     end subroutine reread_system_from_file
     module subroutine duplicate_system(idx)
       integer, intent(in) :: idx
     end subroutine duplicate_system
     module subroutine regenerate_system_pointers()
     end subroutine regenerate_system_pointers
     module function ok_system(isys,level)
       integer, intent(in) :: isys, level
       logical :: ok_system
     end function ok_system
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

end module systems
