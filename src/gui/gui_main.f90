! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
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
  use iso_c_binding, only: c_ptr, c_float, c_int, c_null_ptr
  use gui_interfaces_cimgui, only: ImGuiIO, ImGuiContext, ImVec4
  use gui_window, only: window
  use systemmod, only: system
  use crystalseedmod, only: crystalseed
  implicit none

  private

  ! variables to GUI's structures & data
  real*8, public :: time ! the time
  type(ImGuiIO), pointer, public :: io ! pointer to ImGui's IO object
  type(ImGuiContext), pointer, public :: g ! pointer to ImGui's context
  type(c_ptr), public :: rootwin ! the root window pointer (GLFWwindow*)

  ! GUI control parameters
  real(c_float), public :: tooltip_delay = 0.5 ! tooltip delay, in seconds

  ! GUI colors
  type(ImVec4), parameter, public :: TableCellBg_Mol     = ImVec4(0.43,0.8 ,0.  ,0.3)  ! tree table name cell, molecule
  type(ImVec4), parameter, public :: TableCellBg_MolClus = ImVec4(0.0 ,0.8 ,0.43,0.3)  ! tree table name cell, molecular cluster
  type(ImVec4), parameter, public :: TableCellBg_MolCrys = ImVec4(0.8 ,0.43,0.0 ,0.3)  ! tree table name cell, molecular crystal
  type(ImVec4), parameter, public :: TableCellBg_Crys3d  = ImVec4(0.8 ,0.  ,0.0 ,0.3)  ! tree table name cell, 3d crystal
  type(ImVec4), parameter, public :: TableCellBg_Crys2d  = ImVec4(0.8 ,0.  ,0.43,0.3)  ! tree table name cell, 2d crystal
  type(ImVec4), parameter, public :: TableCellBg_Crys1d  = ImVec4(0.8 ,0.43,0.43,0.3)  ! tree table name cell, 1d crystal
  type(ImVec4), parameter, public :: DialogDir = ImVec4(0.9, 0.9, 0.5, 1.0) ! directories in the dialog
  type(ImVec4), parameter, public :: DialogFile = ImVec4(1.0, 1.0, 1.0, 1.0) ! files in the dialog

  ! systems arrays
  integer, parameter, public :: sys_empty = 0
  integer, parameter, public :: sys_loaded_not_init = 1
  integer, parameter, public :: sys_initializing = 2
  integer, parameter, public :: sys_init = 3
  integer, parameter, public :: sys_loaded_not_init_hidden = 4
  integer, parameter, public :: sys_init_hidden = 5
  type :: sysconf
     integer :: id
     integer :: status = sys_empty
     type(crystalseed) :: seed
     logical :: has_field
     integer :: collapse ! 0 if independent, -1 if master-collapsed, -2 if master-extended, <n> if dependent on n
     type(c_ptr) :: thread_lock = c_null_ptr
     character(len=:), allocatable :: fullname ! full-path name
  end type sysconf
  integer, public :: nsys = 0
  type(system), allocatable, target, public :: sys(:)
  type(sysconf), allocatable, target, public :: sysc(:)

  ! public procedures
  public :: gui_start
  public :: launch_initialization_thread
  public :: system_shorten_names
  public :: add_systems_from_name
  public :: remove_system

  interface
     module subroutine gui_start()
     end subroutine gui_start
     module subroutine launch_initialization_thread()
     end subroutine launch_initialization_thread
     module subroutine system_shorten_names()
     end subroutine system_shorten_names
     module subroutine add_systems_from_name(name,mol,isformat)
       character(len=*), intent(in) :: name
       integer, intent(in) :: mol
       integer, intent(in) :: isformat
     end subroutine add_systems_from_name
     recursive module subroutine remove_system(idx)
       integer, intent(in) :: idx
     end subroutine remove_system
  end interface

contains

end module gui_main
