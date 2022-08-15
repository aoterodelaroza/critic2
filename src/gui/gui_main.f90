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
  use gui_interfaces_cimgui, only: ImGuiIO, ImGuiContext
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

  ! the window stack and named windows
  integer, public :: nwin
  type(window), allocatable, target, public :: win(:)
  integer, public :: iwin_tree
  integer, public :: iwin_console
  integer, public :: iwin_view

  ! systems arrays
  integer, parameter, public :: sys_empty = 0
  integer, parameter, public :: sys_loaded_not_init = 1
  integer, parameter, public :: sys_initializing = 2
  integer, parameter, public :: sys_init = 3
  type :: sysconf
     integer :: id
     integer :: status = sys_empty
     type(crystalseed) :: seed
     logical :: has_field
     type(c_ptr) :: thread_lock = c_null_ptr
  end type sysconf
  integer, public :: nsys = 0
  type(system), allocatable, target, public :: sys(:)
  type(sysconf), allocatable, target, public :: sysc(:)

  ! public procedures
  public :: gui_start

  interface
     module subroutine gui_start()
     end subroutine gui_start
     module subroutine launch_initialization_thread()
     end subroutine launch_initialization_thread
  end interface

contains

end module gui_main
