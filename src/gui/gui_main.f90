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
  use iso_c_binding, only: c_ptr, c_float, c_int
  use gui_interfaces_cimgui, only: ImGuiIO, ImGuiContext
  use gui_window, only: window
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
  integer :: nwin
  type(window), allocatable, target :: win(:)
  integer :: iwin_tree
  integer :: iwin_console
  integer :: iwin_view

  ! the dockspace ID
  integer(c_int) :: iddock = 0

  ! public procedures
  public :: gui_start

  interface
     module subroutine gui_start()
     end subroutine gui_start
  end interface

contains

end module gui_main
