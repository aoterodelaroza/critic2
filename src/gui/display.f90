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

! scene display object: which part of the system is drawn
module display
  use iso_c_binding
  implicit none

  private

  ! Periodicity of an object relative to the Display object of its
  ! view (the radio buttons of iw_periodicity_widget, and the per-object
  ! override consumed by the representations)
  integer(c_int), parameter, public :: pertype_follow = 0 ! draw the cells the Display draws
  integer(c_int), parameter, public :: pertype_none = 1 ! draw the main cell only
  integer(c_int), parameter, public :: pertype_manual = 2 ! draw a count of its own

  !> scene display: which part of the system is drawn
  type scene_display
     integer(c_int) :: ncell(3) = 1 ! number of unit cells drawn along a, b, c
     real*8 :: origin(3) = 0d0 ! translation of the cell contents (fractional; Å in molecules)
     real*8 :: tshift(3) = 0d0 ! origin of the display region (fractional)
     logical :: border = .false. ! atoms at the cell faces are drawn on both faces
     logical :: onemotif = .false. ! translate atoms to display whole molecules
     ! which atoms and molecules are drawn (the masks exist once allocated)
     integer :: atype = 1 ! grouping of the atom Show mask (atlisttype_* in systems; 1 = species)
     logical, allocatable :: ashown(:) ! atom group shown, by group of atype
     logical, allocatable :: mshown(:) ! molecule shown
     real*8 :: timelastreset = 0d0 ! time the show masks were last (re)made
   contains
     procedure :: init => scene_display_init ! defaults for a system
     procedure :: reset_shown => scene_display_reset_shown ! (re)make the Show masks, everything shown
     procedure :: update => scene_display_update ! remake the Show masks if the system changed under them
     procedure :: end => scene_display_end
     procedure :: copy_settings => scene_display_copy_settings ! take another Display's periodic settings (not its masks)
     procedure :: ncells => scene_display_ncells ! cells an object draws, honoring its override
  end type scene_display
  public :: scene_display

  !> Per-object periodicity override of the scene Display (accessed as
  !> r%disp%...): whether the object draws the cells the Display draws,
  !> the main cell only, or a count of its own. ignoresel is for the
  !> transient representations that must draw a fixed atom set
  !> regardless of the Display (the packing preview): the whole cell
  !> contents with a border, no whole molecules, no
  !> display-region shift, no Show masks; the cell count and the origin
  !> translation are kept, so the spheres stay on the atoms.
  type rep_display
     integer(c_int) :: pertype = pertype_follow ! pertype_follow / pertype_none / pertype_manual
     integer(c_int) :: ncell(3) = 1 ! number of unit cells drawn (pertype_manual)
     logical :: ignoresel = .false. ! transients only: draw the whole cell contents (see above)
  end type rep_display
  public :: rep_display

  ! module procedure interfaces
  interface
     module subroutine scene_display_init(disp,isys)
       class(scene_display), intent(inout) :: disp
       integer, intent(in) :: isys
     end subroutine scene_display_init
     module subroutine scene_display_reset_shown(disp,isys)
       class(scene_display), intent(inout) :: disp
       integer, intent(in) :: isys
     end subroutine scene_display_reset_shown
     module subroutine scene_display_update(disp,isys)
       class(scene_display), intent(inout) :: disp
       integer, intent(in) :: isys
     end subroutine scene_display_update
     module subroutine scene_display_end(disp)
       class(scene_display), intent(inout) :: disp
     end subroutine scene_display_end
     module subroutine scene_display_copy_settings(disp,src)
       class(scene_display), intent(inout) :: disp
       type(scene_display), intent(in) :: src
     end subroutine scene_display_copy_settings
     module function scene_display_ncells(disp,ovr) result(n)
       class(scene_display), intent(in) :: disp
       type(rep_display), intent(in) :: ovr
       integer :: n(3)
     end function scene_display_ncells
  end interface

end module display
