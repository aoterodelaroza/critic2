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
submodule (display) proc
  implicit none

contains

  !> Defaults of the Display for system isys: one cell, no shifts;
  !> in a crystal, atoms at the cell faces drawn on both
  !> faces, and whole molecules if there is more than one fragment or a
  !> non-discrete fragment carries reconnection lattice vectors
  !> (dangling pieces split by the cell boundary). Everything shown.
  module subroutine scene_display_init(disp,isys)
    use systems, only: sys, sys_ready, ok_system, atlisttype_species
    class(scene_display), intent(inout) :: disp
    integer, intent(in) :: isys

    integer :: imol, iat

    ! the type's defaults (a passed-object dummy is polymorphic, so the
    ! constructor cannot be assigned to it)
    call disp%end()
    disp%atype = atlisttype_species
    disp%ncell = 1
    disp%origin = 0d0
    disp%tshift = 0d0
    disp%border = .false.
    disp%onemotif = .false.
    if (.not.ok_system(isys,sys_ready)) return
    if (.not.sys(isys)%c%ismolecule) then
       disp%border = .true.
       disp%onemotif = (sys(isys)%c%nmol > 1)
       if (.not.disp%onemotif) then
          moldangler: do imol = 1, sys(isys)%c%nmol
             if (sys(isys)%c%mol(imol)%discrete) cycle
             do iat = 1, sys(isys)%c%mol(imol)%nat
                if (any(sys(isys)%c%mol(imol)%at(iat)%lvec /= 0)) then
                   disp%onemotif = .true.
                   exit moldangler
                end if
             end do
          end do moldangler
       end if
    end if
    call disp%reset_shown(isys)

  end subroutine scene_display_init

  !> (Re)make the atom-group and molecule Show masks of the Display for
  !> system isys, with everything shown, and stamp them.
  module subroutine scene_display_reset_shown(disp,isys)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_ready, ok_system
    class(scene_display), intent(inout) :: disp
    integer, intent(in) :: isys

    if (allocated(disp%ashown)) deallocate(disp%ashown)
    if (allocated(disp%mshown)) deallocate(disp%mshown)
    disp%timelastreset = glfwGetTime()
    if (.not.ok_system(isys,sys_ready)) return

    allocate(disp%ashown(sysc(isys)%attype_number(disp%atype)),disp%mshown(sys(isys)%c%nmol))
    disp%ashown = .true.
    disp%mshown = .true.

  end subroutine scene_display_reset_shown

  !> Remake the Show masks if the system changed under them: the masks
  !> do not exist yet, the geometry or the bonds (which define the
  !> molecules) moved since they were made, or the group and molecule
  !> counts no longer agree with the system.
  module subroutine scene_display_update(disp,isys)
    use systems, only: sys, sysc, sys_ready, ok_system
    class(scene_display), intent(inout) :: disp
    integer, intent(in) :: isys

    logical :: doreset

    if (.not.ok_system(isys,sys_ready)) return
    doreset = .not.allocated(disp%ashown) .or. .not.allocated(disp%mshown)
    doreset = doreset .or. (sysc(isys)%timelastchange_rebond > disp%timelastreset)
    if (.not.doreset) doreset = (size(disp%ashown,1) /= sysc(isys)%attype_number(disp%atype)) .or.&
       (size(disp%mshown,1) /= sys(isys)%c%nmol)
    if (doreset) call disp%reset_shown(isys)

  end subroutine scene_display_update

  !> Deallocate the Display.
  module subroutine scene_display_end(disp)
    class(scene_display), intent(inout) :: disp

    disp%timelastreset = 0d0
    if (allocated(disp%ashown)) deallocate(disp%ashown)
    if (allocated(disp%mshown)) deallocate(disp%mshown)

  end subroutine scene_display_end

  !> Take the periodic settings of another Display src (cells, shifts,
  !> border, whole molecules), leaving the Show masks (sized for this
  !> system) alone. Meant for two crystals.
  module subroutine scene_display_copy_settings(disp,src)
    class(scene_display), intent(inout) :: disp
    type(scene_display), intent(in) :: src

    disp%ncell = src%ncell
    disp%origin = src%origin
    disp%tshift = src%tshift
    disp%border = src%border
    disp%onemotif = src%onemotif

  end subroutine scene_display_copy_settings

  !> Number of unit cells an object is drawn over: the Display's count
  !> (follow), one (main cell only), or the object's own count (manual),
  !> from its override ovr. The single decoder of pertype.
  module function scene_display_ncells(disp,ovr) result(n)
    class(scene_display), intent(in) :: disp
    type(rep_display), intent(in) :: ovr
    integer :: n(3)

    if (ovr%pertype == pertype_none) then
       n = 1
    elseif (ovr%pertype == pertype_manual) then
       n = ovr%ncell
    else
       n = disp%ncell
    end if

  end function scene_display_ncells

end submodule proc
