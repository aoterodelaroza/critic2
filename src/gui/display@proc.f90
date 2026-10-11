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
  !> system isys, with everything shown, and the lattice-vector shifts of
  !> the atom groups, all zero; they are for the bonding revision of the
  !> system now.
  module subroutine scene_display_reset_shown(disp,isys)
    use systems, only: sys, sysc, sys_ready, ok_system
    class(scene_display), intent(inout) :: disp
    integer, intent(in) :: isys

    if (allocated(disp%ashown)) deallocate(disp%ashown)
    if (allocated(disp%ashift)) deallocate(disp%ashift)
    if (allocated(disp%mshown)) deallocate(disp%mshown)
    if (.not.ok_system(isys,sys_ready)) return
    disp%rev = sysc(isys)%rev_rebond

    allocate(disp%ashown(sysc(isys)%attype_number(disp%atype)),disp%mshown(sys(isys)%c%nmol),&
       disp%ashift(3,sysc(isys)%attype_number(disp%atype)))
    disp%ashown = .true.
    disp%ashift = 0
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
    doreset = doreset .or. (disp%rev /= sysc(isys)%rev_rebond)
    if (.not.doreset) doreset = (size(disp%ashown,1) /= sysc(isys)%attype_number(disp%atype)) .or.&
       (size(disp%mshown,1) /= sys(isys)%c%nmol)
    if (doreset) call disp%reset_shown(isys)

  end subroutine scene_display_update

  !> The rows of the atom Show mask and lattice-vector shifts of the
  !> Display (system isys) that hold the cell atoms iat(:), in irow. If
  !> the cell atoms are not whole groups of them, the rows become the
  !> cell atoms first, each shown and shifted as its group was, so that
  !> the rows change only these atoms.
  module subroutine scene_display_rows_for(disp,isys,iat,irow)
    use systems, only: sysc
    class(scene_display), intent(inout) :: disp
    integer, intent(in) :: isys
    integer, intent(in) :: iat(:)
    integer, allocatable, intent(inout) :: irow(:)

    integer, allocatable :: k(:)

    call disp%update(isys)
    if (.not.allocated(disp%ashown)) then
       irow = (/integer ::/)
       return
    end if
    if (.not.sysc(isys)%attype_rows_whole(disp%atype,iat)) then
       k = sysc(isys)%attype_celatom_ids(disp%atype)
       disp%ashown = disp%ashown(k)
       disp%ashift = disp%ashift(:,k)
       disp%atype = sysc(isys)%attype_celatom_type()
    end if
    irow = sysc(isys)%attype_rows(disp%atype,iat)

  end subroutine scene_display_rows_for

  !> Deallocate the Display.
  module subroutine scene_display_end(disp)
    class(scene_display), intent(inout) :: disp

    disp%rev = -1
    if (allocated(disp%ashown)) deallocate(disp%ashown)
    if (allocated(disp%ashift)) deallocate(disp%ashift)
    if (allocated(disp%mshown)) deallocate(disp%mshown)

  end subroutine scene_display_end

  !> Take the periodic settings of another Display src (cells, shifts,
  !> border, whole molecules), leaving the Show masks (sized for this
  !> system) alone. Meant for two crystals.
  module subroutine scene_display_copy_settings(disp,src)
    class(scene_display), intent(inout) :: disp
    type(scene_display), intent(in) :: src

    disp%ncell = src%ncell
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

  !> Write the Display to a new JSON object "display" under p: the
  !> periodic settings and the Show masks, if made.
  module subroutine scene_display_write_json(disp,json,p)
    use json_module, only: json_value, json_core
    class(scene_display), intent(in) :: disp
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(inout) :: p

    type(json_value), pointer :: s

    call json%create_object(s,'display')
    call json%add(p,s)
    call json%add(s,'ncell',int(disp%ncell))
    call json%add(s,'tshift',disp%tshift)
    call json%add(s,'border',disp%border)
    call json%add(s,'onemotif',disp%onemotif)
    call json%add(s,'atype',disp%atype)
    if (allocated(disp%ashown) .and. allocated(disp%ashift) .and. allocated(disp%mshown)) then
       call json%add(s,'ashown',disp%ashown)
       call json%add(s,'ashift',reshape(disp%ashift,(/size(disp%ashift)/)))
       call json%add(s,'mshown',disp%mshown)
    end if

  end subroutine scene_display_write_json

  !> Read the Display from the JSON object p (written by
  !> scene_display_write_json) over the current one: keys that are
  !> missing are left as they are. The masks are taken only if they
  !> are consistent with each other, and are not stamped with a
  !> bonding revision (rev = -1): the caller stamps them once the
  !> system they belong to is built, and update() remakes them if
  !> their sizes do not fit it.
  module subroutine scene_display_read_json(disp,json,p)
    use json_module, only: json_value, json_core
    class(scene_display), intent(inout) :: disp
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(in) :: p

    integer, allocatable :: iv(:)
    real*8, allocatable :: rv(:)
    logical, allocatable :: la(:), lm(:)
    integer :: ival
    logical :: found, lv


    ! (the cell counts and the shifts size what the objects allocate:
    ! they are kept in range)
    call json%get(p,'ncell',iv,found)
    if (found) then
       if (size(iv) == 3) disp%ncell = int(min(max(iv,1),disp_maxncell),c_int)
    end if
    call json%get(p,'tshift',rv,found)
    if (found) then
       if (size(rv) == 3) disp%tshift = min(max(rv,-real(disp_maxshift,8)),real(disp_maxshift,8))
    end if
    call json%get(p,'border',lv,found)
    if (found) disp%border = lv
    call json%get(p,'onemotif',lv,found)
    if (found) disp%onemotif = lv
    call json%get(p,'atype',ival,found)
    if (found) disp%atype = ival

    ! the masks
    disp%rev = -1
    call json%get(p,'ashown',la,found)
    if (.not.found) return
    call json%get(p,'mshown',lm,found)
    if (.not.found) return
    call json%get(p,'ashift',iv,found)
    if (.not.found) return
    if (size(iv) /= 3*size(la)) return
    disp%ashown = la
    disp%mshown = lm
    disp%ashift = reshape(min(max(iv,-disp_maxshift),disp_maxshift),(/3,size(la)/))

  end subroutine scene_display_read_json

end submodule proc
