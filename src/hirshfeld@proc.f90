! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Hirshfeld integration
submodule (hirshfeld) proc
  implicit none

  integer, parameter :: mprops = 1

contains

  !> Set the attractors for Hirshfeld integration. The actual
  !> integration is done using hirsh_weights.
  module subroutine hirsh_grid(s,bas)
    use systemmod, only: system
    use types, only: basindat
    use tools_io, only: ferror, faterr
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    integer :: i

    if (.not.s%isinit) &
       call ferror("hirsh_grid","system not initialized",faterr)
    if (.not.allocated(s%c)) &
       call ferror("hirsh_grid","system does not have crystal",faterr)

    ! Atoms are the attractors in this case
    allocate(bas%xattr(3,s%f(s%iref)%ncpcel))
    bas%xattr = 0d0
    if (bas%atexist) then
       bas%nattr = s%f(s%iref)%ncpcel
       do i = 1, s%f(s%iref)%ncpcel
          bas%xattr(:,i) = s%f(s%iref)%cpcel(i)%x
       end do
    else
       bas%nattr = 0
    end if

  end subroutine hirsh_grid

  !> For system s, calculate the hirshfeld weights for atom idb
  !> (complete list) on a grid and return it in w. The size of w bas%n
  !> determines the size of the grid and bas%f must contain the
  !> promolecular density. The size of w must be consistent with
  !> bas%n.
  module subroutine hirsh_weights(s,bas,idb,w)
    use systemmod, only: system
    use grid3mod, only: grid3
    use fragmentmod, only: fragment
    use types, only: basindat
    use param, only: VSMALL
    type(system), intent(inout) :: s
    type(basindat), intent(in) :: bas
    integer, intent(in) :: idb
    real*8, intent(out) :: w(:,:,:)

    type(fragment) :: fr
    real*8, allocatable :: faux(:,:,:)

    ! initialize
    w = 0d0
    fr%nat = 1
    fr%nspc = 1
    allocate(fr%at(1),fr%spc(1))

    ! Prepare a fragment with this atom
    fr%at(1)%x = s%c%atcel(idb)%x
    fr%at(1)%r = s%c%atcel(idb)%r
    fr%at(1)%is = 1
    fr%at(1)%cidx = idb
    fr%at(1)%idx = s%c%atcel(idb)%idx
    fr%at(1)%lvec = 0
    fr%spc(1) = s%c%spc(s%c%atcel(idb)%is)

    ! Calculate grid with the atomic density. The sum over cells
    ! turns into a sum over periodic copies of the atom.  call
    call s%c%promolecular_array3(faux,bas%n,fr=fr)

    ! hirshfeld weights
    w = faux / max(bas%f,VSMALL)

  end subroutine hirsh_weights

  !> Set the attractors for Voronoi integration and calculate
  !> the assignments of nodes to nuclei (bas%idg). The size
  !> of the grid is given by bas%n.
  module subroutine voronoi_grid(s,bas)
    use systemmod, only: system
    use tools_io, only: ferror, faterr
    use types, only: basindat
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas

    integer :: i

    if (.not.s%isinit) &
       call ferror("voronoi_grid","system not initialized",faterr)
    if (.not.allocated(s%c)) &
       call ferror("voronoi_grid","system does not have crystal",faterr)

    ! Atoms are the attractors in this case
    allocate(bas%xattr(3,s%f(s%iref)%ncpcel))
    bas%xattr = 0d0
    if (bas%atexist) then
       bas%nattr = s%f(s%iref)%ncpcel
       do i = 1, s%f(s%iref)%ncpcel
          bas%xattr(:,i) = s%f(s%iref)%cpcel(i)%x
       end do
    else
       bas%nattr = 0
    end if

    ! assign grid nodes to atoms
    call s%c%nearest_atom_grid(bas%n,bas%idg)

  end subroutine voronoi_grid

  ! Calculate hirshfeld populations and volumes using a mesh.
  module subroutine hirsh_nogrid()
    use fieldmod, only: field
    use meshmod, only: mesh
    use global, only: mesh_type, mesh_level
    use fragmentmod, only: fragment
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_center
    use param, only: VSMALL, im_rho

    type(field) :: fat
    type(mesh) :: m0, mat, mrho
    integer :: i, iat, cidx
    real*8 :: vsum, asum, ntotal, vtotal
    real*8, allocatable :: nat(:), vat(:)
    type(fragment) :: fr
    integer :: prop(1)

    ! header
    write (uout,'("* Hirshfeld atomic electron populations (using mesh integration)")')
    write (uout,'("  Reference field: ",A)') string(sy%iref)

    ! generate the meshes
    call m0%gen(sy%c,mesh_type,mesh_level)
    call mrho%gen(sy%c,mesh_type,mesh_level)
    call mat%gen(sy%c,mesh_type,mesh_level)

    ! report
    write (uout,'("+ Mesh details")')
    call m0%report()

    ! Initialize accumulation and fragment
    allocate(nat(sy%c%nneq),vat(sy%c%nneq))
    ntotal = 0d0
    vtotal = 0d0
    nat = 0d0
    vat = 0d0
    fr%nat = 1
    fr%nspc = 1
    allocate(fr%at(1),fr%spc(1))

    ! calculate the density and promolecular density on the mesh
    prop(1) = im_rho
    call mrho%fill(sy%f(sy%iref),prop,.not.sy%c%ismolecule)
    call m0%fill(sy%f(0),prop,.not.sy%c%ismolecule)

    do iat = 1, sy%c%nneq
       ! Fetch a representative with that iat
       cidx = 0
       do i = 1, sy%c%ncel
          if (sy%c%atcel(i)%idx == iat) then
             cidx = i
             exit
          end if
       end do

       ! Prepare a fragment with this atom
       fr%at(1)%x = sy%c%atcel(cidx)%x
       fr%at(1)%r = sy%c%atcel(cidx)%r
       fr%at(1)%is = 1
       fr%at(1)%idx = iat
       fr%at(1)%cidx = cidx
       fr%at(1)%lvec = 0
       fr%spc(1) = sy%c%spc(sy%c%at(iat)%is)

       ! Load the promolecular field with that fragment
       call fat%load_promolecular(sy%c,-1,"",fr)

       ! Calculate mesh with the atomic density. The sum over cells
       ! turns into a sum over periodic copies of the atom.
       prop(1) = im_rho
       call mat%fill(fat,prop,.not.sy%c%ismolecule)

       ! hirshfeld weights (with mesh weights)
       mat%f(:,1) = mat%f(:,1) / max(m0%f(:,1),VSMALL) * mrho%w

       ! accumulate
       asum = sum(mrho%f(:,1) * mat%f(:,1))
       vsum = sum(mat%f(:,1))
       nat(iat) = asum
       ntotal = ntotal + asum * sy%c%at(iat)%mult
       vat(iat) = vsum
       vtotal = vtotal + vsum * sy%c%at(iat)%mult
    end do

    ! write the results
    write (uout,'("+ Hirshfeld integration results")')
    write (uout,'("# N_hirsh = atomic electron populations")')
    write (uout,'("# V_hirsh = atomic volumes ")')
    write (uout,'("#nneq mult name   N_hirsh          V_hirsh")')
    do iat = 1, sy%c%nneq
       write (uout,'(5(A," "))') string(iat,length=4,justify=ioj_center), &
          string(sy%c%at(iat)%mult,length=4,justify=ioj_center), &
          string(sy%c%at(iat)%name,length=5,justify=ioj_center), &
          string(nat(iat),'f',length=16,decimal=10,justify=3), &
          string(vat(iat),'f',length=16,decimal=10,justify=3)
    end do
    write (uout,'("# total number of electrons: ",A)') string(ntotal,'e',decimal=10)
    write (uout,'("# total volume: ",A)') string(vtotal,'e',decimal=10)
    write (uout,*)

  end subroutine hirsh_nogrid

end submodule proc
