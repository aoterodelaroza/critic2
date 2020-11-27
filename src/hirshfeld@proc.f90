! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (hirshfeld) proc
  implicit none

  integer, parameter :: mprops = 1

contains

  ! calculate hirshfeld charges using a grid
  module subroutine hirsh_props_grid()
    use fragmentmod, only: fragment
    use systemmod, only: sy
    use grid3mod, only: grid3
    use grid1mod, only: grid1, agrid
    use fieldmod, only: type_grid
    use tools_io, only: ferror, faterr, uout, string, ioj_center
    use param, only: VSMALL

    integer :: i, j, k
    real*8 :: dist, rrho, rrho1, rrho2, x(3), xdelta(3,3)
    integer :: n(3)
    logical :: doagain
    integer :: ishl, il, ill, ivec(3), iat, cidx
    real*8 :: lvec(3), vsum, asum, ntotal, vtotal
    real*8, allocatable :: nat(:), vat(:)
    type(grid3) :: hw, hwat
    type(fragment) :: fr

    if (sy%f(sy%iref)%type /= type_grid .or. sy%c%ismolecule) &
       call ferror("hirsh_props_grid","grid hirshfeld only for grid fields in crystals",faterr)

    ! Header
    write (uout,'("* Hirshfeld atomic electron populations")')
    write (uout,'("  Reference field: ",A)') string(sy%iref)
    write (uout,'("  Grid size: ",3(A,X))') (string(sy%f(sy%iref)%grid%n(j)),j=1,3)
    write (uout,*)

    ! Calculate the promolecular density on the grid
    n = sy%f(sy%iref)%grid%n
    call sy%c%promolecular_grid(hw,n)

    ! Initialize accumulation and fragment
    allocate(nat(sy%c%nneq),vat(sy%c%nneq))
    ntotal = 0d0
    vtotal = 0d0
    nat = 0d0
    vat = 0d0
    fr%nat = 1
    fr%nspc = 1
    allocate(fr%at(1),fr%spc(1))

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

       ! Calculate grid with the atomic density. The sum over cells
       ! turns into a sum over periodic copies of the atom.
       call sy%c%promolecular_grid(hwat,n,fr=fr)

       ! hirshfeld weights
       hwat%f = hwat%f / max(hw%f,VSMALL)
       asum = sum(sy%f(sy%iref)%grid%f * hwat%f) * sy%c%omega / real(n(1)*n(2)*n(3),8)
       vsum = sum(hwat%f) * sy%c%omega / real(n(1)*n(2)*n(3),8)
       nat(iat) = asum
       ntotal = ntotal + asum * sy%c%at(iat)%mult
       vat(iat) = vsum
       vtotal = vtotal + vsum * sy%c%at(iat)%mult
    end do

    ! Write the results
    write (uout,'("#nneq mult name   N_hirsh          V_hirsh")')
    do iat = 1, sy%c%nneq
       write (uout,'(5(A,X))') string(iat,length=4,justify=ioj_center), &
          string(sy%c%at(iat)%mult,length=4,justify=ioj_center), &
          string(sy%c%spc(sy%c%at(iat)%is)%name,length=5,justify=ioj_center), &
          string(nat(iat),'f',length=16,decimal=10,justify=3), &
          string(vat(iat),'f',length=16,decimal=10,justify=3)
    end do
    write (uout,'("# total number of electrons: ",A)') string(ntotal,'e',decimal=10)
    write (uout,'("# total volume: ",A)') string(vtotal,'e',decimal=10)
    write (uout,*)

  end subroutine hirsh_props_grid

end submodule proc
