! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> Promolecular and core density calculation using radial grids.
module grd_atomic
  use grid1_tools
  use types
  implicit none

  private

  public :: grda_init
  public :: grda_end
  public :: grda_promolecular

  type(grid1), target, allocatable, public :: agrid(:)
  type(grid1), target, allocatable, public :: cgrid(:,:)

contains

  !> Initialize the atomic density grids for the promolecular densities
  subroutine grda_init(docore,dopro,verbose)
    use struct_basic
    use grid1_tools
    use tools_io
    use param

    logical, intent(in) :: docore, dopro, verbose

    integer :: i, j, iz, iq, nval, nelec, atq
    

    ! allocate memory for all atomic types
    if (docore) then
       if (allocated(cgrid)) deallocate(cgrid)
       allocate(cgrid(maxzat0,maxzat0))
    endif
    if (dopro) then
       if (allocated(agrid)) deallocate(agrid)
       allocate(agrid(maxzat0))
    endif

    ! initialize default agrid and cgrid values
    do i = 1, maxzat0
       if (dopro) then
          agrid(i)%z = 0
          agrid(i)%qat = 0
       end if
       if (docore) then
          do j = 1, maxzat0
             cgrid(i,j)%z = 0
             cgrid(i,j)%qat = 0
          end do
       endif
    end do

    if (verbose) then
       write (uout,'("* Atomic radial density grids")')
       write (uout,'("+ List of atomic charges and atomic numbers")')
       write (uout,'("# ",5(A,2X))') &
          string("nat",length=3,justify=ioj_right), &
          string("name",length=5,justify=ioj_center), &
          string("Z",length=2,justify=ioj_right), &
          string("Q",length=4,justify=ioj_right), &
          string("ZPSP",length=3,justify=ioj_right)
       nelec = 0
       nval = 0
       do i = 1, cr%nneq
          write (uout,'(5(2X,A))') &
             string(i,length=3,justify=ioj_right), &
             string(cr%at(i)%name,length=5,justify=ioj_center), &
             string(cr%at(i)%z,length=2,justify=ioj_right), &
             string(cr%at(i)%qat,'f',length=4,decimal=1,justify=ioj_right), &
             string(cr%at(i)%zpsp,length=3,justify=ioj_right)
          nelec = nelec + cr%at(i)%mult * cr%at(i)%z
          if (cr%at(i)%zpsp >= 0) then
             nval = nval + cr%at(i)%mult * cr%at(i)%zpsp
          else
             nval = nval + cr%at(i)%mult * cr%at(i)%z
          end if
       end do
       write (uout,*)
       write (uout,'("+ Number of electrons: ",A)') string(nelec)
       write (uout,'("+ Number of valence electrons: ",A)') string(nval)
       write (uout,'("+ Number of core electrons: ",A)') string(nelec-nval)
       write (uout,*)
    end if
    
    if (dopro) then
       if (verbose) &
          write (uout,'("* Reading new promolecular density grids")')
       do i = 1, cr%nneq
          if (cr%at(i)%z <= 0 .or. cr%at(i)%z > maxzat) cycle
          if (abs(nint(cr%at(i)%qat) - cr%at(i)%qat) < 1d-10) then
             atq = nint(cr%at(i)%qat)
          else
             atq = 0
          end if
          if (agrid(cr%at(i)%z)%init .and. agrid(cr%at(i)%z)%z == cr%at(i)%z .and.&
             agrid(cr%at(i)%z)%qat == atq) cycle
          call grid1_read_db(agrid(cr%at(i)%z),cr%at(i)%z,atq,verbose)
       end do
       if (verbose) write (uout,*)
    end if

    if (docore) then
       if (verbose) &
          write (uout,'("* Reading new core density grids")')
       do i = 1, cr%nneq
          if (cr%at(i)%zpsp <= 0 .or. cr%at(i)%zpsp > maxzat) cycle
          iz = cr%at(i)%z
          iq = cr%at(i)%zpsp
          if (cgrid(iz,iq)%init .and. cgrid(iz,iq)%z == iz .and.&
             cgrid(iz,iq)%qat == iq) cycle
          call grid1_read_db(cgrid(iz,iq),iz,iq,verbose)
       end do
       if (verbose) write (uout,*)
    endif

  end subroutine grda_init

  !> Deallocate the arrays for the promolecular and core densities
  subroutine grda_end()
    if (allocated(agrid)) deallocate(agrid)
    if (allocated(cgrid)) deallocate(cgrid)
  end subroutine grda_end

  !> Calculate the core or promolecular densities at a point x0
  !> (Cartesian coords) using atomic radial grids. If a fragment is
  !> given, then only the atoms in it contribute.
  subroutine grda_promolecular(x0,f,fp,fpp,nder,iscore,fr,periodic)
    use struct_basic
    use grid1_tools
    use types
    use param

    real*8, intent(in) :: x0(3) !< Point in cryst. coords.
    real*8, intent(out) :: f !< Density
    real*8, intent(out) :: fp(3) !< Density gradient
    real*8, intent(out) :: fpp(3,3) !< Density hessian
    integer, intent(in) :: nder !< Number of derivatives to calculate
    logical, intent(in) :: iscore !< Use atomic cores or promolecular densities?
    type(fragment), intent(in), optional :: fr !< Fragment contributing to the density
    logical, intent(in), optional :: periodic

    integer :: i, j, k, ii
    real*8 :: xc(3), xx(3), r2, r, rinv1, rinv2
    real*8 :: rho, rhop, rhopp, rfac, radd
    integer :: idolist(cr%nenv), nido

    type(grid1), pointer :: g

    ! check that we have an environment
    call cr%checkflags(.true.,init0=.true.,env0=.true.)

    ! initialize 
    xc = x0
    if (present(periodic)) then
       if (periodic) then
          xc = cr%c2x(x0)
          xc = xc - floor(xc)
          xc = cr%x2c(xc)
       end if
    end if
    f = 0d0 
    fp = 0d0
    fpp = 0d0

    ! precompute the list of atoms that contribute
    if (.not.present(fr)) then
       nido = 0
       do i = 1, cr%nenv
          if (cr%at(cr%atenv(i)%idx)%z == 0 .or. cr%at(cr%atenv(i)%idx)%z > maxzat) cycle
          if (iscore .and. (cr%at(cr%atenv(i)%idx)%zpsp <= 0)) cycle
          if (iscore .and. (cr%at(cr%atenv(i)%idx)%z - cr%at(cr%atenv(i)%idx)%zpsp == 0)) cycle
          if (iscore) then
             g => cgrid(cr%at(cr%atenv(i)%idx)%z,cr%at(cr%atenv(i)%idx)%zpsp)
          else
             g => agrid(cr%at(cr%atenv(i)%idx)%z)
          end if
          if (.not.g%init) cycle
          xx = xc - cr%atenv(i)%r
          r2 = dot_product(xx,xx)
          if (r2 > g%rmax2) cycle
          nido = nido + 1
          idolist(nido) = i
       end do
    else
       nido = fr%nat
    end if

    ! do the sum
    do ii = 1, nido
       if (.not.present(fr)) then
          i = idolist(ii)
          if (iscore) then
             g => cgrid(cr%at(cr%atenv(i)%idx)%z,cr%at(cr%atenv(i)%idx)%zpsp)
          else
             g => agrid(cr%at(cr%atenv(i)%idx)%z)
          end if
          xx = xc - cr%atenv(i)%r
       else
          if (iscore) then
             g => cgrid(fr%at(ii)%z,cr%at(fr%at(ii)%idx)%zpsp)
          else
             g => agrid(fr%at(ii)%z)
          end if
          xx = xc - fr%at(ii)%r
       end if

       r2 = dot_product(xx,xx)
       r = max(sqrt(r2),g%r(1))
       r = max(r,1d-14)
       call grid1_interp(g,r,rho,rhop,rhopp)
       rho = max(rho,0d0)

       f = f + rho
       if (nder < 1) cycle
       rinv1 = 1d0 / r
       fp = fp + rhop * xx * rinv1
       if (nder < 2) cycle
       rinv2 = rinv1 * rinv1
       rfac = (rhopp - rhop * rinv1)
       do j = 1, 3
          fpp(j,j) = fpp(j,j) + rhop * rinv1 + rfac * rinv2 * xx(j) * xx(j)
          do k = 1, j-1
             radd = rfac * rinv2 * xx(j) * xx(k)
             fpp(j,k) = fpp(j,k) + radd
             fpp(k,j) = fpp(k,j) + radd
          end do
       end do
    end do

  end subroutine grda_promolecular

end module grd_atomic
