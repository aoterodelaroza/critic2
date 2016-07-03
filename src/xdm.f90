! Copyright (c) 2015 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
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

! This module contains code for the calculation of the dispersion
! energy (and derivatives) using the exchange-hole dipole moment
! (XDM) model.
module xdm
  implicit none

  private

  public :: xdm_driver
  private :: xdm_grid
  private :: xdm_qe
  private :: xdm_mol
  private :: write_cube
  private :: free_volume
  private :: frevol
  private :: edisp_mol

  integer, parameter :: chf_blyp = -1
  integer, parameter :: chf_b3lyp = -2
  integer, parameter :: chf_bhahlyp = -3
  integer, parameter :: chf_camb3lyp = -4
  integer, parameter :: chf_pbe = -5
  integer, parameter :: chf_pbe0 = -6
  integer, parameter :: chf_lcwpbe = -7
  integer, parameter :: chf_pw86 = -8
  integer, parameter :: chf_b971 = -9

contains

  !> Driver for XDM
  subroutine xdm_driver(line)
    use fields
    use global
    use tools_io
    
    character*(*), intent(inout) :: line

    character(len=:), allocatable :: word, chfw
    integer :: lp 
    real*8 :: a1, a2, chf
    logical :: ok

    ! header
    write (uout,'("* Exchange-hole dipole moment (XDM) dispersion calculation ")')
    write (uout,'("  Please cite: ")')
    write (uout,'("  A. D. Becke, E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).")')
    write (uout,'("  A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 136, 174109 (2012).")')
    write (uout,'("  A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).")')
    write (uout,*)

    ! parse the input
    lp = 1
    word = lgetword(line,lp)
    if (equal(word,"grid")) then
       call xdm_grid(line(lp:))
    elseif (equal(word,"qe")) then
       call xdm_qe()
    else
       lp = 1
       ok = eval_next(a1,line,lp)
       ok = ok .and. eval_next(a2,line,lp)
       if (.not.ok) &
          call ferror("xdm_driver","wrong a1 or a2 in XDM",faterr)

       chfw = lgetword(line,lp)
       if (equal(chfw,"blyp")) then
          chf = chf_blyp
       elseif (equal(chfw,"b3lyp")) then
          chf = chf_b3lyp
       elseif (equal(chfw,"bhandhlyp").or.equal(chfw,"bhandh")&
          .or.equal(chfw,"bhah") .or. equal(chfw,"bhahlyp")) then
          chf = chf_bhahlyp
       elseif (equal(chfw,"camb3lyp") .or. equal(chfw,"cam-b3lyp")) then
          chf = chf_camb3lyp
       elseif (equal(chfw,"pbe")) then
          chf = chf_pbe
       elseif (equal(chfw,"pbe0")) then
          chf = chf_pbe0
       elseif (equal(chfw,"lcwpbe").or.equal(chfw,"lc-wpbe")) then
          chf = chf_lcwpbe
       elseif (equal(chfw,"pw86").or.equal(chfw,"pw86pbe")) then
          chf = chf_pw86
       elseif (equal(chfw,"b971").or.equal(chfw,"b97-1")) then
          chf = chf_b971
       elseif (equal(chfw,"hf")) then
          chf = 1.0d0
       else
          lp = 1
          ok = eval_next(chf,chfw,lp)
          if (.not.ok) &
             call ferror("xdm_driver","unknown functional",faterr)
       endif

       call xdm_mol(a1,a2,chf)
    end if

  end subroutine xdm_driver

  !> Calculate XDM using grids.
  subroutine xdm_grid(line)
    use fields
    use struct_basic
    use grid_tools
    use grd_atomic
    use grid1_tools
    use global
    use tools_io
    use param

    character*(*), intent(inout) :: line

    logical :: ok, dopro, docor
    character(len=:), allocatable :: word
    integer :: lp, irho, itau, ielf, ipdens, icor, ilap, igrad, irhoae, ib
    integer :: i, j, k, l, n(3), ntot, ii, jj, nn, m1, m2
    real*8 :: x(3), rhoat, rhocore, rdum1(3), rdum2(3,3), maxc6
    real*8 :: rhos, rhot, grho, lap, taus, ds, qs, rhs, xroot, xshift, xold, expx
    real*8 :: gx, fx, ffx, rmax, rmax2, ri, rhofree, raux1, raux2, wei, db
    real*8 :: mll(3), avoll, a1, a2, c6, c8, c9, c10, rvdw, ri2, db2
    real*8 :: eat, sat(3,3), ri3, rix, fcom, cn0, fll, rvdwx, exx, fxx
    integer :: ll, iat, nvec, imax, jmax, kmax, i3, i3i, l1, l2
    real*8, allocatable :: rc(:,:), alpha(:), ml(:,:), avol(:), afree(:)
    integer, allocatable :: lvec(:,:), ityp(:)
    logical :: isdefa, onlyc
    real*8 :: etotal, sigma(3,3), for(3,cr%ncel), ehadd(6:10)
    integer :: nclean, iclean(8), upto

    real*8, parameter :: ecut = 1d-11

    write (uout,'("+ Using: grid fields.")')
    ! initialization and defaults
    lp = 1
    irho = refden
    itau = -1
    ielf = -1
    ipdens = -1
    icor = -1
    ilap = -1
    igrad = -1
    irhoae = -1
    ib = -1
    a1 = 0.6836 ! defaults for pw86pbe in QE, new version
    a2 = 1.5045
    isdefa = .true.
    nclean = 0
    upto = 10
    onlyc = .false.

    ! parse the input and check sanity
    do while(.true.)
       word = lgetword(line,lp)
       if (equal(word,"rho")) then
          word = getword(line,lp)
          irho = fieldname_to_idx(word)
          if (irho < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(irho)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"tau")) then
          word = getword(line,lp)
          itau = fieldname_to_idx(word)
          if (itau < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(itau)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"elf")) then
          word = getword(line,lp)
          ielf = fieldname_to_idx(word)
          if (ielf < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(ielf)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"pdens")) then
          word = getword(line,lp)
          ipdens = fieldname_to_idx(word)
          if (ipdens < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(ipdens)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"core")) then
          word = getword(line,lp)
          icor = fieldname_to_idx(word)
          if (icor < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(icor)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"lap")) then
          word = getword(line,lp)
          ilap = fieldname_to_idx(word)
          if (ilap < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(ilap)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"grad")) then
          word = getword(line,lp)
          igrad = fieldname_to_idx(word)
          if (igrad < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(igrad)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"rhoae")) then
          word = getword(line,lp)
          irhoae = fieldname_to_idx(word)
          if (irhoae < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(irhoae)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"xb")) then
          word = getword(line,lp)
          ib = fieldname_to_idx(word)
          if (ib < 0) then
             call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
             return
          end if
          if (.not.goodfield(ib)) then
             call ferror("xdm_driver","rho field not allocated",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"xa1")) then
          ok = eval_next(a1,line,lp)
          if (.not.ok) then
             call ferror("xdm_driver","wrong a1",faterr,line,syntax=.true.)
             return
          end if
          isdefa = .false.
       elseif (equal(word,"xa2")) then
          ok = eval_next(a2,line,lp)
          if (.not.ok) then
             call ferror("xdm_driver","wrong a2",faterr,line,syntax=.true.)
             return
          end if
          isdefa = .false.
       elseif (equal(word,"onlyc")) then
          onlyc = .true.
       elseif (equal(word,"upto")) then
          ok = eval_next(upto,line,lp)
          if (.not.ok) then
             call ferror("xdm_driver","wrong upto field",faterr,line,syntax=.true.)
             return
          end if
          if (upto/=6.and.upto/=8.and.upto/=10) then
             call ferror("xdm_driver","upto can only be 6, 8, or 10",faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror("xdm_driver","Unknown extra keyword",faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! Check that the zpsps are defined
    if (any(cr%at(:)%zpsp < 0)) then
       call ferror("xdm_driver","need ZPSP for all atoms for XDM",faterr,line,syntax=.true.)
       return
    end if

    ! Define damping coefficients
    write (uout,'("+ Damping function parameters")')
    write (uout,'("  a1 = ",A)') string(a1,'f',decimal=6)
    write (uout,'("  a2 (ang) = ",A)') string(a2,'f',decimal=6)
    a2 = a2 / bohrtoa
    if (isdefa) &
       call ferror("xdm_driver","Using default a1 and a2 parameters",warning)

    ! we need at least the density and the b or kinetic energy density or elf
    if (.not.goodfield(irho,type_grid)) then
       call ferror("xdm_driver","wrong rho field",faterr,line,syntax=.true.)
       return
    end if
    if (ib > 0) then
       if (.not.goodfield(ib,type_grid)) then
          call ferror("xdm_driver","wrong b field",faterr,line,syntax=.true.)
          return
       end if
    elseif (itau > 0) then
       if (.not.goodfield(itau,type_grid)) then
          call ferror("xdm_driver","wrong tau field",faterr,line,syntax=.true.)
          return
       end if
    elseif (ielf > 0) then
       if (.not.goodfield(ielf,type_grid)) then
          call ferror("xdm_driver","wrong elf field",faterr,line,syntax=.true.)
          return
       end if
       itau = getfieldnum()
       write (uout,'("+ Calculating tau from elf and rho")')
       call taufromelf(ielf,irho,itau)
       nclean = nclean + 1
       iclean(nclean) = itau
       write (uout,'("+ Writing kinetic energy density to: ",A)') trim(fileroot)//"-tau.cube"
       call write_cube(trim(fileroot)//"-tau.cube","Kinetic energy density","Written by critic2 for XDM",f(itau)%n,f(itau)%f)
    else
       call ferror("xdm_driver","no tau or elf field given",faterr,line,syntax=.true.)
       return
    endif

    ! initialize
    n = f(irho)%n
    ntot = n(1)*n(2)*n(3)

    ! allocate the promolecular density
    dopro = .false.
    if (ipdens < 0) then
       ipdens = getfieldnum()
       allocate(f(ipdens)%f(n(1),n(2),n(3)))
       f(ipdens)%n = n
       dopro = .true.
       nclean = nclean + 1
       iclean(nclean) = ipdens
    endif

    ! allocate the core density (only if we don't have b and rhoae)
    docor = .false.
    if (icor < 0 .and..not.(ib > 0 .and. irhoae > 0)) then
       icor = getfieldnum()
       allocate(f(icor)%f(n(1),n(2),n(3)))
       f(icor)%n = n
       docor = .true.
       nclean = nclean + 1
       iclean(nclean) = icor
    endif

    ! calculate the promolecular and core density
    if (dopro .or. docor) then
       if (docor) &
          write (uout,'("+ Calculating core density")')
       if (dopro) &
          write (uout,'("+ Calculating promolecular density")')
       !$omp parallel do private(x,rhoat,rhocore,rdum1,rdum2,k,j,i)
       do k = 1, n(3)
          do j = 1, n(2)
             do i = 1, n(1)
                x = (/real(i-1,8)/n(1), real(j-1,8)/n(2), real(k-1,8)/n(3)/)

                if (dopro) call grda_promolecular(x,rhoat,rdum1,rdum2,0,.false.)
                if (docor) call grda_promolecular(x,rhocore,rdum1,rdum2,0,.true.)
                !$omp critical(write)
                if (dopro) f(ipdens)%f(i,j,k) = rhoat
                if (docor) f(icor)%f(i,j,k) = rhocore
                !$omp end critical(write)
             end do
          end do
       end do
       !$omp end parallel do

       if (dopro) then
          f(ipdens)%init = .true.
          write (uout,'("+ Writing promolecular density to: ",A)') trim(fileroot)//"-pdens.cube"
          call write_cube(trim(fileroot)//"-pdens.cube","Promolecular density","Written by critic2 for XDM",n,f(ipdens)%f)
       end if
       if (docor) then
          f(icor)%init = .true.
          write (uout,'("+ Writing core density to: ",A)') trim(fileroot)//"-core.cube"
          call write_cube(trim(fileroot)//"-core.cube","Core density","Written by critic2 for XDM",n,f(icor)%f)
       endif
    end if

    ! calculate the laplacian and the gradient (if we don't have the b already)
    if (ilap < 0 .and. ib < 0) then
       ilap = getfieldnum()
       write (uout,'("+ Calculating Laplacian of rho")')
       call grid_laplacian(f(irho),f(ilap))
       write (uout,'("+ Writing Laplacian to: ",A)') trim(fileroot)//"-lap.cube"
       call write_cube(trim(fileroot)//"-lap.cube","Laplacian of the electron density","Written by critic2 for XDM",n,f(ilap)%f)
       nclean = nclean + 1
       iclean(nclean) = ilap
    endif
    if (igrad < 0 .and. ib < 0) then
       igrad = getfieldnum()
       write (uout,'("+ Calculating gradient of rho")')
       call grid_gradrho(f(irho),f(igrad))
       write (uout,'("+ Writing gradient to: ",A)') trim(fileroot)//"-grad.cube"
       call write_cube(trim(fileroot)//"-grad.cube","Gradient of the electron density","Written by critic2 for XDM",n,f(igrad)%f)
       nclean = nclean + 1
       iclean(nclean) = igrad
    endif

    ! Calculate the exchange-hole dipole moment on the grid if we don't have it
    if (ib < 0) then
       if (any(f(itau)%n /= n)) &
          call ferror("xdm_driver","incongruent sizes of rho and tau",faterr)
       if (any(f(ilap)%n /= n)) &
          call ferror("xdm_driver","incongruent sizes of rho and lap",faterr)
       if (any(f(igrad)%n /= n)) &
          call ferror("xdm_driver","incongruent sizes of rho and grad",faterr)
       ! allocate a new field for it
       ib = getfieldnum()
       allocate(f(ib)%f(n(1),n(2),n(3)))
       f(ib)%n = n
       write (uout,'("+ Calculating the BR exchange-hole dipole (b)")')

       ! do it
       !$omp parallel do private(rhos,grho,lap,taus,ds,qs,rhs,xroot,xshift,&
       !$omp                     xold,expx,gx,fx,ffx,k,j,i)
       do k = 1, n(3)
          do j = 1, n(2)
             do i = 1, n(1)
                rhos = max(f(irho)%f(i,j,k),1d-14) / 2d0
                grho = f(igrad)%f(i,j,k) / 2d0
                lap = f(ilap)%f(i,j,k) / 2d0
                taus = f(itau)%f(i,j,k) / 2d0

                ds = taus - 0.25d0 * grho**2 / rhos
                qs = 1d0/6d0 * (lap - 2d0 * ds)
                rhs = 2d0/3d0 * pi**(2d0/3d0) * (rhos)**(5d0/3d0) / max(qs,1d-14)

                if (rhs > 0d0) then
                   xroot = 3d0
                   xshift = 1d0
                   do while ( (xroot * exp(-2d0*xroot/3d0)) / (xroot - 2d0) < rhs)
                      xshift = xshift * 0.1d0
                      xroot = 2d0 + xshift
                   end do
                else
                   xroot = 1d0
                   xshift = 1d0
                   do while ( (xroot * exp(-2d0*xroot/3d0)) / (xroot - 2d0) > rhs)
                      xshift = xshift * 0.1d0
                      xroot = 2d0 - xshift
                   end do
                end if

                xold = 2d0
                do while (abs(xroot - xold) > 1d-10)
                   xold = xroot
                   expx = exp(-2d0 * xroot / 3d0)
                   gx = (xroot * expx) / (xroot - 2d0)
                   fx = gx - rhs
                   ffx = gx * (1d0 / xroot - 2d0/3d0 - 1d0 / (xroot - 2d0))
                   xroot = xroot - fx / ffx
                end do

                !$omp critical(write)
                f(ib)%f(i,j,k) = xroot * (exp(-xroot) / (8d0*pi*rhos))**(1d0/3d0)
                !$omp end critical(write)
             end do
          end do
       end do
       !$omp end parallel do

       ! write the cube for b
       f(ib)%init = .true.
       write (uout,'("+ Writing BR exchange-hole dipole (b) to: ",A)') trim(fileroot)//"-b.cube"
       call write_cube(trim(fileroot)//"-b.cube","BR hole dipole moment","Written by critic2 for XDM",n,f(ib)%f)
       nclean = nclean + 1
       iclean(nclean) = ib
    else
       if (any(f(ib)%n /= n)) &
          call ferror("xdm_driver","incongruent sizes of rho and b",faterr)
    end if

    ! calculate the star of lattice vectors
    rmax = 0d0
    do i = 1, cr%nneq
       rmax = max(rmax,cutrad(cr%at(i)%z))
    end do
    call search_lattice(cr%crys2car,rmax,imax,jmax,kmax)
    allocate(lvec(3,(2*imax+1)*(2*jmax+1)*(2*kmax+1)))
    nvec = 0
    do i = -imax,imax
       do j = -jmax,jmax
          do k = -kmax,kmax
             nvec = nvec + 1
             lvec(:,nvec) = (/i,j,k/)
          end do
       end do
    end do
    write (uout,'("+ Star of lattice vector contains: ",I4," (",2(I2,","),I2,")")') nvec,imax,jmax,kmax

    ! allocate space for the atomic properties
    allocate(ml(3,cr%nneq),avol(cr%nneq))
    
    ! check more grid sizes
    if (any(f(ipdens)%n /= n)) &
       call ferror("xdm_driver","incongruent sizes of rho and pdens",faterr)
    if (icor > 0) then
       if (any(f(icor)%n /= n)) &
          call ferror("xdm_driver","incongruent sizes of rho and core",faterr)
    endif
    if (irhoae > 0) then
       if (any(f(irhoae)%n /= n)) &
          call ferror("xdm_driver","incongruent sizes of rho and rhoae",faterr)
    endif

    ! integrate atomic volumes and moments
    avol = 0d0
    ml = 0d0
    write (uout,'("+ Calculating volumes and moments")')
    !$omp parallel do private(x,ri,rhofree,raux1,raux2,rhot,&
    !$omp   wei,db,ri2,db2,i,j,k,l,ll,mll,avoll)
    do iat = 1, cr%nneq
       mll = 0d0
       avoll = 0d0
       do ll = 1, nvec
          do i = 1, n(1)
             do j = 1, n(2)
                do k = 1, n(3)
                   x = (/real(i-1,8)/n(1), real(j-1,8)/n(2), real(k-1,8)/n(3)/) + lvec(:,ll)
                   x = cr%x2c(x - cr%at(iat)%x)
                   if (any(abs(x) > cutrad(cr%at(iat)%z))) cycle
                   ri = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
                   if (ri > cutrad(cr%at(iat)%z)) cycle

                   call grid1_interp(agrid(cr%at(iat)%z),ri,rhofree,raux1,raux2)
                   rhot = f(irho)%f(i,j,k)
                   wei = rhofree * rhot / max(f(ipdens)%f(i,j,k),1d-14)
                   db = max(ri-f(ib)%f(i,j,k),0d0)
                   
                   ri2 = 1d0
                   db2 = 1d0
                   do l = 1, 3
                      ri2 = ri2 * ri
                      db2 = db2 * db
                      mll(l) = mll(l) + wei * (ri2 - db2)**2
                   end do

                   if (irhoae > 0) then
                      wei = rhofree * f(irhoae)%f(i,j,k) / max(f(ipdens)%f(i,j,k),1d-14)
                   else
                      wei = rhofree * (rhot+f(icor)%f(i,j,k)) / max(f(ipdens)%f(i,j,k),1d-14)
                   endif
                   avoll = avoll + wei * ri**3
                end do ! i
             end do ! j 
          end do ! k
       end do ! ll
       !$omp critical(write)
       ml(:,iat) = mll
       avol(iat) = avoll
       !$omp end critical(write)
    end do ! iat
    !$omp end parallel do
    avol = avol * cr%omega / ntot
    ml = ml * cr%omega / ntot
    deallocate(lvec)

    ! free volumes and polarizabilities
    write (uout,'("+ Calculating free volumes and polarizabilities")')
    allocate(ityp(103),afree(cr%nneq),alpha(cr%nneq))
    ityp = 0
    do i = 1, cr%nneq
       if (ityp(cr%at(i)%z) == 0) then
          afree(i) = free_volume(cr%at(i)%z)
          ityp(cr%at(i)%z) = i
       else
          afree(i) = afree(ityp(cr%at(i)%z))
       endif
       alpha(i) = min(avol(i) / afree(i),1d0) * alpha_free(cr%at(i)%z)
    end do
    deallocate(ityp)

    ! Write out some information
    write (uout,*)
    write (uout,'("+ Volumes and moments for non-equivalent atoms")')
    write (uout,'("# All results in atomic units.")')
    write (uout,'("# i          V                 Vfree               M1                 M2                 M3")')
    do ii = 1, cr%ncel
       i = cr%atcel(ii)%idx
       write (uout,'(99(A,X))') string(ii,length=5,justify=ioj_right), &
          string(avol(i),'e',decimal=10,length=18,justify=5), &
          string(afree(i),'e',decimal=10,length=18,justify=5), &
          (string(ml(j,i),'e',decimal=10,length=18,justify=5),j=1,3)
    end do
    write (uout,*)
    deallocate(avol, afree)

    ! Calculate the dispersion coefficents and radii
    write (uout,'("+ Dispersion coefficients")')
    write (uout,'("# All results in atomic units.")')
    write (uout,'("# i   j          C6                 C8                 C10                Rc                Rvdw")')

    allocate(rc(cr%nneq,cr%nneq))
    maxc6 = -1d0
    do ii = 1, cr%ncel
       i = cr%atcel(ii)%idx
       do jj = 1, ii
          j = cr%atcel(jj)%idx
          c6 = alpha(i)*alpha(j)*ml(1,i)*ml(1,j) / (ml(1,i)*alpha(j) + ml(1,j)*alpha(i))
          maxc6 = max(c6,maxc6)
          c8 = 3d0/2d0 * (alpha(i)*alpha(j)*(ml(1,i)*ml(2,j)+ml(2,i)*ml(1,j))) /&
             (ml(1,i)*alpha(j)+ml(1,j)*alpha(i))
          c10 = 2 * alpha(i)*alpha(j) * (ml(1,i)*ml(3,j) + ml(3,i)*ml(1,j)) /&
             (ml(1,i)*alpha(j) + ml(1,j)*alpha(i)) + 21d0/5d0 * alpha(i)*alpha(j)*&
             ml(2,i)*ml(2,j) / (alpha(j)*ml(1,i)+alpha(i)*ml(1,j))
          
          rc(i,j) = (sqrt(c8/c6) + sqrt(c10/c8) + (c10/c6)**(0.25d0)) / 3
          rc(j,i) = rc(i,j)
          rvdw = a1 * rc(i,j) + a2

          write (uout,'(99(A,X))') string(ii,length=5,justify=ioj_right),&
             string(jj,length=5,justify=ioj_right),&
             string(c6,'e',decimal=10,length=18,justify=5),&
             string(c8,'e',decimal=10,length=18,justify=5),&
             string(c10,'e',decimal=10,length=18,justify=5),&
             string(rc(i,j),'e',decimal=10,length=18,justify=5),&
             string(rvdw,'e',decimal=10,length=18,justify=5)
       end do
    end do
    if (.true.) then 
       write (uout,*)
       write (uout,'("# i   j   k             C9")')
       do ii = 1, cr%ncel
          i = cr%atcel(ii)%idx
          do j = 1, i
             do k = 1, j
                c9 = ml(1,i)*ml(1,j)*ml(1,k)*(ml(1,i)/alpha(i)+ml(1,j)/alpha(j)+ml(1,k)/alpha(k)) / &
                   (ml(1,i)/alpha(i) + ml(1,j)/alpha(j)) / (ml(1,i)/alpha(i) + ml(1,k)/alpha(k)) / &
                   (ml(1,j)/alpha(j) + ml(1,k)/alpha(k))
                write (uout,'(99(A,X))') string(ii,length=5,justify=ioj_right), &
                   string(j,length=5,justify=ioj_right),&
                   string(k,length=5,justify=ioj_right),&
                   string(c9,'e',decimal=10,length=18,justify=5)
             end do
          end do
       end do
    end if
    write (uout,*)

    ! calculate energy contributions
    if (onlyc) goto 999
    write (uout,'("+ van der Waals energy")')
    etotal = 0d0
    for = 0d0
    sigma = 0d0
    ehadd = 0d0

    ! set the atomic environment for the sum
    rmax = (maxc6/ecut)**(1d0/6d0)
    rmax2 = rmax*rmax
    do ii = 1, cr%ncel
       i = cr%atcel(ii)%idx
       eat = 0d0
       sat = 0d0
       do jj = 1, cr%nenv
          j = cr%atenv(jj)%idx
          x = cr%atenv(jj)%r - cr%atcel(ii)%r
          ri2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          if (ri2 < 1d-15 .or. ri2>rmax2) cycle
          ri = sqrt(ri2)
          ri3 = ri2 * ri

          fcom = 1.5d0 * alpha(i)*alpha(j) / (ml(1,i)*alpha(j)+ml(1,j)*alpha(i))
          do nn = 6, upto, 2
             ! order R^nn
             i3 = nn / 2 - 1

             ! calculate the Cn using the moments
             cn0 = 0d0
             do i3i = 1, (i3/2)
                l1 = i3i
                l2 = i3-i3i
                fll = fact(2*(l1+l2)) / real((2*l1+1)*(2*l2+1)*&
                   fact(2*l1)*fact(2*l2),8)
                cn0 = cn0 + fll * ml(l1,i) * ml(l2,j)
                if (l1 /= l2) &
                   cn0 = cn0 + fll * ml(l2,i) * ml(l1,j)
             end do
             cn0 = cn0 * fcom

             ! energy contribution
             rvdw = a1 * rc(i,j) + a2
             rvdwx = rvdw**nn
             rix = ri**nn
             exx = cn0 / (rvdwx + rix)
             ehadd(nn) = ehadd(nn) + exx
             eat = eat + exx

             ! force and stress contribution
             fxx = nn * cn0 * rix / ri2 / (rvdwx + rix)**2
             for(:,ii) = for(:,ii) + fxx * x
             do m1 = 1, 3
                sat(m1,m1) = sat(m1,m1) + fxx * x(m1) * x(m1)
                do m2 = 1, m1-1
                   sat(m1,m2) = sat(m1,m2) + fxx * x(m1) * x(m2)
                   sat(m2,m1) = sat(m1,m2)
                end do
             end do
          end do
       end do ! jj
       etotal = etotal + eat 
       sigma = sigma + sat 
    end do ! ii
    etotal= - 0.5d0 * etotal
    sigma = - 0.5d0 * sigma / cr%omega
    ehadd = - 0.5d0 * ehadd
     
    write (uout,'("  Evdw = ",A," Hartree, ",A," Ry")') &
       string(etotal,'e',decimal=10), string(etotal*2,'e',decimal=10)
    do nn = 6, 10, 2
       if (nn > 9) then
          write (uout,'("  Evdw",A," = ",A," Hartree, ",A," Ry")') &
             string(nn), string(ehadd(nn),'e',decimal=10), string(ehadd(nn)*2,'e',decimal=10)
       else
          write (uout,'("  Evdw",A," = ",A," Hartree, ",A," Ry")') &
             string(nn), string(ehadd(nn),'e',decimal=10), string(ehadd(nn)*2,'e',decimal=10)
       end if
    end do

    do i = 1, cr%ncel
       write (uout,'("  Fvdw (",A,") = ",3(A,X)," Hartree/bohr")') &
          string(i), (string(for(j,i),'e',decimal=10),j=1,3)
       write (uout,'("            = ",3(A,X)," Ry/bohr")') &
          (string(for(j,i)*2,'e',decimal=10),j=1,3)
    end do
    write (uout,'("  sigma_vdw (Hy/bohr**3) = ",3(A,X)," ")') (string(sigma(1,j),'e',decimal=10),j=1,3)
    write (uout,'("                           ",3(A,X)," ")') (string(sigma(2,j),'e',decimal=10),j=1,3)
    write (uout,'("                           ",3(A,X)," ")') (string(sigma(3,j),'e',decimal=10),j=1,3)
    write (uout,'("  sigma_vdw (Ry/bohr**3) = ",3(A,X)," ")') (string(sigma(1,j)*2,'e',decimal=10),j=1,3)
    write (uout,'("                           ",3(A,X)," ")') (string(sigma(2,j)*2,'e',decimal=10),j=1,3)
    write (uout,'("                           ",3(A,X)," ")') (string(sigma(3,j)*2,'e',decimal=10),j=1,3)
    write (uout,'("  sigma_vdw (GPa) = ",3(A,X)," ")') (string(sigma(1,j)*autogpa,'e',decimal=10),j=1,3)
    write (uout,'("                    ",3(A,X)," ")') (string(sigma(2,j)*autogpa,'e',decimal=10),j=1,3)
    write (uout,'("                    ",3(A,X)," ")') (string(sigma(3,j)*autogpa,'e',decimal=10),j=1,3)
    write (uout,*)
    
999 continue

    ! cleanup
    deallocate(rc, alpha, ml)
    do i = 1, nclean
       call fields_unload(iclean(i))
    end do

  end subroutine xdm_grid

  !> Calculate XDM from the information in a QE output
  subroutine xdm_qe()
    use tools_io
    use struct_basic
    
    integer :: i, j, jj, nn, i3
    integer :: lu, idx, idx1, idx2
    character(len=:), allocatable :: line, str
    logical :: ok
    real*8 :: a1, a2, maxc6, rmax, rmax2, x(3), ri, ri2, ri3
    real*8 :: cn0, rvdwx, rix, exx, eat, etotal
    real*8 :: c6(cr%ncel,cr%ncel), c8(cr%ncel,cr%ncel), c10(cr%ncel,cr%ncel)
    real*8 :: rc(cr%ncel,cr%ncel), rvdw(cr%ncel,cr%ncel)
    
    real*8, parameter :: ecut = 1d-11

    write (uout,'("* Sum the XDM dispersion energy using a QE output")')
    write (uout,'("+ Reading coefficients from the file: ",A)') string(cr%file)

    lu = fopen_read(string(cr%file))
    call cr%build_env(150d0)

    main: do while (getline(lu,line))
       ! read the parameters
       if (trim(line) == "* XDM dispersion") then
          ok = getline(lu,line)
          idx = index(line,"=")
          str = trim(line(idx+1:))
          read(str,*) a1
          ok = getline(lu,line)
          ok = getline(lu,line)
          idx = index(line,"=")
          str = trim(line(idx+1:))
          read(str,*) a2
       end if
       ! read the dispersion coefficients and the radii
       if (trim(line) == "+ Dispersion coefficients") then
          do i = 1, cr%ncel
             do j = 1, i
                ok = getline(lu,line)
                read(line,*) idx1, idx2, c6(i,j), c8(i,j), c10(i,j), rc(i,j), rvdw(i,j)
                if (idx1 /= i .or. idx2 /= j) then
                   write (*,*) idx1, i
                   write (*,*) idx2, j
                   call ferror("trick","read indices do not match",faterr)
                end if
                c6(j,i) = c6(i,j)
                c8(j,i) = c8(i,j)
                c10(j,i) = c10(i,j)
                rc(j,i) = rc(i,j)
                rvdw(j,i) = rvdw(i,j)
             end do
          end do
       end if
    end do main
    call fclose(lu)

    ! calculate the energy
    maxc6 = maxval(c6)
    ! set the atomic environment for the sum
    rmax = (maxc6/ecut)**(1d0/6d0)
    rmax2 = rmax*rmax
    etotal = 0d0
    do i = 1, cr%ncel
       eat = 0d0
       do jj = 1, cr%nenv
          j = cr%atenv(jj)%cidx
          x = cr%atenv(jj)%r - cr%atcel(i)%r
          ri2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          if (ri2 < 1d-15 .or. ri2>rmax2) cycle
          ri = sqrt(ri2)
          ri3 = ri2 * ri
    
          do nn = 6, 10, 2
             ! order R^nn
             i3 = nn / 2 - 1
    
             ! check the ii and the jj, also how the env is generated
             cn0 = 0d0
             if (nn == 6) then
                cn0 = c6(i,j)
             elseif (nn == 8) then
                cn0 = c8(i,j)
             elseif (nn == 10) then
                cn0 = c10(i,j)
             end if
    
             ! energy contribution
             rvdwx = rvdw(i,j)**nn
             rix = ri**nn
             exx = cn0 / (rvdwx + rix)
             eat = eat + exx
          end do
       end do ! jj
       etotal = etotal + eat 
    end do ! ii
    etotal= -0.5d0 * etotal
    
    write (uout,'("  Evdw = ",A," Hartree, ",A," Ry")') &
       string(etotal,'e',decimal=10), string(etotal*2,'e',decimal=10)
    
    write (uout,*)

  end subroutine xdm_qe

  !> Calculate XDM in molecules
  subroutine xdm_mol(a1o,a2o,chf)
    use meshmod
    use fields
    use global
    use struct_basic
    use grd_atomic
    use grid1_tools
    use tools_math
    use tools_io
    use types
    use param
    
    real*8, intent(in) :: a1o, a2o
    real*8, intent(in) :: chf

    type(molmesh) :: m
    integer :: nelec
    integer :: id(7), prop(7), i, j, lu, luh, iz
    real*8 :: rho, rhop, rhopp, x(3), r, a1, a2, nn, rb
    real*8 :: mm(3,cr%ncel), v(cr%ncel)
    
    ! only for molecules
    if (.not.cr%ismolecule) &
       call ferror("xdm_mol","molecular XDM only for molecules",faterr)

    ! only for wfn or dftb
    if (f(refden)%type /= type_wfn .and. f(refden)%type /= type_dftb) &
       call ferror("xdm_mol","molecular XDM only for wfn and dftb fields",faterr)

    ! only for closed shells
    if (f(refden)%type == type_wfn) then
       if (f(refden)%wfntyp /= 0) &
          call ferror("xdm_mol","open shell wavefunctions not supported",faterr)
    end if

    ! write some info to the output
    write (uout,'("a1             ",A)') string(a1o,'f',12,6)
    write (uout,'("a2(ang)        ",A)') string(a2o,'f',12,6)
    a1 = a1o
    a2 = a2o / bohrtoa

    if (chf < 0d0) then
       select case (nint(chf))
       case (chf_blyp)
          write(uout,'("a_hf           blyp")')
       case (chf_b3lyp)
          write(uout,'("a_hf           b3lyp")')
       case (chf_bhahlyp)
          write(uout,'("a_hf           bhandhlyp")')
       case (chf_camb3lyp)
          write(uout,'("a_hf           camb3lyp")')
       case (chf_pbe)
          write(uout,'("a_hf           pbe")')
       case (chf_pbe0)
          write(uout,'("a_hf           pbe0")')
       case (chf_lcwpbe)
          write(uout,'("a_hf           lcwpbe")')
       case (chf_pw86)
          write(uout,'("a_hf           pw86pbe")')
       case (chf_b971)
          write(uout,'("a_hf           b971")')
       end select
    elseif (abs(chf-1d0) < 1d-9) then
       write(uout,'("a_hf           hf")')
    else
       write(uout,'("a_hf           ",A)') string(chf,'f',12,6)
    endif

    ! prepare the mesh
    m = genmesh(cr)
    write (uout,'("mesh size      ",A)') string(m%n)

    ! properties to calculate 
    do i = 1, 5
       id(i) = i
    end do
    prop(1) = im_rho
    prop(2) = im_gradrho
    prop(3) = im_gkin
    prop(4) = im_null ! for promolecular / hirshfeld weights
    prop(5) = im_null ! for the atomic density contribution
    prop(6) = im_b    

    ! fill the mesh with those properties
    call fillmesh(m,f(refden),id,prop)

    ! fill the promolecular and the atomic densities
    m%f(:,4:5) = 0d0
    lu = fopen_scratch()
    nelec = 0
    do i = 1, cr%ncel
       iz = cr%at(cr%atcel(i)%idx)%z
       if (iz < 1) cycle
       nelec = nelec + iz

       do j = 1, m%n
          x = m%x(:,j) - cr%atcel(i)%r
          r = norm(x)
          call grid1_interp(agrid(iz),r,rho,rhop,rhopp)
          m%f(j,4) = m%f(j,4) + rho
          m%f(j,5) = rho
       enddo
       write (lu) (m%f(j,5),j=1,m%n)
    enddo

    ! create the temporary file with the weights
    luh = fopen_scratch()
    rewind(lu)
    do i = 1, cr%ncel
       iz = cr%at(cr%atcel(i)%idx)%z
       if (iz < 1) cycle
       read (lu) (m%f(j,5),j=1,m%n)
       write (luh) (max(m%f(j,5),1d-40)/max(m%f(j,4),1d-40),j=1,m%n)
    enddo
    call fclose(lu)

    ! integrate the number of electrons
    nn = sum(m%f(:,1) * m%w)
    write (uout,'("nelec          ",A)') string(nelec)
    write (uout,'("nelec (promol) ",A)') string(sum(m%f(:,4) * m%w),'f',12,6)
    write (uout,'("nelec, total   ",A)') string(nn,'f',12,6)
    if (abs(nn - nelec) > 0.1d0) &
       call ferror("xdm_mol","inconsistent nelec. I hope you know what you are doing",warning)

    ! calculate moments and volumes
    mm = 0d0
    v = 0d0
    rewind(luh)
    do i = 1, cr%ncel
       iz = cr%at(cr%atcel(i)%idx)%z
       if (iz < 1) cycle
       read (luh) (m%f(j,4),j=1,m%n)

       ! calculate hole dipole and moments
       do j = 1, m%n
          r = norm(m%x(:,j)-cr%atcel(i)%r)
          rb = max(0.d0,r-m%f(j,6))

          mm(1,i) = mm(1,i) + m%w(j) * m%f(j,4) * m%f(j,1) * (r-rb)**2
          mm(2,i) = mm(2,i) + m%w(j) * m%f(j,4) * m%f(j,1) * (r**2-rb**2)**2
          mm(3,i) = mm(3,i) + m%w(j) * m%f(j,4) * m%f(j,1) * (r**3-rb**3)**2
          v(i) = v(i) + m%w(j) * m%f(j,4) * m%f(j,1) * r**3
       enddo
    enddo
    call fclose(luh)

    ! calculate and output energy and derivatives
    call edisp_mol(a1,a2,chf,v,mm)

  end subroutine xdm_mol

  !> Write the header of a cube file
  subroutine write_cube(file,line1,line2,n,c)
    use global
    use tools_io
    use struct_basic

    character*(*), intent(in) :: file, line1, line2
    integer, intent(in) :: n(3)
    real*8, intent(in) :: c(0:n(1)-1,0:n(2)-1,0:n(3)-1)

    integer :: i, j, k, lu
    real*8 :: xx(3)

    ! open the file
    lu = fopen_write(file)

    ! cube title and grid dimensions
    write(lu,'(A)') trim(line1)
    write(lu,'(A)') trim(line2)
    write(lu,'(I5,3(F12.6))') cr%ncel, 0d0, 0d0, 0d0
    xx = cr%x2c((/1d0,0d0,0d0/))
    write(lu,'(I5,3(F12.6))') n(1), xx / real(n(1),8)
    xx = cr%x2c((/0d0,1d0,0d0/))
    write(lu,'(I5,3(F12.6))') n(2), xx / real(n(2),8)
    xx = cr%x2c((/0d0,0d0,1d0/))
    write(lu,'(I5,3(F12.6))') n(3), xx / real(n(3),8)

    ! write the atoms
    do i = 1, cr%ncel
       write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') cr%at(cr%atcel(i)%idx)%z, 0d0, cr%atcel(i)%r
    end do

    ! write the field
    do i = 0, n(1)-1
       do j = 0, n(2)-1
             if (precisecube) then
                write (lu,'(6(1x,e22.14))') (c(i,j,k),k=0,n(3)-1)
             else
                write (lu,'(1p,6(1x,e12.5))') (c(i,j,k),k=0,n(3)-1)
             end if
       enddo
    enddo

    ! close the file
    call fclose(lu)

  end subroutine write_cube

  function free_volume(iz) result(afree)
    use tools_io
    use grd_atomic
    use grid1_tools
    use param

    integer, intent(in) :: iz
    real*8 :: afree

    integer, parameter :: ngau = 251
    integer :: k
    real*8 :: rmid, h, q, r, rwei, rhofree, raux1, raux2

    afree = 0d0
    rmid = 1d0 / real(iz,8)**(1d0/3d0)
    h = 1d0 / real(ngau+1,8)
    do k = 1, ngau
       q = h * k
       r = rmid * q / (1d0-q)
       rwei = 4d0*pi*h * r**2 * rmid/(1d0-q)**2
       call grid1_interp(agrid(iz),r,rhofree,raux1,raux2)
       afree = afree + rhofree * rwei * r**3
    end do
  end function free_volume

  function frevol(z,chf)
    use tools_io

    integer, intent(in) :: z
    real*8, intent(in) :: chf
    real*8 :: frevol

    ! DKH LSDA/UGBS Free Atomic Volumes
    real*8, parameter :: frevol0(0:103) = (/&
       0.d0,&
       9.194D0,   4.481D0,  91.957D0,  61.357D0,  49.813D0,  36.728D0,&
       27.633D0,  23.517D0,  19.322D0,  15.950D0, 109.359D0, 103.064D0,&
       120.419D0, 104.229D0,  86.782D0,  77.133D0,  66.372D0,  57.336D0,&
       203.093D0, 212.202D0, 183.101D0, 162.278D0, 143.250D0, 108.209D0,&
       123.098D0, 105.735D0,  92.944D0,  83.794D0,  75.750D0,  81.177D0,&
       118.371D0, 116.334D0, 107.474D0, 103.221D0,  95.111D0,  87.605D0,&
       248.772D0, 273.748D0, 249.211D0, 223.801D0, 175.809D0, 156.831D0,&
       160.042D0, 136.654D0, 127.754D0,  97.024D0, 112.778D0, 121.627D0,&
       167.906D0, 172.030D0, 165.500D0, 163.038D0, 153.972D0, 146.069D0,&
       341.992D0, 385.767D0, 343.377D0, 350.338D0, 334.905D0, 322.164D0,&
       310.337D0, 299.537D0, 289.567D0, 216.147D0, 268.910D0, 259.838D0,&
       251.293D0, 243.174D0, 235.453D0, 228.284D0, 229.617D0, 209.971D0,&
       197.541D0, 183.236D0, 174.685D0, 164.139D0, 150.441D0, 135.765D0,&
       125.297D0, 131.258D0, 185.769D0, 195.671D0, 193.036D0, 189.142D0,&
       185.919D0, 181.089D0, 357.787D0, 407.283D0, 383.053D0, 362.099D0,&
       346.565D0, 332.462D0, 319.591D0, 308.095D0, 297.358D0, 300.572D0,&
       275.792D0, 266.317D0, 257.429D0, 209.687D0, 203.250D0, 230.248D0,&
       236.878D0/)

    real*8, parameter :: frevol_blyp(0:36) = (/0.0d0,&
       8.6751280810827840d0,  4.3645522950717863d0,  89.719664495180297d0,  61.278566735200307d0,&
       49.519428604126382d0,  36.855287686102294d0,  27.995151669578650d0,  23.764998306645893d0,&
       19.585085705364346d0,  16.198873185715303d0,  110.72206109235110d0,  104.99418885629647d0,&
       124.86364864987824d0,  105.67572783021383d0,  87.717325899499158d0,  77.848560257666222d0,&
       67.289916054772661d0,  58.134006480572936d0,                0.00d0,                0.00d0,&
       189.14940377667708d0,  171.31844286164150d0,  149.09546424400813d0,  114.14072711487664d0,&
       128.27850058556785d0,  108.89501751641639d0,  100.22570628848071d0,  88.525764030777225d0,&
       79.682125997948589d0,  86.320056678306727d0,  125.54483633211143d0,  121.00196086145128d0,&
       110.76811978260324d0,  106.26389443012715d0,  98.394240990660435d0,  90.352795412830417d0/)

    real*8, parameter :: frevol_b3lyp(0:36) = (/0.0d0,&
       8.3489695345990036d0,   4.2304165146306838d0,   88.786443805705602d0,   60.659572325858221d0,&
       48.180704821245364d0,   35.739786376347816d0,   27.113474668623077d0,   22.954305454378311d0,&
       18.919035719008878d0,   15.664743495672592d0,   110.10781957063872d0,   104.58688460317792d0,&
       121.82866689851529d0,   103.64265930616017d0,   86.233229483159292d0,   76.534831202674539d0,&
       66.225085507692555d0,   57.258468294377920d0,                 0.00d0,                 0.00d0,&
       191.24458275864666d0,   171.33960226014415d0,   154.31690117090523d0,   111.84266101395626d0,&
       128.53304731251799d0,   116.41853334783917d0,   108.02725265724374d0,   84.734127226732753d0,&
       78.821628181721010d0,   86.690015512660722d0,   122.56838049698425d0,   118.57117842661357d0,&
       108.88205342457658d0,   104.55842043337798d0,   96.968864870029904d0,   89.143550757884739d0/)

    real*8, parameter :: frevol_bhahlyp(0:36) = (/0.0d0,&
       8.0162383356457330d0,   4.0855044443453972d0,   89.521965961222222d0,   60.766620056468859d0,&
       47.163104244372654d0,   34.732519196898075d0,   26.261248810982419d0,   22.163686460916928d0,&
       18.257837007647932d0,   15.123544024501175d0,   112.24236986192338d0,   105.88927909724653d0,&
       120.53914101629064d0,   102.42697418515390d0,   85.163324240312747d0,   75.599161377353852d0,&
       65.445419322876759d0,   56.592207744291457d0,                 0.00d0,                 0.00d0,&
       196.29478033696057d0,   175.14878735114101d0,   157.82044436715444d0,   112.66706636630363d0,&
       131.50432127575891d0,   120.17186617063507d0,   110.69809115434559d0,   102.41457915920718d0,&
       79.806952905523545d0,   88.942968328671427d0,   121.67064092786063d0,   117.25525753962938d0,&
       107.64004292512693d0,   103.49115496709287d0,   96.074020468035101d0,   88.360967365230920d0/)

    real*8, parameter :: frevol_camb3lyp(0:36) = (/0.0d0,&
       8.4912418798953802d0,   4.3046720809866850d0,   88.495476866973974d0,   60.708687086264462d0,&
       48.017293545760708d0,   35.742021952558758d0,   27.186861831258817d0,   23.059619196158092d0,&
       19.022193061849432d0,   15.756088275099161d0,   109.34523733398444d0,   104.23708524052515d0,&
       119.90584483487682d0,   102.80245055238821d0,   85.960831224028368d0,   76.454433439613283d0,&
       66.279667949988067d0,   57.375114690853394d0,                 0.00d0,                 0.00d0,&
       192.25186364675062d0,   171.60008159670065d0,   154.54149488373386d0,   111.51614466860299d0,&
       105.82575187101907d0,   97.244815362081795d0,   90.518519096216096d0,   84.233565981799387d0,&
       78.777041862773274d0,   86.830269580599790d0,   119.93150277142236d0,   117.12848860140892d0,&
       108.22197759832967d0,   104.14881532532733d0,   96.808752555226278d0,   89.137842527944699d0/)

    real*8, parameter :: frevol_pbe(0:36) = (/0.0d0,&
       8.7017290470668396d0,   4.3452296013273557d0,   90.345016391875745d0,   60.583880555189161d0,&
       49.405610809309636d0,   36.699238724745918d0,   27.804307905805462d0,   23.504820577931341d0,&
       19.366722293791256d0,   16.015015785319534d0,   114.11568983081936d0,   105.34255389919765d0,&
       122.27570960010875d0,   104.11673088180522d0,   86.659855994802655d0,   76.696371875869517d0,&
       66.330422698625398d0,   57.338280708600728d0,                 0.00d0,                 0.00d0,&
       187.56578096056526d0,   171.16842233126320d0,   148.57152896051596d0,   113.40726679106993d0,&
       128.41370162442536d0,   107.47276564091520d0,   98.987585702398064d0,   87.019526141682050d0,&
       79.950858621796073d0,   86.179295605944361d0,   123.32536942981970d0,   119.49220381812802d0,&
       109.53573383239170d0,   104.78641436128312d0,   97.036678428323626d0,   89.131992369076343d0/)

    real*8, parameter :: frevol_pbe0(0:36) = (/0.0d0,&
       8.2794385587230224d0,   4.1765568461439084d0,   89.990651247904850d0,   60.150518670639258d0,&
       47.925007649610208d0,   35.403450375407488d0,   26.774856262986901d0,   22.577665436425793d0,&
       18.604506038051770d0,   15.403636840460491d0,   114.72060124074002d0,   105.56153583642701d0,&
       119.94555610703479d0,   102.21839478488732d0,   85.124898059747593d0,   75.344227406670839d0,&
       65.219744182377752d0,   56.417325659383977d0,                 0.00d0,                 0.00d0,&
       191.49634784581932d0,   172.37802114577727d0,   155.38093772078227d0,   111.48097347813784d0,&
       129.94413257650865d0,   118.24195915775856d0,   108.60919379569880d0,   84.698903894685841d0,&
       79.451081251238904d0,   87.077874282230383d0,   120.92043418719176d0,   117.15417864570686d0,&
       107.60403177731271d0,   103.07063942096924d0,   95.595468315759462d0,   87.905274745875047d0/)

    real*8, parameter :: frevol_lcwpbe(0:36) = (/0.0d0,&
       8.2370522934694321d0,   4.3223022069392556d0,   88.889190621676747d0,   59.167955275706255d0,&
       46.644645536860530d0,   35.000149688018325d0,   26.874237675680991d0,   22.830136179756057d0,&
       18.961662156609091d0,   15.787003768893198d0,   115.52305481049332d0,   105.20996979820804d0,&
       115.58130215151486d0,   99.382061772188109d0,   83.591029347688959d0,   74.254876662748089d0,&
       64.635098656590586d0,   56.191387396031978d0,                 0.00d0,                 0.00d0,&
       191.95897600535326d0,   172.74995562085155d0,   155.27567955463721d0,   110.25556775344671d0,&
       129.26382250173134d0,   117.09487172302605d0,   107.19200454429361d0,   83.698218802709803d0,&
       78.438819467992630d0,   85.426026058281309d0,   114.73764666891768d0,   113.36022260917842d0,&
       105.50701995208520d0,   101.36644049619710d0,   94.506386358808882d0,   87.304323578968550d0/)

    real*8, parameter :: frevol_pw86(0:36) = (/0.0d0,&
       8.5848597505149957d0,   4.3045044427345758d0,   89.105017427126995d0,   60.116125959883448d0,&
       48.948725444378958d0,   36.639492704666679d0,   27.890608359059900d0,   23.482178111604007d0,&
       19.392669723193279d0,   16.068196774925887d0,   110.77812603458771d0,   103.32122266742829d0,&
       121.98526328640968d0,   104.32311632196750d0,   86.866106676280893d0,   76.813580132172405d0,&
       66.479441471513695d0,   57.487390010861731d0,                 0.00d0,                 0.00d0,&
       185.73638891654375d0,   168.40972266496894d0,   146.57709921981228d0,   114.30230293458976d0,&
       127.03123778878130d0,   106.92921677417986d0,   79.779742298798681d0,   87.455617557953332d0,&
       79.471676692123495d0,   85.082884646536542d0,   123.88868121486288d0,   120.63195902241362d0,&
       110.68925500617975d0,   105.58567728669152d0,   97.798158973601787d0,   89.844087877827207d0/)

    real*8, parameter :: frevol_b971(0:36) = (/0.0d0,&
       8.2753044447676203d0,   4.1538172892410463d0,   92.548449604186061d0,   60.079247546006165d0,&
       47.862796246981773d0,   35.305970222823603d0,   26.699111937155749d0,   22.644191719940501d0,&
       18.675265336738804d0,   15.464799224948079d0,   121.31651007378332d0,   105.68066738843501d0,&
       120.26430745377735d0,   101.96967180988391d0,   85.008007953828127d0,   75.622028007253363d0,&
       65.525564639136547d0,   56.692633782240122d0,                 0.00d0,                 0.00d0,&
       193.61923913260267d0,   171.80907721991764d0,   154.23278648930369d0,   114.13224357718498d0,&
       107.38843940948567d0,   98.695645988131020d0,   92.056435539203278d0,   86.394941236711290d0,&
       80.059364904746786d0,   87.457830726555457d0,   121.47718849519879d0,   116.93338080131248d0,&
       107.07435287066555d0,   103.20705595068554d0,   95.879065901986593d0,   88.146215884084469d0/)

    real*8 :: rchf

    ! gaussian basis set, %HF=(0,25,50,75,100) for periods 1 and 2.
    real*8 frevol1(5,0:10)
    save frevol1
    data frevol1/&
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
       8.7017290470668396d0, 8.2794385587230224d0, 7.9275669093485428d0, 7.6272753143131924d0, 7.3655545253236543d0,&
       4.3452296013273557d0, 4.1765568461439084d0, 4.0288815322803657d0, 3.8977496795733644d0, 3.7798775917054876d0,&
       90.345016391875745d0, 89.990651247904850d0, 89.843447659370099d0, 89.862545631289919d0, 90.018126381136568d0,&
       60.583880555189161d0, 60.150518670639258d0, 59.789689626639294d0, 59.490664095467629d0, 59.244839599587813d0,&
       49.405610809309636d0, 47.925007649610208d0, 46.767078359393196d0, 45.827336436528093d0, 45.040918655916336d0,&
       36.699238724745918d0, 35.403450375407488d0, 34.355914497332748d0, 33.486552145008638d0, 32.747999094847692d0,&
       27.804307905805462d0, 26.774856262986901d0, 25.931460984272665d0, 25.224788496044660d0, 24.620417659921344d0,&
       23.504820577931341d0, 22.577665436425793d0, 21.811051051552955d0, 21.165993022839672d0, 20.614016545011395d0,&
       19.366722293791256d0, 18.604506038051770d0, 17.967854713848929d0, 17.427825507512161d0, 16.962774503019283d0,&
       16.015015785319534d0, 15.403636840460491d0, 14.886731288556328d0, 14.443878154969715d0, 14.059464430879629d0/

    ! pure GGA
    if (abs(chf) < 1d-10) then
       if (z > 10) then
          frevol = frevol0(z)
       else
          frevol = frevol1(1,z)
       endif
       return
    endif

    ! special cases
    if (chf < 0d0 .and. (z<19.or.z>20.and.z<37)) then ! up to Kr except K and Ca
       select case(nint(chf))
       case(chf_blyp) 
          frevol = frevol_blyp(z)
       case(chf_b3lyp) 
          frevol = frevol_b3lyp(z)
       case(chf_bhahlyp) 
          frevol = frevol_bhahlyp(z)
       case(chf_camb3lyp) 
          frevol = frevol_camb3lyp(z)
       case(chf_pbe) 
          frevol = frevol_pbe(z)
       case(chf_pbe0) 
          frevol = frevol_pbe0(z)
       case(chf_lcwpbe) 
          frevol = frevol_lcwpbe(z)
       case(chf_pw86) 
          frevol = frevol_pw86(z)
       case(chf_b971) 
          frevol = frevol_b971(z)
       case default
          call ferror("frevol","unknown functional",faterr)
       end select
    else
       ! general hybrid
       if (chf < 0d0) then
          select case(nint(chf))
          case(chf_blyp) 
             rchf = 0d0
          case(chf_b3lyp) 
             rchf = 0.2d0
          case(chf_bhahlyp) 
             rchf = 0.5d0
          case(chf_camb3lyp) 
             rchf = 0.2d0
          case(chf_pbe) 
             rchf = 0.0d0
          case(chf_pbe0) 
             rchf = 0.25d0
          case(chf_lcwpbe) 
             rchf = 0.25d0
          case(chf_pw86) 
             rchf = 0.0d0
          case(chf_b971) 
             rchf = 0.21d0
          case default
             call ferror("frevol","unknown functional",faterr)
          end select
       else
          rchf = chf
       endif

       if (z > 10) then
          frevol = frevol0(z)
       else
          if (rchf < 0.25d0) then
             frevol = frevol1(1,z) + (frevol1(2,z)-frevol1(1,z)) / 0.25d0 * (rchf-0d0)
          elseif (rchf < 0.50d0) then
             frevol = frevol1(2,z) + (frevol1(3,z)-frevol1(2,z)) / 0.25d0 * (rchf-0.25d0)
          elseif (rchf < 0.75d0) then
             frevol = frevol1(3,z) + (frevol1(4,z)-frevol1(3,z)) / 0.25d0 * (rchf-0.50d0)
          else
             frevol = frevol1(4,z) + (frevol1(5,z)-frevol1(4,z)) / 0.25d0 * (rchf-0.75d0)
          endif
       endif
    endif

  endfunction frevol

  !> Calculate the dispersion energy in a molecule
  subroutine edisp_mol(a1,a2,chf,v,mm)
    use struct_basic
    use tools_io
    use param

    real*8, intent(in) :: a1, a2, chf
    real*8, intent(in) :: v(cr%ncel), mm(3,cr%ncel)

    integer :: i, j, k, k1, k2, iz
    real*8 :: d, atpol(cr%ncel), fac, rvdw, c6, c8, c10, rc
    real*8 :: c6com, c8com, c10com, xij(3), ifac
    real*8 :: e, f(3,cr%ncel), q(3,cr%ncel,3,cr%ncel), qfac

    ! write moments and volumes
    write (uout,'("moments and volumes ")')
    write (uout,'("# i At        <M1^2>             <M2^2>              <M3^2>           Volume              Vfree")')
    do i = 1, cr%ncel
       iz = cr%at(cr%atcel(i)%idx)%z
       write (uout,'(99(A,X),5(E18.10,X))') string(i,3), string(cr%at(cr%atcel(i)%idx)%name,2),&
          (string(mm(j,i),'e',18,10),j=1,3), string(v(i),'e',18,10), string(frevol(iz,chf),'e',18,10)
    enddo
    write (uout,'("#")')
    
    do i = 1, cr%ncel
       iz = cr%at(cr%atcel(i)%idx)%z
       if (iz < 1) cycle
       atpol(i) = v(i) * alpha_free(iz) / frevol(iz,chf)
    enddo

    write (uout,'("coefficients and distances (a.u.)")')
    write (uout,'("# i  j       dij            C6               C8               C10              Rc           Rvdw")') 
    e = 0d0
    f = 0d0
    q = 0d0
    do i = 1, cr%ncel
       iz = cr%at(cr%atcel(i)%idx)%z
       if (iz < 1) cycle
       do j = i, cr%ncel
          iz = cr%at(cr%atcel(j)%idx)%z
          if (iz < 1) cycle
          xij = cr%atcel(j)%r-cr%atcel(i)%r
          d = sqrt(dot_product(xij,xij))
          fac = atpol(i)*atpol(j)/(mm(1,i)*atpol(j)+mm(1,j)*atpol(i))
          c6 = fac*mm(1,i)*mm(1,j)
          c8 = 1.5d0*fac*(mm(1,i)*mm(2,j)+mm(2,i)*mm(1,j))
          c10 = 2.d0*fac*(mm(1,i)*mm(3,j)+mm(3,i)*mm(1,j))&
             +4.2d0*fac*mm(2,i)*mm(2,j)
          rc = (sqrt(c8/c6) + sqrt(sqrt(c10/c6)) +&
             sqrt(c10/c8)) / 3.D0
          rvdw = a1 * rc + a2
          if (d > 1d-5) then
             e = e - c6 / (rvdw**6 + d**6) - c8 / (rvdw**8+d**8) - &
                c10 / (rvdw**10 + d**10)
             c6com = 6.d0*c6*d**4/(rvdw**6+d**6)**2
             c8com = 8.d0*c8*d**6/(rvdw**8+d**8)**2
             c10com = 10.d0*c10*d**8/(rvdw**10+d**10)**2
             f(:,i) = f(:,i) + (c6com+c8com+c10com) * xij
             f(:,j) = f(:,j) - (c6com+c8com+c10com) * xij
             do k1 = 1, 3
                do k2 = 1, 3
                   if (k1 == k2) then
                      ifac = 1d0
                   else
                      ifac = 0d0
                   endif
                   qfac = &
                      c6com  * (-ifac - 4*xij(k1)*xij(k2)/d**2 + 12*xij(k1)*xij(k2)*d**4/(rvdw**6+d**6)) + &
                      c8com  * (-ifac - 6*xij(k1)*xij(k2)/d**2 + 16*xij(k1)*xij(k2)*d**6/(rvdw**8+d**8)) + &
                      c10com * (-ifac - 8*xij(k1)*xij(k2)/d**2 + 20*xij(k1)*xij(k2)*d**8/(rvdw**10+d**10)) 
                   q(k1,i,k2,j) = qfac
                   q(k2,j,k1,i) = qfac
                enddo
             enddo
          endif
          write (uout,'(I3,X,I3,1p,E14.6,X,3(E16.9,X),2(E13.6,X))') &
             i, j, d, c6, c8, c10, rc, rvdw
       end do
    end do
    write (uout,'("#")')

    ! sum rules for the second derivatives
    do i = 1, cr%ncel
       do k1 = 1, 3
          do k2 = 1, 3
             q(k1,i,k2,i) = 0d0
             do j = 1, cr%ncel
                if (j == i) cycle
                q(k1,i,k2,i) = q(k1,i,k2,i) - q(k1,i,k2,j)
             enddo
          enddo
       enddo
    enddo

    write (uout,'("dispersion energy ",1p,E20.12)') e
    write (uout,'("dispersion forces ")')
    write (uout,'("# i          Fx                   Fy                   Fz")')
    do i = 1, cr%ncel
       write (uout,'(I3,X,1p,3(E20.12,X))') i, f(:,i)
    enddo
    write (uout,'("#")')
    write (uout,'("dispersion force constant matrix ")')
    write (uout,'("# i  xyz   j   xyz    Exixj ")')
    do i = 1, cr%ncel
       do k1 = 1, 3
          do j = 1, i-1
             do k2 = 1, 3
                write (uout,'(4(I3,X),1p,E20.12)') i, k1, j, k2, q(k1,i,k2,j)
             enddo
          enddo
          do k2 = 1, k1
             write (uout,'(4(I3,X),1p,E20.12)') i, k1, i, k2, q(k1,i,k2,i)
          enddo
       enddo
    enddo
    write (uout,'("#"/)')

  end subroutine edisp_mol

end module xdm
