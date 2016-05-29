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

contains

  !> Driver for XDM
  subroutine xdm_driver(line)
    use tools_io
    
    character*(*), intent(inout) :: line

    character(len=:), allocatable :: word
    integer :: lp 

    ! header
    write (uout,'("* Exchange-hole dipole moment (XDM) dispersion calculation ")')
    write (uout,'("  Please cite: ")')
    write (uout,'("  A. D. Becke, E. R. Johnson, J. Chem. Phys. 127, 154108 (2007).")')
    write (uout,'("  A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 136, 174109 (2012).")')
    write (uout,'("  A. Otero de la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).")')
    write (uout,*)

    lp = 1
    word = lgetword(line,lp)
    if (equal(word,"grid")) then
       call xdm_grid(line(lp:))
    elseif (equal(word,"qe")) then
       call xdm_qe()
    else
       call ferror('xdm_driver',"Unknown keyword in XDM",faterr,syntax=.true.)
       return
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
    
    integer :: i, j, ii, jj, nn, i3
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

end module xdm
