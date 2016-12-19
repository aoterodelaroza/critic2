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

!>  The following applies to the lebedev subroutines of this file:
!>  This subroutine is part of a set of subroutines that generate
!>  Lebedev grids [1-6] for integration on a sphere. The original
!>  C-code [1] was kindly provided by Dr. Dmitri N. Laikov and
!>  translated into fortran by Dr. Christoph van Wuellen.
!>  This subroutine was translated from C to fortran77 by hand.
!>  Users of this code are asked to include reference [1] in their
!>  publications, and in the user- and programmers-manuals
!>  describing their codes.
!>  This code was distributed through CCL (http://www.ccl.net/).
!>  [1] V.I. Lebedev, and D.N. Laikov
!>      "A quadrature formula for the sphere of the 131st
!>       algebraic order of accuracy"
!>      Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!>
!>  [2] V.I. Lebedev
!>      "A quadrature formula for the sphere of 59th algebraic
!>       order of accuracy"
!>      Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
!>
!>  [3] V.I. Lebedev, and A.L. Skorokhodov
!>      "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!>      Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
!>
!>  [4] V.I. Lebedev
!>      "Spherical quadrature formulas exact to orders 25-29"
!>      Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
!>
!>  [5] V.I. Lebedev
!>      "Quadratures on a sphere"
!>      Computational Mathematics and Mathematical Physics, Vol. 16,
!>      1976, pp. 10-24.
!>
!>  [6] V.I. Lebedev
!>      "Values of the nodes and weights of ninth to seventeenth
!>       order Gauss-Markov quadrature formulae invariant under the
!>       octahedron group with inversion"
!>      Computational Mathematics and Mathematical Physics, Vol. 15,
!>      1975, pp. 44-51.

!> Quadrature schemes and integration-related tools.
module integration
  implicit none

  private

  ! integration properties
  public :: intgrid_driver
  private :: intgrid_multipoles
  private :: intgrid_deloc_wfn
  private :: intgrid_deloc_wannier
  public :: int_radialquad
  public :: gauleg_msetnodes
  public :: gauleg_mquad
  public :: lebedev_msetnodes
  public :: lebedev_mquad
  private :: quadpack_f
  public :: int_output
  private :: int_output_multipoles
  private :: int_output_deloc_wfn
  private :: int_output_deloc_wannier
  private :: assign_strings
  public :: int_reorder_gridout
  private :: int_gridbasins
  private :: unpackidx
  private :: packidx

  ! grid integration types
  integer, parameter :: imtype_bader = 1
  integer, parameter :: imtype_yt = 2

contains

  subroutine intgrid_driver(line)
    use bader, only: bader_integrate
    use yt, only: yt_integrate, yt_weights, ytdata, ytdata_clean
    use fields, only: f, type_grid, type_wfn, integ_prop, itype_v, itype_deloc, itype_mpoles,&
       itype_fval, itype_f, itype_lapval, itype_lap, itype_gmod, goodfield, nprops,&
       writegrid_cube, grdall
    use grid_tools, only: grid_rhoat, grid_laplacian, grid_gradrho
    use struct_basic, only: cr
    use global, only: refden, eval_next, dunit, iunit, iunitname0, fileroot
    use tools_io, only: ferror, faterr, lgetword, equal, isexpression_or_word, uout,&
       string, fclose
    use types, only: field

    character*(*), intent(in) :: line

    real*8, parameter :: ratom_def0 = 1d0

    character(len=:), allocatable :: word, file, expr
    character*3 :: basinfmt
    integer :: i, k, n(3), nn, ntot
    integer :: lp, lp2, imtype, i1, i2, i3
    logical :: ok, nonnm, noatoms, atexist, pmask(nprops), dowcube, dobasins
    logical :: plmask(nprops)
    real*8 :: ratom, ratom_def, padd, lprop(nprops)
    real*8 :: x(3), x2(3)
    integer, allocatable :: idg(:,:,:), icp(:), idprop(:)
    real*8, allocatable :: psum(:,:), xgatt(:,:)
    real*8, allocatable :: w(:,:,:)
    integer :: fid, nattr, luw
    character*60 :: reason(nprops)
    real*8, allocatable :: sij(:,:,:,:,:), mpole(:,:,:)
    logical :: dodelocwfn, dodelocwan, dompole, fillgrd
    real*8, allocatable :: fbasin(:,:,:), fint(:,:,:,:)
    type(field) :: faux
    type(ytdata) :: dat
    
    ! only grids
    if (f(refden)%type /= type_grid) then
       call ferror("intgrid_driver","BADER/YT can only be used with grids",faterr,line,syntax=.true.)
       return
    end if

    ! method and header
    lp = 1
    word = lgetword(line,lp)
    if (equal(word,"yt")) then
       imtype = imtype_yt
    elseif (equal(word,"bader")) then
       imtype = imtype_bader
    else
       call ferror("intgrid_driver","wrong method",faterr,line,syntax=.true.)
       return
    endif

    ! parse the input
    ratom_def = ratom_def0
    nonnm = .true.
    noatoms = .false.
    dowcube = .false.
    dobasins = .false.
    basinfmt = "obj"
    expr = ""
    do while(.true.)
       word = lgetword(line,lp)
       if (equal(word,"nnm")) then
          nonnm = .false.
       elseif (equal(word,"noatoms")) then
          noatoms = .true.
       elseif (equal(word,"ratom")) then
          nonnm = .false.
          ok = eval_next(ratom_def,line,lp)
          if (.not.ok) then
             call ferror("intgrid_driver","wrong RATOM keyword",faterr,line,syntax=.true.)
             return
          end if
          ratom_def = ratom_def / dunit
       elseif (equal(word,"wcube")) then
          dowcube = .true.
       elseif (equal(word,"basins")) then
          dobasins = .true.
          lp2 = lp
          word = lgetword(line,lp)
          if (equal(word,"obj")) then
             basinfmt = "obj"
          elseif (equal(word,"ply")) then
             basinfmt = "ply"
          elseif (equal(word,"off")) then
             basinfmt = "off"
          else
             lp = lp2
             basinfmt = "obj"
          end if
       elseif (equal(word,"discard")) then
          ok = isexpression_or_word(expr,line,lp)
          if (.not. ok) then
             call ferror("intgrid_driver","wrong DISCARD keyword",faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror("intgrid_driver","Unknown extra keyword",faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! field number of grid points
    n = f(refden)%n
    ntot = n(1)*n(2)*n(3)

    ! distance for atom assignment
    if (noatoms) then
       atexist = .false.
       ratom = ratom_def
    elseif (nonnm) then
       atexist = .true.
       ratom = 1d40
    else
       atexist = .true.
       ratom = ratom_def
    end if

    ! prepare the array for the basin field
    allocate(fbasin(n(1),n(2),n(3)))
    if (f(refden)%usecore) then
       faux = f(refden)
       call grid_rhoat(faux,faux,3)
       fbasin = f(refden)%f + faux%f
       if(allocated(faux%f)) deallocate(faux%f)
    else
       fbasin = f(refden)%f
    end if

    ! call the integration method 
    luw = 0
    if (imtype == imtype_bader) then
       write (uout,'("* Henkelman et al. integration ")')
       write (uout,'("  Please cite: ")')
       write (uout,'("    G. Henkelman, A. Arnaldsson, and H. Jonsson, Comput. Mater. Sci. 36, 254-360 (2006).")')
       write (uout,'("    E. Sanville, S. Kenny, R. Smith, and G. Henkelman, J. Comput. Chem. 28, 899-908 (2007).")')
       write (uout,'("    W. Tang, E. Sanville, and G. Henkelman, J. Phys.: Condens. Matter 21, 084204 (2009).")')
       write (uout,'("+ Distance atomic assignment (",A,"): ",A)') iunitname0(iunit),&
          string(max(ratom,0d0),'e',decimal=4)
       if (len_trim(expr) > 0) &
          write (uout,'("+ Discard attractor expression: ",A)') trim(expr)
       call bader_integrate(cr,fbasin,expr,atexist,ratom,nattr,xgatt,idg)
       write (uout,'("+ Attractors in BADER: ",A)') string(nattr)
    elseif (imtype == imtype_yt) then
       write (uout,'("* Yu-Trinkle integration ")')
       write (uout,'("  Please cite: ")')
       write (uout,'("  Min Yu, Dallas Trinkle, J. Chem. Phys. 134 (2011) 064111.")')
       write (uout,'("+ Distance atomic assignment (",A,"): ",A)') iunitname0(iunit),&
          string(max(ratom,0d0),'e',decimal=4)
       if (len_trim(expr) > 0) &
          write (uout,'("+ Discard attractor expression: ",A)') trim(expr)
       call yt_integrate(cr,fbasin,expr,atexist,ratom,nattr,xgatt,idg,luw)
       write (uout,'("+ Attractors in YT: ",A)') string(nattr)
    endif
    deallocate(fbasin)

    ! reorder the attractors
    call int_reorder_gridout(cr,f(refden),nattr,xgatt,idg,atexist,ratom,luw,icp)
    write (uout,'("+ Attractors after reordering: ",A)') string(nattr)
    write (uout,*)

    ! set the properties mask
    write (uout,'("+ Calculating properties"/)') 
    dodelocwfn = .false.
    dodelocwan = .false.
    dompole = .false.
    pmask = .false.
    do k = 1, nprops
       if (.not.integ_prop(k)%used) then
          reason(k) = "not used (?)"
          cycle
       end if
       if (integ_prop(k)%itype == itype_v) then
          pmask(k) = .true.
          cycle
       end if
       fid = integ_prop(k)%fid
       if (.not.goodfield(fid)) then
          reason(k) = "unknown or invalid field"
          cycle
       elseif (integ_prop(k)%itype == itype_deloc) then
          if (f(fid)%type == type_wfn .and. cr%ismolecule) then
             dodelocwfn = .true.
             reason(k) = "DIs integrated separately (see table below)"
          elseif (f(fid)%type == type_grid .and. allocated(f(fid)%fwan)) then
             dodelocwan = .true.
             reason(k) = "DIs integrated separately (see table below)"
          else
             reason(k) = "Integrable field not a molecular wavefunction/wannier set."
          end if
          cycle
       elseif (integ_prop(k)%itype == itype_mpoles) then
          dompole = .true.
          reason(k) = "multipoles integrated separately (see table below)"
          cycle
       end if
       pmask(k) = .true.
    end do

    ! allocate and fill the integrable properties
    allocate(idprop(nprops),fint(n(1),n(2),n(3),count(pmask)))
    nn = 0
    do k = 1, nprops
       fillgrd = .false.
       if (integ_prop(k)%itype == itype_v) cycle
       if (.not.pmask(k).and..not.integ_prop(k)%itype == itype_mpoles) cycle
       fid = integ_prop(k)%fid
       nn = nn + 1
       idprop(k) = nn

       ! Use the grid values if available. Use FFT when possible.
       if (f(fid)%type == type_grid) then
          if (integ_prop(k)%itype == itype_fval .or.&
              integ_prop(k)%itype == itype_f.and..not.f(fid)%usecore.or.&
              integ_prop(k)%itype == itype_mpoles.and..not.f(fid)%usecore) then
             fint(:,:,:,nn) = f(fid)%f
          elseif (integ_prop(k)%itype == itype_lapval .or.&
                  integ_prop(k)%itype == itype_lap.and..not.f(fid)%usecore) then
             faux = f(fid)
             call grid_laplacian(f(fid),faux)
             fint(:,:,:,nn) = faux%f
          elseif (integ_prop(k)%itype == itype_gmod.and..not.f(fid)%usecore) then
             faux = f(fid)
             call grid_gradrho(f(fid),faux)
             fint(:,:,:,nn) = faux%f
          else
             fillgrd = .true.
          end if
          if(allocated(faux%f)) deallocate(faux%f)
       else
          fillgrd = .true.
       endif

       ! Calculate the grid values using the general routine (grd).
       if (fillgrd) then
          plmask = .false.
          plmask(k) = .true.
          !$omp parallel do private (x,x2,lprop) schedule(dynamic)
          do i1 = 1, n(1)
             x(1) = real(i1-1,8) / n(1)
             do i2 = 1, n(2)
                x(2) = real(i2-1,8) / n(2)
                do i3 = 1, n(3)
                   x(3) = real(i3-1,8) / n(3)
                   x2 = cr%x2c(x)
                   call grdall(x2,lprop,plmask)
                   !$omp critical (write)
                   fint(i1,i2,i3,nn) = lprop(k)
                   !$omp end critical (write)
                end do
             end do
          end do
          !$omp end parallel do
       end if
    end do

    ! compute weights and integrate the scalar field properties
    allocate(psum(nprops,nattr))
    psum = 0d0
    if (imtype == imtype_yt) then
       allocate(w(n(1),n(2),n(3)))
       w = 0d0
       call yt_weights(luw=luw,dout=dat)
    end if
    !$omp parallel do private(padd) firstprivate(w) schedule(dynamic)
    do i = 1, nattr
       if (imtype == imtype_yt) then
          call yt_weights(din=dat,idb=i,w=w)
       end if
       do k = 1, nprops
          if (.not.integ_prop(k)%used) cycle
          if (integ_prop(k)%itype == itype_v) then
             if (imtype == imtype_bader) then
                padd = count(idg == i) * cr%omega / ntot
             else
                padd = sum(w) * cr%omega / ntot
             endif
          elseif (pmask(k)) then
             if (imtype == imtype_bader) then
                padd = sum(fint(:,:,:,idprop(k)),idg==i) * cr%omega / ntot
             else
                padd = sum(w * fint(:,:,:,idprop(k))) * cr%omega / ntot
             endif
          endif
          !$omp critical (accum)
          psum(k,i) = psum(k,i) + padd
          !$omp end critical (accum)
       end do
    end do
    !$omp end parallel do

    ! write weight cubes
    if (dowcube) then
       do i = 1, nattr
          if (imtype == imtype_yt) then
             call yt_weights(din=dat,idb=i,w=w)
          else
             w = 0d0
             where (idg == i)
                w = 1d0
             end where
          endif
          file = trim(fileroot) // "_wcube_" // string(i,2,pad0=.true.) // ".cube"
          call writegrid_cube(cr,w,file,.false.)
       end do
       write (uout,'("* Weights written to ",A,"_wcube_*.cube")') trim(fileroot)
    end if

    ! clean up
    if (imtype == imtype_yt) then
       if (allocated(w)) deallocate(w)
       call ytdata_clean(dat)
    end if

    ! calculate atomic basin multipoles
    if (dompole) &
       call intgrid_multipoles(fint,idprop,nattr,xgatt,idg,imtype,luw,mpole)

    ! localization and delocalization indices - molecular wfn case
    if (dodelocwfn) then
       call intgrid_deloc_wfn(nattr,xgatt,idg,imtype,luw,sij)
    elseif (dodelocwan) then
       call intgrid_deloc_wannier(nattr,xgatt,idg,imtype,luw,sij)
    end if

    ! output the results for the attractor and the scalar properties
    call int_output(pmask,reason,nattr,icp,xgatt,psum,.false.,sij,mpole)

    ! multipole output
    if (dompole) &
       call int_output_multipoles(nattr,icp,mpole)

    ! localization and delocalization indices output
    if (dodelocwfn) then
       call int_output_deloc_wfn(nattr,icp,sij)
    else
       call int_output_deloc_wannier(nattr,icp,xgatt,sij)
    end if

    ! bains plotting
    if (dobasins) &
       call int_gridbasins(basinfmt,nattr,icp,xgatt,idg,imtype,luw)

    ! clean up YT weight file
    if (imtype == imtype_yt) then
       call fclose(luw)
    endif
    deallocate(psum,idg,icp,xgatt)

  end subroutine intgrid_driver

  !> Calculate the multipole moments using integration on a grid. The
  !> input is: the array containing the integrable field information
  !> on the basin field grid (fint), the index array relating
  !> 1->nprops to the fint array, the calculated number of attractors
  !> (natt), their positions in crystallographic coordinates (xgatt),
  !> the attractor assignment for each of the grid nodes (idg), the
  !> integration type (imtype: bader or yt), the weight file for YT,
  !> and the output multipolar moments.
  subroutine intgrid_multipoles(fint,idprop,natt,xgatt,idg,imtype,luw,mpole)
    use yt, only: yt_weights, ytdata, ytdata_clean
    use fields, only: nprops, integ_prop, itype_mpoles, f
    use struct_basic, only: cr
    use global, only: refden
    use tools_math, only: tosphere, genrlm_real

    real*8, intent(in) :: fint(:,:,:,:)
    integer, intent(in) :: idprop(:)
    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: imtype
    integer, intent(in) :: luw
    real*8, allocatable, intent(inout) :: mpole(:,:,:)

    integer :: i, j, k, l, m, n(3), ntot, np
    integer :: ix
    integer :: fid, lmax
    real*8, allocatable :: w(:,:,:)
    real*8 :: dv(3), r, tp(2), p(3)
    real*8, allocatable :: rrlm(:)
    type(ytdata) :: dat

    ! allocate space for the calculated multipoles
    lmax = -1
    np = 0
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_mpoles) cycle
       lmax = max(integ_prop(l)%lmax,lmax)
       np = np + 1
    end do
    if (np == 0) return

    if (allocated(mpole)) deallocate(mpole)
    allocate(mpole((lmax+1)*(lmax+1),natt,np))

    ! size of the grid and YT weights
    n = f(refden)%n
    ntot = n(1)*n(2)*n(3)
    if (imtype == imtype_yt) then
       allocate(w(n(1),n(2),n(3)))
       call yt_weights(luw=luw,dout=dat)
    endif

    ! run over all properties for which multipole calculation is active
    np = 0
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_mpoles) cycle
       fid = integ_prop(l)%fid
       np = np + 1

       ! allocate space
       lmax = integ_prop(l)%lmax
       allocate(rrlm((lmax+1)*(lmax+1)))

       ! calcualate the multipoles
       mpole(:,:,np) = 0d0
       if (imtype == imtype_bader) then
          !$omp parallel do private(p,ix,dv,r,tp,rrlm) schedule(dynamic)
          do i = 1, n(1)
             p(1) = real(i-1,8)
             do j = 1, n(2)
                p(2) = real(j-1,8)
                do k = 1, n(3)
                   p(3) = real(k-1,8)
                   ix = idg(i,j,k)
                   dv = p/real(n,8) - xgatt(:,ix)
                   call cr%shortest(dv,r)
                   call tosphere(dv,r,tp)
                   call genrlm_real(lmax,r,tp,rrlm)

                   !$omp critical (accum)
                   mpole(:,ix,np) = mpole(:,ix,np) + rrlm * fint(i,j,k,idprop(l))
                   !$omp end critical (accum)
                end do
             end do
          end do
          !$omp end parallel do
       else
          !$omp parallel do private(p,dv,r,tp,rrlm) firstprivate(w) schedule(dynamic)
          do m = 1, natt
             call yt_weights(din=dat,idb=m,w=w)
             do i = 1, n(1)
                do j = 1, n(2)
                   do k = 1, n(3)
                      if (abs(w(i,j,k)) < 1d-15) cycle
                      p = real((/i-1,j-1,k-1/),8)
                      dv = p/real(n,8) - xgatt(:,m)
                      call cr%shortest(dv,r)
                      call tosphere(dv,r,tp)
                      call genrlm_real(lmax,r,tp,rrlm)

                      !$omp critical (accum)
                      mpole(:,m,np) = mpole(:,m,np) + rrlm * fint(i,j,k,idprop(l)) * w(i,j,k)
                      !$omp end critical (accum)
                   end do
                end do
             end do
          end do
          !$omp end parallel do
       endif
       mpole = mpole * cr%omega / ntot
       deallocate(rrlm)
    end do
    if (imtype == imtype_yt) then
       deallocate(w)
       call ytdata_clean(dat)
    endif

  end subroutine intgrid_multipoles

  !> Calculate localization and delocalization indices using the
  !> basin assignment found by YT or BADER and a wfn scalar field.
  subroutine intgrid_deloc_wfn(natt,xgatt,idg,imtype,luw,sij)
    use yt, only: yt_weights
    use wfn_private, only: wfn_rho2
    use fields, only: nprops, integ_prop, itype_deloc, f, type_wfn
    use struct_basic, only: cr
    use global, only: refden
    use tools_io, only: fopen_scratch, fclose

    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: imtype
    integer, intent(in) :: luw
    real*8, intent(inout), allocatable :: sij(:,:,:,:,:)

    integer :: i, j, k, l, m, i1, i2, ndeloc
    integer :: fid, n(3), ntot, lumo, ix
    real*8, allocatable :: w(:,:,:)
    real*8, allocatable :: xmo(:)
    real*8 :: x(3), rho, auxg(3), auxh(3,3), gkin, vir
    real*8 :: stress(3,3)

    ndeloc = 0
    do l = 1, nprops
       ! only wfn fields for now
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_wfn) cycle
       ndeloc = ndeloc + 1
    end do
    if (ndeloc == 0) return

    ! allocate output sij
    if (allocated(sij)) deallocate(sij)
    allocate(sij(f(fid)%nmo,f(fid)%nmo,natt,1,ndeloc))

    ! size of the grid
    n = f(refden)%n
    ntot = n(1)*n(2)*n(3)

    if (imtype == imtype_yt) then
       allocate(w(n(1),n(2),n(3)))
    endif

    ! run over all properties for which multipole calculation is active
    ndeloc = 0
    do l = 1, nprops
       ! only wfn fields for now
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_wfn) cycle
       ndeloc = ndeloc + 1

       ! calculate the overlap matrix
       allocate(xmo(f(fid)%nmo))
       sij(:,:,:,:,ndeloc) = 0d0
       xmo = 0d0
       if (imtype == imtype_bader) then
          !$omp parallel do private (i,j,k,x,rho,auxg,auxh,ix,i1,i2,gkin,vir) firstprivate(xmo) schedule(dynamic)
          do i = 1, n(1)
             do j = 1, n(2)
                do k = 1, n(3)
                   x = cr%x2c(real((/i-1,j-1,k-1/),8) / n)
                   call wfn_rho2(f(fid),x,0,rho,auxg,auxh,gkin,vir,stress,xmo)
                   ix = idg(i,j,k)
                   do i1 = 1, f(fid)%nmo
                      do i2 = i1, f(fid)%nmo
                         !$omp critical (sijwrite)
                         sij(i1,i2,ix,1,ndeloc) = sij(i1,i2,ix,1,ndeloc) + xmo(i1) * xmo(i2)
                         !$omp end critical (sijwrite)
                      end do
                   end do
                end do
             end do
          end do
          !$omp end parallel do
       else
          ! write the MO vectors to a file
          lumo = fopen_scratch()
          do i = 1, n(1)
             do j = 1, n(2)
                do k = 1, n(3)
                   x = cr%x2c(real((/i-1,j-1,k-1/),8) / n)
                   call wfn_rho2(f(fid),x,0,rho,auxg,auxh,gkin,vir,stress,xmo)
                   ix = idg(i,j,k)
                   write (lumo) xmo
                end do
             end do
          end do

          ! calculate the atomic overlap matrix and kill the file
          do m = 1, natt
             call yt_weights(luw=luw,idb=m,w=w)
             rewind(lumo)
             do i = 1, n(1)
                do j = 1, n(2)
                   do k = 1, n(3)
                      read(lumo) xmo
                      do i1 = 1, f(fid)%nmo
                         do i2 = i1, f(fid)%nmo
                            sij(i1,i2,m,1,ndeloc) = sij(i1,i2,m,1,ndeloc) + xmo(i1) * xmo(i2) * w(i,j,k)
                         end do
                      end do
                   end do
                end do
             end do
          end do
          call fclose(lumo)
       end if
       deallocate(xmo)

       ! symmetrize and scale
       sij(:,:,:,:,ndeloc) = sij(:,:,:,:,ndeloc) * cr%omega / ntot
       do ix = 1, natt
          do i1 = 1, f(fid)%nmo
             do i2 = i1, f(fid)%nmo
                sij(i2,i1,ix,1,ndeloc) = sij(i1,i2,ix,1,ndeloc)
             end do
          end do
       end do
    end do
    if (allocated(w)) deallocate(w)

  end subroutine intgrid_deloc_wfn

  !> Calculate localization and delocalization indices using the
  !> basin assignment found by YT or BADER and a grid-related dat file.
  subroutine intgrid_deloc_wannier(natt,xgatt,idg,imtype,luw,sij)
    use yt, only: yt_weights, ytdata, ytdata_clean
    use fields, only: f, integ_prop, nprops, itype_deloc, type_grid
    use struct_basic, only: cr
    use global, only: refden
    use types, only: realloc
    use tools_io, only: uout, string

    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: imtype
    integer, intent(in) :: luw
    real*8, allocatable, intent(inout) :: sij(:,:,:,:,:)

    integer :: is, nspin, ndeloc, natt1
    integer :: ia, ja, ka, iba, ib, jb, kb, ibb
    integer :: i, j, l
    integer :: nwan(3), m1, m2, m3
    integer :: fid, n(3), p(3)
    integer :: nbnd, nlat, nmo, imo, jmo, imo1, jmo1
    real*8, allocatable :: w(:,:,:)
    real*8, allocatable :: psic(:,:,:), psic2(:,:,:)
    real*8 :: x(3), xs(3), d2, padd
    logical :: found
    integer, allocatable :: idg1(:,:,:), iatt(:), ilvec(:,:)
    logical, allocatable :: wmask(:,:,:)
    type(ytdata) :: dat

    ! run over properties with wannier delocalization indices; allocate sij
    ndeloc = 0
    nspin = 1
    nmo = 0
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_grid .or..not.allocated(f(fid)%fwan)) cycle
       ndeloc = ndeloc + 1
       nspin = max(nspin,size(f(fid)%fwan,5))
       nbnd = size(f(fid)%fwan,4)
       nlat = f(fid)%nwan(1)*f(fid)%nwan(2)*f(fid)%nwan(3)
       nmo = max(nmo,nlat * nbnd)
    end do
    if (ndeloc == 0) return

    ! size of the grid and header
    n = f(refden)%n
    write (uout,'("+ Calculating atomic overlap matrices")')

    ! yt data
    if (imtype == imtype_yt) then
       call yt_weights(luw=luw,dout=dat)
    end if

    ! Recalculate the number of attractors without cell translation symmetry.
    allocate(iatt(natt))
    natt1 = natt
    do i = 1, natt
       iatt(i) = i
    enddo
    allocate(ilvec(3,natt))
    ilvec = 0
    write (uout,'(99(A,X))') "  Attractors before remapping =", string(natt)
    if (imtype == imtype_bader) then
       allocate(idg1(n(1),n(2),n(3)))
       do m3 = 1, n(3)
          do m2 = 1, n(2)
             do m1 = 1, n(1)
                idg1(m1,m2,m3) = idg(m1,m2,m3)
                p = (/m1,m2,m3/)
                x = real(p-1,8) / n - xgatt(:,idg(m1,m2,m3))
                xs = x
                call cr%shortest(xs,d2)
                p = nint(x - cr%c2x(xs))
                if (any(p /= 0)) then
                   found = .false.
                   do i = natt+1, natt1
                      if (iatt(i) == idg(m1,m2,m3) .and. all(p == ilvec(:,i))) then
                         found = .true.
                         idg1(m1,m2,m3) = i
                         exit
                      end if
                   end do
                   if (.not.found) then
                      natt1 = natt1 + 1
                      if (natt1 > size(ilvec,2)) then
                         call realloc(ilvec,3,2*natt1)
                         call realloc(iatt,2*natt1)
                      end if
                      ilvec(:,natt1) = p
                      idg1(m1,m2,m3) = natt1
                      iatt(natt1) = idg(m1,m2,m3)
                   end if
                end if
             end do
          end do
       end do
    else
       allocate(w(n(1),n(2),n(3)))
       do i = 1, natt
          call yt_weights(din=dat,idb=i,w=w)
          do m3 = 1, n(3)
             do m2 = 1, n(2)
                do m1 = 1, n(1)
                   if (abs(w(m1,m2,m3)) < 1d-15) cycle
                   p = (/m1,m2,m3/)
                   x = real(p-1,8) / n - xgatt(:,i)
                   xs = x
                   call cr%shortest(xs,d2)
                   p = nint(x - cr%c2x(xs))
                   if (any(p /= 0)) then
                      found = .false.
                      do j = natt+1, natt1
                         if (iatt(j) == i .and. all(p == ilvec(:,j))) then
                            found = .true.
                            exit
                         end if
                      end do
                      if (.not.found) then
                         natt1 = natt1 + 1
                         if (natt1 > size(ilvec,2)) then
                            call realloc(ilvec,3,2*natt1)
                            call realloc(iatt,2*natt1)
                         end if
                         ilvec(:,natt1) = p
                         iatt(natt1) = i
                      end if
                   end if
                end do
             end do
          end do
       end do
       deallocate(w)
    end if
    call realloc(ilvec,3,natt1)
    call realloc(iatt,natt1)
    write (uout,'(99(A,X))') "  Attractors after remapping =", string(natt1)

    ! header
    if (allocated(sij)) deallocate(sij)
    allocate(sij(nmo,nmo,natt,nspin,ndeloc))
    sij = 0d0

    ! run over all properties
    ndeloc = 0
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_grid .or..not.allocated(f(fid)%fwan)) cycle
       ndeloc = ndeloc + 1

       ! assign values to some integers
       nbnd = size(f(fid)%fwan,4)
       nwan = f(fid)%nwan
       nlat = nwan(1)*nwan(2)*nwan(3)
       nmo = nlat * nbnd
       nspin = size(f(fid)%fwan,5)

       ! write out some info
       write (uout,'("# Integrated property (number ",A,"): ",A)') string(l), string(integ_prop(l)%prop_name)
       write (uout,'(99(A,X))') "  Number of bands =", string(nbnd)
       write (uout,'(99(A,X))') "  ... lattice translations =", (string(f(fid)%nwan(j)),j=1,3)
       write (uout,'(99(A,X))') "  ... Wannier functions =", string(nmo)
       write (uout,'(99(A,X))') "  ... spin channels =", string(nspin)
     
       write (uout,'(99(A,X))') "  Calculating overlaps..."
       allocate(psic(n(1),n(2),n(3)))
       sij(:,:,:,:,ndeloc) = 0d0

       if (imtype == imtype_bader) then
          ! bader integration
          do is = 1, nspin
             do imo = 1, nmo
                call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nwan)
                do jmo = imo, nmo
                   call unpackidx(jmo,ib,jb,kb,ibb,nmo,nbnd,nwan)
                   psic = f(fid)%fwan(ia*n(1)+1:(ia+1)*n(1),ja*n(2)+1:(ja+1)*n(2),ka*n(3)+1:(ka+1)*n(3),iba,is) * &
                      f(fid)%fwan(ib*n(1)+1:(ib+1)*n(1),jb*n(2)+1:(jb+1)*n(2),kb*n(3)+1:(kb+1)*n(3),ibb,is)
                   do i = 1, natt1
                      call packidx(ia+ilvec(1,i),ja+ilvec(2,i),ka+ilvec(3,i),iba,imo1,nmo,nbnd,nwan)
                      call packidx(ib+ilvec(1,i),jb+ilvec(2,i),kb+ilvec(3,i),ibb,jmo1,nmo,nbnd,nwan)
                      sij(imo1,jmo1,iatt(i),is,ndeloc) = sij(imo1,jmo1,iatt(i),is,ndeloc) + sum(psic,idg1==i)
                      sij(jmo1,imo1,iatt(i),is,ndeloc) = sij(imo1,jmo1,iatt(i),is,ndeloc)
                   end do
                end do
             end do
          end do
       else
          ! yt integration
          allocate(w(n(1),n(2),n(3)),wmask(n(1),n(2),n(3)),psic2(n(1),n(2),n(3)))
          w = 0d0
          wmask = .false.
          psic2 = 0d0
          psic = 0d0
          !$omp parallel do private(p,x,xs,d2,ia,ja,ka,iba,ib,jb,kb,ibb,padd,imo1,jmo1)&
          !$omp firstprivate(w,wmask,psic2,psic) schedule(dynamic)
          do i = 1, natt1
             call yt_weights(din=dat,idb=iatt(i),w=w)
             wmask = .false.
             do m3 = 1, n(3)
                do m2 = 1, n(2)
                   do m1 = 1, n(1)
                      if (abs(w(m1,m2,m3)) < 1d-15) cycle
                      p = (/m1,m2,m3/)
                      x = real(p-1,8) / n - xgatt(:,iatt(i))
                      xs = x
                      call cr%shortest(xs,d2)
                      p = nint(x - cr%c2x(xs))
                      wmask(m1,m2,m3) = all(p == ilvec(:,i))
                   end do
                end do
             end do
             psic2 = 0d0
             do is = 1, nspin
                do imo = 1, nmo
                   call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nwan)
                   where (wmask)
                      psic2 = f(fid)%fwan(ia*n(1)+1:(ia+1)*n(1),ja*n(2)+1:(ja+1)*n(2),ka*n(3)+1:(ka+1)*n(3),iba,is) * w
                   end where
                   do jmo = imo, nmo
                      call unpackidx(jmo,ib,jb,kb,ibb,nmo,nbnd,nwan)
                      where (wmask)
                         psic =  psic2 * &
                            f(fid)%fwan(ib*n(1)+1:(ib+1)*n(1),jb*n(2)+1:(jb+1)*n(2),kb*n(3)+1:(kb+1)*n(3),ibb,is)
                      end where
                      padd = sum(psic,wmask)
                      call packidx(ia+ilvec(1,i),ja+ilvec(2,i),ka+ilvec(3,i),iba,imo1,nmo,nbnd,nwan)
                      call packidx(ib+ilvec(1,i),jb+ilvec(2,i),kb+ilvec(3,i),ibb,jmo1,nmo,nbnd,nwan)
                      !$omp critical (sum)
                      sij(imo1,jmo1,iatt(i),is,ndeloc) = sij(imo1,jmo1,iatt(i),is,ndeloc) + padd
                      sij(jmo1,imo1,iatt(i),is,ndeloc) = sij(imo1,jmo1,iatt(i),is,ndeloc)
                      !$omp end critical (sum)
                   end do
                end do
             end do
          end do
          !$omp end parallel do
          deallocate(w,wmask,psic2)
       end if
       ! scale (the omega comes from wannier)
       sij(:,:,:,:,ndeloc) = sij(:,:,:,:,ndeloc) / (n(1)*n(2)*n(3))
       write (uout,'("  Done."/)')

       ! wrap up
       deallocate(psic)
    end do

    ! clean up
    if (allocated(idg1)) deallocate(idg1)
    deallocate(iatt,ilvec)

  end subroutine intgrid_deloc_wannier

  !> Do a radial numerical quadrature on the given ray, between the
  !> selected radii and return the properties. The r^2 factor is
  !> included.
  subroutine int_radialquad(x,theta,phi,r0,rend,lprop,abserr,neval,iaserr,ierr)
    use quadpack, only: dqags, dqng, dqag
    use global, only: int_radquad_type, int_gauleg, int_radquad_nr, int_qags, &
       int_radquad_relerr, int_qng, int_radquad_abserr, int_qag, int_iasprec
    use fields, only: nprops, grdall
    use tools_math, only: gauleg

    real*8, intent(in) :: x(3)  !< The center of the basin (cartesian coords)
    real*8, intent(in) :: theta !< Polar angle of the ray
    real*8, intent(in) :: phi   !< Azimuthal angle or the ray
    real*8, intent(in) :: r0    !< Left limit of the radial interval
    real*8, intent(in) :: rend  !< Right limit of the radial interval
    real*8, intent(out) :: lprop(Nprops) !< The integrated properties
    real*8, intent(out) :: abserr !< Estimated absolute error
    integer, intent(out) :: neval !< Number of evaluations of grdall
    real*8, intent(out) :: iaserr(Nprops) !< Estimated IAS precision error
    integer, intent(out) :: ierr

    integer :: k, isign
    real*8 :: unit(3), xaux(3), r1, r2
    real*8, dimension(Nprops) :: atprop
    integer :: ier
    ! quadpack internals
    integer, parameter :: limit = 10000
    integer, parameter :: lenw = 4*limit
    integer :: iwork(limit)
    real*8 :: work(lenw)
    integer :: last, key
    ! radial integration
    real*8, allocatable :: rpoints(:), rweights(:)

    unit = (/ sin(theta) * cos(phi),&
              sin(theta) * sin(phi),&
              cos(theta) /)
    ierr = 0

    if (rend > r0) then
       r1 = r0
       r2 = rend
       isign = 1
    else
       r1 = rend
       r2 = r0
       isign = -1
    end if

    select case(INT_radquad_type)
    case (INT_gauleg)
       ! allocate and set radial nodes
       allocate(rpoints(INT_radquad_nr),rweights(INT_radquad_nr))
       call gauleg(r1,r2,rpoints,rweights,INT_radquad_nr)

       ! do the quadrature
       lprop = 0d0
       do k = 1, INT_radquad_nr
          xaux = x + rpoints(k) * unit
          call grdall(xaux,atprop)
          lprop = lprop + rpoints(k)**2 * rweights(k) * atprop
       end do
       lprop = lprop * isign
       lprop(1) = isign * (r2**3 - r1**3) / 3d0 ! hardwire volume
       abserr = 0d0
       neval = INT_radquad_nr

       ! deallocate
       deallocate(rpoints,rweights)

    case (INT_qags)
       call dqags(quadpack_f,unit,x,r0,rend,INT_radquad_abserr,&
          INT_radquad_relerr,lprop,abserr,neval,&
          ier,limit,lenw,last,iwork,work)
       ierr = max(ierr,ier)

    case (INT_qng)
       call dqng(quadpack_f,unit,x,r0,rend,INT_radquad_abserr,&
          INT_radquad_relerr,lprop,abserr,neval,ier)
       ierr = max(ierr,ier)

    case (INT_qag)
       if (INT_radquad_nr >= 30) then
          key = 6
       else if (INT_radquad_nr >= 25) then
          key = 5
       else if (INT_radquad_nr >= 20) then
          key = 4
       else if (INT_radquad_nr >= 15) then
          key = 3
       else if (INT_radquad_nr >= 10) then
          key = 2
       else
          key = 1
       end if
       call dqag(quadpack_f,unit,x,r0,rend,INT_radquad_abserr,&
          INT_radquad_relerr,key,lprop,abserr,neval,&
          ier,limit,lenw,last,iwork,work)
       ierr = max(ierr,ier)
    end select

    ! IAS precision error
    xaux = x + rend * unit
    call grdall(xaux,atprop)
    iaserr = abs(atprop * rend**2 * INT_iasprec)
    neval = neval + 1

  end subroutine int_radialquad

  !> Sets the angular nodes for the quadrature of a sphere by treating
  !> ntheta and nphi as 1d Gauss-Legendre quadratures. For a given
  !> theta (polar angle), the number of phi (azimuthal angle) points
  !> is int(nphi*abs(sin(theta)))+1, where nphi is provided in the input.
  !> minisurface version
  subroutine gauleg_msetnodes(srf,ntheta,nphi)
    use types, only: minisurf
    use param, only: pi
    use tools_math, only: gauleg

    type(minisurf), intent(inout) :: srf !< Surface where the node info is placed
    integer, intent(in) :: ntheta !< Number of phi points
    integer, intent(in) :: nphi   !< Number of azimuthal points

    integer :: j, k
    integer :: realnphi
    ! Nodes and weights
    real*8, allocatable, dimension(:) :: tpoints
    real*8, allocatable, dimension(:) :: tweights
    real*8, allocatable, dimension(:) :: ppoints
    real*8, allocatable, dimension(:) :: pweights

    allocate(tpoints(ntheta),tweights(ntheta))
    allocate(ppoints(nphi+1),pweights(nphi+1))
    ! place the gauss-legendre nodes on the surface
    call gauleg(0d0,pi,tpoints,tweights,ntheta)
    srf%nv = 0
    do j = 1, ntheta
       realnphi = int(nphi*abs(sin(tpoints(j))))+1
       call gauleg(0d0,2d0*pi,ppoints,pweights,realnphi)
       do k = 1, realnphi
          srf%nv = srf%nv + 1
          srf%r(srf%nv)  = -1d0
          srf%th(srf%nv) = tpoints(j)
          srf%ph(srf%nv) = ppoints(k)
       end do
    end do

    deallocate(tpoints,tweights,ppoints,pweights)

  end subroutine gauleg_msetnodes

  !> Sum the gauss-legendre 2d (ntheta * nphi) quadrature. The srf
  !> minisurface contains the rays generated by gauleg_msetnodes, with
  !> radii corresponding to the IAS. ntheta and nphi. The resulting
  !> integrated properties are returned in the lprop vector. abserr is
  !> the integrated error of the radial quadrature. iaserr is the error
  !> caused by IAS inaccuracies.  neval is the number of evaluations of
  !> grdall used. rbeta is a reference radius (beta-sphere).
  subroutine gauleg_mquad(srf,ntheta,nphi,rbeta,lprop,abserr,neval,iaserr)
    use tools_math, only: gauleg
    use tools_io, only: ferror, warning, uout
    use fields, only: nprops
    use types, only: minisurf
    use param, only: pi

    type(minisurf), intent(inout) :: srf !< Surface representing the basin
    integer, intent(in) :: ntheta !< Number of polar points
    integer, intent(in) :: nphi   !< Number of azimuthal points
    real*8, intent(in) :: rbeta   !< Beta-spehre radius
    real*8, intent(out) :: lprop(Nprops) !< Properties vector
    real*8, intent(out) :: abserr !< Integrated radial quad. error
    integer, intent(out) :: neval !< Number of evaluations
    real*8, intent(out) :: iaserr(Nprops) !< Integrated IAS precision error

    integer :: i, j, v, leval
    real*8 :: rprop(Nprops), pprop(Nprops)
    real*8 :: rerr, perr, riaserr(Nprops), piaserr(Nprops)
    integer :: realnphi, ierr, err, vidx(ntheta,nphi+1)
    ! Nodes and weights
    real*8, allocatable, dimension(:) :: tpoints
    real*8, allocatable, dimension(:) :: tweights
    real*8, allocatable, dimension(:) :: ppoints
    real*8, allocatable, dimension(:) :: pweights

    allocate(tpoints(ntheta),tweights(ntheta))
    allocate(ppoints(nphi+1),pweights(nphi+1))
    ppoints = 0d0
    pweights = 0d0

    ! initialize
    neval = 0
    lprop = 0d0
    abserr = 0d0
    iaserr = 0d0
    err = 0

    ! set theta-nodes
    call gauleg(0d0,pi,tpoints,tweights,ntheta)

    ! fill idx array
    v = 0
    do i = 1, ntheta
       do j = 1, int(nphi*abs(sin(tpoints(i))))+1
          v = v + 1
          vidx(i,j) = v
       end do
    end do

    !$omp parallel do reduction(+:neval,abserr) reduction(max:err) private(&
    !$omp realnphi,pprop,rprop,perr,piaserr,leval,ierr,rerr,&
    !$omp riaserr,v) firstprivate(ppoints,pweights) schedule(guided)
    do i = 1, ntheta
       realnphi = int(nphi*abs(sin(tpoints(i))))+1
       call gauleg(0d0,2d0*pi,ppoints,pweights,realnphi)
       pprop = 0d0
       perr = 0d0
       piaserr = 0d0
       do j = 1, realnphi
          v = vidx(i,j)
          call int_radialquad(srf%n,srf%th(v),srf%ph(v),rbeta,srf%r(v),&
             rprop,rerr,leval,riaserr,ierr)
          pprop = pprop + rprop * pweights(j)
          perr = perr + rerr * pweights(j)
          piaserr = piaserr + riaserr * pweights(j)
          neval = neval + leval
          err = max(err,ierr)
       end do
       !$omp critical (sum_gauleg_mquad)
       lprop = lprop + sin(srf%th(v)) * tweights(i) * pprop
       iaserr = iaserr + sin(srf%th(v)) * tweights(i) * piaserr
       !$omp end critical (sum_gauleg_mquad)
       abserr = abserr + sin(srf%th(v)) * tweights(i) * perr
    end do
    !$omp end parallel do

    deallocate(tpoints,tweights,ppoints,pweights)

    if (err /= 0) then
       call ferror('gauleg_mquad','Radial integration had non-zero error code',warning)
       write (uout,'(a,I2)') " ier = ", err
       write (uout,'(a)') " Check the routine documentation for more info."
    end if

  end subroutine gauleg_mquad

  !> Sets the angular nodes for the quadrature of a sphere using
  !> Lebedev method (see above for references), where npts is the
  !> number of points, that must be one for which a cubature
  !> exists (6, 14, 26, 38, 50, 74,
  !> 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974,
  !> 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,
  !> 5294, 5810). The (theta,phi) information is written to the surface.
  !> minisurface version.
  subroutine lebedev_msetnodes(srf,npts)
    use tools_math, only: select_lebedev
    use types, only: minisurf

    type(minisurf), intent(inout) :: srf !< Surface where the node info is placed
    integer, intent(in) :: npts !< Number of quadrature points

    integer :: j
    ! Nodes and weights
    real*8 :: xleb(5810)
    real*8 :: yleb(5810)
    real*8 :: zleb(5810)
    real*8 :: wleb(5810)

    call select_lebedev(npts,xleb,yleb,zleb,wleb)

    srf%nv = 0
    do j = 1, npts
       srf%nv = srf%nv + 1
       srf%r(srf%nv)  = -1d0
       srf%th(srf%nv) = acos(zleb(j))
       srf%ph(srf%nv) = atan2(yleb(j),xleb(j))
    end do

  end subroutine lebedev_msetnodes

  !> Sum the Lebedev 2d quadrature with npts points. The srf
  !> minisurface contains the rays generated by lebedev_msetnodes, with
  !> radii corresponding to the IAS. ntheta and nphi. The resulting
  !> integrated properties are returned in the lprop vector. abserr is
  !> the integrated error of the radial quadrature. iaserr is the error
  !> caused by IAS inaccuracies.  neval is the number of evaluations of
  !> grdall used. rbeta is a reference radius (beta-sphere).
  subroutine lebedev_mquad(srf,npts,rbeta,lprop,abserr,neval,iaserr)
    use tools_math, only: select_lebedev
    use tools_io, only: uout, ferror, warning
    use fields, only: nprops
    use types, only: minisurf

    type(minisurf), intent(inout) :: srf !< Surface representing the basin
    integer, intent(in) :: npts   !< Number of points
    real*8, intent(in) :: rbeta   !< Beta-spehre radius
    real*8, intent(out) :: lprop(Nprops) !< Properties vector
    real*8, intent(out) :: abserr !< Integrated radial quad. error
    integer, intent(out) :: neval !< Number of evaluations
    real*8, intent(out) :: iaserr(Nprops) !< Integrated IAS precision error

    integer :: i, v, leval
    real*8 :: rprop(Nprops)
    real*8 :: rerr, riaserr(Nprops)
    integer :: ierr, err
    ! Nodes and weights
    real*8 :: xleb(5810)
    real*8 :: yleb(5810)
    real*8 :: zleb(5810)
    real*8 :: wleb(5810)

    call select_lebedev(npts,xleb,yleb,zleb,wleb)

    neval = 0
    lprop = 0d0
    abserr = 0d0
    iaserr = 0d0
    v = 0
    err = 0

    !$omp parallel do private(rprop,rerr,leval,riaserr,ierr) reduction(+:abserr,neval)&
    !$omp reduction(max:err) schedule(guided)
    do i = 1, npts
       call int_radialquad(srf%n,srf%th(i),srf%ph(i),rbeta,srf%r(i),&
          rprop,rerr,leval,riaserr,ierr)
       !$omp critical (sum_lebedev_mquad)
       lprop = lprop + rprop * wleb(i)
       iaserr = iaserr + riaserr * wleb(i)
       !$omp end critical (sum_lebedev_mquad)
       abserr = abserr + rerr  * wleb(i)
       neval = neval + leval
       err = max(err,ierr)
    end do
    !$omp end parallel do

    if (err /= 0) then
       call ferror('lebedev_mquad','Radial integration had non-zero error code',warning)
       write (uout,'(a,I2)') " ier = ", err
       write (uout,'(a)') " Check the routine documentation for more info."
    end if

  end subroutine lebedev_mquad

  !> Dummy function for quadpack integration
  function quadpack_f(x,unit,xnuc)
    use fields, only: Nprops, grdall

    real*8, intent(in) :: x
    real*8, intent(in) :: unit(3)
    real*8, intent(in) :: xnuc(3)
    real*8 :: quadpack_f(Nprops)

    real*8 :: xaux(3)

    xaux = xnuc + x * unit
    call grdall(xaux,quadpack_f)
    quadpack_f = quadpack_f * x * x

  end function quadpack_f

  !> Output routine for all integration methods. All input
  !> arguments. pmask: properties mask (skips inactive
  !> properties). reason: explanation why a property has been
  !> deactivated. nattr: number of attractors. icp: identifier
  !> relating the attractors to known critical points. xattr:
  !> attractor positions (crysatllographic coords.). aprop(i,j):
  !> integrated scalar property j for attractor i. The j index runs
  !> over active properties only (pmask = .true.). Usesym: deactivate
  !> writing the attractor multiplicities in non-symmetric cases.  di:
  !> delocalization indices (optional). mpole: multipoles (optional).
  subroutine int_output(pmask,reason,nattr,icp,xattr,aprop,usesym,sij,mpole)
    use fields, only: integ_prop, itype_v, itype_expr, itype_mpoles, itype_names,&
       itype_deloc, f, type_wfn, nprops
    use struct_basic, only: cr
    use global, only: iunitname0, iunit, dunit
    use varbas, only: cp, cpcel
    use tools_io, only: uout, string, ioj_left, ioj_center, ioj_right

    logical, intent(in) :: pmask(nprops)
    character*(*), intent(in) :: reason(nprops)
    integer, intent(in) :: nattr
    integer, intent(in) :: icp(nattr)
    real*8, intent(in) :: xattr(3,nattr)
    real*8, intent(in) :: aprop(nprops,nattr)
    logical, intent(in) :: usesym
    real*8, intent(in), allocatable, optional :: sij(:,:,:,:,:)
    real*8, intent(in), allocatable, optional :: mpole(:,:,:)

    integer :: i, j, k, l, n, ip, ipmax, iplast, nn, ndeloc
    integer :: fid, nacprop(5), lmax
    real*8 :: x(3), sump(nprops), xmult
    character(len=:), allocatable :: saux, itaux, label, cini
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    character(len=:), allocatable :: sout1, sout2

    ! List of integrable properties accepted/rejected and why some
    ! were rejected
    write (uout,'("* List of properties integrated in the attractor basins")')
    write (uout,'("+ The ""Label"" entries will be used to identify the integrable")')
    write (uout,'("  in the tables of integrated atomic properties. The entries with")')
    write (uout,'("  an ""x"" will not be integrated.")')
    write (uout,'("# id    Label      fid  Field       Additional")')
    do i = 1, nprops
       fid = integ_prop(i)%fid
       if (.not.pmask(i)) then
          label = "--inactive--"
          cini = "x "
          saux = "Reason: " // string(reason(i))
       else
          label = integ_prop(i)%prop_name
          cini = "  "
       end if
       if (integ_prop(i)%itype == itype_v) then
          if (pmask(i)) saux = ""
          itaux = "--"
       elseif (integ_prop(i)%itype == itype_expr) then
          if (pmask(i)) saux = "expr = " // string(integ_prop(i)%expr)
          itaux = "--"
       elseif (integ_prop(i)%itype == itype_mpoles) then
          if (pmask(i)) saux = "Lmax = " // string(integ_prop(i)%lmax)
          itaux = string(fid)
       else
          if (pmask(i)) saux = ""
          itaux = string(fid)
       end if

       write (uout,'(A2,99(A,X))') &
          cini,&
          string(i,4,ioj_left), &
          string(label,12,ioj_center), &
          string(itaux,2,ioj_right), &
          string(itype_names(integ_prop(i)%itype),12,ioj_left),&
          string(saux)
    end do
    write (uout,*)

    ! List of attractors and positions
    write (uout,'("* List of attractors integrated")')
    if (.not.cr%ismolecule) then
       write (uout,'("# Id   cp   ncp   Name  Z   mult           Position (cryst.) ")')
    else
       write (uout,'("# Id   cp   ncp   Name  Z   mult           Position (",A,") ")') iunitname0(iunit)
    endif
    do i = 1, nattr
       call assign_strings(i,icp(i),usesym,scp,sncp,sname,smult,sz)
       if (.not.cr%ismolecule) then
          x = xattr(:,i)
       else
          x = (cr%x2c(xattr(:,i)) + cr%molx0) * dunit
       endif
       write (uout,'(2X,99(A,X))') & 
          string(i,4,ioj_left), scp, sncp, sname, sz, &
          smult, (string(x(j),'f',12,7,4),j=1,3)
    end do
    write (uout,*)

    ! Integrated scalar atomic properties
    write (uout,'("* Integrated atomic properties")')
    iplast = 0
    do ip = 0, (count(pmask)-1)/5
       ! show only the properties that are active
       nacprop = 0
       ipmax = 0
       do i = iplast+1, nprops
          if (pmask(i)) then
             ipmax = ipmax + 1
             nacprop(ipmax) = i
          end if
          if (ipmax == 5) exit
       end do
       if (ipmax == 0) exit
       iplast = nacprop(ipmax)

       ! Table header for this set of properties
       write (uout,'("# Integrable properties ",A," to ",A)') string(nacprop(1)), string(nacprop(ipmax))
       write (uout,'("# Id   cp   ncp   Name  Z   mult ",5(A,X))') &
          (string(integ_prop(nacprop(j))%prop_name,15,ioj_center),j=1,ipmax)

       ! Table rows
       sump = 0d0
       do i = 1, nattr
          call assign_strings(i,icp(i),usesym,scp,sncp,sname,smult,sz)
          ! add to the sum
          if (icp(i) > 0 .and. usesym) then
             xmult = cp(cpcel(icp(i))%idx)%mult
          else
             xmult = 1
          endif
          do j = 1, ipmax
             sump(nacprop(j)) = sump(nacprop(j)) + aprop(nacprop(j),i) * xmult
          end do
          ! table entry
          write (uout,'(2X,99(A,X))') &
             string(i,4,ioj_left), scp, sncp, sname, sz, smult, &
             (string(aprop(nacprop(j),i),'e',15,8,4),j=1,ipmax)
       end do
       write (uout,'(32("-"),99(A))') ("----------------",j=1,ipmax)
       write (uout,'(2X,"Sum                         ",99(A,X))') &
          (string(sump(nacprop(j)),'e',15,8,4),j=1,ipmax)
       write (uout,*)
    end do
    
  end subroutine int_output

  !> Output calculate multipole moments. Inputs: number of attractors
  !> (nattr), CP identifiers for the attractors (icp), and the
  !> multipoles themselves (mpole).
  subroutine int_output_multipoles(nattr,icp,mpole)
    use fields, only: integ_prop, itype_mpoles, nprops, f
    use tools_io, only: uout, string, ioj_left, ioj_center
    integer, intent(in) :: nattr
    integer, intent(in) :: icp(nattr)
    real*8, intent(in), allocatable, optional :: mpole(:,:,:)
    
    integer :: i, j, k, nn, n, l, lmax, fid
    integer, allocatable :: l_(:), m_(:)
    character(len=:), allocatable :: lbl
    character*3 :: ls, ms
    character(len=:), allocatable :: sncp, scp, sname, sz, smult

    if (.not.present(mpole)) return
    if (.not.allocated(mpole)) return
    n = 0
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_mpoles) cycle
       fid = integ_prop(l)%fid
       write (uout,'("* Basin multipole moments (using real solid harmonics)")')
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(integ_prop(l)%prop_name)
       write (uout,'("+ The calculated multipoles are: ")')
       write (uout,'("    Q_lm^A = int_A rho(r) * Rlm(r) dr ")')
       write (uout,'("  where the integral is over the basin of A, and Rlm is a real solid harmonic.")')
       write (uout,'("  The coordinates are referred to the attractor of the A basin. All quantities")')
       write (uout,'("  in atomic units.")')
       write (uout,*)
       lmax = integ_prop(l)%lmax
       n = n + 1

       ! figure out the indices for the labels
       allocate(l_((lmax+1)*(lmax+1)),m_((lmax+1)*(lmax+1)))
       nn = 0
       do i = 0, lmax
          do j = -i, i
             nn = nn + 1
             l_(nn) = i
             m_(nn) = j
          end do
       end do

       ! pretty output
       nn = (lmax+1)**2/5
       if (mod((lmax+1)**2,5) /= 0) nn = nn + 1
       do i = 1, nn
          ! header
          lbl = "# Id   cp   ncp   Name  Z   "
          do j = (i-1)*5+1,min(5*i,(lmax+1)**2)
             if (j == 1) then
                lbl = lbl // " " // string("1",15,justify=ioj_center)
             elseif (j == 2) then
                lbl = lbl // " " // string("x",15,justify=ioj_center)
             elseif (j == 3) then
                lbl = lbl // " " // string("z",15,justify=ioj_center)
             elseif (j == 4) then
                lbl = lbl // " " // string("y",15,justify=ioj_center)
             elseif (j == 5) then
                lbl = lbl // " " // string("sq(3)/2(x^2-y^2)",15,justify=ioj_center)
             elseif (j == 6) then
                lbl = lbl // " " // string("sq(3)xz",15,justify=ioj_center)
             elseif (j == 7) then
                lbl = lbl // " " // string("(3z^2-r^2)/2",15,justify=ioj_center)
             elseif (j == 8) then
                lbl = lbl // " " // string("sq(3)yz",15,justify=ioj_center)
             elseif (j == 9) then
                lbl = lbl // " " // string("sq(3)xy",15,justify=ioj_center)
             else
                write (ls,'(I3)') l_(j)
                write (ms,'(I3)') abs(m_(j))
                if (m_(j) <= 0) then
                   lbl = lbl // " " // string("C("//string(ls)//","//string(ms)//")",15,justify=ioj_center)
                else
                   lbl = lbl // " " // string("S("//string(ls)//","//string(ms)//")",15,justify=ioj_center)
                end if
             endif
          end do
          write (uout,'(A)') lbl

          ! body
          do j = 1, nattr
             call assign_strings(j,icp(j),.false.,scp,sncp,sname,smult,sz)
             write (uout,'(2X,99(A,X))') &
                string(j,4,ioj_left), scp, sncp, sname, sz, &
                (string(mpole((i-1)*5+1+k,j,n),'e',15,8,4),k=0,min(4,size(mpole,1)-(i-1)*5-1))
          enddo
          if (i < nn) then
             write (uout,*)
          else
             write (uout,'(32("-"),99(A))') ("----------------",j=1,5)
          end if
       end do
       deallocate(l_,m_)
       write (uout,*)
    end do

  end subroutine int_output_multipoles

  !> Output delocalization indices in a molecular/wfn calculation.
  !> Inputs: number of attractors (nattr), CP identifiers for the
  !> attractors (icp), and the atomic overlap matrix (sij).
  subroutine int_output_deloc_wfn(nattr,icp,sij)
    use fields, only: integ_prop, itype_deloc, f, type_wfn, nprops
    use wfn_private, only: wfn_rhf
    use struct_basic, only: cr
    use tools_io, only: uout, string, ioj_left
    integer, intent(in) :: nattr
    integer, intent(in) :: icp(nattr)
    real*8, intent(in), allocatable, optional :: sij(:,:,:,:,:)

    integer :: i, j, l, fid
    integer :: ndeloc
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    character(len=:), allocatable :: sout1, sout2
    real*8, allocatable :: di(:,:)
    real*8 :: fspin

    if (.not.present(sij)) return
    if (.not.allocated(sij)) return
    if (.not.cr%ismolecule) return

    ! space for the DIs
    allocate(di(nattr,nattr))

    ndeloc = 0
    do l = 1, nprops
       if (.not.cr%ismolecule) cycle
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_wfn) cycle
       ndeloc = ndeloc + 1

       ! calculate localization and delocalization indices using
       ! the atomic overlap matrix
       di = 0d0
       if (f(fid)%wfntyp == wfn_rhf) then
          fspin = 4d0
       else
          fspin = 1d0
       end if
       do i = 1, nattr
          do j = i, nattr
             di(i,j) = fspin * sum(sij(:,:,i,1,ndeloc)*sij(:,:,j,1,ndeloc))
             di(j,i) = di(i,j)
          end do
       end do

       ! Header
       write (uout,'("* Localization and delocalization indices")')
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(integ_prop(l)%prop_name)

       ! output the lambdas
       write (uout,'("+ Localization indices (lambda(A))")')
       write (uout,'("# Id   cp   ncp   Name  Z      lambda(A)  ")')
       do j = 1, nattr
          call assign_strings(j,icp(j),.false.,scp,sncp,sname,smult,sz)
          write (uout,'(2X,99(A,X))') &
             string(j,4,ioj_left), scp, sncp, sname, sz, &
             string(0.5d0*di(j,j),'e',15,8,4)
       enddo
       write (uout,*)

       write (uout,'("+ Delocalization indices (delta(A,B))")')
       write (uout,'("#     ----- atom A -----          ----- atom B -----                       ")')
       write (uout,'("#   Id   cp   ncp   Name  Z     Id   cp   ncp   Name  Z        delta(A,B)  ")')
       do i = 1, nattr
          call assign_strings(i,icp(i),.false.,scp,sncp,sname,smult,sz)
          sout1 = "| " // string(i,4,ioj_left) // " " // scp // " " // sncp // " " // sname // " " // sz // " " // " |"
          do j = i+1, nattr
             call assign_strings(j,icp(j),.false.,scp,sncp,sname,smult,sz)
             sout2 = string(j,4,ioj_left) // " " // scp // " " // sncp // " " // sname // " " // sz // " " // " |"
             write (uout,'(2X,99(A,X))') string(sout1), string(sout2), &
                string(di(i,j),'e',15,8,4)
          end do
       end do
       write (uout,*)
    end do
          
    ! clean up
    deallocate(di)

  end subroutine int_output_deloc_wfn

  !> Output delocalization indices in a wannier set calculation.
  !> Inputs: number of attractors (nattr), CP identifiers for the
  !> attractors (icp), and the atomic overlap matrix (sij).
  subroutine int_output_deloc_wannier(natt,icp,xgatt,sij)
    use fields, only: integ_prop, itype_deloc, f, type_grid, nprops
    use global, only: refden, iunit, iunitname0, dunit
    use struct_basic, only: cr, crystal
    use tools, only: qcksort
    use tools_io, only: uout, string, ioj_left, ioj_right
    integer, intent(in) :: natt
    integer, intent(in) :: icp(natt)
    real*8, intent(in) :: xgatt(3,natt)
    real*8, intent(in), allocatable, optional :: sij(:,:,:,:,:)

    integer :: l, m, fid, ndeloc, nspin, nbnd, nlat, nmo, nwan(3)
    real*8 :: fspin, xnn, xli, r1(3), d2, asum
    integer, allocatable :: imap(:,:), io(:), ilvec(:,:), idat(:)
    real*8, allocatable :: fa(:,:,:,:), diout(:), dist(:)
    integer :: imo, jmo, ia, ja, ka, iba
    integer :: ic, jc, kc, i, j, k, idx(3), is
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    type(crystal) :: cr1

    if (.not.present(sij)) return
    if (.not.allocated(sij)) return

    write (uout,'("* Localization and delocalization indices")')

    ndeloc = 0
    do l = 1, nprops
       ! skip the incorrect properties
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_grid .or..not.allocated(f(fid)%fwan)) cycle
       ndeloc = ndeloc + 1

       ! header
       write (uout,'("+ Integrated property (number ",A,"): ",A/)') string(l), string(integ_prop(l)%prop_name)

       ! some integers for the run
       nwan = f(fid)%nwan
       nspin = size(f(fid)%fwan,5)
       nbnd = size(f(fid)%fwan,4)
       nlat = nwan(1)*nwan(2)*nwan(3)
       nmo = nlat * nbnd

       ! spin factor
       if (nspin == 1) then
          fspin = 2d0
       else
          fspin = 1d0
       end if
       
       ! calculate the index mapping in order to build the fa matrix
       ! kmo = imap(imo,jlat) gives the Wannier function index kmo
       ! resulting from taking the Wannier function imo and
       ! translating by lattice vector jlat
       allocate(imap(nmo,nlat))
       do imo = 1, nmo
          call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nwan)
          k = 0
          do ic = 0, nwan(1)-1
             do jc = 0, nwan(2)-1
                do kc = 0, nwan(3)-1
                   k = k + 1
                   idx = (/ic+ia, jc+ja, kc+ka/)
                   idx = modulo(idx,f(fid)%nwan)
                   call packidx(idx(1),idx(2),idx(3),iba,imap(imo,k),nmo,nbnd,nwan)
                end do
             end do
          end do
       end do

       ! calculate the values of the fa matrix
       allocate(fa(natt,natt,nlat,nspin))
       fa = 0d0
       do is = 1, nspin
          do k = 1, nlat
             do i = 1, natt
                do j = 1, natt
                   fa(i,j,k,is) = 0d0
                   do imo = 1, nmo
                      do jmo = 1, nmo
                         fa(i,j,k,is) = fa(i,j,k,is) + sij(imo,jmo,i,is,ndeloc) * sij(imap(imo,k),imap(jmo,k),j,is,ndeloc)
                      end do
                   end do
                end do
             end do
          end do
       end do
       deallocate(imap)
       
       ! localization indices
       write (uout,'("+ Localization indices")')
       write (uout,'("# Id   cp   ncp   Name  Z       LI(A)           N(A)")')
       do i = 1, natt
          call assign_strings(i,icp(i),.false.,scp,sncp,sname,smult,sz)
          xli = sum(abs(fa(i,i,1,:))) * fspin
          xnn = sum(abs(fa(i,:,:,:))) * fspin
          write (uout,'(2X,99(A,X))') & 
             string(i,4,ioj_left), scp, sncp, sname, sz, &
             string(xli,'f',15,8,4), string(xnn,'f',12,8,4)
       end do
       write (uout,*)

       ! Prepare the supercell crystal. For c%shortest, we need
       ! c%isortho, c%car2crys, c%crys2car, c%nws, and c%ivws. The
       ! last two come from a call to c%wigner.
       call cr1%init()
       do i = 1, 3
          cr1%crys2car(:,i) = cr%crys2car(:,i) * nwan(i)
          cr1%car2crys(i,:) = cr%car2crys(i,:) / nwan(i)
       end do
       call cr1%wigner((/0d0,0d0,0d0/),nvec=cr1%nws,vec=cr1%ivws)
       cr1%isortho = (cr1%nws <= 6)
       if (cr1%isortho) then
          do i = 1, cr1%nws
             cr1%isortho = cr1%isortho .and. (count(abs(cr1%ivws(:,i)) == 1) == 1)
          end do
       endif

       write (uout,'("+ Delocalization indices Each block gives information about a single atom in")')
       write (uout,'("  the main cell. First line: localization index, next lines: delocaliazation index")')
       write (uout,'("  with all atoms in the environment, last line: sum of LI + 0.5 * DIs, equal to the")')
       write (uout,'("  atomic population. Distances are in ",A,".")') iunitname0(iunit)
       allocate(dist(natt*nlat),io(natt*nlat),diout(natt*nlat),ilvec(3,natt*nlat),idat(natt*nlat))
       do i = 1, natt
          call assign_strings(i,icp(i),.false.,scp,sncp,sname,smult,sz)
          write (uout,'("# Attractor ",A," (cp=",A,", ncp=",A,", name=",A,", Z=",A,") at: ",3(A,2X))') & 
             string(i), trim(scp), trim(sncp), trim(adjustl(sname)), trim(sz), (trim(string(xgatt(j,i),'f',12,7)),j=1,3)
          write (uout,'("# Id   cp   ncp   Name  Z    Latt. vec.     ----  Cryst. coordinates ----       Distance        LI/DI")')

          ! precompute the localization/delocalization indices for this atom and
          ! location and distance information
          dist = 0d0
          k = 0
          m = 0
          do ic = 0, nwan(1)-1
             do jc = 0, nwan(2)-1
                do kc = 0, nwan(3)-1
                   k = k + 1
                   do j = 1, natt
                      m = m + 1
                      io(m) = m
                      r1 = (xgatt(:,j) + (/ic,jc,kc/) - xgatt(:,i)) / real(nwan,8)
                      call cr1%shortest(r1,d2)
                      dist(m) = sqrt(d2) * dunit
                      diout(m) = 2d0 * sum(abs(fa(i,j,k,:))) * fspin
                      if (dist(m) < 1d-5) diout(m) = diout(m) / 2d0
                      idat(m) = j
                      ilvec(:,m) = nint(xgatt(:,i) + cr1%c2x(r1) * nwan - xgatt(:,j))
                   end do
                end do
             end do
          end do

          ! sort by increasing distance and output for this atom
          call qcksort(dist,io,1,natt*nlat)
          asum = 0d0
          do m = 1, natt*nlat
             j = io(m)
             if (dist(j) < 1d-5) then
                write (uout,'(2X,"Localization index",71("."),A)') string(diout(j),'f',12,8,4)
                asum = asum + diout(j)
             else
                call assign_strings(j,icp(idat(j)),.false.,scp,sncp,sname,smult,sz)
                r1 = xgatt(:,idat(j)) + ilvec(:,j)
                write (uout,'(2X,99(A,X))') string(j,4,ioj_left), scp, sncp, sname, sz,&
                   (string(ilvec(k,j),3,ioj_right),k=1,3), (string(r1(k),'f',12,7,4),k=1,3),&
                   string(dist(j),'f',12,7,4), string(diout(j),'f',12,8,4)
                asum = asum + 0.5d0 * diout(j)
             end if
          end do
          write (uout,'(2X,"Total (atomic population)",64("."),A)') string(asum,'f',12,8,4)
          write (uout,*)
       end do

       ! clean up
       call cr1%end()
       deallocate(fa,dist,io,diout,ilvec,idat)
    end do

  contains

  end subroutine int_output_deloc_wannier

  !> Assign strings for attractor i. icp is the CP identifier for the
  !> attractor. usesym controls whether the multiplicit will be
  !> written to the output. The output strings are scp (the cell CP
  !> identifier, or -- if not a CP attractor), sncp (the
  !> non-equivalent CP identifier, or -- if not a CP attractor), sname
  !> (the name of the CP), smult (CP multiplicity), and sz (atomic
  !> number, if the CP is an atom).
  subroutine assign_strings(i,icp,usesym,scp,sncp,sname,smult,sz)
    use varbas, only: cp, cpcel
    use struct_basic, only: cr
    use tools_io, only: string, ioj_left, ioj_center
    integer, intent(in) :: i
    integer, intent(in) :: icp
    logical, intent(in) :: usesym
    character(len=:), allocatable, intent(out) :: sncp, scp, sname, sz, smult

    integer :: idx

    if (icp > 0) then
       ! this is a cp
       idx = cpcel(icp)%idx
       scp = string(icp,4,ioj_left)
       sncp = string(idx,4,ioj_left)
       sname = string(cp(idx)%name,6,ioj_center)
       smult = string(cp(idx)%mult,4,ioj_center)
       if (cp(idx)%isnuc) then
          sz = string(cr%at(idx)%z,2,ioj_left)
       else
          sz = "--"
       endif
    else
       ! this is an unknown nnm
       scp = " -- "
       sncp = " -- "
       sname = "  ??  "
       smult = string(1,4,ioj_center)
       sz = "--"
    end if
    if (.not.usesym) smult = " -- "
    
  end subroutine assign_strings

  !> The attractors coming out of YT and BADER are not in order
  !> compatible with the crystal structure. Reorder them, including
  !> the weights of the YT. On output, give the identity 
  !> of the attractors (icp) in the complete CP list.
  subroutine int_reorder_gridout(c,ff,nattr,xgatt,idg,atexist,ratom,luw,icp)
    use autocp, only: addcp
    use varbas, only: ncpcel, cpcel, nearest_cp
    use struct_basic, only: crystal
    use tools_io, only: ferror, faterr, fopen_scratch, fclose
    use types, only: field, realloc

    type(crystal), intent(in) :: c
    type(field), intent(in) :: ff
    integer, intent(inout) :: nattr
    real*8, intent(inout), allocatable :: xgatt(:,:)
    integer, intent(inout), allocatable :: idg(:,:,:)
    logical, intent(in) :: atexist
    real*8, intent(in) :: ratom
    integer, intent(inout) :: luw
    integer, intent(inout), allocatable :: icp(:)

    integer :: i, j
    integer, allocatable :: idgaux(:,:,:)
    real*8, allocatable :: xattr(:,:)
    integer, allocatable :: assigned(:)
    integer :: nn, nid, nattr0, luw2, n(3), lvec(3), nbasin, nvec
    real*8 :: dist
    integer, allocatable :: nlo(:), ibasin(:), ibasin2(:), iio(:), inear(:,:)
    real*8, allocatable :: fnear(:,:)

    n = ff%n

    ! reorder the maxima and assign maxima to atoms according to ratom
    if (allocated(icp)) deallocate(icp)
    allocate(icp(nattr),xattr(3,nattr),assigned(nattr))
    icp = 0
    assigned = 0
    nattr0 = nattr
    ! assign attractors to atoms
    if (atexist) then
       do i = 1, nattr0
          nid = 0
          call c%nearest_atom(xgatt(:,i),nid,dist,lvec)
          if (dist < ratom) then
             assigned(i) = nid
          else
             ! maybe the closest point is a known nnm
             call nearest_cp(xgatt(:,i),nid,dist,ff%typnuc)
             if (dist < ratom) then
                assigned(i) = nid
             end if
          end if
       end do
    end if
    ! create the new known attractors in the correct order
    nattr = 0
    do i = 1, ncpcel
       if (any(assigned == i)) then
          nattr = nattr + 1
          icp(nattr) = i
          xattr(:,nattr) = cpcel(i)%x
       endif
    end do
    ! the rest are their own nnm, add them to the CP list, accumulate
    ! using a radius equal to ratom.
    do i = 1, nattr0
       if (assigned(i) > 0) cycle
       nattr = nattr + 1
       assigned(i) = nattr
       icp(nattr) = 0

       nn = 1
       xattr(:,nattr) = xgatt(:,i)
       do j = i+1, nattr0
          if (assigned(j) > 0) cycle
          if (c%are_lclose(xgatt(:,i),xgatt(:,j),ratom)) then
             nn = nn + 1
             assigned(j) = nattr
          end if
       end do
       call addcp(c%x2c(xattr(:,nattr)),"",ff%typnuc)
    end do
    deallocate(xgatt)

    ! update the idg
    allocate(idgaux(size(idg,1),size(idg,2),size(idg,3)))
    idgaux = idg
    do i = 1, nattr0
       if (assigned(i) /= i) then
          where (idg == i)
             idgaux = assigned(i)
          end where
       end if
    end do
    call move_alloc(idgaux,idg)

    ! update the weights the YT file
    if (luw /= 0) then
       ! read all the info from the scratch file
       rewind(luw)
       read (luw) nbasin, nn, nvec
       if (nattr0 /= nbasin) &
          call ferror('int_reorder_gridout','inconsistent number of attractors in yt checkpoint',faterr)
       allocate(ibasin(nn),ibasin2(nn),nlo(nn),inear(nvec,nn),fnear(nvec,nn),iio(nn))
       read (luw) nlo
       read (luw) ibasin
       read (luw) iio
       read (luw) inear
       read (luw) fnear

       ! build a new ibasin array
       ibasin2 = ibasin
       do i = 1, nattr0
          where (ibasin == i)
             ibasin2 = assigned(i)
          end where
       end do

       ! write the new data to a new scratch file
       luw2 = fopen_scratch()
       write (luw2) nattr0, nn, nvec
       write (luw2) nlo 
       write (luw2) ibasin2
       write (luw2) iio 
       write (luw2) inear
       write (luw2) fnear
       call fclose(luw)
       luw = luw2
    end if
    deallocate(assigned)
    
    if (allocated(xgatt)) deallocate(xgatt)
    call realloc(xattr,3,nattr)
    call move_alloc(xattr,xgatt)

  end subroutine int_reorder_gridout

  !> Plot the atomic basins found by YT or BADER to graphical files.
  !> Use format fmt (obj, ply, off). nattr: number of attractors.
  !> icp: identification of the attractors as CPs. xgatt: crsyt.
  !> coords. of the attractors. idg: basin assignment of each
  !> point in the grid. imtype: type of grid integration carried out.
  !> luw: logical unit for YT weights.
  subroutine int_gridbasins(fmt,nattr,icp,xgatt,idg,imtype,luw)
    use global, only: fileroot
    use struct_basic, only: cr, crystal
    use struct_writers, only: struct_write_3dmodel
    use graphics, only: graphics_open, graphics_close, graphics_polygon
    use tools_io, only: string, uout
    use yt, only: yt_weights, ytdata_clean, ytdata
    character*3, intent(in) :: fmt
    integer, intent(in) :: nattr
    integer, intent(in), allocatable :: icp(:)
    real*8, intent(in) :: xgatt(3,nattr)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: imtype
    integer, intent(in) :: luw

    integer :: lu, lumtl
    character(len=:), allocatable :: str
    integer :: i1, i2, i3, n(3), i, j, k, p(3), q(3)
    real*8 :: x(3), xd(3), d2
    type(crystal) :: caux
    real*8, allocatable :: xface(:,:), w(:,:,:)
    integer, allocatable :: idg0(:,:,:)
    type(ytdata) :: dat

    integer, parameter :: rgb1(3) = (/128,128,128/)

    ! prepare wigner-seitz tetrahedra
    do i = 1, 3
       n(i) = size(idg,i)
    end do
    caux%isinit = .true.
    caux%aa = cr%aa / real(n,8)
    caux%bb = cr%bb
    call caux%set_cryscar()
    call caux%wigner((/0d0,0d0,0d0/),nvec=caux%nws,vec=caux%ivws,&
       nvert_ws=caux%nvert_ws,nside_ws=caux%nside_ws,iside_ws=caux%iside_ws,&
       vws=caux%vws)
    allocate(xface(3,size(caux%iside_ws,1)))

    ! output
    write (uout,'("* Basins written to ",A,"_basins-*.",A/)') trim(fileroot), fmt

    ! write the unit cell
    str = trim(fileroot) // "_basins-cell." // fmt
    call struct_write_3dmodel(cr,str,fmt,(/1,1,1/),.true.,.false.,.true.,&
       .true.,.true.,-1d0,(/0d0,0d0,0d0/),-1d0,(/0d0,0d0,0d0/))

    ! prepare the idg array
    allocate(idg0(n(1),n(2),n(3)))
    if (imtype == imtype_bader) then
       idg0 = idg
    elseif (imtype == imtype_yt) then
       allocate(w(n(1),n(2),n(3)))
       call yt_weights(luw=luw,dout=dat)
       do i = 1, nattr
          call yt_weights(din=dat,idb=i,w=w)
          where(w >= 0.5d0)
             idg0 = i
          end where
       end do
       call ytdata_clean(dat)
       deallocate(w)
    endif

    ! write the basins
    do i = 1, nattr
       ! name this file
       str = trim(fileroot) // "_basins-" // string(i) // "." // fmt
       call graphics_open(fmt,str,lu,lumtl)

       do i1 = 1, n(1)
          do i2 = 1, n(2)
             do i3 = 1, n(3)
                if (idg0(i1,i2,i3) == i) then
                   p = (/i1,i2,i3/)
                   x = real(p-1,8) / n

                   ! move to the ws of the attractor
                   xd = x - xgatt(:,i)
                   call cr%shortest(xd,d2)

                   ! convert to Cartesian
                   x = cr%x2c(xgatt(:,i)) + xd

                   ! plot, if on the border of the cell
                   do j = 1, caux%nws
                      q = modulo(p + caux%ivws(:,j) - 1,n) + 1
                      if (idg0(q(1),q(2),q(3)) /= i) then
                         do k = 1, caux%nside_ws(j)
                            xface(:,k) = x + caux%x2c(caux%vws(:,caux%iside_ws(k,j)))
                         end do
                         call graphics_polygon(fmt,lu,xface,rgb1)
                      end if
                   end do
                end if
             end do
          end do
       end do

       call graphics_close(fmt,lu,lumtl)
    end do

    deallocate(xface)

  end subroutine int_gridbasins

  !> Unpacking routine for use in Wannier delocalization index calculations.
  subroutine unpackidx(idx,io,jo,ko,bo,nmo,nbnd,nwan)
    integer, intent(in) :: idx, nmo, nbnd, nwan(3)
    integer, intent(out) :: io, jo, ko, bo

    integer :: iaux

    ! unpack
    iaux = modulo(idx-1,nmo)
    bo = modulo(iaux,nbnd) + 1
    iaux = (idx-1 - (bo-1)) / nbnd
    ko = modulo(iaux,nwan(3)) + 1
    iaux = (iaux - (ko-1)) / nwan(3)
    jo = modulo(iaux,nwan(2)) + 1
    iaux = (iaux - (jo-1)) / nwan(2)
    io = modulo(iaux,nwan(1)) + 1

    ! translate
    io = io - 1
    jo = jo - 1
    ko = ko - 1

  end subroutine unpackidx

  !> Unpacking routine for use in Wannier delocalization index calculations.
  subroutine packidx(io,jo,ko,bo,idx,nmo,nbnd,nwan)
    integer, intent(in) :: io, jo, ko, bo, nmo, nbnd, nwan(3)
    integer, intent(out) :: idx

    integer :: zio, zjo, zko

    ! transformed indices
    zio = modulo(io,nwan(1)) + 1
    zjo = modulo(jo,nwan(2)) + 1
    zko = modulo(ko,nwan(3)) + 1

    ! translate and pack
    idx = (bo-1) + nbnd * ((zko-1) + nwan(3) * ((zjo-1) + nwan(2) * (zio-1))) + 1

  end subroutine packidx

end module integration
