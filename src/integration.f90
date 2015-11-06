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
  use fields, only: nprops
  implicit none

  private

  ! integration properties
  public :: intgrid_driver
  private :: intgrid_output
  private :: intgrid_multipoles
  private :: intgrid_deloc_wfn
  private :: intgrid_deloc_brf
  public :: int_radialquad
  public :: gauleg_msetnodes
  public :: gauleg_mquad
  public :: lebedev_msetnodes
  public :: lebedev_mquad
  private :: quadpack_f
  public :: ppty_output
  public :: int_output

  ! grid integration types
  integer, parameter :: itype_bader = 1
  integer, parameter :: itype_yt = 2

contains

  subroutine intgrid_driver(line)
    use bader
    use yt
    use surface
    use fields
    use struct_basic
    use global
    use types, only: gk
    use tools_math
    use tools_io
    use param

    character*(*), intent(in) :: line

    real*8, parameter :: ratom_def0 = 1d0

    character(len=:), allocatable :: word
    integer :: i, j, k, n(3), nn, ntot
    integer :: lp, itype, nid, lvec(3), p(3)
    logical :: ok, nonnm, noatoms, doflist, doescher, pmask(nprops)
    real*8 :: ratom, dv(3), dist, r, tp(2), ratom_def
    integer, allocatable :: idg(:,:,:), idatt(:), idx(:)
    real*8, allocatable :: psum(:,:), xgatt(:,:)
    real(kind=gk), allocatable :: w(:,:,:)
    integer :: fid, natt, ix, luw
    character*3, allocatable :: sidx(:)

    ! only grids
    if (f(refden)%type /= type_grid) &
       call ferror("intgrid_driver","bader can only be used with grids",faterr,line)

    ! method and header
    lp = 1
    word = lgetword(line,lp)
    if (equal(word,"yt")) then
       itype = itype_yt
    elseif (equal(word,"bader")) then
       itype = itype_bader
    else
       call ferror("intgrid_driver","wrong method",faterr,line)
    endif

    ! parse the input
    ratom_def = ratom_def0
    nonnm = .true.
    noatoms = .false.
    doflist = .false.
    doescher = .false.
    do while(.true.)
       word = lgetword(line,lp)
       if (equal(word,"nnm")) then
          nonnm = .false.
       elseif (equal(word,"noatoms")) then
          noatoms = .true.
       elseif (equal(word,"ratom")) then
          nonnm = .false.
          ok = eval_next(ratom_def,line,lp)
          if (.not.ok) &
             call ferror("intgrid_driver","wrong RATOM keyword",faterr,line)
       elseif (equal(word,"escher")) then
          doescher = .true.
       elseif (equal(word,"fieldlist")) then
          doflist = .true.
       elseif (len_trim(word) > 0) then
          call ferror("intgrid_driver","Unknown extra keyword",faterr,line)
       else
          exit
       end if
    end do

    ! field number of grid points
    n = f(refden)%n
    ntot = n(1)*n(2)*n(3)

    ! distance for atom assignment
    if (noatoms) then
       ratom = -1d0
    elseif (nonnm) then
       ratom = 1d40
    else
       ratom = ratom_def
    end if

    ! call the integration method proper
    luw = 0
    if (itype == itype_bader) then
       write (uout,'("* Henkelman et al. integration ")')
       write (uout,'("  Please cite: ")')
       write (uout,'("    G. Henkelman, A. Arnaldsson, and H. Jonsson, Comput. Mater. Sci. 36, 254-360 (2006).")')
       write (uout,'("    E. Sanville, S. Kenny, R. Smith, and G. Henkelman, J. Comput. Chem. 28, 899-908 (2007).")')
       write (uout,'("    W. Tang, E. Sanville, and G. Henkelman, J. Phys.: Condens. Matter 21, 084204 (2009).")')
       write (uout,'("+ Distance atomic assignment (bohr): ",A)') string(max(ratom,0d0),'e',decimal=4)
       call bader_integrate(cr,f(refden),natt,xgatt,idatt,idg,ratom)
    elseif (itype == itype_yt) then
       write (uout,'("* Yu-Trinkle integration ")')
       write (uout,'("  Please cite: ")')
       write (uout,'("  Min Yu, Dallas Trinkle, J. Chem. Phys. 134 (2011) 064111.")')
       write (uout,'("+ Distance atomic assignment (bohr): ",A)') string(max(ratom,0d0),'e',decimal=4)
       call yt_integrate(natt,xgatt,idatt,idg,ratom,luw)
    endif
    write (uout,*)

    ! reorder the maxima
    allocate(idx(natt),sidx(natt))
    ! first the atoms
    idx = 0
    nn = 0
    do i = 1, cr%ncel
       if (all(idatt /= i)) cycle
       do k = 1, natt
          if (idatt(k) == i) exit
       end do
       nn = nn + 1
       idx(nn) = k
       sidx(nn) = trim(adjustl(cr%at(cr%atcel(i)%idx)%name))
    end do
    !! then the rest
    do i = 1, natt
       if (any(idx == i)) cycle
       nn = nn + 1
       idx(nn) = i
       write (sidx(nn),'("n",I2.2)') nn
    end do

    ! find the properties mask and output the rejected, summary of integrable properties
    pmask = .false.
    write (uout,'("+ List of integrable properties used in this run")')
    do k = 1, nprops
       if (integ_prop(k)%itype == itype_v) then
          pmask(k) = .true.
       else
          fid = integ_prop(k)%fid
          if (fid < 0 .or. fid > mf) then
             write (uout,'("  Integrable ",A," is rejected: unknown field")') string(k)
             cycle
          endif
          if (f(fid)%type /= type_grid) then
             write (uout,'("  Integrable ",A," is rejected: not a grid")') string(k)
             cycle
          endif
          if (any(f(fid)%n /= f(refden)%n)) then
             write (uout,'("  Integrable ",A," is rejected: grid size not the same as reference field")') &
                string(k)
             cycle
          endif
          if (.not.integ_prop(k)%itype == itype_f .and..not.integ_prop(k)%itype == itype_fval) then
             write (uout,'("  Integrable ",A," is rejected: only field values (f,fval) can be integrated on a grid")') &
                string(k)
             cycle
          endif
          pmask(k) = .true.
       end if
       if (pmask(k)) then
          write (uout,'("  Integrable ",A," will be used under header label ''",A,"''")') &
             string(k), string(integ_prop(k)%prop_name)
       endif
    end do
    write (uout,*)

    ! compute weights and integrate the scalar field properties
    if (itype == itype_yt) then
       rewind(luw)
    endif
    allocate(psum(nprops,natt))
    allocate(w(n(1),n(2),n(3)))
    psum = 0d0
    do i = 1, natt
       if (itype == itype_yt) then
          read (luw) w
       end if
       do k = 1, nprops
          if (.not.integ_prop(k)%used) cycle
          if (integ_prop(k)%itype == itype_v) then
             if (itype == itype_bader) then
                psum(k,i) = count(idg == i) * cr%omega / ntot
             else
                psum(k,i) = sum(w) * cr%omega / ntot
             endif
          elseif (pmask(k)) then
             fid = integ_prop(k)%fid
             if (itype == itype_bader) then
                w = 0d0
                where (idg == i)
                   w = f(fid)%f
                end where
                psum(k,i) = sum(w) * cr%omega / ntot
             else
                psum(k,i) = sum(w * f(fid)%f) * cr%omega / ntot
             endif
          endif
       end do
    end do
    deallocate(w)

    ! output all normal properties
    call intgrid_output(noatoms,doescher,doflist,natt,idatt,psum,pmask,xgatt,idx,sidx)

    ! compute and output multipoles
    call intgrid_multipoles(natt,xgatt,idatt,idg,itype,luw,idx,sidx)

    ! localization and delocalization indices
    call intgrid_deloc_wfn(natt,xgatt,idatt,idg,itype,luw,idx,sidx)
    call intgrid_deloc_brf(natt,xgatt,idatt,idg,itype,luw,idx,sidx)

    ! clean up YT checkpoint files
    if (itype == itype_yt) then
       call fclose(luw)
    endif
    deallocate(idx,sidx,psum,idatt,idg,xgatt)

  end subroutine intgrid_driver

  !> Output the results of an integration on a grid (YT or bader).
  subroutine intgrid_output(noatoms,doescher,doflist,nbasin,idatt,psum,pmask,xcoord,idx,sidx)
    use fields
    use global
    use struct
    use struct_basic
    use tools_io

    logical, intent(in) :: noatoms, doescher, doflist
    integer, intent(in) :: nbasin
    integer, intent(in) :: idatt(nbasin)
    real*8, intent(in) :: psum(nprops,nbasin)
    logical, intent(in) :: pmask(nprops)
    real*8, intent(in) :: xcoord(3,nbasin)
    integer, intent(in) :: idx(nbasin)
    character*3, intent(in) :: sidx(nbasin)

    character*3 :: sat
    real*8 :: auxs, tsum(2), asum(2), rpmask(nprops), rpmask2(nprops), dv(3)
    integer :: i, j, k, iaux, fid
    integer :: niprops, ntot, ntyp(100), nn, lu
    character(len=:), allocatable :: lbl
    character(len=:), allocatable :: oline, order
    real*8, allocatable :: flist(:,:)

    ! charge and volume summary
    if (.not.noatoms) then
       write (uout,'("   at   name     Z  mult    Volume(a.u.)        Num. elec.           Charge")')
       auxs = 0d0
       tsum = 0d0
       do i = 1, cr%nneq
          do k = 1, cr%ncel
             if (cr%atcel(k)%idx == i) exit
          end do
          if (cr%at(i)%zpsp > 0) then
             iaux = cr%at(i)%zpsp
          else
             iaux = cr%at(i)%z
          end if

          asum = 0d0
          do j = 1, nbasin
             if (idatt(j) == k) then
                asum(:) = asum(:) + psum(1:2,j)
             end if
          end do

          write (uout,'(99(A,X))') &
             string(i,length=4,justify=ioj_right), &
             string(cr%at(i)%name,length=10,justify=ioj_center), &
             string(iaux,length=3,justify=ioj_right), &
             string(cr%at(i)%mult,length=3,justify=ioj_right), &
             string(asum(1),'f',decimal=10,length=18,justify=6), &
             string(asum(2),'f',decimal=10,length=18,justify=6), &
             string(real(iaux,8)-asum(2),'f',decimal=10,length=18,justify=6)
          tsum(1) = tsum(1) + cr%at(i)%mult * asum(1)
          tsum(2) = tsum(2) + cr%at(i)%mult * asum(2)
          auxs = auxs + cr%at(i)%mult * (real(iaux,8) - asum(2))
       end do
       write (uout,'(80("-"))')
       write (uout,'(" Total    ",14X,3(A,1X)/)') &
          string(tsum(1),'f',decimal=10,length=18,justify=6), &
          string(tsum(2),'f',decimal=10,length=18,justify=6), &
          string(auxs,'f',decimal=10,length=18,justify=6)
    endif

    ! per-basin summary
    ! header
    write (uout,'("* List of basins and local properties (atoms not necessarily in input order)")')
    oline = " bas atom               Approx. pos              "
    niprops = 0
    do k = 1, nprops
       if (pmask(k)) then
          oline = oline // " " // string(integ_prop(k)%prop_name,length=19,justify=ioj_center)
       end if
    end do
    write (uout,'(A)') oline

    ! body
    do i = 1, nbasin
       j = idx(i)
       niprops = 0
       do k = 1, nprops
          if (pmask(k)) then
             niprops = niprops + 1
             rpmask(niprops) = psum(k,j)
          end if
       end do
       ! write (uout,'(I3,X,A3,X,3(F11.7,X),1p,999(E18.11,X))') i, sidx(i), xcoord(:,j), (rpmask(k),k=1,niprops)
       write (uout,'(999(A,X))') &
          string(i,length=4,justify=ioj_right), &
          string(sidx(i),length=5,justify=ioj_center), &
          (string(xcoord(k,j),'f',decimal=6,length=12,justify=4),k=1,3), &
          (string(rpmask(k),'e',decimal=10,length=18,justify=6),k=1,niprops)
    end do
    write (uout,'(80("-"))')

    ! total and cell
    rpmask = 0d0
    rpmask2 = 0d0
    niprops = 0
    do k = 1, nprops
       if (pmask(k)) then
          niprops = niprops + 1
          do i = 1, nbasin
             rpmask(niprops) = rpmask(niprops) + psum(k,i)
          end do
          if (integ_prop(k)%itype == itype_v) then
             rpmask2(niprops) = cr%omega
          else
             ntot = f(integ_prop(k)%fid)%n(1)*f(integ_prop(k)%fid)%n(2)*f(integ_prop(k)%fid)%n(3)
             rpmask2(niprops) = sum(f(integ_prop(k)%fid)%f) * cr%omega / ntot
          end if
       end if
    end do
    write (uout,'("Total integration",33X,1p,999(A,X))') &
          (string(rpmask(k),'e',decimal=10,length=18,justify=6),k=1,niprops)
    write (uout,'("Cell ",45X,1p,999(A,X))') &
          (string(rpmask2(k),'e',decimal=10,length=18,justify=6),k=1,niprops)
    write (uout,*)

    if (doescher) then
       ! write an escher file
       write(uout,'(" Writing the maxima escher file: ",A/)') trim(fileroot)//"_bader.m"
       order = trim(fileroot) // "_bader.m"
       call struct_write(order)

       ! count number of atoms per type
       ntyp = 0
       do i = 1, cr%ncel
          ntyp(cr%at(cr%atcel(i)%idx)%z) = ntyp(cr%at(cr%atcel(i)%idx)%z) + 1
       end do

       nn = count(idatt < 0)
       lu = fopen_append(order)
       write (lu,'("cr.nat += ",I6,";")') nn
       write (lu,'("cr.ntyp += 1;")')
       write (lu,'("cr.ztyp = [cr.ztyp 105];")')
       write (lu,'("cr.attyp = [cr.attyp,""n@""];")')
       lbl = "cr.typ = [cr.typ "
       do j = 1, nn
          lbl = lbl // " " // string(count(ntyp>0)+1)
       end do
       lbl = lbl // "];"
       write (lu,'(A)') lbl
       write (lu,'("cr.x = [cr.x")')
       do i = 1, nbasin
          if (idatt(i) > 0) cycle
          dv = xcoord(:,i)
          write (lu,'(2X,1p,3(E22.14,X))') dv
       end do
       write (lu,'("  ];")')
       call fclose(lu)
    end if

    ! write the value of all integrable fields on the maxima
    if (doflist) then
       ! header
       write (uout,'("* Integrable field values at the maxima")')
       oline = string("bas",3) // " " // string("at",3) // string("Approx. pos.",35,justify=ioj_center)
       niprops = 0
       do k = 1, nprops
          if (pmask(k).and.integ_prop(k)%itype/=itype_v) &
             oline = oline // " " // string(integ_prop(k)%prop_name,17,justify=ioj_center)
       end do
       write (uout,'(A)') oline

       ! allocate and calculate the field value list
       allocate(flist(nbasin,nprops))
       do k = 1, nprops
          if (integ_prop(k)%itype == itype_v .or..not.pmask(k)) cycle
          fid = integ_prop(k)%fid
          do i = 1, nbasin
             if (idatt(i) > 0) then
                flist(i,k) = grd0(f(fid),cr%atcel(idatt(i))%r)
             else
                dv = cr%x2c(xcoord(:,i))
                flist(i,k) = grd0(f(fid),dv)
             endif
          end do
       end do

       ! body
       do i = 1, nbasin
          j = idx(i)
          niprops = 0
          do k = 1, nprops
             if (pmask(k) .and. integ_prop(k)%itype /= itype_v) then
                niprops = niprops + 1
                rpmask(niprops) = flist(j,k)
             end if
          end do
          write (uout,'(I3,X,A3,X,3(F11.7,X),1p,999(E18.11,X))') i, sidx(i), xcoord(:,j), (rpmask(k),k=1,niprops)
       end do
       write (uout,'(X,80("-")/)')
       deallocate(flist)
    end if

  end subroutine intgrid_output

  !> Calculate and output the multipole moments using integration on a grid.
  subroutine intgrid_multipoles(natt,xgatt,idatt,idg,itype,luw,idx,sidx)
    use fields
    use struct_basic
    use global
    use types, only: gk
    use tools_math
    use tools_io

    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idatt(natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: itype
    integer, intent(in) :: luw
    integer, intent(in) :: idx(natt)
    character*3, intent(in) :: sidx(natt)

    character(len=:), allocatable :: lbl
    character*3 :: ls, ms
    integer :: i, j, k, l, m, n(3), nn, ntot, nl
    integer :: ix
    integer :: fid, lmax
    real*8, allocatable :: mpole(:,:)
    real(kind=gk), allocatable :: w(:,:,:)
    integer, allocatable :: l_(:), m_(:)
    real*8 :: dv(3), r, tp(2), p(3)
    real*8, allocatable :: rrlm(:)

    ! size of the grid
    n = f(refden)%n
    ntot = n(1)*n(2)*n(3)

    if (itype == itype_yt) then
       allocate(w(n(1),n(2),n(3)))
    endif
    ! run over all properties for which multipole calculation is active
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_mpoles) cycle
       fid = integ_prop(l)%fid
       write (uout,'("* Basin multipole moments (using real solid harmonics)")')
       write (uout,'("+ Reference field (number ",A,"): ",A)') string(refden), string(integ_prop(refden)%prop_name)
       write (uout,'("+ Integrated field (number ",A,"): ",A)') string(l), string(integ_prop(l)%prop_name)
       write (uout,'("+ The calculated multipoles are: ")')
       write (uout,'("    Q_lm^A = int_A rho(r) * Rlm(r) dr ")')
       write (uout,'("  where the integral is over the basin of A, and Rlm is a real solid harmonic.")')
       write (uout,'("  The coordinates are referred to the attractor of the A basin.")')
       write (uout,*)

       ! allocate space
       lmax = integ_prop(l)%lmax
       allocate(rrlm((lmax+1)*(lmax+1)))
       allocate(mpole((lmax+1)*(lmax+1),natt))

       ! calcualate the multipoles
       mpole = 0d0
       if (itype == itype_bader) then
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
                   mpole(:,ix) = mpole(:,ix) + rrlm * f(fid)%f(i,j,k)
                end do
             end do
          end do
       else
          rewind(luw)
          do m = 1, natt
             read(luw) w
             do i = 1, n(1)
                do j = 1, n(2)
                   do k = 1, n(3)
                      if (abs(w(i,j,k)) < 1d-15) cycle
                      p = real((/i-1,j-1,k-1/),8)
                      dv = p/real(n,8) - xgatt(:,m)
                      call cr%shortest(dv,r)
                      call tosphere(dv,r,tp)
                      call genrlm_real(lmax,r,tp,rrlm)
                      mpole(:,m) = mpole(:,m) + rrlm * f(fid)%f(i,j,k) * w(i,j,k)
                   end do
                end do
             end do
          end do
       endif
       mpole = mpole * cr%omega / ntot

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
          lbl = " n  at   "
          do j = (i-1)*5+1,min(5*i,(lmax+1)**2)
             if (j == 1) then
                lbl = lbl // " " // string("1",18,justify=ioj_center)
             elseif (j == 2) then
                lbl = lbl // " " // string("x",18,justify=ioj_center)
             elseif (j == 3) then
                lbl = lbl // " " // string("z",18,justify=ioj_center)
             elseif (j == 4) then
                lbl = lbl // " " // string("y",18,justify=ioj_center)
             elseif (j == 5) then
                lbl = lbl // " " // string("sq(3)/2(x^2-y^2)",18,justify=ioj_center)
             elseif (j == 6) then
                lbl = lbl // " " // string("sq(3)xz",18,justify=ioj_center)
             elseif (j == 7) then
                lbl = lbl // " " // string("(3z^2-r^2)/2",18,justify=ioj_center)
             elseif (j == 8) then
                lbl = lbl // " " // string("sq(3)yz",18,justify=ioj_center)
             elseif (j == 9) then
                lbl = lbl // " " // string("sq(3)xy",18,justify=ioj_center)
             else
                write (ls,'(I3)') l_(j)
                write (ms,'(I3)') abs(m_(j))
                if (m_(j) <= 0) then
                   lbl = lbl // " " // string("C("//string(ls)//","//string(ms)//")",18,justify=ioj_center)
                else
                   lbl = lbl // " " // string("S("//string(ls)//","//string(ms)//")",18,justify=ioj_center)
                end if
             endif
          end do
          write (uout,'(A)') lbl

          ! body
          do j = 1, natt
             write (uout,'(99(A,X))') string(j,length=4,justify=ioj_right),&
                string(sidx(j),length=5,justify=ioj_center),&
                (string(mpole((i-1)*5+1+k,idx(j)),'f',decimal=10,length=18,justify=6),k=0,min(4,size(mpole,1)-(i-1)*5-1))
          enddo
          if (i < nn) then
             write (uout,*)
          else
             write (uout,'(X,80("-")/)')
          end if
       end do
       deallocate(mpole,l_,m_,rrlm)
    end do
    if (itype == itype_yt) then
       deallocate(w)
    endif

  end subroutine intgrid_multipoles

  !> Calculate localization and delocalization indices using the
  !> basin assignment found by YT or BADER and a wfn scalar field.
  subroutine intgrid_deloc_wfn(natt,xgatt,idatt,idg,itype,luw,idx,sidx)
    use wfn_private
    use fields
    use struct_basic
    use types, only: gk
    use global
    use tools_io

    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idatt(natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: itype
    integer, intent(in) :: luw
    integer, intent(in) :: idx(natt)
    character*3, intent(in) :: sidx(natt)

    integer :: i, j, k, l, m, i1, i2
    integer :: fid, n(3), ntot, lumo, ix, nn
    real(kind=gk), allocatable :: w(:,:,:)
    real*8, allocatable :: xmo(:), sij(:,:,:)
    real*8, allocatable :: li(:), di(:,:)
    real*8 :: x(3), rho, auxg(3), auxh(3,3), gkin, vir
    character*3 :: sat1, sat2

    ! size of the grid
    n = f(refden)%n
    ntot = n(1)*n(2)*n(3)

    if (itype == itype_yt) then
       allocate(w(n(1),n(2),n(3)))
    endif
    ! run over all properties for which multipole calculation is active
    do l = 1, nprops
       ! only wfn fields for now
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_wfn) cycle
       write (uout,'("* Localization and delocalization indices")')
       write (uout,'("+ Reference field (number ",I2,"): ",A)') refden, trim(integ_prop(refden)%prop_name)
       write (uout,'("+ Integrated field (number ",I2,"): ",A)') l, trim(integ_prop(l)%prop_name)

       ! calculate the overlap matrix
       allocate(sij(f(fid)%nmo,f(fid)%nmo,natt),xmo(f(fid)%nmo))
       sij = 0d0
       xmo = 0d0
       if (itype == itype_bader) then
          !$omp parallel do private (i,j,k,x,rho,auxg,auxh,ix,i1,i2,gkin,vir) firstprivate(xmo) schedule(dynamic)
          do i = 1, n(1)
             do j = 1, n(2)
                do k = 1, n(3)
                   x = cr%x2c(real((/i-1,j-1,k-1/),8) / n)
                   call wfn_rho2(f(fid),x,0,rho,auxg,auxh,gkin,vir,xmo)
                   ix = idg(i,j,k)
                   do i1 = 1, f(fid)%nmo
                      do i2 = i1, f(fid)%nmo
                         !$omp critical (sijwrite)
                         sij(i1,i2,ix) = sij(i1,i2,ix) + xmo(i1) * xmo(i2)
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
                   call wfn_rho2(f(fid),x,0,rho,auxg,auxh,gkin,vir,xmo)
                   ix = idg(i,j,k)
                   write (lumo) xmo
                end do
             end do
          end do

          ! calculate the atomic overlap matrix and kill the file
          rewind(luw)
          do m = 1, natt
             read(luw) w
             rewind(lumo)
             do i = 1, n(1)
                do j = 1, n(2)
                   do k = 1, n(3)
                      read(lumo) xmo
                      do i1 = 1, f(fid)%nmo
                         do i2 = i1, f(fid)%nmo
                            sij(i1,i2,m) = sij(i1,i2,m) + xmo(i1) * xmo(i2) * w(i,j,k)
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
       sij = sij * cr%omega / ntot
       do ix = 1, natt
          do i1 = 1, f(fid)%nmo
             do i2 = i1, f(fid)%nmo
                sij(i2,i1,ix) = sij(i1,i2,ix)
             end do
          end do
       end do

       ! calculate localization and delocalization indices
       allocate(li(natt),di(natt,natt))
       if (f(fid)%wfntyp == wfn_rhf) then
          do i = 1, natt
             li(i) = 2d0 * sum(sij(:,:,i)*sij(:,:,i))
             do j = i+1, natt
                di(i,j) = 4d0 * sum(sij(:,:,i)*sij(:,:,j))
                di(j,i) = di(i,j)
             end do
          end do
       else
          do i = 1, natt
             li(i) = 0.5d0 * sum(sij(:,:,i)*sij(:,:,i))
             do j = i, natt
                di(i,j) = sum(sij(:,:,i)*sij(:,:,j))
                di(j,i) = di(i,j)
             end do
          end do
       end if

       ! output the lambdas
       write (uout,'("+ Localization indices (lambda)")')
       do j = 1, natt
          write (uout,'(I2,X,A3,5(F16.10,X))') j, sidx(j), li(idx(j))
       enddo
       write (uout,*)

       write (uout,'("+ Delocalization indices (lambda)")')
       do i = 1, natt
          do j = i+1, natt
             write (uout,'(I2,X,A3,X,I2,X,A3,5(F16.10,X))') i, sidx(i), j, sidx(j), di(idx(i),idx(j))
          end do
       end do
       write (uout,*)

       ! wrap up
       deallocate(sij,di,li)
    end do

    if (itype == itype_yt) then
       deallocate(w)
    endif

  end subroutine intgrid_deloc_wfn

  !> Calculate localization and delocalization indices using the
  !> basin assignment found by YT or BADER and a grid-related dat file.
  subroutine intgrid_deloc_brf(natt,xgatt,idatt,idg,itype,luw,idx,sidx)
    use bader
    use wfn_private
    use fields
    use struct_basic
    use global
    use types
    use tools_io
    use tools_math
    use tools

    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idatt(natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: itype
    integer, intent(in) :: luw
    integer, intent(in) :: idx(natt)
    character*3, intent(in) :: sidx(natt)
    character(len=:), allocatable :: aux

    integer :: iaa, jaa, kaa
    integer :: ia, ja, ka, iba, ib, jb, kb, ibb
    integer :: i, j, k, l, m, in, jn, kn, ia1, ia2, ia3
    integer :: i0, i1, j0, j1, k0, k1, m1, m2, m3, n0(3)
    integer :: fid, n(3), lumo, ix, nn, kmo, lmo
    integer :: rat(3), idx1(3), idx2(3), idxw(3)
    real(kind=gk), allocatable :: w(:,:,:)
    integer :: nbnd, nnr, nlat, ilat, jlat, nmo, imo, jmo
    real*8, allocatable :: psic(:,:,:), psic2(:,:,:)
    real*8, allocatable :: sij(:,:,:)
    real*8, allocatable :: fa(:,:)
    real*8 :: li, asum, asum2, x(3), fval1, fval2
    logical :: found
    integer :: natt1
    real*8 :: x0(3,3), x0inv(3,3), r1(3), r2(3)
    real*8, allocatable :: xgatt1(:,:), dist(:)
    integer, allocatable :: idg1(:,:,:), idatt1(:), lvec(:,:)
    integer, allocatable :: io(:)
    type(crystal) :: cr1
    type(field) :: f1

    ! size of the grid
    n = f(refden)%n

    if (itype == itype_yt) then
       allocate(w(n(1),n(2),n(3)))
    endif
    ! run over all properties for which multipole calculation is active
    do l = 1, nprops
       ! check that it is a brf file
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_deloc) cycle
       fid = integ_prop(l)%fid
       if (f(fid)%type /= type_grid .or..not.allocated(f(fid)%fwan)) cycle
       write (uout,'("* Localization and delocalization indices")')
       write (uout,'("+ Reference field (number ",A,"): ",A)') string(refden), string(integ_prop(refden)%prop_name)
       write (uout,'("+ Integrated field (number ",A,"): ",A)') string(l), string(integ_prop(l)%prop_name)

       ! check -- is the density correct?
       if (refden == 1 .and. fid == 2) then
          write (uout,'("+ Density check +")')
          write (uout,'("Maximum rho diff: ",A)') string(maxval(f(fid)%f - f(refden)%f),'e')
          write (uout,'("Minimum rho diff: ",A)') string(minval(f(fid)%f - f(refden)%f),'e')
          write (uout,'("Average rho diff: ",A)') string(sum(f(fid)%f - f(refden)%f) / (f(fid)%n(1)*f(fid)%n(2)*f(fid)%n(3)),'e')
          write (uout,'("Integral 1: ",A)') string(sum(f(fid)%f) * cr%omega / (f(fid)%n(1)*f(fid)%n(2)*f(fid)%n(3)),'e')
          write (uout,'("Integral 2: ",A)') string(sum(f(refden)%f) * cr%omega / (f(refden)%n(1)*f(refden)%n(2)*f(refden)%n(3)),'e')
       else
          write (uout,*) "rho needs to be field 1 and xsf field 2"
          stop 1
       endif

       ! assign values to some integers
       nbnd = size(f(fid)%fwan,4)
       nlat = f(fid)%nwan(1)*f(fid)%nwan(2)*f(fid)%nwan(3)
       nmo = nlat * nbnd
       i0 = 0
       i1 = f(fid)%nwan(1)-1
       j0 = 0
       j1 = f(fid)%nwan(2)-1
       k0 = 0
       k1 = f(fid)%nwan(3)-1

       ! ! check -- are the wannier functions orthonormal?
       ! if (allocated(f(fid)%fwan)) then
       !    do imo = 1, nmo
       !       call unpackidx(imo,ia,ja,ka,iba)
       !       write (uout,'(99(A))') "+ Overlaps of: b=", string(iba), " Rx=", string(ia), &
       !          " Ry=", string(ja), " Rz=", string(ka)
       !       do jmo = 1, nmo
       !          call unpackidx(jmo,ib,jb,kb,ibb)
       !          asum = 0d0
       !          do m1 = 0, f(fid)%nwan(1)-1
       !             do m2 = 0, f(fid)%nwan(2)-1
       !                do m3 = 0, f(fid)%nwan(3)-1
       !                   idx1 = (/m1,m2,m3/) + (/ia,ja,ka/)
       !                   idx1 = modulo(idx1,f(fid)%nwan)
       !                   idx2 = (/m1,m2,m3/) + (/ib,jb,kb/)
       !                   idx2 = modulo(idx2,f(fid)%nwan)
       ! 
       !                   asum = asum + sum(f(fid)%fwan(idx1(1)*n(1)+1:idx1(1)*n(1)+n(1),idx1(2)*n(2)+1:idx1(2)*n(2)+n(2),idx1(3)*n(3)+1:idx1(3)*n(3)+n(3),iba) * &
       !                      f(fid)%fwan(idx2(1)*n(1)+1:idx2(1)*n(1)+n(1),idx2(2)*n(2)+1:idx2(2)*n(2)+n(2),idx2(3)*n(3)+1:idx2(3)*n(3)+n(3),ibb))
       !                end do
       !             end do
       !          end do
       !          asum = asum / (n(1)*n(2)*n(3))
       !          write (uout,'(99(A))') "b=", string(ibb), " Rx=", string(ib), &
       !             " Ry=", string(jb), " Rz=", string(kb), " ", string(asum,'e')
       !       end do
       !    end do
       ! else
       !    write (uout,*) "field 2 does not have wannier functions"
       !    stop 1
       ! end if

       ! build the supercell
       cr1 = cr
       x0 = 0d0
       do i = 1, 3
          x0(i,i) = real(f(fid)%nwan(i),8)
       end do
       call cr1%newcell(x0,.false.)

       ! build the field for the supercell
       f1%type = type_grid
       f1%usecore = .false.
       n0 = f(fid)%n * f(fid)%nwan
       f1%n = n0
       f1%mode = f(refden)%mode
       allocate(f1%f(n0(1),n0(2),n0(3)))
       do i = 0, f(fid)%nwan(1)-1
          do j = 0, f(fid)%nwan(2)-1
             do k = 0, f(fid)%nwan(3)-1
                f1%f(i*n(1)+1:(i+1)*n(1),j*n(2)+1:(j+1)*n(2),k*n(3)+1:(k+1)*n(3)) = f(refden)%f
             end do
          end do
       end do
       f1%c2x = cr1%car2crys
       f1%x2c = cr1%crys2car
       f1%init = .true.

       ! run bader integration on the supercell !xxxx! ratom
       call bader_integrate(cr1,f1,natt1,xgatt1,idatt1,idg1,1d40)

       ! check -- list of attractors with identity, volumes, and charges
       ! natt1 -> number of attractors
       ! xgatt1 -> cryst coords of the attractor
       ! idatt1 -> id of the basin (>0 atom, <0 nnm)
       ! idg1 -> array the size of the lattice with the id of the basin (1 to natt1)
       allocate(w(n0(1),n0(2),n0(3)))
       write (uout,'("+ List of basins and charges")')
       do i = 1, natt1
          r1 = xgatt1(:,i) * n0 / real(n,8)
          r2 = floor(r1)
          r1 = r1 - r2
          if (idatt1(i) > 0) then
             aux = cr1%at(cr1%atcel(idatt1(i))%idx)%name
             r1 = cr1%atcel(idatt1(i))%x * n0 / real(n,8)
             r1 = r1 - floor(r1)
          else
             aux = "nnm" // string(idatt1(i))
          endif
          w = 0d0
          where(idg1 == i)
             w = f1%f
          end where
          asum = sum(w) * cr1%omega / (n0(1)*n0(2)*n0(3))
          asum2 = count(idg1 == i) * cr1%omega / (n0(1)*n0(2)*n0(3))
          write (uout,'(99A)') string(i), " (id=", string(idatt1(i)), ") x= ", &
             (string(r1(j),'e',12,4),j=1,3), "atom = ", string(aux), " ", (string(nint(r2(j)),3),j=1,3), &
             " q = ", string(asum,'e',15,6), " v = ", string(asum2,'e',15,6)
       end do

       ! calculate the overlap matrix
       allocate(sij(nmo,nmo,natt1))
       sij = 0d0
       
       if (itype == itype_bader) then
          write (*,*) "+ Calculating overlap matrices..."
          allocate(psic(n0(1),n0(2),n0(3)),psic2(n0(1),n0(2),n0(3)))
          do imo = 1, nmo
             write (*,*) "imo", imo
             call unpackidx(imo,ia,ja,ka,iba)
             ! prepare the supergrid
             do m1 = 0, f(fid)%nwan(1)-1
                do m2 = 0, f(fid)%nwan(2)-1
                   do m3 = 0, f(fid)%nwan(3)-1
                      idx1 = (/m1,m2,m3/) + (/ia,ja,ka/)
                      idx1 = modulo(idx1,f(fid)%nwan)
                      psic(m1*n(1)+1:m1*n(1)+n(1),m2*n(2)+1:m2*n(2)+n(2),m3*n(3)+1:m3*n(3)+n(3)) = &
                         f(fid)%fwan(idx1(1)*n(1)+1:idx1(1)*n(1)+n(1),idx1(2)*n(2)+1:idx1(2)*n(2)+n(2),idx1(3)*n(3)+1:idx1(3)*n(3)+n(3),iba)
                   end do
                end do
             end do
             do jmo = imo, nmo
                write (*,*) "jmo", jmo
                call unpackidx(jmo,ib,jb,kb,ibb)
                ! prepare the supergrid
                do m1 = 0, f(fid)%nwan(1)-1
                   do m2 = 0, f(fid)%nwan(2)-1
                      do m3 = 0, f(fid)%nwan(3)-1
                         idx2 = (/m1,m2,m3/) + (/ib,jb,kb/)
                         idx2 = modulo(idx2,f(fid)%nwan)
                         psic2(m1*n(1)+1:m1*n(1)+n(1),m2*n(2)+1:m2*n(2)+n(2),m3*n(3)+1:m3*n(3)+n(3)) = &
                            psic(m1*n(1)+1:m1*n(1)+n(1),m2*n(2)+1:m2*n(2)+n(2),m3*n(3)+1:m3*n(3)+n(3)) * &
                            f(fid)%fwan(idx2(1)*n(1)+1:idx2(1)*n(1)+n(1),idx2(2)*n(2)+1:idx2(2)*n(2)+n(2),idx2(3)*n(3)+1:idx2(3)*n(3)+n(3),ibb)
                      end do
                   end do
                end do
                do i = 1, natt1
                   w = 0d0
                   where(idg1 == i)
                      w = psic2
                   end where
                   sij(imo,jmo,i) = sum(w)
                   sij(jmo,imo,i) = sij(imo,jmo,i)
                end do
             end do
          end do
       else
          write (uout,*) "Not implemented with YT"
          stop 1
       end if

       ! scale (the omega comes from wannier)
       sij = sij / (n(1)*n(2)*n(3))

       ! non-spin-polarized
       sij = 2d0 * sij

       write (uout,'("+ Charge check using the overlap matrices")')
       asum2 = 0d0
       do i = 1, natt1
          asum = 0d0
          do imo = 1, nmo
             asum = asum + sij(imo,imo,i)
          end do
          write (uout,'(" N(A) -- sum_Rn wRn^2 for atom ",I2,X,F12.6)') i, asum
          asum2 = asum2 + asum
       end do
       write (uout,'(" N(total) -- sum_A sum_Rn wRn^2 ",F12.6)') asum2
       
       write (uout,'("+ Orthonormality check")')
       do imo = 1, nmo
          do jmo = 1, nmo
             asum = sum(sij(imo,jmo,:))
             write (uout,'("Orb. ",I2,X,I2,X,F12.6)') imo, jmo, asum
          end do
       end do

       ! calculate localization and delocalization indices
       allocate(fa(natt1,natt1))
       fa = 0d0
       do i = 1, natt1
          do j = i, natt1
             fa(i,j) = sum(sij(:,:,i) * sij(:,:,j))
             fa(j,i) = fa(i,j)
          end do
       end do

       ! localization indices
       write (uout,'("+ Localization indices (lambda)")')
       do i = 1, natt1
          r1 = xgatt1(:,i) * n0 / real(n,8)
          r2 = floor(r1)
          r1 = r1 - r2
          if (idatt1(i) > 0) then
             aux = cr1%at(cr1%atcel(idatt1(i))%idx)%name
             r1 = cr1%atcel(idatt1(i))%x * n0 / real(n,8)
             r1 = r1 - floor(r1)
          else
             aux = "nnm" // string(idatt1(i))
          endif

          asum = 0d0
          do imo = 1, nmo
             asum = asum + sij(imo,imo,i)
          end do
          write (uout,'(99A)') string(i), " (id=", string(idatt1(i)), ") x= ", &
             (string(r1(j),'e',12,4),j=1,3), "atom = ", string(aux), " ", (string(nint(r2(j)),3),j=1,3), &
             " li = ", string(0.5d0*fa(i,i),'e',15,6), " q(li) = ", string(0.5d0*sum(fa(i,:)),'e',15,6), &
             " q = ", string(asum,'e',15,6)
       end do
       write (uout,*)
       
       write (uout,'("+ Delocalization indices (delta)")')
       allocate(dist(natt1),io(natt1))
       do i = 1, natt1
          r1 = xgatt1(:,i) * n0 / real(n,8)
          r2 = floor(r1)
          if (.not.all(nint(r2) == 0)) cycle
          r1 = r1 - r2
          if (idatt1(i) > 0) then
             aux = cr1%at(cr1%atcel(idatt1(i))%idx)%name
             r1 = cr1%atcel(idatt1(i))%x * n0 / real(n,8)
             r1 = r1 - floor(r1)
          else
             aux = "nnm" // string(idatt1(i))
          endif

          write (uout,'(99A)') string(i), " (id=", string(idatt1(i)), ") x= ", &
             (string(r1(j),'e',12,4),j=1,3), "atom = ", string(aux), " ", (string(nint(r2(j)),3),j=1,3)

          dist = 0d0
          do j = 1, natt1
             io(j) = j
             if (j == i) cycle
             dist(j) = cr1%eql_distance(xgatt1(:,i),xgatt1(:,j))
          end do
          
          call qcksort(dist,io,1,natt1)
          do m = 2, natt1
             j = io(m)
             r1 = xgatt1(:,j) * n0 / real(n,8)
             r2 = floor(r1)
             r1 = r1 - r2
             if (idatt1(j) > 0) then
                aux = cr1%at(cr1%atcel(idatt1(j))%idx)%name
                r1 = cr1%atcel(idatt1(j))%x * n0 / real(n,8)
                r1 = r1 - floor(r1)
             else
                aux = "nnm" // string(idatt1(j))
             endif
             write (uout,'(99A)') "+ ", string(j), " (id=", string(idatt1(j)), ") x= ", &
                (string(r1(k),'e',12,4),k=1,3), "atom = ", string(aux), " ", (string(nint(r2(k)),3),k=1,3),&
                "dist = ", string(dist(j),'f'), " di = ", string(fa(i,j),'f')
          end do
       end do
       write (uout,*)

       ! wrap up
       deallocate(sij,fa)
    end do

    if (itype == itype_yt) then
       deallocate(w)
    endif

  contains

    subroutine unpackidx(idx,io,jo,ko,bo)
      integer, intent(in) :: idx
      integer, intent(out) :: io, jo, ko, bo

      integer :: iaux
      integer :: inn, jnn, knn

      inn = i1-i0+1
      jnn = j1-j0+1
      knn = k1-k0+1

      ! unpack
      bo = modulo(idx-1,nbnd) + 1
      iaux = (idx-1 - (bo-1)) / nbnd
      ko = modulo(iaux,knn) + 1
      iaux = (iaux - (ko-1)) / knn
      jo = modulo(iaux,jnn) + 1
      iaux = (iaux - (jo-1)) / jnn
      io = modulo(iaux,inn) + 1

      ! translate
      io = io - 1 + i0
      jo = jo - 1 + j0
      ko = ko - 1 + k0

    end subroutine unpackidx

    subroutine packidx(io,jo,ko,bo,idx)
      integer, intent(in) :: io, jo, ko, bo
      integer, intent(out) :: idx

      integer :: zio, zjo, zko
      integer :: inn, jnn, knn

      inn = i1-i0+1
      jnn = j1-j0+1
      knn = k1-k0+1
      zio = modulo(io-i0,inn) + 1
      zjo = modulo(jo-j0,jnn) + 1
      zko = modulo(ko-k0,knn) + 1

      ! translate and pack
      idx = (bo-1) + nbnd * ((zko-1) + knn * ((zjo-1) + jnn * (zio-1))) + 1

    end subroutine packidx

  end subroutine intgrid_deloc_brf

  !> Do a radial numerical quadrature on the given ray, between the
  !> selected radii and return the properties. The r^2 factor is
  !> included.
  subroutine int_radialquad(x,theta,phi,r0,rend,lprop,abserr,neval,iaserr,ierr)
    use quadpack
    use varbas
    use global
    use fields
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
    use surface
    use tools_math, only: gauleg
    use types
    use param

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
    use surface
    use fields
    use tools_math, only: gauleg
    use tools_io, only: ferror, warning, uout
    use types
    use param

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
    use surface
    use tools_math, only: select_lebedev
    use types

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
    use surface
    use fields
    use tools_math, only: select_lebedev
    use tools_io, only: uout, ferror, warning
    use types

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

  !> Output the integrated atomic properties. id is the attractor id
  !> from the non-equivalent CP list. If id==0, then the properties
  !> correspond to the unit cell. atprop is the properties
  !> vector. Only valid for the appropriate scalar field.
  subroutine ppty_output(id,atprop,maskprop)
    use fields
    use varbas
    use global
    use struct_basic
    use tools_io

    integer, intent(in) :: id
    real*8, intent(in) :: atprop(Nprops)
    logical, intent(in) :: maskprop(Nprops)

    character(len=:), allocatable :: saux
    real*8 :: totaln
    integer :: i, j, zaux

    totaln = 0d0
    do i = 1, cr%nneq
       if (cr%at(i)%zpsp > 0) then
          zaux = cr%at(i)%zpsp
       else
          zaux = cr%at(i)%z
       end if
       totaln = totaln + zaux * cr%at(i)%mult
    end do

    if (id == 0) then
       write (uout,'("* Cell properties ")')
    else
       if (id <= cr%nneq) then
          saux = cr%at(id)%name
          if (cr%at(id)%zpsp > 0) then
             zaux = cr%at(id)%zpsp
          else
             zaux = cr%at(id)%z
          end if
       else
          saux = "NNM"
          zaux = 0
       end if
    end if

    if (id > 0) then
       write (uout,'("+ Integrated properties for attractor: ",A," (",A,")")') &
          string(id), string(saux)
       write (uout,'("  Pos.: ",3(A,X))') (string(cp(id)%x(j),'f',decimal=6),j=1,3)
       write (uout,'("  Z: ",A)'), string(zaux)
       write (uout,'("  Cell multiplicity: ",A)') string(cp(id)%mult)
    end if
    !
    if (maskprop(1)) then
       write (uout,'("  Volume (<1>): ",A)') string(atprop(1),'f',decimal=10)
       if (id == 0) then
          write (uout,'("  Cell volume: ",A)') string(cr%omega,'f',decimal=10)
          write (uout,'("  Cell volume error: ",A)') &
             string(cr%omega - atprop(1),'e',decimal=10)
          write (uout,'("  Rel. cell volume error: ",A)') &
             string((cr%omega - atprop(1)) / cr%omega,'e',decimal=10)
       end if
    end if
    !
    if (maskprop(2)) then
       write (uout,'("  Average population (<2>): ",A)') string(atprop(2),'f',decimal=10)
       if (id > 0) then
          write (uout,'("  Charge (Z - <2>): ",A)') string(zaux - atprop(2),'f',decimal=10)
       else
          write (uout,'("  Cell number of electrons: ",A)') string(totaln,'f',decimal=10)
          write (uout,'("  Cell charge error: ",A)') string(totaln - atprop(2),'e',decimal=10)
          write (uout,'("  Rel. cell charge error: ",A)') string((totaln - atprop(2)) / totaln,'e',decimal=10)
       end if
    end if
    do i = 3, nprops
       if (.not.integ_prop(i)%used) cycle
       if (id > 0) then
          write (uout,'("  Property (<",A,":",A,">): ",A)') string(i), &
             string(integ_prop(i)%prop_name), string(atprop(i),'e',decimal=10)
       else
          write (uout,'("  Cell property (<",A,":",A,">) : ",A)') string(i), &
             string(integ_prop(i)%prop_name), string(atprop(i),'e',decimal=10)
       end if
    end do
    write (uout,*)

  end subroutine ppty_output

  !> Output routine for all integration methods
  subroutine int_output(pmask,reason,nattr,icp,xattr,aprop,usesym,di,mpole)
    use fields
    use struct_basic
    use global
    use varbas
    use types
    use tools_io

    logical, intent(in) :: pmask(nprops)
    character*(*), intent(in) :: reason(nprops)
    integer, intent(in) :: nattr
    integer, intent(in) :: icp(nattr)
    real*8, intent(in) :: xattr(3,nattr)
    real*8, intent(in) :: aprop(nprops,nattr)
    logical, intent(in) :: usesym
    real*8, intent(in), optional :: di(nattr,nattr)
    real*8, intent(in), optional :: mpole(:,:)

    integer :: i, j, ip, ipmax, iplast
    integer :: fid, idx, nacprop(5)
    real*8 :: x(3), sump(nprops), xmult
    character(len=:), allocatable :: saux, itaux, label, cini
    character(len=:), allocatable :: sncp, scp, sname, sz, smult

    ! List of integrable properties and why they were rejected
    write (uout,'("* List of properties integrated in the attractor basins")')
    write (uout,'("+ The ""Label"" entries will be used to identify the integrable")')
    write (uout,'("  in the tables of integrated atomic properties. The entries with")')
    write (uout,'("  an ""x"" will not be integrated.")')
    write (uout,'("# id    Label    fid  Field       Additional")')
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
          string(i,2,ioj_left), &
          string(label,12,ioj_center), &
          string(itaux,2,ioj_right), &
          string(itype_names(integ_prop(i)%itype),12,ioj_left),&
          string(saux)
    end do
    write (uout,*)

    ! List of attractors
    write (uout,'("* List of attractors integrated")')
    if (.not.cr%ismolecule) then
       write (uout,'("# Id cp   ncp   Name  Z   mult           Position (cryst.) ")')
    else
       write (uout,'("# Id cp   ncp   Name  Z   mult           Position (",A,") ")') iunitname0(iunit)
    endif
    do i = 1, nattr
       if (icp(i) > 0) then
          ! this is a cp
          idx = cpcel(icp(i))%idx
          scp = string(icp(i),4,ioj_left)
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
       if (.not.cr%ismolecule) then
          x = xattr(:,i)
       else
          x = (cr%x2c(xattr(:,i)) + cr%molx0) * dunit
       endif
       write (uout,'(2X,99(A,X))') & 
          string(i,2,ioj_left), scp, sncp, sname, sz, &
          smult, (string(x(j),'f',12,7,4),j=1,3)
    end do
    write (uout,*)

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
       write (uout,'("# Id cp   ncp   Name  Z   mult ",5(A,X))') &
          (string(integ_prop(nacprop(j))%prop_name,15,ioj_center),j=1,ipmax)

       ! Table rows
       sump = 0d0
       do i = 1, nattr
          if (icp(i) > 0) then
             ! this is a cp
             idx = cpcel(icp(i))%idx
             scp = string(icp(i),4,ioj_left)
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
             string(i,2,ioj_left), scp, sncp, sname, sz, smult, &
             (string(aprop(nacprop(j),i),'e',15,8,4),j=1,ipmax)
       end do
       write (uout,'(24("-",99(A,X)))') ("---------------",j=1,ipmax)
       write (uout,'(2X,"Sum                         ",99(A,X))') &
          (string(sump(nacprop(j)),'e',15,8,4),j=1,ipmax)
       write (uout,*)
    end do
    
  end subroutine int_output

end module integration
