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

!> Basic plotting capabilities: contour diagrams, 1d, 2d and 3d representations.
submodule (rhoplot) proc
  implicit none

  !xx! private procedures
  ! subroutine contour(ff,r0,r1,r2,nx,ny,niso,ziso,rootname,dognu,dolabels)
  ! subroutine relief(rootname,outfile,zmin,zmax)
  ! subroutine colormap(rootname,outfile,cmopt,r0,r1,r2,dolabels)
  ! subroutine hallarpuntos(ff,zc,x,y,nx,ny)
  ! subroutine ordenarpuntos (luw,calpha,ziso)
  ! subroutine linea (x,y,n,luw,ziso)
  ! subroutine plotvec(r0,r1,r2,autocheck,udat)
  ! subroutine plotvec(r0,r1,r2,autocheck,udat)
  ! subroutine autochk(rp0)
  ! subroutine write_fichlabel(rootname)
  ! subroutine write_fichgnu(rootname,dolabels,docontour,dograds)

  ! contour lines common variables.
  integer, parameter :: nptcorte = 20000
  integer :: npuntos
  real*8  :: xc(nptcorte), yc(nptcorte), pic(nptcorte), pjc(nptcorte)

  ! gradient field walk paths database.
  integer, parameter :: MORIG = 1000
  integer, parameter :: mncritp = 1000
  integer :: norig, grpproj
  integer :: grpatr(MORIG), grpup(MORIG), grpdwn(MORIG)
  real*8  :: grphcutoff
  real*8  :: grpx(3,MORIG)
  real*8  :: grpcpeps
  real*8  :: r01, r02, a012, sinalfa, cosalfa, cotgalfa
  real*8  :: rp0(3), rp1(3), rp2(3), rp01(3), rp02(3), rpn(4)
  real*8  :: amat(3,3), bmat(3,3)
  real*8  :: outcpx, outcpy
  integer :: newncriticp
  integer :: newtypcrit(mncritp)
  real*8  :: newcriticp(3,mncritp)
  integer :: cpup(mncritp), cpdn(mncritp)
  integer :: indmax
  logical :: isneg

  ! label height
  real*8 :: RHOP_Hmax = 1d-1

  real*8, parameter  :: epsdis = 1d-4 !< distance cutoff (bohr) for in-plane
  real*8, parameter  :: epsf = 1d-5 !< distance cutoff (cryst) for strict tests
  integer, parameter :: RHOP_Mstep = 4000 !< max. number of gp steps

contains

  ! Calculate properties at a point or at a set of points given by the user.
  module subroutine rhoplot_point(line)
    use systemmod, only: sy
    use global, only: eval_next, dunit0, iunit
    use arithmetic, only: eval
    use tools_io, only: ferror, faterr, lgetword, equal, getword, &
       isexpression_or_word, uout, string, fopen_read, fclose, getline
    use types, only: scalar_value, realloc
    use param, only: bohrtoa
    character*(*), intent(in) :: line

    type(scalar_value) :: res
    logical :: ok, doall
    integer :: lp, lp2, i, j, ifi, lu
    real*8 :: x0(3), xx(3), rdum, x1, y1, z1
    character(len=:), allocatable :: word, expr, str, errmsg
    integer :: np
    real*8, allocatable :: xp(:,:)

    ! check if the next field is a file
    lp = 1
    word = getword(line,lp)
    if (len_trim(word) == 0) goto 999
    inquire(file=word,exist=ok)

    if (.not.ok) then
       ! read the point
       np = 1
       allocate(xp(3,1))
       lp = 1
       ok = eval_next(xp(1,1),line,lp)
       ok = ok .and. eval_next(xp(2,1),line,lp)
       ok = ok .and. eval_next(xp(3,1),line,lp)
       if (.not. ok) goto 999
    else
       np = 0
       allocate(xp(3,10))
       lu = fopen_read(word)
       do while (getline(lu,str))
          read (str,*,err=998) x1, y1, z1
          np = np + 1
          if (np > size(xp,2)) call realloc(xp,3,2*np)
          xp(1,np) = x1
          xp(2,np) = y1
          xp(3,np) = z1
       end do
       call realloc(xp,3,np)
       call fclose(lu)
    end if

    ! read additional options
    ifi = sy%iref
    doall = .false.
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'all')) then
          doall = .true.
       elseif (equal(word,'field')) then
          lp2 = lp
          word = getword(line,lp)
          ifi = sy%fieldname_to_idx(word)
          if (ifi < 0) then
             lp = lp2
             ok = isexpression_or_word(expr,line,lp)
             if (.not.ok) then
                call ferror('rhoplot_point','wrong FIELD in POINT',faterr,line,syntax=.true.)
                return
             end if
          else
             if (.not.sy%goodfield(ifi)) then
                call ferror('rhoplot_point','field not allocated',faterr,line,syntax=.true.)
                return
             end if
          end if
       elseif(len_trim(word) < 1) then
          exit
       else
          call ferror('rhoplot_point','Unknown keyword in POINT',faterr,line,syntax=.true.)
          return
       end if
    end do

    do i = 1, np
       x0 = xp(:,i)
       write (uout,'("* POINT ",3(A,"  "))') (string(x0(j),'f',decimal=7),j=1,3)
       if (.not.sy%c%ismolecule) then
          xx = sy%c%x2c(x0)
          write (uout,'("  Coordinates (bohr): ",3(A,"  "))') (string(xx(j),'f',decimal=7),j=1,3)
          write (uout,'("  Coordinates (ang): ",3(A,"  "))') (string(xx(j)*bohrtoa,'f',decimal=7),j=1,3)
       else
          write (uout,'("  Coordinates (ang): ",3(A,"  "))') (string(x0(j),'f',decimal=7),j=1,3)
          xx = x0 / dunit0(iunit) - sy%c%molx0
          x0 = sy%c%c2x(xx)
       endif
       if (ifi > -1) then
          write (uout,'("+ Field: ",A)') string(ifi)
          call sy%propty(ifi,x0,res,.false.,.true.,doall)
       else
          rdum = sy%eval(expr,errmsg,xx)
          if (len_trim(errmsg) > 0) then
             call ferror('critic','error in expression: '//trim(expr),faterr,line,syntax=.true.)
             return
          end if
          write (uout,'("  Expression (",A,"): ",A)') string(expr), string(rdum,'e',decimal=9)
       endif
       write (uout,*)
    end do

    return

998 continue
    call ferror('critic','error reading file in POINT keyword',faterr,line,syntax=.true.)
    return
999 continue
    call ferror('critic','wrong syntax in POINT keyword',faterr,line,syntax=.true.)
    return

  end subroutine rhoplot_point

  ! Calculate properties on a line
  module subroutine rhoplot_line(line)
    use systemmod, only: sy
    use global, only: eval_next, dunit0, iunit
    use arithmetic, only: eval
    use tools_io, only: ferror, faterr, lgetword, equal, getword, equal,&
       isexpression_or_word, fopen_write, uout, string, fclose
    use types, only: scalar_value, vstring
    character*(*), intent(in) :: line

    integer :: lp, lp2, nti, id, luout, np
    real*8 :: x0(3), x1(3), xp(3), dist, rhopt, lappt, xout(3)
    character(len=:), allocatable :: word, outfile, prop, expr
    type(scalar_value) :: res
    logical :: ok, iok
    integer :: i, j
    real*8, allocatable :: rhoout(:), lapout(:)
    type(vstring) :: lerrmsg

    ! read the points
    lp = 1
    ok = eval_next(x0(1),line,lp)
    ok = ok .and. eval_next(x0(2),line,lp)
    ok = ok .and. eval_next(x0(3),line,lp)
    ok = ok .and. eval_next(x1(1),line,lp)
    ok = ok .and. eval_next(x1(2),line,lp)
    ok = ok .and. eval_next(x1(3),line,lp)
    ok = ok .and. eval_next(np,line,lp)
    if (.not. ok) then
       call ferror('rhoplot_line','wrong syntax in LINE',faterr,line,syntax=.true.)
       return
    end if

    ! at least two points
    np = max(np,2)

    ! read additional options
    nti = 0
    prop = ""
    id = sy%iref
    outfile = ""
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'file')) then
          outfile = getword(line,lp)
          if (len_trim(outfile) < 1) then
             call ferror('rhoplot_line','file name not found',faterr,line,syntax=.true.)
             return
          end if
       else if (equal(word,'field')) then
          lp2 = lp
          word = getword(line,lp)
          id = sy%fieldname_to_idx(word)
          if (id < 0) then
             lp = lp2
             ok = isexpression_or_word(expr,line,lp)
             if (.not.ok) then
                call ferror('rhoplot_line','wrong FIELD in LINE',faterr,line,syntax=.true.)
                return
             end if
          else
             if (.not.sy%goodfield(id)) then
                call ferror('rhoplot_line','field not allocated',faterr,line,syntax=.true.)
                return
             end if
          end if
       elseif (equal(word,'gx')) then
          nti = 1
          prop = word
       else if (equal(word,'gy')) then
          nti = 2
          prop = word
       else if (equal(word,'gz')) then
          nti = 3
          prop = word
       else if (equal(word,'gmod')) then
          nti = 4
          prop = word
       else if (equal(word,'hxx')) then
          nti = 5
          prop = word
       else if (equal(word,'hxy') .or. equal(word,'hyx')) then
          nti = 6
          prop = word
       else if (equal(word,'hxz') .or. equal(word,'hzx')) then
          nti = 7
          prop = word
       else if (equal(word,'hyy')) then
          nti = 8
          prop = word
       else if (equal(word,'hyz') .or. equal(word,'hzy')) then
          nti = 9
          prop = word
       else if (equal(word,'hzz')) then
          nti = 10
          prop = word
       else if (equal(word,'lap')) then
          nti = 11
          prop = word
       else if (len_trim(word) > 0) then
          call ferror('rhoplot_line','Unknown keyword in LINE',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! open the output
    if (len_trim(outfile) > 0) then
       luout = fopen_write(outfile)
    else
       luout = uout
       write (uout,'("* LINE ")')
    end if

    ! header
    write (luout,'("# Field values (f) along a line (d = distance).")')
    write (luout,'("# ",4a15,1p,2a20,0p)') "x","y","z","d","f",string(prop)

    ! calculate the line
    if (.not.sy%c%ismolecule) then
       x0 = sy%c%x2c(x0)
       x1 = sy%c%x2c(x1)
    else
       x0 = x0 / dunit0(iunit) - sy%c%molx0
       x1 = x1 / dunit0(iunit) - sy%c%molx0
    endif
    allocate(rhoout(np),lapout(np))
    lapout = 0d0

    lerrmsg%s = ""
    !$omp parallel do private(xp,dist,res,iok,rhopt,lappt) firstprivate(lerrmsg)
    do i=1,np
       xp = x0 + (x1 - x0) * real(i-1,8) / real(np-1,8)
       if (id >= 0) then
          if (nti == 0) then
             call sy%f(id)%grd(xp,0,res,periodic=.not.sy%c%ismolecule)
          elseif (nti >= 1 .and. nti <= 4) then
             call sy%f(id)%grd(xp,1,res,periodic=.not.sy%c%ismolecule)
          else
             call sy%f(id)%grd(xp,2,res,periodic=.not.sy%c%ismolecule)
          end if

          rhopt = res%f
          select case(nti)
          case (0)
             lappt = rhopt
          case (1)
             lappt = res%gf(1)
          case (2)
             lappt = res%gf(2)
          case (3)
             lappt = res%gf(3)
          case (4)
             lappt = res%gfmod
          case (5)
             lappt = res%hf(1,1)
          case (6)
             lappt = res%hf(1,2)
          case (7)
             lappt = res%hf(1,3)
          case (8)
             lappt = res%hf(2,2)
          case (9)
             lappt = res%hf(2,3)
          case (10)
             lappt = res%hf(3,3)
          case (11)
             lappt = res%del2f
          end select
       else
          rhopt = sy%eval(expr,lerrmsg%s,xp)
          lappt = rhopt
       end if

       !$omp critical (iowrite)
       rhoout(i) = rhopt
       lapout(i) = lappt
       if (len_trim(lerrmsg%s) > 0) &
          call ferror('rhoplot_line',lerrmsg%s,faterr)
       !$omp end critical (iowrite)
    enddo
    !$omp end parallel do

    ! write the line to output
    do i = 1, np
       xp = x0 + (x1 - x0) * real(i-1,8) / real(np-1,8)
       dist = norm2(xp-x0) * dunit0(iunit)
       if (.not.sy%c%ismolecule) then
          xout = sy%c%c2x(xp)
       else
          xout = (xp + sy%c%molx0) * dunit0(iunit)
       end if
       if (nti == 0) then
          write (luout,'(99(" ",A))') (string(xout(j),'f',18,13),j=1,3),&
             string(dist,'f',18,13), string(rhoout(i),'e',20,13)
       else
          write (luout,'(99(" ",A))') (string(xout(j),'f',18,13),j=1,3),&
             string(dist,'f',18,13), string(rhoout(i),'e',20,13),&
             string(lapout(i),'e',20,13)
       end if
    end do
    write (luout,*)
    if (len_trim(outfile) > 0) then
       call fclose(luout)
       write (uout,'("* LINE written to file: ",A/)') string(outfile)
    end if

    ! clean up
    deallocate(rhoout,lapout)

  end subroutine rhoplot_line

  ! Calculate properties on a 3d cube
  module subroutine rhoplot_cube(line)
    use systemmod, only: sy
    use fieldmod, only: type_grid
    use grid3mod, only: grid3
    use global, only: eval_next, dunit0, iunit, fileroot
    use arithmetic, only: eval
    use tools_io, only: lgetword, faterr, ferror, equal, getword, &
       isexpression_or_word, uout, string, isinteger
    use types, only: scalar_value, vstring
    use param, only: eye
    use iso_c_binding, only: c_loc
    character*(*), intent(in) :: line

    integer :: lp, nti, id, nn(3)
    real*8 :: x0(3), x1(3), xp(3), lappt
    real*8 :: rgr, dd(3), xd(3,3)
    integer :: lp2, nder
    character(len=:), allocatable :: word, outfile, expr, wext1
    type(scalar_value) :: res
    logical :: ok, doortho
    integer :: ix, iy, iz, i, ibnd, ik, inr(3), ispin
    real*8, allocatable :: lf(:,:,:), fake(:,:,:), ff(:,:,:)
    logical :: dogrid, useexpr, doheader
    integer :: outform, ishift(3), dopsi
    type(grid3) :: faux
    complex*16, allocatable :: caux(:,:,:)
    type(vstring) :: lerrmsg

    integer, parameter :: outform_cube = 1
    integer, parameter :: outform_bincube = 2
    integer, parameter :: outform_vasp = 3
    integer, parameter :: outform_xsf = 4

    integer, parameter :: psi_none = 0
    integer, parameter :: psi_mlwf = 1
    integer, parameter :: psi_wannier = 2
    integer, parameter :: psi_unk = 3
    integer, parameter :: psi_psink = 4

    integer, parameter :: nti_none = 0
    integer, parameter :: nti_f = 1
    integer, parameter :: nti_gx = 2
    integer, parameter :: nti_gy = 3
    integer, parameter :: nti_gz = 4
    integer, parameter :: nti_gmod = 5
    integer, parameter :: nti_hxx = 6
    integer, parameter :: nti_hxy = 7
    integer, parameter :: nti_hxz = 8
    integer, parameter :: nti_hyy = 9
    integer, parameter :: nti_hyz = 10
    integer, parameter :: nti_hzz = 11
    integer, parameter :: nti_lap = 12
    integer, parameter :: nti_real = 13
    integer, parameter :: nti_imag = 14
    integer, parameter :: nti_abs = 15

    ! read the points
    doortho = .false.
    lp = 1
    lp2 = 1
    dogrid = .false.
    dopsi = psi_none
    doheader = .false.
    word = lgetword(line,lp)
    ok = eval_next(x0(1),word,lp2)
    if (ok) then
       ! read initial and final points
       ok = ok .and. eval_next(x0(2),line,lp)
       ok = ok .and. eval_next(x0(3),line,lp)
       ok = ok .and. eval_next(x1(1),line,lp)
       ok = ok .and. eval_next(x1(2),line,lp)
       ok = ok .and. eval_next(x1(3),line,lp)
       if (.not. ok) then
          call ferror('rhoplot_cube','wrong CUBE syntax',faterr,line,syntax=.true.)
          return
       end if

       ! If it is a molecule, that was Cartesian
       if (sy%c%ismolecule) then
          x0 = sy%c%c2x(x0 / dunit0(iunit) - sy%c%molx0)
          x1 = sy%c%c2x(x1 / dunit0(iunit) - sy%c%molx0)
       endif

       ! cubic cube
       xd = 0d0
       do i = 1, 3
          xd(i,i) = x1(i) - x0(i)
       end do
    elseif (equal(word,'cell')) then
       ! whole cell, maybe non-cubic
       xd = eye
       x0 = 0d0
    elseif (equal(word,'grid')) then
       dogrid = .true.
       xd = eye
       x0 = 0d0
    elseif (equal(word,'mlwf').or.equal(word,'wannier').or.equal(word,'unk').or.equal(word,'psink')) then
       dogrid = .true.
       inr = 0
       if (equal(word,'mlwf')) then
          dopsi = psi_mlwf
       elseif (equal(word,'unk')) then
          dopsi = psi_unk
       elseif (equal(word,'psink')) then
          dopsi = psi_psink
       else
          dopsi = psi_wannier
       end if
       ok = isinteger(ibnd,line,lp)
       if (dopsi == psi_unk.or.dopsi == psi_psink) then
          ok = ok .and. isinteger(ik,line,lp)
       end if
       if (dopsi /= psi_unk) then
          ok = ok .and. isinteger(inr(1),line,lp)
          ok = ok .and. isinteger(inr(2),line,lp)
          ok = ok .and. isinteger(inr(3),line,lp)
       end if
       if (.not. ok) then
          call ferror('rhoplot_cube','wrong MLWF/WANNIER/UNK/PSINK syntax',faterr,line,syntax=.true.)
          return
       end if
       xd = eye
       x0 = 0d0
       ispin = 1
    endif

    ! calculate the distances
    x0 = sy%c%x2c(x0)
    do i = 1, 3
       xd(:,i) = sy%c%x2c(xd(:,i))
       dd(i) = norm2(xd(:,i))
    end do

    if (.not.dogrid) then
       ! read number of points or grid resolution
       lp2 = lp
       ok = eval_next(nn(1),line,lp)
       ok = ok .and. eval_next(nn(2),line,lp)
       ok = ok .and. eval_next(nn(3),line,lp)
       if (.not. ok) then
          lp = lp2
          ok = eval_next(rgr,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_cube','wrong CUBE syntax',faterr,line,syntax=.true.)
             return
          end if
          rgr = rgr / dunit0(iunit)
          nn = nint(dd / rgr) + 1
       else
          do i = 1, 3
             nn(i) = max(nn(i),2)
          end do
       end if
    end if

    ! read additional options
    ishift = 0
    nti = nti_none
    nder = 0
    id = sy%iref
    useexpr = .false.
    outform = outform_cube
    outfile = trim(fileroot) // ".cube"
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'file')) then
          outfile = getword(line,lp)
          if (len_trim(outfile) < 1) then
             call ferror('rhoplot_cube','file name not found',faterr,line,syntax=.true.)
             return
          end if
          wext1 = outfile(index(outfile,'.',.true.)+1:)
          if (equal(wext1,'cube')) then
             outform = outform_cube
          elseif (equal(wext1,'bincube')) then
             outform = outform_bincube
          elseif (equal(wext1,'xsf')) then
             outform = outform_xsf
          else
             outform = outform_vasp
          end if
       else if (equal(word,'field')) then
          lp2 = lp
          word = getword(line,lp)
          id = sy%fieldname_to_idx(word)
          if (id < 0) then
             lp = lp2
             ok = isexpression_or_word(expr,line,lp)
             useexpr = .true.
             if (.not.ok) then
                call ferror('rhoplot_cube','wrong FIELD in CUBE',faterr,line,syntax=.true.)
                return
             end if
          else
             if (.not.sy%goodfield(id)) then
                call ferror('rhoplot_cube','field not allocated',faterr,line,syntax=.true.)
                return
             end if
          end if
       else if (equal(word,'header')) then
          doheader = .true.
       else if (equal(word,'spin')) then
          ok = isinteger(ispin,line,lp)
          if (.not.ok) then
             call ferror('rhoplot_cube','wrong SPIN keyword in CUBE',faterr,line,syntax=.true.)
             return
          end if
          if (ispin < 1 .or. ispin > 2) then
             call ferror('rhoplot_cube','SPIN keyword in CUBE must be either 1 or 2',faterr,line,syntax=.true.)
             return
          end if
       else if (equal(word,'f')) then
          nti = nti_f
          nder = 0
       elseif (equal(word,'gx')) then
          nti = nti_gx
          nder = 1
       else if (equal(word,'gy')) then
          nti = nti_gy
          nder = 1
       else if (equal(word,'gz')) then
          nti = nti_gz
          nder = 1
       else if (equal(word,'gmod')) then
          nti = nti_gmod
          nder = 1
       else if (equal(word,'hxx')) then
          nti = nti_hxx
          nder = 2
       else if (equal(word,'hxy') .or. equal(word,'hyx')) then
          nti = nti_hxy
          nder = 2
       else if (equal(word,'hxz') .or. equal(word,'hzx')) then
          nti = nti_hxz
          nder = 2
       else if (equal(word,'hyy')) then
          nti = nti_hyy
          nder = 2
       else if (equal(word,'hyz') .or. equal(word,'hzy')) then
          nti = nti_hyz
          nder = 2
       else if (equal(word,'hzz')) then
          nti = nti_hzz
          nder = 2
       else if (equal(word,'lap')) then
          nti = nti_lap
          nder = 2
       else if (equal(word,'real')) then
          nti = nti_real
       else if (equal(word,'imag')) then
          nti = nti_imag
       else if (equal(word,'abs')) then
          nti = nti_abs
       else if (equal(word,'shift')) then
          if (.not.dogrid) then
             call ferror('rhoplot_cube','SHIFT can only be used with the GRID option in CUBE',faterr,line,syntax=.true.)
             return
          end if
          ok = isinteger(ishift(1),line,lp)
          ok = ok .and. isinteger(ishift(2),line,lp)
          ok = ok .and. isinteger(ishift(3),line,lp)
       else if (equal(word,'ortho')) then
          doortho = .true.
       else if (len_trim(word) > 0) then
          call ferror('rhoplot_cube','Unknown keyword in CUBE',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! re-do the lattice vectors if orthogonal vectors were requested
    if (doortho) then
       x1 = x0 + xd(:,1) + xd(:,2) + xd(:,3)
       xd = 0d0
       do i = 1, 3
          xd(i,i) = x1(i) - x0(i)
       end do
    end if

    ! step sizes and various checks
    if (dogrid) then
       ok = (id > -1)
       if (ok) ok = ok .and. sy%f(id)%type == type_grid .and. allocated(sy%f(id)%grid)
       if (.not.ok) then
          call ferror('rhoplot_cube','CUBE GRID/MLWF/... can only be used with grid fields',faterr,syntax=.true.)
          return
       end if
       nn = sy%f(id)%grid%n
       if (dopsi /= psi_none.and..not.sy%f(id)%grid%isqe) then
          call ferror('rhoplot_cube','CUBE MLWF/WANNIER/... requires a QE wavefunction file (pwc)',faterr,syntax=.true.)
          return
       end if
       if (dopsi == psi_mlwf.and..not.sy%f(id)%grid%iswan) then
          call ferror('rhoplot_cube','CUBE MLWF requires a wannier90 checkpoint file (chk)',faterr,syntax=.true.)
          return
       end if
       if (dopsi /= psi_none.and.(ibnd < 1 .or. ibnd > sy%f(id)%grid%qe%nbnd)) then
          call ferror('rhoplot_cube','CUBE MLWF/WANNIER/...: incorrect band number',faterr,syntax=.true.)
          return
       end if
       if (sy%f(id)%grid%partial) then
          call ferror('rhoplot_cube','CUBE MLWF/WANNIER/...: cannot be used with non-periodic grids',faterr,syntax=.true.)
          return
       end if
    end if

    ! convert to one-step vectors in each direction
    do i = 1, 3
       xd(:,i) = xd(:,i) / real(nn(i),8)
    end do

    ! write cube header
    write (uout,'("* CUBE written to file: ",A/)') string(outfile)
    if (doheader) then
       if (outform == outform_bincube) then
          call ferror("rhoplot_cube","BINCUBE format is incompatible with HEADER",faterr)
       elseif (outform == outform_cube) then
          call sy%c%writegrid_cube(fake,outfile,.true.,.false.,xd,x0+sy%c%molx0)
       elseif (outform == outform_vasp) then
          call sy%c%writegrid_vasp(fake,outfile,.true.)
       elseif (outform == outform_xsf) then
          call sy%c%writegrid_xsf(fake,outfile,.true.)
       endif
       return
    end if

    ! calculate properties
    if (dogrid) then
       if (dopsi == psi_none) then
          ! GRID keyword
          if (sy%f(id)%usecore) then
             call sy%c%promolecular_array3(ff,sy%f(id)%grid%n,sy%f(id)%zpsp)
             call faux%from_array3(ff,sy%c%m_x2c,c_loc(sy%c))
             deallocate(ff)
             faux%f = faux%f + sy%f(id)%grid%f
          else
             faux = sy%f(id)%grid
          end if
       else
          allocate(caux(nn(1),nn(2),nn(3)))
          if (dopsi == psi_mlwf .or. dopsi == psi_wannier) then
             ! MLWF and WANNIER keywords
             call sy%f(id)%grid%get_qe_wnr_standalone(sy%f(id)%c%omega,ibnd,ispin,inr,(dopsi == psi_mlwf),caux)
          elseif (dopsi == psi_unk .or. dopsi == psi_psink) then
             ! UNK and PSINK keyword
             call sy%f(id)%grid%get_qe_psink_standalone(sy%f(id)%c%omega,ibnd,ik,ispin,(dopsi == psi_psink),inr,caux)
          end if
          if (nti == nti_none .or. nti == nti_abs) then
             faux%f = abs(caux)
          elseif (nti == nti_real) then
             faux%f = real(caux)
          elseif (nti == nti_imag) then
             faux%f = aimag(caux)
          end if
       end if

       if (outform == outform_bincube) then
          call sy%c%writegrid_cube(faux%f,outfile,.false.,.true.,xd0=xd,x00=x0+sy%c%molx0)
       elseif (outform == outform_cube) then
          call sy%c%writegrid_cube(faux%f,outfile,.false.,.false.,xd0=xd,x00=x0+sy%c%molx0,ishift0=ishift)
       elseif (outform == outform_vasp) then
          call sy%c%writegrid_vasp(faux%f,outfile,.false.)
       elseif (outform == outform_xsf) then
          call sy%c%writegrid_xsf(faux%f,outfile,.false.)
       endif
    else
       ! no special keywords used
       ok = .false.
       if (useexpr) then
          call faux%new_eval(c_loc(sy),c_loc(sy%c),nn,expr,sy%c%m_x2c)
          ok = faux%isinit
       end if
       allocate(lf(nn(1),nn(2),nn(3)))
       if (ok) then
          lf = faux%f
          call faux%end()
       else
          lerrmsg%s = ""
          !$omp parallel do private(xp,res,lappt) firstprivate(lerrmsg)
          do iz = 0, nn(3)-1
             do iy = 0, nn(2)-1
                do ix = 0, nn(1)-1
                   xp = x0 + real(ix,8) * xd(:,1) + real(iy,8) * xd(:,2) &
                      + real(iz,8) * xd(:,3)

                   if (.not.useexpr) then
                      call sy%f(id)%grd(xp,nder,res,periodic=.not.sy%c%ismolecule)
                      select case(nti)
                      case (nti_none,nti_f)
                         lappt = res%f
                      case (nti_gx)
                         lappt = res%gf(1)
                      case (nti_gy)
                         lappt = res%gf(2)
                      case (nti_gz)
                         lappt = res%gf(3)
                      case (nti_gmod)
                         lappt = res%gfmod
                      case (nti_hxx)
                         lappt = res%hf(1,1)
                      case (nti_hxy)
                         lappt = res%hf(1,2)
                      case (nti_hxz)
                         lappt = res%hf(1,3)
                      case (nti_hyy)
                         lappt = res%hf(2,2)
                      case (nti_hyz)
                         lappt = res%hf(2,3)
                      case (nti_hzz)
                         lappt = res%hf(3,3)
                      case (nti_lap)
                         lappt = res%del2f
                      end select
                   else
                      lappt = sy%eval(expr,lerrmsg%s,xp)
                   end if
                   !$omp critical (fieldwrite)
                   lf(ix+1,iy+1,iz+1) = lappt
                   if (len_trim(lerrmsg%s) > 0) &
                      call ferror('rhoplot_cube',lerrmsg%s,faterr)
                   !$omp end critical (fieldwrite)
                end do
             end do
          end do
          !$omp end parallel do
       end if
       ! cube body
       if (outform == outform_bincube) then
          call sy%c%writegrid_cube(lf,outfile,.false.,.true.,xd,x0+sy%c%molx0)
       elseif (outform == outform_cube) then
          call sy%c%writegrid_cube(lf,outfile,.false.,.false.,xd,x0+sy%c%molx0)
       elseif (outform == outform_vasp) then
          call sy%c%writegrid_vasp(lf,outfile,.false.)
       elseif (outform == outform_xsf) then
          call sy%c%writegrid_xsf(lf,outfile,.false.)
       endif
       if (allocated(lf)) deallocate(lf)
    end if

  end subroutine rhoplot_cube

  !> Calculate properties on a plane.
  module subroutine rhoplot_plane(line)
    use systemmod, only: sy
    use global, only: eval_next, dunit0, iunit, fileroot, iunitname0, iunit, dunit0
    use arithmetic, only: eval
    use tools_io, only: ferror, faterr, lgetword, equal, getword, &
       isexpression_or_word, fopen_write, uout, string, fclose
    use tools_math, only: plane_scale_extend, assign_ziso, niso_manual,&
       niso_lin, niso_log, niso_atan, niso_bader
    use types, only: scalar_value, realloc, vstring
    character*(*), intent(in) :: line

    integer :: lp2, lp, nti, id, luout, nx, ny, niso_type, niso, nn
    real*8 :: x0(3), x1(3), x2(3), xp(3), du, dv, rhopt
    real*8 :: uu(3), vv(3), fmin, fmax
    real*8 :: sx0, sy0, zx0, zx1, zy0, zy1, zmin, zmax, rdum
    logical :: docontour, dorelief, docolormap, fset
    character(len=:), allocatable :: word, outfile, root0, expr
    type(scalar_value) :: res
    logical :: ok
    integer :: ix, iy, cmopt, nder
    real*8, allocatable :: ff(:,:), ziso(:)
    type(vstring) :: lerrmsg

    ! read the points
    lp = 1
    ok = eval_next(x0(1),line,lp)
    ok = ok .and. eval_next(x0(2),line,lp)
    ok = ok .and. eval_next(x0(3),line,lp)
    ok = ok .and. eval_next(x1(1),line,lp)
    ok = ok .and. eval_next(x1(2),line,lp)
    ok = ok .and. eval_next(x1(3),line,lp)
    ok = ok .and. eval_next(x2(1),line,lp)
    ok = ok .and. eval_next(x2(2),line,lp)
    ok = ok .and. eval_next(x2(3),line,lp)
    if (.not. ok) then
       call ferror('rhoplot_plane','Wrong PLANE command: x0, x1, x2',faterr,line,syntax=.true.)
       return
    end if
    if (sy%c%ismolecule) then
       x0 = sy%c%c2x(x0 / dunit0(iunit) - sy%c%molx0)
       x1 = sy%c%c2x(x1 / dunit0(iunit) - sy%c%molx0)
       x2 = sy%c%c2x(x2 / dunit0(iunit) - sy%c%molx0)
    endif

    ok = eval_next(nx,line,lp)
    ok = ok .and. eval_next(ny,line,lp)
    if (.not. ok) then
       call ferror('rhoplot_plane','Wrong PLANE command: nx and ny',faterr,line,syntax=.true.)
       return
    end if

    ! at least two points
    nx = max(nx,2)
    ny = max(ny,2)

    ! read additional options
    RHOP_Hmax = 1d-1
    fset = .false.
    sx0 = 1d0
    sy0 = 1d0
    zx0 = 0d0
    zx1 = 0d0
    zy0 = 0d0
    zy1 = 0d0
    nti = 0
    docontour = .false.
    dorelief = .false.
    docolormap = .false.
    id = sy%iref
    outfile = trim(fileroot) // "_plane.dat"
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'file')) then
          outfile = getword(line,lp)
          if (len_trim(outfile) < 1) then
             call ferror('rhoplot_plane','file name not found',faterr,line,syntax=.true.)
             return
          end if
       else if (equal(word,'field')) then
          lp2 = lp
          word = getword(line,lp)
          id = sy%fieldname_to_idx(word)
          if (id < 0) then
             lp = lp2
             ok = isexpression_or_word(expr,line,lp)
             if (.not.ok) then
                call ferror('rhoplot_plane','wrong FIELD in PLANE',faterr,line,syntax=.true.)
                return
             end if
          else
             if (.not.sy%goodfield(id)) then
                call ferror('rhoplot_plane','field not allocated',faterr,line,syntax=.true.)
                return
             end if
          end if
       elseif (equal(word,'scale')) then
          ok = eval_next(sx0,line,lp)
          ok = ok .and. eval_next(sy0,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_plane','wrong SCALE keyword in PLANE',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'extendx')) then
          ok = eval_next(zx0,line,lp)
          ok = ok .and. eval_next(zx1,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_plane','wrong EXTENDX keyword in PLANE',faterr,line,syntax=.true.)
             return
          end if
          zx0 = zx0 / dunit0(iunit)
          zx1 = zx1 / dunit0(iunit)
       elseif (equal(word,'extendy')) then
          ok = eval_next(zy0,line,lp)
          ok = ok .and. eval_next(zy1,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_plane','wrong EXTENDY keyword in PLANE',faterr,line,syntax=.true.)
             return
          end if
          zy0 = zy0 / dunit0(iunit)
          zy1 = zy1 / dunit0(iunit)
       elseif (equal(word,'labelz')) then
          ok = eval_next(RHOP_Hmax,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_plane','wrong LABELZ keyword in PLANE',faterr,line,syntax=.true.)
             return
          end if
       else if (equal(word,'relief')) then
          dorelief = .true.
          zmin = -1d0
          zmax = 1d0
          ok = eval_next(zmin,line,lp)
          ok = ok .and. eval_next(zmax,line,lp)
          if (.not.ok) then
             call ferror('rhoplot_plane','wrong levels in RELIEF',faterr,line,syntax=.true.)
             return
          end if
       else if (equal(word,'contour')) then
          lp2 = lp
          word = lgetword(line,lp)
          docontour = .true.
          if (equal(word,'lin')) then
             niso_type = niso_lin
          elseif (equal(word,'log')) then
             niso_type = niso_log
          elseif (equal(word,'atan')) then
             niso_type = niso_atan
          elseif (equal(word,'bader')) then
             niso_type = niso_bader
          else
             niso_type = niso_manual
             lp = lp2
             ok = eval_next(rdum,line,lp)
             if (ok) then
                niso = 1
                if (allocated(ziso)) deallocate(ziso)
                allocate(ziso(1))
                ziso(1) = rdum
                do while (.true.)
                   lp2 = lp
                   ok = eval_next(rdum,line,lp)
                   if (.not.ok) exit
                   niso = niso + 1
                   if (niso > size(ziso,1)) call realloc(ziso,2*niso)
                   ziso(niso) = rdum
                end do
                call realloc(ziso,niso)
             else
                call ferror("rhoplot_plane","Unknown contour keyword",faterr,line,syntax=.true.)
             end if
          end if
          if (niso_type == niso_lin .or. niso_type == niso_log .or. niso_type == niso_atan) then
             ok = eval_next(niso,line,lp)
             if (.not.ok) then
                call ferror("rhoplot_plane","number of isovalues not found",faterr,line,syntax=.true.)
                return
             end if
             lp2 = lp
             ok = eval_next(fmin,line,lp)
             ok = ok .and. eval_next(fmax,line,lp)
             if (.not.ok) then
                lp = lp2
             else
                fset = .true.
             end if
          end if
       else if (equal(word,'colormap')) then
          docolormap = .true.
          lp2 = lp
          word = lgetword(line,lp)
          cmopt = 0
          if (equal(word,'log')) then
             cmopt = 1
          elseif (equal(word,'atan')) then
             cmopt = 2
          else
             lp = lp2
          end if
       elseif (equal(word,'f')) then
          nti = 0
       else if (equal(word,'gx')) then
          nti = 1
       else if (equal(word,'gy')) then
          nti = 2
       else if (equal(word,'gz')) then
          nti = 3
       else if (equal(word,'gmod')) then
          nti = 4
       else if (equal(word,'hxx')) then
          nti = 5
       else if (equal(word,'hxy') .or. equal(word,'hyx')) then
          nti = 6
       else if (equal(word,'hxz') .or. equal(word,'hzx')) then
          nti = 7
       else if (equal(word,'hyy')) then
          nti = 8
       else if (equal(word,'hyz') .or. equal(word,'hzy')) then
          nti = 9
       else if (equal(word,'hzz')) then
          nti = 10
       else if (equal(word,'lap')) then
          nti = 11
       else if (len_trim(word) > 0) then
          call ferror('rhoplot_plane','Unknown keyword in PLANE',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! root file name
    nn = index(outfile,".",.true.)
    if (nn == 0) nn = len(trim(outfile)) + 1
    root0 = outfile(1:nn-1)

    ! Transform to Cartesian, extend and scale, and set up the plane
    ! vectors.
    x0 = sy%c%x2c(x0)
    x1 = sy%c%x2c(x1)
    x2 = sy%c%x2c(x2)
    call plane_scale_extend(x0,x1,x2,sx0,sy0,zx0,zx1,zy0,zy1)
    uu = (x1-x0) / real(nx-1,8)
    vv = (x2-x0) / real(ny-1,8)
    du = norm2(uu) * dunit0(iunit)
    dv = norm2(vv) * dunit0(iunit)

    ! allocate space for field values on the plane
    allocate(ff(nx,ny))
    if (nti == 0) then
       nder = 0
    elseif (nti >= 1 .and. nti <= 4) then
       nder = 1
    else
       nder = 2
    end if

    lerrmsg%s = ""
    !$omp parallel do private (xp,res,rhopt) firstprivate(lerrmsg)
    do ix = 1, nx
       do iy = 1, ny
          xp = x0 + real(ix-1,8) * uu + real(iy-1,8) * vv

          if (id >= 0) then
             call sy%f(id)%grd(xp,nder,res,periodic=.not.sy%c%ismolecule)
             select case(nti)
             case (0)
                rhopt = res%f
             case (1)
                rhopt = res%gf(1)
             case (2)
                rhopt = res%gf(2)
             case (3)
                rhopt = res%gf(3)
             case (4)
                rhopt = res%gfmod
             case (5)
                rhopt = res%hf(1,1)
             case (6)
                rhopt = res%hf(1,2)
             case (7)
                rhopt = res%hf(1,3)
             case (8)
                rhopt = res%hf(2,2)
             case (9)
                rhopt = res%hf(2,3)
             case (10)
                rhopt = res%hf(3,3)
             case (11)
                rhopt = res%del2f
             end select
          else
             rhopt = sy%eval(expr,lerrmsg%s,xp)
          endif
          !$omp critical (write)
          ff(ix,iy) = rhopt
          if (len_trim(lerrmsg%s) > 0) &
             call ferror('rhoplot_plane',lerrmsg%s,faterr)
          !$omp end critical (write)
       end do
    end do
    !$omp end parallel do

    ! open the output
    if (len_trim(outfile) > 0) then
       luout = fopen_write(outfile)
    else
       luout = uout
       write (uout,'("* PLANE ")')
    end if

    ! header
    write (luout,'("# Field values (and derivatives) on a plane (units=",A,").")') iunitname0(iunit)
    write (luout,'("# x y z u v f ")')

    x0 = sy%c%c2x(x0)
    x1 = sy%c%c2x(x1)
    x2 = sy%c%c2x(x2)
    uu = sy%c%c2x(uu)
    vv = sy%c%c2x(vv)
    do ix = 1, nx
       do iy = 1, ny
          xp = x0 + real(ix-1,8) * uu + real(iy-1,8) * vv
          xp = (sy%c%x2c(xp) + sy%c%molx0) * dunit0(iunit)
          write (luout,'(" ",5(f15.10," "),1p,1(e18.10," "),0p)') xp, real(ix-1,8)*du, real(iy-1,8)*dv, ff(ix,iy)
       end do
       write (luout,*)
    end do

    ! contour/relief/colormap plots
    if (docontour) then
       if (.not.fset) then
          fmin = minval(ff)
          fmax = maxval(ff)
       end if
       call assign_ziso(niso_type,niso,ziso,fmin,fmax)
       call contour(ff,x0,x1,x2,nx,ny,niso,ziso,root0,.true.,.true.)
    end if
    if (dorelief) call relief(root0,string(outfile),zmin,zmax)
    if (docolormap) call colormap(root0,string(outfile),cmopt,x0,x1,x2,.true.)

    if (len_trim(outfile) > 0) then
       call fclose(luout)
       write (uout,'("* PLANE written to file: ",A)') string(outfile)
       write (uout,*)
    end if

    if (allocated(ziso)) deallocate(ziso)
    deallocate(ff)

  end subroutine rhoplot_plane

  !> Plot of gradient paths and contours in the style of aimpac's grdvec.
  module subroutine rhoplot_grdvec()
    use systemmod, only: sy
    use global, only: fileroot, eval_next, dunit0, iunit
    use tools_io, only: uout, uin, ucopy, getline, lgetword, equal,&
       faterr, ferror, string, ioj_right, fopen_write, getword, fclose
    use tools_math, only: plane_scale_extend, assign_ziso, &
       niso_manual, niso_atan, niso_lin, niso_log, niso_bader
    use types, only: scalar_value, realloc
    character(len=:), allocatable :: line, word, datafile, rootname
    integer :: lpold, lp, udat, ll, i, j
    integer :: updum, dndum, updum1, dndum1
    real*8  :: xp(3), rhopt
    logical :: doagain, ok, autocheck
    real*8  :: r0(3), r1(3), r2(3), xdum
    real*8  :: r0out(3), r1out(3), r2out(3)
    real*8  :: q0(3), xo0(3), xo1(3), xo2(3)
    integer :: cpid
    integer :: niso_type, nfi, ix, iy
    real*8 :: sx0, sy0, zx0, zx1, zy0, zy1, rdum
    real*8 :: x0(3), uu(3), vv(3), fmin, fmax
    logical :: docontour, dograds, goodplane, fset
    integer :: n1, n2, niso, nder
    type(scalar_value) :: res
    real*8, allocatable :: ff(:,:), ziso(:)

    ! Header
    write (uout,'("* GRDVEC: gradient paths and contours in 2d")')

    ! Initialization
    grpcpeps = 1d-2
    grphcutoff = 1.0d-3
    grpproj = 1
    rootname = trim(fileroot)
    datafile = trim(fileroot) // ".dat"
    autocheck = .false.
    norig = 0
    newncriticp = 0
    cpup = 0
    cpdn = 0
    outcpx = 1d0
    outcpy = 1d0
    docontour = .false.
    dograds = .false.
    nfi = 0
    sx0 = 1d0
    sy0 = 1d0
    zx0 = 0d0
    zx1 = 0d0
    zy0 = 0d0
    zy1 = 0d0
    RHOP_Hmax = 1d-1

    !.Read user options:
    ll = len(line)
    doagain = getline(uin,line,ucopy=ucopy)
    goodplane = .false.
    do while (doagain)
       lp=1
       word = lgetword(line,lp)
       ok = (len_trim(word)>0)  .and. lp.le.ll

       if (equal(word,'files').or.equal(word,'root').or.equal(word,'oname')) then
          ! skip spaces
          rootname = adjustl(line(lp:))
          datafile = trim(rootname) // ".dat"

       else if (equal(word,'plane')) then
          ok = eval_next (r0(1), line, lp)
          ok = ok .and. eval_next (r0(2), line, lp)
          ok = ok .and. eval_next (r0(3), line, lp)
          ok = ok .and. eval_next (r1(1), line, lp)
          ok = ok .and. eval_next (r1(2), line, lp)
          ok = ok .and. eval_next (r1(3), line, lp)
          ok = ok .and. eval_next (r2(1), line, lp)
          ok = ok .and. eval_next (r2(2), line, lp)
          ok = ok .and. eval_next (r2(3), line, lp)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','Bad limits for crystal',faterr,line,syntax=.true.)
             return
          end if
          if (sy%c%ismolecule) then
             r0 = sy%c%c2x(r0 / dunit0(iunit) - sy%c%molx0)
             r1 = sy%c%c2x(r1 / dunit0(iunit) - sy%c%molx0)
             r2 = sy%c%c2x(r2 / dunit0(iunit) - sy%c%molx0)
          endif
          goodplane = .true.

       elseif (equal(word,'scale')) then
          ok = eval_next(sx0,line,lp)
          ok = ok .and. eval_next(sy0,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','wrong SCALE keyword in PLANE',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'extendx')) then
          ok = eval_next(zx0,line,lp)
          ok = ok .and. eval_next(zx1,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','wrong EXTENDX keyword in PLANE',faterr,line,syntax=.true.)
             return
          end if
          zx0 = zx0 / dunit0(iunit)
          zx1 = zx1 / dunit0(iunit)
       elseif (equal(word,'extendy')) then
          ok = eval_next(zy0,line,lp)
          ok = ok .and. eval_next(zy1,line,lp)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','wrong EXTENDY keyword in PLANE',faterr,line,syntax=.true.)
             return
          end if
          zy0 = zy0 / dunit0(iunit)
          zy1 = zy1 / dunit0(iunit)
       else if (equal(word,'outcp')) then
          ok = eval_next(outcpx, line, lp)
          ok = ok .and. eval_next(outcpy, line, lp)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','Bad outcp options',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'hmax')) then
          ok = eval_next (xdum, line, lp)
          if (ok) RHOP_Hmax = xdum / dunit0(iunit)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','Wrong hmax line',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'cp')) then
          newncriticp = newncriticp + 1
          if (newncriticp .gt. mncritp) then
             call ferror('rhoplot_grdvec','too many points in a check order. Increase MNCRITP',faterr,syntax=.true.)
             return
          end if
          ok = eval_next (cpid, line, lp)
          ok = ok .and. eval_next (cpup(newncriticp), line, lp)
          ok = ok .and. eval_next (cpdn(newncriticp), line, lp)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','bad cp option',faterr,line,syntax=.true.)
             return
          end if
          if (cpid <= 0 .or. cpid > sy%f(sy%iref)%ncpcel) then
             call ferror('rhoplot_grdvec','cp not recognized',faterr,line,syntax=.true.)
             return
          end if
          newcriticp(:,newncriticp) = sy%f(sy%iref)%cpcel(cpid)%x
          newtypcrit(newncriticp) = sy%f(sy%iref)%cpcel(cpid)%typ
          dograds = .true.
          autocheck = .true.
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'cpall')) then
          ! copy cps
          do i = 1, sy%f(sy%iref)%ncpcel
             newncriticp = newncriticp + 1
             newcriticp(:,newncriticp) = sy%f(sy%iref)%cpcel(i)%x
             newtypcrit(newncriticp) = sy%f(sy%iref)%cpcel(i)%typ
          end do
          dograds = .true.
          autocheck = .true.
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'bcpall')) then
          ok = eval_next (updum, line, lp)
          if (.not.ok) updum = 2
          ok = eval_next (dndum, line, lp)
          if (.not.ok) dndum = 0

          ! copy bcps
          do i = 1, sy%f(sy%iref)%ncpcel
             if (sy%f(sy%iref)%cpcel(i)%typ == -1) then
                newncriticp = newncriticp + 1
                newcriticp(:,newncriticp) = sy%f(sy%iref)%cpcel(i)%x
                newtypcrit(newncriticp) = sy%f(sy%iref)%cpcel(i)%typ
                cpup(newncriticp) = updum
                cpdn(newncriticp) = dndum
             end if
          end do
          dograds = .true.
          autocheck = .true.
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'rbcpall')) then
          ok = eval_next (updum, line, lp)
          if (.not.ok) updum = 2
          ok = eval_next (dndum, line, lp)
          if (.not.ok) dndum = 0
          ok = eval_next (updum1, line, lp)
          if (.not.ok) updum1 = 0
          ok = eval_next (dndum1, line, lp)
          if (.not.ok) dndum1 = 2

          ! copy rcps and bcps
          do i = 1, sy%f(sy%iref)%ncpcel
             if (sy%f(sy%iref)%cpcel(i)%typ == -1) then
                newncriticp = newncriticp + 1
                newcriticp(:,newncriticp) = sy%f(sy%iref)%cpcel(i)%x
                newtypcrit(newncriticp) = sy%f(sy%iref)%cpcel(i)%typ
                cpup(newncriticp) = updum
                cpdn(newncriticp) = dndum
             else if (sy%f(sy%iref)%cpcel(i)%typ == 1) then
                newncriticp = newncriticp + 1
                newcriticp(:,newncriticp) = sy%f(sy%iref)%cpcel(i)%x
                newtypcrit(newncriticp) = sy%f(sy%iref)%cpcel(i)%typ
                cpup(newncriticp) = updum1
                cpdn(newncriticp) = dndum1
             end if
          end do
          dograds = .true.
          autocheck = .true.
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'orig')) then
          dograds = .true.
          norig = norig + 1
          if (norig .gt. MORIG) then
             call ferror('rhoplot_grdvec','Too many ORIGIN points. Increase MORIG',faterr,syntax=.true.)
             return
          end if
          ok = eval_next (grpx(1,norig), line, lp)
          ok = ok .and. eval_next (grpx(2,norig), line, lp)
          ok = ok .and. eval_next (grpx(3,norig), line, lp)
          ok = ok .and. eval_next (grpatr(norig), line, lp)
          ok = ok .and. eval_next (grpup(norig), line, lp)
          ok = ok .and. eval_next (grpdwn(norig), line, lp)
          if (.not. ok) then
             call ferror('rhoplot_grdvec','Bad limits for 3Dc plot',faterr,line,syntax=.true.)
             return
          end if
          if (sy%c%ismolecule) &
             grpx(:,norig) = sy%c%c2x(grpx(:,norig) / dunit0(iunit) - sy%c%molx0)
          ok = check_no_extra_word()
          if (.not.ok) return

       elseif (equal (word,'check')) then
          ok = check_no_extra_word()
          if (.not.ok) return
          ! read the user-entered points:
          ok = getline(uin,line,.true.,ucopy)
          lp = 1
          word = lgetword (line,lp)
          do while (ok.and..not.equal(word, 'endcheck').and..not.equal(word, 'end'))
             newncriticp = newncriticp + 1
             if (newncriticp .gt. mncritp) then
                call ferror('rhoplot_grdvec','too many points in a check order. Increase MNCRITP',faterr,syntax=.true.)
                return
             end if
             lp = 1
             ok = eval_next (newcriticp(1,newncriticp), line, lp)
             ok = ok .and. eval_next (newcriticp(2,newncriticp), line, lp)
             ok = ok .and. eval_next (newcriticp(3,newncriticp), line, lp)
             q0 = sy%c%x2c(newcriticp(:,newncriticp))
             call sy%f(sy%iref)%grd(q0,2,res,periodic=.not.sy%c%ismolecule)
             newtypcrit(newncriticp) = res%s
             q0 = sy%c%c2x(q0)

             ok = ok .and. getline(uin,line,.true.,ucopy)
             lp = 1
             word = lgetword(line,lp)
          enddo
          dograds = .true.
          autocheck = .true.
          ok = check_no_extra_word()
          if (.not.ok) return

       elseif (equal(word,'contour')) then

          ! pass contours to the gnu file
          docontour = .true.

          ! read the field
          word = lgetword(line,lp)
          if (equal(word,'f') .or. equal(word,'rho')) then
             nfi = 0
          else if (equal(word,'gx')) then
             nfi = 1
          else if (equal(word,'gy')) then
             nfi = 2
          else if (equal(word,'gz')) then
             nfi = 3
          else if (equal(word,'gmod')) then
             nfi = 4
          else if (equal(word,'hxx')) then
             nfi = 5
          else if (equal(word,'hxy') .or. equal(word,'hyx')) then
             nfi = 6
          else if (equal(word,'hxz') .or. equal(word,'hzx')) then
             nfi = 7
          else if (equal(word,'hyy')) then
             nfi = 8
          else if (equal(word,'hyz') .or. equal(word,'hzy')) then
             nfi = 9
          else if (equal(word,'hzz')) then
             nfi = 10
          else if (equal(word,'lap')) then
             nfi = 11
          else
             call ferror('rhoplot_grdvec','contour field keyword needed',faterr,line,syntax=.true.)
             return
          end if

          ! read the number of points
          ok = eval_next (n1, line, lp)
          ok = ok .and. eval_next (n2, line, lp)
          if (.not.ok) then
             call ferror('rhoplot_grdvec','contour number of points needed',faterr,line,syntax=.true.)
             return
          end if
          n1 = max(n1,2)
          n2 = max(n2,2)

          ok = .true.
          fset = .false.
          lpold = lp
          word = lgetword(line,lp)
          if (equal(word,'atan') .or. equal(word,'log') .or. equal(word,'lin')) then
             if (equal(word,'log')) then
                niso_type = niso_log
             elseif (equal(word,'lin')) then
                niso_type = niso_lin
             else
                niso_type = niso_atan
             end if
             ok = eval_next (niso, line, lp)
             if (ok) then
                lpold = lp
                ok = eval_next (fmin, line, lp)
                ok = ok .and. eval_next (fmax, line, lp)
                if (.not.ok) then
                   lp = lpold
                   ok = .true.
                else
                   fset = .true.
                end if
             end if
          else if (equal(word,'bader')) then
             niso_type = niso_bader
          else
             lp = lpold
             niso_type = niso_manual
             if (allocated(ziso)) deallocate(ziso)
             allocate(ziso(1))
             niso = 0
             do while (.true.)
                ok = eval_next(rdum,line,lp)
                if (.not.ok) exit
                niso = niso + 1
                if (niso > size(ziso,1)) call realloc(ziso,2*niso)
                ziso(niso) = rdum
             end do
             if (niso == 0) then
                call ferror("rhoplot_grdvec","wrong contour values",faterr,line,syntax=.true.)
                return
             end if
             call realloc(ziso,niso)
             ok = .true.
          end if

          if (.not.ok) then
             call ferror('rhoplot_grdvec','wrong contour values',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return

       else if (equal(word,'endgrdvec').or.equal(word,'end')) then
          ok = check_no_extra_word()
          if (.not.ok) return
          goto 999
       else
          call ferror('rhoplot_grdvec','Unkown keyword in GRDVEC',faterr,line,syntax=.true.)
          return
       endif
       doagain = getline(uin,line,ucopy=ucopy)
    enddo
    call ferror('rhoplot_grdvec','Unexpected end of input',faterr,line,syntax=.true.)
    return
999 continue
    if (.not.goodplane) then
       call ferror('rhoplot_grdvec','No PLANE given in GRDVEC',faterr,syntax=.true.)
       return
    end if

    ! extend and scale
    r0 = sy%c%x2c(r0)
    r1 = sy%c%x2c(r1)
    r2 = sy%c%x2c(r2)
    r0out = r0
    r1out = r1
    r2out = r2
    call plane_scale_extend(r0,r1,r2,sx0,sy0,zx0,zx1,zy0,zy1)
    call plane_scale_extend(r0out,r1out,r2out,sx0*outcpx,sy0*outcpy,zx0,zx1,zy0,zy1)
    r0 = sy%c%c2x(r0)
    r1 = sy%c%c2x(r1)
    r2 = sy%c%c2x(r2)
    r0out = sy%c%c2x(r0out)
    r1out = sy%c%c2x(r1out)
    r2out = sy%c%c2x(r2out)

    ! output the plane
    write (uout,'("* Name of the output data file: ",a)') string(datafile)
    if (.not.sy%c%ismolecule) then
       xo0 = r0
       xo1 = r1
       xo2 = r2
    else
       xo0 = (sy%c%x2c(r0) + sy%c%molx0) * dunit0(iunit)
       xo1 = (sy%c%x2c(r1) + sy%c%molx0) * dunit0(iunit)
       xo2 = (sy%c%x2c(r2) + sy%c%molx0) * dunit0(iunit)
    end if
    write (uout,'("  Plane origin: ",3(A," "))') (string(xo0(j),'f',12,6,ioj_right),j=1,3)
    write (uout,'("  Plane x-end:  ",3(A," "))') (string(xo1(j),'f',12,6,ioj_right),j=1,3)
    write (uout,'("  Plane y-end:  ",3(A," "))') (string(xo2(j),'f',12,6,ioj_right),j=1,3)

    ! calculate the contour plot
    if (docontour) then
       allocate(ff(n1,n2))
       x0 = sy%c%x2c(r0)
       uu = sy%c%x2c((r1-r0) / real(n1-1,8))
       vv = sy%c%x2c((r2-r0) / real(n2-1,8))
       if (nfi == 0) then
          nder = 0
       else if (nfi >= 1 .and. nfi <= 4) then
          nder = 1
       else
          nder = 2
       end if
       !$omp parallel do private(xp,res,rhopt)
       do ix = 1, n1
          do iy = 1, n2
             xp = x0 + real(ix-1,8) * uu + real(iy-1,8) * vv
             call sy%f(sy%iref)%grd(xp,nder,res,periodic=.not.sy%c%ismolecule)
             select case(nfi)
             case (0)
                rhopt = res%f
             case (1)
                rhopt = res%gf(1)
             case (2)
                rhopt = res%gf(2)
             case (3)
                rhopt = res%gf(3)
             case (4)
                rhopt = res%gfmod
             case (5)
                rhopt = res%hf(1,1)
             case (6)
                rhopt = res%hf(1,2)
             case (7)
                rhopt = res%hf(1,3)
             case (8)
                rhopt = res%hf(2,2)
             case (9)
                rhopt = res%hf(2,3)
             case (10)
                rhopt = res%hf(3,3)
             case (11)
                rhopt = res%del2f
             end select
             !$omp critical (write)
             ff(ix,iy) = rhopt
             !$omp end critical (write)
          end do
       end do
       !$omp end parallel do
       if (.not.fset) then
          fmin = minval(ff)
          fmax = maxval(ff)
       end if
       call assign_ziso(niso_type,niso,ziso,fmin,fmax)
       call contour(ff,r0,r1,r2,n1,n2,niso,ziso,rootname,.false.,.false.)
       if (allocated(ziso)) deallocate(ziso)
       deallocate(ff)
    end if

    udat = fopen_write(datafile)
    indmax = nint(max(maxval(abs(r0out)),maxval(abs(r1out)),maxval(abs(r2out))))
    call plotvec(r0,r1,r2,autocheck,udat)
    call fclose(udat)

    ! print labels (plane info is common)
    call write_fichlabel(rootname)

    ! print gnuplot (plane info is common)
    call write_fichgnu(rootname,.true.,docontour,dograds)
    write (uout,*)

  contains

    function check_no_extra_word()
      character(len=:), allocatable :: aux2
      logical :: check_no_extra_word
      aux2 = getword(line,lp)
      check_no_extra_word = .true.
      if (len_trim(aux2) > 0) then
         call ferror('rhoplot_grdvec','Unknown extra keyword',faterr,line,syntax=.true.)
         check_no_extra_word = .false.
      end if
    end function check_no_extra_word

  end subroutine rhoplot_grdvec

  !xx! private procedures

  !> Contour plots using the 2d-field ff, defined on a plane
  !> determined by poitns r0, r1, r2 (crystallographic coords.). The
  !> number of points in each direction is nx and ny. ziso(1:niso) is the
  !> array contaiing the contour levels. rootname is the root for all
  !> the files generated (.iso, .neg.iso, -grd.dat, -label.gnu,
  !> .gnu). If dolabels, write the labels file. If dognu, write the
  !> gnu file.
  subroutine contour(ff,r0,r1,r2,nx,ny,niso,ziso,rootname,dognu,dolabels)
    use systemmod, only: sy
    use tools_io, only: fopen_write, uout, string, faterr, ferror, fclose
    use tools_math, only: cross, det3, matinv
    integer, intent(in) :: nx, ny
    real*8, intent(in) :: ff(nx,ny)
    real*8, intent(in) :: r0(3), r1(3), r2(3)
    integer, intent(in) :: niso
    real*8, intent(in) :: ziso(niso)
    character*(*), intent(in) :: rootname
    logical, intent(in) :: dognu, dolabels

    character(len=:), allocatable :: root0, fichiso, fichiso1, fichgnu
    integer :: lud, lud1
    integer :: i, j
    real*8 :: r012, ua, va, ub, vc
    real*8, allocatable :: x(:), y(:)

    ! set rootname
    root0 = rootname

    ! name files
    fichgnu = trim(root0) // '-contour.gnu'
    fichiso = trim(root0) // '.iso'
    fichiso1 = trim(root0) // '.neg.iso'

    ! connect units for writing.
    lud = fopen_write(fichiso)
    lud1 = fopen_write(fichiso1)

    write (uout,'("* Name of the contour lines file: ",a)') string(fichiso)
    write (uout,'("* Name of the negative contour lines file: ",a)') string(fichiso1)

    ! geometry
    rp0 = sy%c%x2c(r0)
    rp1 = sy%c%x2c(r1)
    rp2 = sy%c%x2c(r2)
    rp01 = rp1 - rp0
    rp02 = rp2 - rp0
    r01 = norm2(rp01)
    r02 = norm2(rp02)
    r012 = dot_product(rp01,rp02)
    cosalfa = r012/r01/r02
    sinalfa = sqrt(max(1-cosalfa**2,0d0))

    indmax = nint(max(maxval(abs(r0)),maxval(abs(r1)),maxval(abs(r2))))

    ! normal vector and plane equation
    rpn(1:3) = cross(rp01,rp02)
    rpn(4) = -dot_product(rpn(1:3),rp0)
    rpn = rpn / (r01*r02)

    ! plane to cartesian
    amat(:,1) = rp01
    amat(:,2) = rp02
    amat(:,3) = rpn(1:3)

    ! cartesian to plane
    if (abs(det3(amat)) < 1d-15) &
       call ferror('contour','Error in the input plane: singular matrix',faterr)
    bmat = amat
    call matinv(bmat,3)

    ! define plane limits
    ua = 0d0
    ub = r01
    va = 0d0
    vc = r02

    ! calculate grid in each direction
    isneg = minval(ff) < 0d0
    allocate(x(nx))
    allocate(y(ny))
    do i = 1, nx
       do j = 1, ny
          x(i) = ua + (i-1) * (ub-ua) / real(nx-1,8)
          y(j) = va + (j-1) * (vc-va) / real(ny-1,8)
       end do
    end do

    do i = 1, niso
       call hallarpuntos(ff,ziso(i),x,y,nx,ny)
       if (ziso(i).gt.0) then
          call ordenarpuntos(lud,cosalfa,ziso(i))
       else
          call ordenarpuntos(lud1,cosalfa,ziso(i))
       endif
    enddo
    call fclose(lud)
    call fclose(lud1)

    if (dolabels) call write_fichlabel(root0)
    if (dognu) call write_fichgnu(root0,dolabels,.true.,.false.)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)

  end subroutine contour

  !> Write a gnuplot template for the relief plot
  subroutine relief(rootname,outfile,zmin,zmax)
    use tools_io, only: fopen_write, uout, string, fclose
    real*8, intent(in) :: zmin, zmax
    character*(*), intent(in) :: rootname, outfile

    character(len=:), allocatable :: file
    integer :: lu

    ! file name
    file = trim(rootname) // '-relief.gnu'

    ! connect unit
    lu = fopen_write(file)
    write (uout,'("* Gnuplot file (relief): ",a)') string(file)

    write (lu,'("set terminal pdfcairo")')
    write (lu,'("set output """,A,"-relief.pdf""")') rootname
    write (lu,'("set encoding iso_8859_1")')
    write (lu,'("")')
    write (lu,'("set style line 1 lt 1 lc rgb ""#000000"" ")')
    write (lu,'("")')
    write (lu,'("# Define the zrange and the capping functions")')
    write (lu,'("zmin = ",A)') string(zmin,'e',12,5)
    write (lu,'("zmax = ",A)') string(zmax,'e',12,5)
    write (lu,'("stats """,A,""" u 6 nooutput")') outfile
    write (lu,'("min(x) = (x<zmin)?min=x:zmin")')
    write (lu,'("max(x) = (x>zmax)?max=zmax:x")')
    write (lu,'("set zrange [(zmin<STATS_min)?STATS_min:zmin:(zmax>STATS_max)?STATS_max:zmax]")')
    write (lu,'("")')
    write (lu,'("# tics, etc")')
    write (lu,'("unset colorbox")')
    write (lu,'("unset title")')
    write (lu,'("set format x ""%.1f""")')
    write (lu,'("set format y ""%.1f""")')
    write (lu,'("")')
    write (lu,'("# Surface definition")')
    write (lu,'("set pm3d depthorder hidden3d 1")')
    write (lu,'("set hidden3d")')
    write (lu,'("set style fill transparent solid 0.7")')
    write (lu,'("set palette rgb 9,9,3")')
    write (lu,'("set view 60,45")')
    write (lu,'("set size ratio -1")')
    write (lu,'("")')
    write (lu,'("splot """,A,""" u 4:5:(max($6)) ls 1 w pm3d notitle")') outfile

    ! wrap up
    call fclose(lu)

  end subroutine relief

  !> Write a gnuplot template for the color map plot
  subroutine colormap(rootname,outfile,cmopt,r0,r1,r2,dolabels)
    use systemmod, only: sy
    use tools_math, only: cross, det3, matinv
    use tools_io, only: fopen_write, uout, string, fclose, ferror ,faterr
    character*(*), intent(in) :: rootname, outfile
    integer, intent(in) :: cmopt
    real*8, intent(in) :: r0(3), r1(3), r2(3)
    logical :: dolabels

    character(len=:), allocatable :: file, fichlabel
    integer :: lu
    real*8 :: r012

    ! file name
    file = trim(rootname) // '-colormap.gnu'

    ! geometry
    rp0 = sy%c%x2c(r0)
    rp1 = sy%c%x2c(r1)
    rp2 = sy%c%x2c(r2)
    rp01 = rp1 - rp0
    rp02 = rp2 - rp0
    r01 = norm2(rp01)
    r02 = norm2(rp02)
    r012 = dot_product(rp01,rp02)
    cosalfa = r012/r01/r02
    sinalfa = sqrt(max(1-cosalfa**2,0d0))
    indmax = nint(max(maxval(abs(r0)),maxval(abs(r1)),maxval(abs(r2))))

    ! normal vector and plane equation
    rpn(1:3) = cross(rp01,rp02)
    rpn(4) = -dot_product(rpn(1:3),rp0)
    rpn = rpn / (r01*r02)

    ! plane to cartesian
    amat(:,1) = rp01
    amat(:,2) = rp02
    amat(:,3) = rpn(1:3)

    ! cartesian to plane
    if (abs(det3(amat)) < 1d-15) &
       call ferror('colormap','Error in the input plane: singular matrix',faterr)
    bmat = amat
    call matinv(bmat,3)

    ! connect unit
    lu = fopen_write(file)
    write (uout,'("* Gnuplot file (colormap): ",a)') string(file)

    write (lu,'("set encoding iso_8859_1")')
    write (lu,'("set terminal postscript eps color enhanced ""Helvetica""")')
    write (lu,'("set output """,A,"-colormap.eps""")') rootname
    write (lu,'("")')
    write (lu,'("# line styles")')
    write (lu,'("set style line 1 lt 1 lw 1 lc rgb ""#000000""")')
    write (lu,'("set style line 2 lt 1 lw 1 lc rgb ""#000000""")')
    write (lu,'("")')
    write (lu,'("# title, key, size")')
    write (lu,'("unset title")')
    write (lu,'("unset key")')
    write (lu,'("set size ratio -1")')
    write (lu,'("")')
    write (lu,'("# set pm3d at b map interpolate 5,5")')
    write (lu,'("set pm3d at b map")')
    write (lu,'("")')
    write (lu,'("# tics")')
    write (lu,'("set cbtics")')
    write (lu,'("")')
    write (lu,'("# color schemes")')
    write (lu,'("set palette defined ( 0 ""red"", 1 ""white"", 2 ""green"" ) ")')
    write (lu,'("")')
    write (lu,'("# set contours")')
    write (lu,'("unset clabel")')
    write (lu,'("set contour base")')
    write (lu,'("set cntrparam bspline")')
    write (lu,'("# set cntrparam levels incremental -min,step,max")')
    write (lu,'("")')

    if (dolabels) then
       fichlabel = trim(rootname) // "-label.gnu"
       write (lu,'("load ''",A,"''")') string(fichlabel)
       call write_fichlabel(rootname)
    end if

    if (cmopt == 1) then
       write (lu,'("splot """,A,""" u 4:5:(log(abs($6))) ls 1 w pm3d notitle")') outfile
    elseif (cmopt == 2) then
       write (lu,'("splot """,A,""" u 4:5:(2/pi*atan($6)) ls 1 w pm3d notitle")') outfile
    else
       write (lu,'("splot """,A,""" u 4:5:6 ls 1 w pm3d notitle")') outfile
    end if

    ! wrap up
    call fclose(lu)

  end subroutine colormap

  !> Find contour with value = zc on a surface given by a grid.
  !> uses linear interpolation.
  subroutine hallarpuntos(ff,zc,x,y,nx,ny)
    use param, only: zero
    integer, intent(in) :: nx, ny
    real*8, intent(in) :: ff(nx,ny)
    real*8, intent(in) :: zc
    real*8, intent(in) :: x(:), y(:)

    integer :: i, j
    real*8 :: xa, ya, za, xb, yb, zb
    real*8 :: zazc, zbzc

    ! initialize
    npuntos = 0

    ! run over columns
    do i = 1, nx
       xa = x(i)
       ya = y(1)
       za = ff(i,1)
       zazc = za-zc

       ! check if it is an intersection
       if (zazc.eq.zero) then
          npuntos = npuntos+1
          if (npuntos.gt.size(xc)) then
             return
          end if
          xc(npuntos) = xa
          yc(npuntos) = ya
          pic(npuntos) = i
          pjc(npuntos) = 1d0
       endif

       do j = 2,ny
          xb = x(i)
          yb = y(j)
          zb = ff(i,j)
          zbzc = zb-zc

          ! sign changed, interpolate and write a point
          if ((zazc*zbzc).lt.zero) then
             npuntos = npuntos+1
             if (npuntos.gt.size(xc)) then
                return
             end if
             xc(npuntos) = xb
             yc(npuntos) = ya-zazc*(yb-ya)/(zb-za)
             pic(npuntos) = i
             pjc(npuntos) = j-0.5d0
          else if (zbzc.eq.zero) then
             npuntos = npuntos+1
             if (npuntos.gt.size(xc)) then
                return
             end if
             xc(npuntos) = xb
             yc(npuntos) = yb
             pic(npuntos) = i
             pjc(npuntos) = j
          endif
          ! reassign zazc
          ya = yb
          za = zb
          zazc = zbzc
       enddo
    enddo

    ! run over rows
    do j = 1,ny
       xa = x(1)
       ya = y(j)
       za = ff(1,j)
       zazc = za-zc

       do i = 2, nx
          xb = x(i)
          yb = y(j)
          zb = ff(i,j)
          zbzc = zb-zc
          ! sign changed, interpolate and write a point
          if ((zazc*zbzc).lt.zero) then
             npuntos = npuntos+1
             if (npuntos.gt.size(xc)) then
                return
             end if
             xc(npuntos) = xa-zazc*(xb-xa)/(zb-za)
             yc(npuntos) = yb
             pic(npuntos) = i-0.5d0
             pjc(npuntos) = j
          endif
          ! reassign zazc
          xa = xb
          za = zb
          zazc = zbzc
       enddo
    enddo

  end subroutine hallarpuntos

  !> Determines the connectivity of the set of contour points.
  subroutine ordenarpuntos(luw,calpha,ziso)
    use param, only: one, half, zero
    real*8, parameter :: eps = 0.10d0

    integer, intent(in) :: luw
    real*8, intent(in) :: calpha, ziso

    real*8 :: salpha

    logical :: lc(size(xc)), cerrada, primerabusqueda, hallado
    logical :: malla0, malla1, malla2, incompleta
    integer :: punto(-1*size(xc)-1:size(xc)+1), puntoinicial
    integer :: puntofinal, puntoactual, salto, nptoscurva
    integer :: nptosestudiados
    real*8  :: x(nptcorte+1), y(nptcorte+1)
    integer :: k, n1
    real*8  :: pi0, pj0, pi1, pj1, pi2, pj2, si1, sj1

    ! inline functions
    logical :: adyacentes
    adyacentes(pi0,pj0,pi1,pj1,si1,sj1)=(abs(abs(pi0-pi1)-si1) < eps).and.&
       (abs(abs(pj0-pj1)-sj1) < eps)

    ! initialize
    do k = 1, npuntos
       lc(k) = .false.
    enddo
    nptosestudiados = 0

    ! run over all points
    x = 0d0
    y = 0d0
    n1 = 1
    do while (nptosestudiados.lt.npuntos)

       ! search for a non-used point and start a curve
       k = n1
       do while ((lc(k)) .and. (k.le.npuntos))
          k = k+1
       enddo
       if (k.gt.npuntos) then
          return
       end if
       ! the isoline starts at k
       n1 = k

       lc(n1) = .true.
       pi1 = pic(n1)
       pj1 = pjc(n1)

       ! initialize search vars.
       puntoinicial = 0
       puntofinal = 0
       puntoactual = 0
       punto(puntoactual) = n1
       salto = 1
       primerabusqueda = .true.
       incompleta = .true.
       cerrada = .false.

       ! search
       do while (incompleta)
          hallado = .true.
          do while (hallado)
             hallado = .false.

             ! (pi1,pj1) is the last found point. Search for connection
             ! with the rest of the points in the isoline
             malla1 = (dabs(pi1-dint(pi1)).le.eps) .and.              &
                &                (dabs(pj1-dint(pj1)).le.eps)
             k = n1
             do while ((.not.hallado) .and. (k.lt.npuntos))
                k = k+1
                if (.not.lc(k)) then

                   !.search in non-connected points
                   pi2 = pic(k)
                   pj2 = pjc(k)
                   malla2 = (dabs(pi2-dint(pi2)).le.eps) .and.        &
                      &                      (dabs(pj2-dint(pj2)).le.eps)
                   if (malla1 .or. malla2) then

                      ! one of (pi1,pj1) and (pi2,pj2) is in the grid
                      if (adyacentes(pi1,pj1,pi2,pj2,one,half)        &
                         .or. adyacentes(pi1,pj1,pi2,pj2,half,one)  &
                         .or. adyacentes(pi1,pj1,pi2,pj2,one,zero)  &
                         .or. adyacentes(pi1,pj1,pi2,pj2,zero,one)  &
                         .or. adyacentes(pi1,pj1,pi2,pj2,half,zero) &
                         .or. adyacentes(pi1,pj1,pi2,pj2,zero,half) &
                         .or. adyacentes(pi1,pj1,pi2,pj2,one,one)) then
                         hallado = .true.
                         puntoactual = puntoactual+salto
                         punto(puntoactual) = k
                         puntoinicial = min(puntoinicial,puntoactual)
                         puntofinal = max(puntofinal,puntoactual)
                         lc(k) = .true.
                         pi1 = pi2
                         pj1 = pj2
                      endif
                   else

                      ! none of the points is in the grid. points are
                      si1 = zero
                      sj1 = zero
                      if (dabs(pi1-dint(pi1)).le.eps) si1 = one
                      if (dabs(pj1-dint(pj1)).le.eps) sj1 = one
                      if (adyacentes(pi1,pj1,pi2,pj2,half,half)       &
                         .or. adyacentes(pi1,pj1,pi2,pj2,si1,sj1)) then
                         hallado = .true.
                         puntoactual = puntoactual+salto
                         punto(puntoactual) = k
                         puntoinicial = min(puntoinicial,puntoactual)
                         puntofinal = max(puntofinal,puntoactual)
                         lc(k) = .true.
                         pi1 = pi2
                         pj1 = pj2
                      endif
                   endif
                endif
             enddo
          enddo

          ! no adyacent point found. Probe if it is a closed or open curve
          if ((puntofinal-puntoinicial+1).gt.2) then
             pi0 = pic(punto(puntoinicial))
             pj0 = pjc(punto(puntoinicial))
             malla0 = (dabs(pi0-dint(pi0)).le.eps) .and.              &
                &                (dabs(pj0-dint(pj0)).le.eps)
             pi1 = pic(punto(puntofinal))
             pj1 = pjc(punto(puntofinal))
             malla1 = (dabs(pi1-dint(pi1)).le.eps) .and.              &
                &                (dabs(pj1-dint(pj1)).le.eps)
             if (malla0 .or. malla1) then
                if (adyacentes(pi0,pj0,pi1,pj1,one,half)              &
                   .or. adyacentes(pi0,pj0,pi1,pj1,half,one)           &
                   .or. adyacentes(pi0,pj0,pi1,pj1,one,zero)           &
                   .or. adyacentes(pi0,pj0,pi1,pj1,zero,one)           &
                   .or. adyacentes(pi0,pj0,pi1,pj1,half,zero)          &
                   .or. adyacentes(pi0,pj0,pi1,pj1,zero,half)          &
                   .or. adyacentes(pi0,pj0,pi1,pj1,one,one)) then
                   cerrada = .true.
                endif
             else
                si1 = zero
                sj1 = zero
                if (dabs(pi1-dint(pi1)).le.eps) si1 = one
                if (dabs(pj1-dint(pj1)).le.eps) sj1 = one
                if (adyacentes(pi0,pj0,pi1,pj1,half,half)             &
                   .or. adyacentes(pi0,pj0,pi1,pj1,si1,sj1)) then
                   cerrada = .true.
                endif
             endif
          endif

          ! closed
          if (cerrada) then
             incompleta = .false.
             puntofinal = puntofinal+1
             punto(puntofinal) = punto(puntoinicial)
             ! open -> back to the first point to explore the other branch
          else if (primerabusqueda) then
             primerabusqueda = .false.
             puntoactual = 0
             salto = -1
             pi1 = pic(punto(puntoactual))
             pj1 = pjc(punto(puntoactual))
             ! curve is complete
          else
             incompleta = .false.
          endif
       enddo

       ! write
       nptoscurva = puntofinal-puntoinicial+1
       do k = puntoinicial, puntofinal
          x(k-puntoinicial+1) = xc(punto(k))
          y(k-puntoinicial+1) = yc(punto(k))
       enddo

       ! transform to non-orthogonal coordinates
       salpha = sqrt(1d0-calpha**2)
       x = x + y * calpha
       y = y * salpha

       call linea(x,y,nptoscurva,luw,ziso)

       ! update number of points
       if (cerrada) nptoscurva = nptoscurva-1
       nptosestudiados = nptosestudiados+nptoscurva
    enddo

  end subroutine ordenarpuntos

  !> Write (x(n),y(n)) curve in luw.
  subroutine linea(x,y,n,luw,ziso)
    use global, only: dunit0, iunit
    use tools_io, only: string
    integer, intent(in) :: n
    real*8, dimension(n), intent(in) :: x, y
    integer, intent(in) :: luw
    real*8, intent(in) :: ziso

    integer :: i

    write (luw,*)
    write (luw,'("# z = ",A)') string(ziso,'e',20,14)
    do i = 1, n
       write (luw,20) x(i) * dunit0(iunit), y(i) * dunit0(iunit)
    enddo
20  format (1p, 2(" ",e15.8))

  end subroutine linea

  !> Plot of the gradient vector field in the plane defined
  !> by the vectors (r1-r0) & (r2-r0).
  subroutine plotvec(r0,r1,r2,autocheck,udat)
    use systemmod, only: sy
    use global, only: dunit0, iunit, iunitname0, prunedist, gcpchange
    use tools_math, only: cross, matinv
    use tools_io, only: uout, string, ioj_right, ioj_left
    use param, only: pi
    use types, only: scalar_value, gpathp
    integer, intent(in) :: udat
    logical, intent(in) :: autocheck
    real*8, dimension(3), intent(in) :: r0, r1, r2

    integer :: nptf, i, j, iorig, up1d, up2d, ntotpts, nindex, ntype
    real*8  :: xstart(3), phii, u1, v1, u, v
    real*8  :: r012, v1d(3), v2da(3), v2db(3), xo0(3), xo1(3), xo2(3)
    real*8  :: xtemp(3), c1coef, c2coef, plen
    integer :: ier
    type(scalar_value) :: res
    type(gpathp), allocatable :: xpath(:)

    ! plane metrics
    rp0 = sy%c%x2c(r0)
    rp1 = sy%c%x2c(r1)
    rp2 = sy%c%x2c(r2)
    r01 = 0d0
    r02 = 0d0
    r012 = 0d0
    do i = 1, 3
       rp01(i) = rp1(i) - rp0(i)
       rp02(i) = rp2(i) - rp0(i)
       r01 = r01 + rp01(i) * rp01(i)
       r02 = r02 + rp02(i) * rp02(i)
       r012 = r012 + rp01(i) * rp02(i)
    enddo
    r01 = sqrt(r01)
    r02 = sqrt(r02)
    cosalfa = r012 / (r01*r02)
    a012 = acos(cosalfa)
    sinalfa = sin(a012)
    cotgalfa = cosalfa / sinalfa

    ! normal vector
    rpn(1:3) = cross(rp01,rp02)
    rpn(4) = -dot_product(rpn(1:3),rp0)
    rpn = rpn / (r01*r02)

    ! plane to cartesian
    amat(:,1) = rp01
    amat(:,2) = rp02
    amat(:,3) = rpn(1:3)

    ! cartesian to plane
    bmat = amat
    call matinv(bmat,3)

    write (uout,'("* Plot of the gradient vector field in the plane:")')
    write (uout,'("    r = r0 + u * (r1 - r0) + v * (r2 - r0)")')
    write (uout,'("  where the parametric coordinates u and v go from 0 to 1.")')
    if (.not.sy%c%ismolecule) then
       write (uout,'("+ Crystal coordinates of r0: ",3(A,"  "))') (string(r0(j),'f',12,6,ioj_right),j=1,3)
       write (uout,'("+ Crystal coordinates of r1: ",3(A,"  "))') (string(r1(j),'f',12,6,ioj_right),j=1,3)
       write (uout,'("+ Crystal coordinates of r2: ",3(A,"  "))') (string(r2(j),'f',12,6,ioj_right),j=1,3)
       write (uout,'("+ Cartesian coordinates of r0: ",3(A," "))') (string(rp0(j),'f',12,6,ioj_right),j=1,3)
       write (uout,'("+ Cartesian coordinates of r1: ",3(A," "))') (string(rp1(j),'f',12,6,ioj_right),j=1,3)
       write (uout,'("+ Cartesian coordinates of r2: ",3(A," "))') (string(rp2(j),'f',12,6,ioj_right),j=1,3)
    else
       xo0 = (sy%c%x2c(r0) + sy%c%molx0) * dunit0(iunit)
       xo1 = (sy%c%x2c(r1) + sy%c%molx0) * dunit0(iunit)
       xo2 = (sy%c%x2c(r2) + sy%c%molx0) * dunit0(iunit)
       write (uout,'("+ Coordinates of r0: ",3(A," "))') (string(xo0(j),'f',12,6,ioj_right),j=1,3)
       write (uout,'("+ Coordinates of r1: ",3(A," "))') (string(xo1(j),'f',12,6,ioj_right),j=1,3)
       write (uout,'("+ Coordinates of r2: ",3(A," "))') (string(xo2(j),'f',12,6,ioj_right),j=1,3)
    endif
    ! Check the in-plane CPs
    if (autocheck) call autochk(rp0)

    !.Run over the points defined as origins:
    write (uout,'("+ List of critical points that act as in-plane gradient path generators")')
    write (uout,'("# i       xcrys        ycrys        zcrys    iatr   up down     xplane       yplane")')
    do iorig = 1, norig
       ! calculate in-plane coordinates
       xtemp = sy%c%x2c(grpx(:,iorig))
       xtemp = xtemp - rp0
       xtemp = matmul(bmat,xtemp)
       u = xtemp(1)*r01 + xtemp(2)*r02*cosalfa
       v = xtemp(2)*r02*sinalfa

       xtemp = grpx(:,iorig)
       if (sy%c%ismolecule) &
          xtemp = (sy%c%x2c(xtemp) + sy%c%molx0) * dunit0(iunit)
       write (uout,'(99(A,"  "))') string(iorig,length=5,justify=ioj_left), &
          (string(xtemp(j),'f',decimal=6,length=11,justify=4),j=1,3),&
          string(grpatr(iorig),length=3,justify=ioj_right), &
          string(grpup(iorig),length=3,justify=ioj_right), &
          string(grpdwn(iorig),length=3,justify=ioj_right), &
          string(u,'f',decimal=6,length=11,justify=4), string(v,'f',decimal=6,length=11,justify=4)
    enddo
    write (uout,*)

    ! write the data file header
    write (udat,'("# Gradient path information")')
    write (udat,'("# u(",A,") v(",A,") x(",A,") y(",A,") z(",A,") color")') (trim(iunitname0(iunit)),j=1,5)

    write (uout,'("+ List of gradient paths traced")')
    write (uout,'("# i       xcrys        ycrys        zcrys        type    up down    pts")')
    do iorig = 1, norig
       ntotpts = 0
       grpx(:,iorig) = sy%c%x2c(grpx(:,iorig))

       if (grpatr(iorig) .eq. 1) then
          ! 3D attraction or repulsion critical points:
          do i = 1, grpup(iorig)
             phii = (i-1) * 2d0 * pi / max(grpup(iorig)-1,1)
             u1 = gcpchange * sin(phii)
             v1 = gcpchange * cos(phii)
             u = (u1 - v1 * cotgalfa) / r01
             v = v1 / (r02 * sinalfa)
             xstart = grpx(:,iorig) + u * rp01 + v * rp02
             call sy%f(sy%iref)%gradient(xstart,+1,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
             nptf = size(xpath,1)
             call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
             ntotpts = ntotpts + nptf
          enddo
          do i = 1, grpdwn(iorig)
             phii = (i-1) * 2d0 * pi / max(grpdwn(iorig)-1,1)
             u1 = gcpchange * sin(phii)
             v1 = gcpchange * cos(phii)
             u = (u1 - v1 * cotgalfa) / r01
             v = v1 / (r02 * sinalfa)
             xstart = grpx(:,iorig) + u * rp01 + v * rp02
             call sy%f(sy%iref)%gradient(xstart,-1,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
             nptf = size(xpath,1)
             call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
             ntotpts = ntotpts + nptf
          enddo
          nindex = 3
          if (grpup(iorig) > 0) then
             ntype = 3
          else
             ntype = -3
          end if

          xtemp = (grpx(:,iorig) + sy%c%molx0) * dunit0(iunit)

          write (uout,'(4(A,"  "),"(",A,","A,") ",3(A,"  "))') &
             string(iorig,length=5,justify=ioj_left), &
             (string(xtemp(j),'f',decimal=6,length=11,justify=4),j=1,3),&
             string(nindex,length=3,justify=ioj_right),&
             string(ntype,length=3,justify=ioj_right),&
             string(grpup(iorig),length=3,justify=ioj_right),&
             string(grpdwn(iorig),length=3,justify=ioj_right),&
             string(ntotpts,length=5,justify=ioj_right)
       else
          ! A (3,-1) or (3,+3) critical point:
          xstart = grpx(:,iorig)
          call sy%f(sy%iref)%grd(xstart,2,res,periodic=.not.sy%c%ismolecule)
          if (res%r .eq. 3) then
             if (res%s .eq. -1) then
                up1d = +1
                up2d = -1
                if (grpup(iorig)  .eq. 0) up1d = 0
                if (grpdwn(iorig) .eq. 0) up2d = 0
                v2da(1) = res%hfevec(1,1)
                v2da(2) = res%hfevec(2,1)
                v2da(3) = res%hfevec(3,1)
                v2db(1) = res%hfevec(1,2)
                v2db(2) = res%hfevec(2,2)
                v2db(3) = res%hfevec(3,2)
                v1d(1)  = res%hfevec(1,3)
                v1d(2)  = res%hfevec(2,3)
                v1d(3)  = res%hfevec(3,3)
             else if (res%s .eq. +1) then
                up1d = -1
                up2d = +1
                if (grpup(iorig)  .eq. 0) up2d = 0
                if (grpdwn(iorig) .eq. 0) up1d = 0
                v1d(1)  = res%hfevec(1,1)
                v1d(2)  = res%hfevec(2,1)
                v1d(3)  = res%hfevec(3,1)
                v2da(1) = res%hfevec(1,2)
                v2da(2) = res%hfevec(2,2)
                v2da(3) = res%hfevec(3,2)
                v2db(1) = res%hfevec(1,3)
                v2db(2) = res%hfevec(2,3)
                v2db(3) = res%hfevec(3,3)
             else
                up1d = 0
                up2d = 0
             endif
             if (up2d .ne. 0) then
                !.2D walk:
                c1coef = dot_product(v2db,rpn(1:3))
                c2coef = -dot_product(v2da,rpn(1:3))

                if (abs(c1coef) < 1d-10 .and. abs(c2coef) < 1d-10) then
                   ! the plot plane is coincident with the ias tangent plane
                   ! use v2da.
                   xstart = grpx(:,iorig) + gcpchange * v2da
                   call sy%f(sy%iref)%gradient(xstart,-1,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
                   nptf = size(xpath,1)
                   call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
                   ntotpts = ntotpts + nptf
                   xstart = grpx(:,iorig) - gcpchange * v2da
                   call sy%f(sy%iref)%gradient(xstart,up2d,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
                   nptf = size(xpath,1)
                   call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
                   ntotpts = ntotpts + nptf
                else
                   xtemp = c1coef * v2da + c2coef * v2db
                   xstart = grpx(:,iorig) + gcpchange * xtemp
                   call sy%f(sy%iref)%gradient(xstart,up2d,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
                   nptf = size(xpath,1)
                   call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
                   ntotpts = ntotpts + nptf
                   xstart = grpx(:,iorig) - gcpchange * xtemp
                   call sy%f(sy%iref)%gradient(xstart,up2d,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
                   nptf = size(xpath,1)
                   call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
                   ntotpts = ntotpts + nptf
                end if

             endif
             if (up1d .ne. 0) then
                !
                !.................1D walk:
                !
                xstart = grpx(:,iorig) + gcpchange * v1d
                call sy%f(sy%iref)%gradient(xstart,up1d,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
                nptf = size(xpath,1)
                call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
                xstart = grpx(:,iorig) - gcpchange * v1d
                call sy%f(sy%iref)%gradient(xstart,up1d,nptf,ier,.false.,plen,xpath,prunedist,pathini=grpx(:,iorig))
                nptf = size(xpath,1)
                call wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
                ntotpts = ntotpts + nptf
             endif
          endif

          xtemp = (grpx(:,iorig) + sy%c%molx0) * dunit0(iunit)

          write (uout,'(4(A,"  "),"(",A,","A,") ",3(A,"  "))') &
             string(iorig,length=5,justify=ioj_left), &
             (string(xtemp(j),'f',decimal=6,length=11,justify=4),j=1,3),&
             string(res%r,length=3,justify=ioj_right),&
             string(res%s,length=3,justify=ioj_right),&
             string(grpup(iorig),length=3,justify=ioj_right),&
             string(grpdwn(iorig),length=3,justify=ioj_right),&
             string(ntotpts,length=5,justify=ioj_right)
       endif
    enddo
    write (uout,*)

  end subroutine plotvec

  !> Write the gradient path to the udat logical unit.
  subroutine wrtpath(xpath,nptf,udat,rp0,r01,r02,cosalfa,sinalfa)
    use systemmod, only: sy
    use global, only: prunedist, dunit0, iunit
    use types, only: gpathp
    use param, only: jmlcol, icrd_crys
    type(gpathp) :: xpath(*)
    integer, intent(in) :: nptf
    integer, intent(in) :: udat
    real*8, intent(in) :: rp0(3), r01, r02, cosalfa, sinalfa

    integer :: i, j, nid1, nid2, iz, rgb(3)
    real*8 :: xxx, yyy, zzz, u, v, h, uort, vort
    real*8 :: dist1, dist2, dd
    logical :: wasblank

    ! identify the endpoints
    nid1 = sy%c%identify_atom(xpath(1)%x,icrd_crys,dist=dist1,distmax=1.1d0*prunedist)
    nid2 = sy%c%identify_atom(xpath(nptf)%x,icrd_crys,dist=dist2,distmax=1.1d0*prunedist)
    rgb = (/0,0,0/)
    if (nid1 > 0 .and. (dist1 < dist2 .or. nid2 == 0)) then
       iz = sy%c%spc(sy%c%atcel(nid1)%is)%z
       if (iz /= 1) rgb = jmlcol(:,iz)
    elseif (nid2 > 0 .and. (dist2 < dist1 .or. nid1 == 0)) then
       iz = sy%c%spc(sy%c%atcel(nid2)%is)%z
       if (iz /= 1) rgb = jmlcol(:,iz)
    endif

    write (udat,*)
    write (udat,*)
    wasblank = .true.
    dd = dunit0(iunit)
    do i = 1, nptf
       !.transform the point to the plotting plane coordinates:
       xxx = xpath(i)%r(1) - rp0(1)
       yyy = xpath(i)%r(2) - rp0(2)
       zzz = xpath(i)%r(3) - rp0(3)
       u = bmat(1,1)*xxx + bmat(1,2)*yyy + bmat(1,3)*zzz
       v = bmat(2,1)*xxx + bmat(2,2)*yyy + bmat(2,3)*zzz
       h = bmat(3,1)*xxx + bmat(3,2)*yyy + bmat(3,3)*zzz

       !.clip the line if it is out of the plotting area:
       if (u < 0d0 .or. u > 1d0 .or. v < 0d0 .or. v > 1d0) then
          ! print the last point so that gnuplot interpolates to the edge
          if (i-1 >= 1) then
             uort = u*r01 + v*r02*cosalfa
             vort = v*r02*sinalfa
             if (abs(h) < grphcutoff .or. grpproj > 0) then
                write (udat,200) uort*dd, vort*dd, (xpath(i)%r(j)*dd, j = 1, 3), (rgb(1) * 256 + rgb(2)) * 256 + rgb(3)
             end if
          end if
          if (.not.wasblank) write (udat,*)
          wasblank = .true.
       else
          wasblank = .false.
          uort = u*r01 + v*r02*cosalfa
          vort = v*r02*sinalfa
          if (abs(h) < grphcutoff .or. grpproj > 0) then
             write (udat,200) uort*dd, vort*dd, (xpath(i)%r(j)*dd, j = 1, 3), (rgb(1) * 256 + rgb(2)) * 256 + rgb(3)
          end if
       endif
    enddo

200 format (2f15.9, 1p, 3e15.6, I10)

  end subroutine wrtpath

  !> Check a user-entered collection of points to see if they
  !> are critical points, remove repeated points, and check which
  !> equivalents (within the main cell) lie on the plotting plane.
  subroutine autochk(rp0)
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_left, ioj_right, faterr, ferror
    use param, only: one
    use types, only: scalar_value
    integer :: i, j, k, l, ncopies
    integer :: iorde(2*indmax+1), indcell(3,(2*indmax+1)**3), iii, inum
    real*8  :: xp(3)
    real*8  :: x0(3), x1(3), uu, vv, hh, hmin
    type(scalar_value) :: res
    real*8 :: rp0(3)

    !
    write (uout,'("+ List of candidate in-plane CPs")')
    write (uout,'("# cp       x           y           z")')
    do i = 1, newncriticp
       newcriticp(:,i) = newcriticp(:,i) - floor(newcriticp(:,i))
       write (uout,'("  ",4(A,"  "))') string(i,length=3,justify=ioj_left), &
          (string(newcriticp(j,i),'f',10,6,ioj_right),j=1,3)
    enddo
    write (uout,*)

    inum = 2*indmax+1
    do i=1,inum
       iorde(i) = i-1-indmax
    end do

    iii = 0
    do i = 1, inum
       do j = 1, inum
          do k = 1, inum
             iii = iii + 1
             indcell(1,iii) = iorde(i)
             indcell(2,iii) = iorde(j)
             indcell(3,iii) = iorde(k)
          enddo
       enddo
    enddo

    write (uout,'("+ Pruning actions on the in-plane CP list")')
    do i = 1, newncriticp
       xp = newcriticp(:,i)

       !.discard repeated points
       do j = i+1, newncriticp
          if (sy%c%are_close(newcriticp(:,j),xp,epsf)) then
             write (uout,'("  CP ",A," is equivalent to ",A," -> Rejected!")') string(j), string(i)
             cycle
          endif
       enddo

       !.Get the properties
       xp = sy%c%x2c(newcriticp(:,i))
       call sy%f(sy%iref)%grd(xp,2,res,periodic=.not.sy%c%ismolecule)
       if (.not.res%valid) then
          write (uout,'("  CP ",A," has an invalid value -> Rejected!")') string(i)
          cycle
       elseif (res%gfmod > grpcpeps) then
          write (uout,'("  CP ",A," has a large gradient: ",A," -> Rejected!")') &
             string(i), string(res%gfmod,'e',decimal=6)
          cycle
       else
          newtypcrit(i) = res%s
       endif

       !.Determine the points to be added as gradient path origins:
       ncopies = 0
       hmin = 1d30
       do l = 1, (inum)**3
          xp = sy%c%x2c(newcriticp(:,i) + indcell(:,l))

          !.transform the point to the plotting plane coordinates:
          x0 = xp - rp0
          x1 = matmul(bmat,x0)
          uu = x1(1)
          vv = x1(2)
          hh = x1(3)

          !.clip if out of the plotting area
          hmin = min(hmin,abs(hh))
          if (uu >= -epsdis+(one-outcpx) .and. uu <= outcpx+epsdis .and.&
             vv >= -epsdis+(one-outcpy) .and. vv <= outcpy+epsdis .and.&
             abs(hh) <= RHOP_Hmax) then
             ncopies = ncopies + 1
             if (norig .ge. MORIG) then
                call ferror('autochk', 'Too many origins. Increase MORIG',faterr,syntax=.true.)
                return
             endif
             norig = norig + 1
             grpx(:,norig) = newcriticp(:,i) + indcell(:,l)
             if (cpup(i) > 0 .or. cpdn(i) > 0) then
                if (newtypcrit(i) == -3 .or. newtypcrit(i) == +3) then
                   grpatr(norig) = 1
                else
                   grpatr(norig) = 0
                endif
                grpup(norig) = cpup(i)
                grpdwn(norig) = cpdn(i)
             else
                if (newtypcrit(i) == -3) then
                   grpatr(norig) = 1
                   grpup(norig)  = 0
                   grpdwn(norig) = 36
                else if (newtypcrit(i) == +3) then
                   grpatr(norig) = 1
                   grpup(norig)  = 36
                   grpdwn(norig) = 0
                else
                   grpatr(norig) = 0
                   grpup(norig)  = 2
                   grpdwn(norig) = 2
                endif
             end if
          endif
       enddo
       write (uout,'("  CP ",A," created ",A," copies in plane (hmin = ",A,")")') &
          string(i), string(ncopies), string(hmin,'f',decimal=6)
    end do
    write (uout,*)

  end subroutine autochk

  !> Writes the "label" file to rootname-label.gnu. The labels file
  !> contains the list of critical points contained in the plot
  !> plane, ready to be read in gnuplot.
  subroutine write_fichlabel(rootname)
    use global, only: dunit0, iunit
    use systemmod, only: sy
    use tools_math, only: det3
    use tools_io, only: uout, string, fopen_write, fclose, nameguess
    use param, only: one
    character*(*), intent(in) :: rootname

    character(len=:), allocatable :: fichlabel
    integer :: i, j, k, iii
    integer :: inum, iorde(2*indmax+1), indcell(3,(2*indmax+1)**3)
    integer :: lul
    real*8 :: xp(3), xxx, yyy, zzz, uu, vv, hh, u, v
    character*2 :: cpletter

    fichlabel = trim(rootname) // "-label.gnu"

    write (uout,'("* Name of the labels file: ",a)') string(fichlabel)
    write (uout,*)

    ! Create labels file
    lul = fopen_write(fichlabel)
    inum = 2*indmax+1
    do i=1,inum
       iorde(i) = i-1-indmax
    end do
    iii = 0
    do i = 1, inum
       do j = 1, inum
          do k = 1, inum
             iii = iii + 1
             indcell(1,iii) = iorde(i)
             indcell(2,iii) = iorde(j)
             indcell(3,iii) = iorde(k)
          enddo
       enddo
    enddo

    do i = 1, sy%f(sy%iref)%ncpcel
       do j = 1, (inum)**3
          xp = sy%c%x2c(sy%f(sy%iref)%cpcel(i)%x + indcell(:,j))
          xxx = xp(1) - rp0(1)
          yyy = xp(2) - rp0(2)
          zzz = xp(3) - rp0(3)
          uu = bmat(1,1)*xxx + bmat(1,2)*yyy + bmat(1,3)*zzz
          vv = bmat(2,1)*xxx + bmat(2,2)*yyy + bmat(2,3)*zzz
          hh = bmat(3,1)*xxx + bmat(3,2)*yyy + bmat(3,3)*zzz
          u = uu*r01+vv*r02*cosalfa
          v = vv*r02*sinalfa

          if (uu.ge.-epsdis .and. uu.le.one+epsdis .and. &
             vv.ge.-epsdis .and. vv.le.one+epsdis .and. &
             abs(hh).le.RHOP_Hmax) then
             ! assign cp letter, check if it is a nucleus
             select case (sy%f(sy%iref)%cpcel(i)%typ)
             case (3)
                cpletter = "c"
             case (1)
                cpletter = "r"
             case (-1)
                cpletter = "b"
             case (-3)
                if (sy%f(sy%iref)%cpcel(i)%isnuc) then
                   cpletter = nameguess(sy%c%spc(sy%c%at(sy%f(sy%iref)%cpcel(i)%idx)%is)%z,.true.)
                else
                   cpletter = "n"
                endif
             end select
             write (lul,'(3a,f12.6,a,f12.6,a)') 'set label "',trim(cpletter),'" at ',u * dunit0(iunit),&
                ',',v * dunit0(iunit),' center front'
          end if
       end do
    end do
    call fclose(lul)

  end subroutine write_fichlabel

  !> Writes the gnuplot file rootname.gnu. If dolabels, write the
  !> labels file too. If docontour, write the .iso (postivie iso
  !> -values) and the .neg.iso (negative iso-values). If dograds,
  !> write the file containing the gradient paths.
  subroutine write_fichgnu(rootname,dolabels,docontour,dograds)
    use global, only: iunitname0, iunit, dunit0
    use tools_io, only: uout, string, fopen_write, string, fclose
    character*(*), intent(in) :: rootname
    logical, intent(in) :: dolabels, docontour, dograds

    integer :: lun
    character(len=:), allocatable :: fichgnu, fichlabel, fichiso, fichiso1, fichgrd
    character(len=:), allocatable :: swri
    real*8 :: dd

    dd = dunit0(iunit)
    fichgnu = trim(rootname) // '.gnu'
    fichlabel = trim(rootname) // "-label.gnu"
    fichiso = trim(rootname) // ".iso"
    fichiso1 = trim(rootname) // ".neg.iso"
    fichgrd = trim(rootname) // ".dat"

    write (uout,'("* Name of the gnuplot script file: ",a/)') string(fichgnu)

    lun=fopen_write(fichgnu)

    write (lun,*) "set terminal postscript eps color enhanced size 10,7 'Helvetica' 40"
    write (lun,*) "set style line 1 lt 1 lc rgbcolor variable lw 4"
    write (lun,*) 'set style line 2 lt 1 lc rgbcolor "#007700"'
    write (lun,*) 'set style line 3 lt 4 lc rgbcolor "#0000ff"'
    write (lun,*) 'set output "' // trim(rootname) // '.eps"'
    write (lun,*) 'set size ratio -1'
    write (lun,*) 'unset key'
    ! title
    write (lun,*) 'set xlabel "x (',trim(iunitname0(iunit)),' )"'
    write (lun,*) 'set ylabel "y (',trim(iunitname0(iunit)),' )"'
    ! ranges
    write (lun,18) min(0d0,r02*cosalfa)*dd, max(r01,r01+r02*cosalfa)*dd, 0d0, r02*sinalfa*dd
    ! labels
    if (dolabels) &
       write (lun,*) 'load "',string(fichlabel),'"'
    ! contours
    swri = ""
    if (docontour) then
       swri = trim(swri) // """" // string(fichiso) // """ with lines ls 2"
       if (isneg) &
          swri = trim(swri) // ", """ // string(fichiso1) // """ with lines ls 3"
       if (dograds) swri = trim(swri) // ", "
    end if
    if (dograds) then
       swri = swri // """" // string(fichgrd) // """ u 1:2:6 with lines ls 1"
    end if

    write (lun,'("plot ",A)') swri
    call fclose(lun)

18  format (' set xrange [',f7.3,':',f7.3,'] '/' set yrange [',f7.3,':',f7.3,']')

  end subroutine write_fichgnu

end submodule proc
