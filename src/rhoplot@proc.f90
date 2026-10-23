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
  ! subroutine hallarpuntos(ff,zc,x,y,nx,ny,np,xc,yc,pic,pjc)
  ! subroutine ordenarpuntos(np,xc,yc,pic,pjc,calpha,iso,nline,npt,lptr,xy,liso)
  ! subroutine linea (x,y,n,luw,ziso)
  ! subroutine plotvec(r0,r1,r2,autocheck,udat)
  ! subroutine plotvec(r0,r1,r2,autocheck,udat)
  ! subroutine autochk(rp0)
  ! subroutine write_fichlabel(rootname)
  ! subroutine write_fichgnu(rootname,dolabels,docontour,dograds)

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
          call sy%propty(ifi,x0,res,.false.,doall)
       else
          rdum = sy%eval(expr,errmsg,xx)
          if (len_trim(errmsg) > 0) then
             call ferror('critic','error in expression ' // trim(expr) // ': ' // errmsg,faterr,line,syntax=.true.)
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
    use tools_io, only: ferror, faterr, lgetword, equal, getword, equal,&
       isexpression_or_word, fopen_write, uout, string, fclose
    character*(*), intent(in) :: line

    integer :: lp, lp2, nti, id, luout, np
    real*8 :: x0(3), x1(3), xp(3), dist, xout(3)
    character(len=:), allocatable :: word, outfile, prop, expr, errmsg
    logical :: ok
    integer :: i, j
    real*8, allocatable :: rhoout(:), lapout(:)

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
    if (.not.allocated(expr)) expr = ""
    call rhoplot_line_values(sy,x0,x1,np,id,expr,nti,rhoout,lapout,errmsg)
    if (len_trim(errmsg) > 0) &
       call ferror('rhoplot_line',errmsg,faterr)

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

  !> The field request that evaluates property nti of a field: the
  !> value (0), a first derivative (1 to 4: the gradient components and
  !> its norm) or a second derivative (the rest; as in LINE, PLANE and
  !> CUBE).
  pure module function rhoplot_nti_request(nti) result(request)
    use types, only: field_evaluation_avail
    integer, intent(in) :: nti
    type(field_evaluation_avail) :: request

    if (nti == 0) then
       call request%field_only()
    elseif (nti >= 1 .and. nti <= 4) then
       call request%field_nder1()
    else
       call request%field_nder2()
    end if

  end function rhoplot_nti_request

  !> Calculate the values of a field along the line from x0 to x1
  !> (Cartesian, bohr, cell frame of system s) at np points. If id
  !> >= 0, use field id and return its value in f and the property
  !> nti in g (0 = f, 1-3 = gradient components, 4 = gradient norm,
  !> 5-10 = Hessian xx, xy, xz, yy, yz, zz, 11 = Laplacian);
  !> otherwise, evaluate the expression expr and return its value in
  !> both f and g. Returns a non-empty errmsg if the evaluation
  !> failed.
  module subroutine rhoplot_line_values(s,x0,x1,np,id,expr,nti,f,g,errmsg)
    use systemmod, only: system
    use types, only: scalar_value, vstring, field_evaluation_avail
    type(system), intent(inout) :: s
    real*8, intent(in) :: x0(3), x1(3)
    integer, intent(in) :: np, id, nti
    character*(*), intent(in) :: expr
    real*8, intent(out) :: f(np), g(np)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i
    real*8 :: xp(3), rhopt, lappt
    type(scalar_value) :: res
    type(vstring) :: lerrs
    type(field_evaluation_avail) :: request
    character(len=:), allocatable :: lerrmsg

    request = rhoplot_nti_request(nti)

    ! calculate the line; the private error string is in a vstring:
    ! gfortran mishandles private deferred-length character variables
    lerrmsg = ""
    f = 0d0
    g = 0d0
    !$omp parallel do private(xp,res,rhopt,lappt,lerrs) schedule(dynamic)
    do i=1,np
       xp = x0 + (x1 - x0) * real(i-1,8) / real(max(np-1,1),8)
       lerrs%s = ""
       if (id >= 0) then
          call s%f(id)%grd(xp,request,res,periodic=.not.s%c%ismolecule)

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
          case default
             lappt = 0d0
          end select
       else
          rhopt = s%eval(expr,lerrs%s,xp)
          lappt = rhopt
       end if

       f(i) = rhopt
       g(i) = lappt
       if (len_trim(lerrs%s) > 0) then
          !$omp critical (iowrite)
          lerrmsg = lerrs%s
          !$omp end critical (iowrite)
       end if
    enddo
    !$omp end parallel do
    errmsg = lerrmsg

  end subroutine rhoplot_line_values

  ! Calculate properties on a 3d cube
  module subroutine rhoplot_cube(line)
    use systemmod, only: sy
    use crystalmod, only: crystal
    use global, only: eval_next, dunit0, iunit, fileroot, abort_requested
    use tools_io, only: lgetword, faterr, ferror, equal, getword, &
       isexpression_or_word, uout, string, isinteger
    use param, only: eye
    character*(*), intent(in) :: line

    integer :: lp, nti, id, nn(3)
    real*8 :: x0(3), x1(3), x00(3,3)
    real*8 :: rgr, xd(3,3)
    integer :: lp2
    character(len=:), allocatable :: word, outfile, expr, errmsg
    logical :: ok, doortho
    integer :: i, ibnd, ik, inr(3), ispin, ipart
    real*8, allocatable :: lf(:,:,:)
    logical :: dogrid, useexpr, doheader, iscell
    integer :: outform, ishift(3), dopsi
    integer :: ia(3), na(3)
    logical :: iused(3)
    type(crystal), pointer :: cr
    type(crystal), target :: nc

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
    dopsi = cubepsi_none
    doheader = .false.
    iscell = .false.
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
       iscell = .true.
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
          dopsi = cubepsi_mlwf
       elseif (equal(word,'unk')) then
          dopsi = cubepsi_unk
       elseif (equal(word,'psink')) then
          dopsi = cubepsi_psink
       else
          dopsi = cubepsi_wannier
       end if
       ok = isinteger(ibnd,line,lp)
       if (dopsi == cubepsi_unk.or.dopsi == cubepsi_psink) then
          ok = ok .and. isinteger(ik,line,lp)
       end if
       if (dopsi /= cubepsi_unk) then
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
    else
       call ferror('rhoplot_cube','wrong CUBE syntax',faterr,line,syntax=.true.)
       return
    endif

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
          ! calculate the axis lengths to convert resolution into number of points
          do i = 1, 3
             nn(i) = nint(norm2(sy%c%x2c(xd(:,i))) / rgr) + 1
          end do
       else
          do i = 1, 3
             nn(i) = max(nn(i),2)
          end do
       end if
    end if

    ! the crystal pointer defaults to the current structure
    cr => sy%c

    ! read additional options
    ishift = 0
    nti = nti_none
    id = sy%iref
    useexpr = .false.
    expr = ""
    outform = cubefmt_cube
    outfile = trim(fileroot) // ".cube"
    ia = -1
    na = -1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'file')) then
          outfile = getword(line,lp)
          if (len_trim(outfile) < 1) then
             call ferror('rhoplot_cube','file name not found',faterr,line,syntax=.true.)
             return
          end if
          outform = rhoplot_cube_outform(outfile)
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
       elseif (equal(word,'gx')) then
          nti = nti_gx
       else if (equal(word,'gy')) then
          nti = nti_gy
       else if (equal(word,'gz')) then
          nti = nti_gz
       else if (equal(word,'gmod')) then
          nti = nti_gmod
       else if (equal(word,'hxx')) then
          nti = nti_hxx
       else if (equal(word,'hxy') .or. equal(word,'hyx')) then
          nti = nti_hxy
       else if (equal(word,'hxz') .or. equal(word,'hzx')) then
          nti = nti_hxz
       else if (equal(word,'hyy')) then
          nti = nti_hyy
       else if (equal(word,'hyz') .or. equal(word,'hzy')) then
          nti = nti_hyz
       else if (equal(word,'hzz')) then
          nti = nti_hzz
       else if (equal(word,'lap')) then
          nti = nti_lap
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
       else if (equal(word,'newcell')) then
          iused = .false.
          do i = 1, 3
             word = lgetword(line,lp)
             if (equal(word,'a')) then
                ia(i) = 1
             elseif (equal(word,'b')) then
                ia(i) = 2
             elseif (equal(word,'c')) then
                ia(i) = 3
             else
                call ferror('rhoplot_cube','error in axes passed to NEWCELL option',faterr,line,syntax=.true.)
                return
             end if
             if (iused(ia(i))) then
                call ferror('rhoplot_cube','repeated axes in NEWCELL option',faterr,line,syntax=.true.)
                return
             end if
             iused(ia(i)) = .true.

             ok = isinteger(na(i),line,lp)
             if (.not.ok) then
                call ferror('rhoplot_cube','error reading axis repetition in NEWCELL',faterr,line,syntax=.true.)
                return
             end if
             if (na(i) < 1) then
                call ferror('rhoplot_cube','error reading axis repetition in NEWCELL',faterr,line,syntax=.true.)
                return
             end if
          end do

          ! build the corresponding supercell
          x00 = 0d0
          do i = 1, 3
             x00(ia(i),i) = na(i)
          end do
          nc = sy%c
          call nc%newcell(x00,noenv=.true.,errmsg=errmsg)
          if (len_trim(errmsg) > 0) then
             call ferror('rhoplot_cube',errmsg,faterr,line,syntax=.true.)
             return
          end if
          cr => nc
       else if (len_trim(word) > 0) then
          call ferror('rhoplot_cube','Unknown keyword in CUBE',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! the wavefunction keywords take only REAL, IMAG, and ABS, and
    ! only they do
    if (dopsi /= cubepsi_none) then
       if (nti /= nti_none .and. nti < nti_real) then
          call ferror('rhoplot_cube','MLWF/WANNIER/UNK/PSINK in CUBE accept only REAL, IMAG, or ABS',&
             faterr,line,syntax=.true.)
          return
       end if
    elseif (nti >= nti_real) then
       call ferror('rhoplot_cube','REAL, IMAG, and ABS in CUBE require MLWF, WANNIER, UNK, or PSINK',&
          faterr,line,syntax=.true.)
       return
    end if

    ! calculate the distances
    x0 = cr%x2c(x0)
    do i = 1, 3
       xd(:,i) = cr%x2c(xd(:,i))
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
       errmsg = rhoplot_cube_gridcheck(sy,id,dopsi,ibnd,any(na > 0))
       if (len(errmsg) > 0) then
          call ferror('rhoplot_cube',errmsg,faterr,syntax=.true.)
          return
       end if
       nn = sy%f(id)%grid%n
    else
       if (any(na > 0)) then
          call ferror('rhoplot_cube','NEWCELL can only be used with GRID',faterr,syntax=.true.)
          return
       end if
    end if

    ! convert to one-step vectors in each direction
    if (any(na > 0)) then
       do i = 1, 3
          xd(:,i) = xd(:,i) / real(nn(ia(i)) * na(i),8)
       end do
    else
       do i = 1, 3
          xd(:,i) = xd(:,i) / real(nn(i),8)
       end do
    end if

    ! write cube header
    write (uout,'("* CUBE written to file: ",A/)') string(outfile)
    if (doheader) then
       call rhoplot_cube_write(cr,outform,lf,outfile,.true.,xd,x0,ishift,errmsg)
       if (len(errmsg) > 0) call ferror("rhoplot_cube",errmsg,faterr)
       return
    end if

    ! calculate properties
    if (dogrid) then
       if (nti == nti_real) then
          ipart = cubepart_real
       elseif (nti == nti_imag) then
          ipart = cubepart_imag
       else
          ipart = cubepart_abs
       end if
       call rhoplot_cube_gridvalues(sy,id,ia,na,dopsi,ibnd,ik,inr,ispin,ipart,lf)
    else
       ! the kernel numbers the properties from 0 (f): nti_f - 1
       if (useexpr) id = -1
       call rhoplot_cube_values(sy,x0,xd,nn,id,expr,max(nti-1,0),iscell.and..not.doortho,lf,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('rhoplot_cube',errmsg,faterr)
    end if
    if (abort_requested) then
       write (uout,'("+ The calculation was cancelled: no file was written."/)')
       return
    end if

    ! cube body
    call rhoplot_cube_write(cr,outform,lf,outfile,.false.,xd,x0,ishift,errmsg)

  end subroutine rhoplot_cube

  !> Calculate the values of a field on the points x0 + ix*xd(:,1) +
  !> iy*xd(:,2) + iz*xd(:,3), ix = 0..nn(1)-1, ..., of a cube
  !> (Cartesian, bohr, cell frame of system s) in lf(nn(1),nn(2),nn(3)).
  !> If id >= 0, use field id and its property nti (0 = f, 1-3 =
  !> gradient components, 4 = gradient norm, 5-10 = Hessian xx, xy,
  !> xz, yy, yz, zz, 11 = Laplacian); otherwise, evaluate the
  !> expression expr. If oncell, the points are the nodes of a grid
  !> on the whole cell, and an expression of grid fields is evaluated
  !> on the grids directly. Returns a non-empty errmsg if the
  !> evaluation failed. Reports the progress by lines and stops when
  !> an abort is requested (lf is then incomplete).
  module subroutine rhoplot_cube_values(s,x0,xd,nn,id,expr,nti,oncell,lf,errmsg)
    use systemmod, only: system
    use grid3mod, only: grid3
    use global, only: abort_poll, progress_start, progress_step
    use types, only: scalar_value, vstring, field_evaluation_avail
    use iso_c_binding, only: c_loc
    type(system), intent(inout), target :: s
    real*8, intent(in) :: x0(3), xd(3,3)
    integer, intent(in) :: nn(3), id, nti
    character*(*), intent(in) :: expr
    logical, intent(in) :: oncell
    real*8, allocatable, intent(inout) :: lf(:,:,:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: ix, iy, iz
    real*8 :: xp(3), lappt
    type(scalar_value) :: res
    type(vstring) :: lerrs
    type(field_evaluation_avail) :: request
    type(grid3) :: faux
    character(len=:), allocatable :: lerrmsg

    errmsg = ""
    if (allocated(lf)) deallocate(lf)
    allocate(lf(nn(1),nn(2),nn(3)))

    ! an expression of grid fields on the whole cell is evaluated on
    ! the grids directly
    if (id < 0 .and. oncell) then
       call faux%new_eval(c_loc(s),c_loc(s%c),nn,expr,s%c%m_x2c)
       if (faux%isinit) then
          lf = faux%f
          call faux%end()
          return
       end if
    end if

    request = rhoplot_nti_request(nti)

    ! calculate the cube; the private error string is in a vstring:
    ! gfortran mishandles private deferred-length character variables
    lerrmsg = ""
    call progress_start(nn(3)*nn(2),"lines")
    !$omp parallel do private(xp,res,lappt,lerrs,ix,iy) schedule(dynamic)
    do iz = 0, nn(3)-1
       do iy = 0, nn(2)-1
          if (abort_poll()) cycle
          do ix = 0, nn(1)-1
             xp = x0 + real(ix,8) * xd(:,1) + real(iy,8) * xd(:,2) &
                + real(iz,8) * xd(:,3)
             lerrs%s = ""
             if (id >= 0) then
                call s%f(id)%grd(xp,request,res,periodic=.not.s%c%ismolecule)
                select case(nti)
                case (0)
                   lappt = res%f
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
                case default
                   lappt = 0d0
                end select
             else
                lappt = s%eval(expr,lerrs%s,xp)
             end if
             lf(ix+1,iy+1,iz+1) = lappt
             if (len_trim(lerrs%s) > 0) then
                !$omp critical (fieldwrite)
                lerrmsg = lerrs%s
                !$omp end critical (fieldwrite)
             end if
          end do
          call progress_step()
       end do
    end do
    !$omp end parallel do
    errmsg = lerrmsg

  end subroutine rhoplot_cube_values

  !> Check that the nodes of field id of system s can be written as a
  !> cube: the field is a grid and, for a wavefunction (ipsi /=
  !> cubepsi_none, band ibnd), a periodic QE grid with the right
  !> files. supercell: the grid is to be repeated in a supercell.
  !> Returns the reason if it cannot, an empty string otherwise.
  module function rhoplot_cube_gridcheck(s,id,ipsi,ibnd,supercell) result(errmsg)
    use systemmod, only: system
    use fieldmod, only: type_grid
    type(system), intent(in) :: s
    integer, intent(in) :: id, ipsi, ibnd
    logical, intent(in) :: supercell
    character(len=:), allocatable :: errmsg

    logical :: ok

    errmsg = ""
    ok = (id > -1)
    if (ok) ok = ok .and. s%f(id)%type == type_grid .and. allocated(s%f(id)%grid)
    if (.not.ok) then
       errmsg = 'CUBE GRID/MLWF/... can only be used with grid fields'
       return
    end if
    if (ipsi /= cubepsi_none) then
       if (.not.s%f(id)%grid%isqe) then
          errmsg = 'CUBE MLWF/WANNIER/... requires a QE wavefunction file (pwc)'
       elseif (ipsi == cubepsi_mlwf.and..not.s%f(id)%grid%iswan) then
          errmsg = 'CUBE MLWF requires a wannier90 checkpoint file (chk)'
       elseif (ibnd < 1 .or. ibnd > s%f(id)%grid%qe%nbnd) then
          errmsg = 'CUBE MLWF/WANNIER/...: incorrect band number'
       end if
       if (len(errmsg) > 0) return
    end if
    if (s%f(id)%grid%partial) then
       errmsg = 'CUBE MLWF/WANNIER/...: cannot be used with non-periodic grids'
    elseif (supercell .and. ipsi /= cubepsi_none) then
       errmsg = 'MLWF/WANNIER/...: incompatible with NEWCELL'
    end if

  end function rhoplot_cube_gridcheck

  !> Return in f the nodes of grid field id of system s (checked with
  !> rhoplot_cube_gridcheck). If na(1) > 0, the grid is repeated
  !> na(i) times along the cell axis ia(i) (i = 1..3, a permutation),
  !> as in CUBE GRID NEWCELL. If the field has a core augmentation,
  !> the core density is added (an abort stops it: f is then
  !> incomplete). If ipsi /= cubepsi_none, f is instead the part ipart
  !> (cubepart_*) of the wavefunction: band ibnd, k-point ik (UNK and
  !> PSINK), spin ispin, and lattice translation inr (all but UNK).
  module subroutine rhoplot_cube_gridvalues(s,id,ia,na,ipsi,ibnd,ik,inr,ispin,ipart,f)
    use systemmod, only: system
    use global, only: abort_requested
    type(system), intent(inout) :: s
    integer, intent(in) :: id, ia(3), na(3), ipsi, ibnd, ik, inr(3), ispin, ipart
    real*8, allocatable, intent(inout) :: f(:,:,:)

    integer :: i1, i2, i3, j, k, nn(3)
    real*8, allocatable :: ff(:,:,:), ff2(:,:,:)
    complex*16, allocatable :: caux(:,:,:)

    if (allocated(f)) deallocate(f)
    nn = s%f(id)%grid%n
    if (ipsi == cubepsi_none) then
       if (na(1) > 0) then
          ! GRID keyword, with NEWCELL -> rearrange grid
          allocate(ff2(nn(ia(1)),nn(ia(2)),nn(ia(3))))
          do k = 1, nn(3)
             do j = 1, nn(2)
                if (ia(2) == 1 .and. ia(3) == 2) then     ! ia(1) == 3
                   ff2(k,:,j) = s%f(id)%grid%f(:,j,k)
                elseif (ia(2) == 1 .and. ia(3) == 3) then ! ia(1) == 2
                   ff2(j,:,k) = s%f(id)%grid%f(:,j,k)
                elseif (ia(2) == 2 .and. ia(3) == 1) then ! ia(1) == 3
                   ff2(k,j,:) = s%f(id)%grid%f(:,j,k)
                elseif (ia(2) == 2 .and. ia(3) == 3) then ! ia(1) == 1
                   ff2(:,j,k) = s%f(id)%grid%f(:,j,k)
                elseif (ia(2) == 3 .and. ia(3) == 1) then ! ia(1) == 2
                   ff2(j,k,:) = s%f(id)%grid%f(:,j,k)
                elseif (ia(2) == 3 .and. ia(3) == 2) then ! ia(1) == 1
                   ff2(:,k,j) = s%f(id)%grid%f(:,j,k)
                end if
             end do
          end do

          ! multiply rearranged grid
          allocate(f(na(1)*nn(ia(1)),na(2)*nn(ia(2)),na(3)*nn(ia(3))))
          do i1 = 1, na(1)
             do i2 = 1, na(2)
                do i3 = 1, na(3)
                   f((i1-1)*nn(ia(1))+1:i1*nn(ia(1)),(i2-1)*nn(ia(2))+1:i2*nn(ia(2)),&
                      (i3-1)*nn(ia(3))+1:i3*nn(ia(3))) = ff2
                end do
             end do
          end do
       else
          ! GRID keyword, current cell
          f = s%f(id)%grid%f
          if (s%f(id)%usecore) then
             call s%c%promolecular_array3(ff,nn,s%f(id)%zpsp)
             if (abort_requested) return
             f = f + ff
          end if
       end if
    else
       allocate(caux(nn(1),nn(2),nn(3)))
       if (ipsi == cubepsi_mlwf .or. ipsi == cubepsi_wannier) then
          ! MLWF and WANNIER keywords
          call s%f(id)%grid%get_qe_wnr_standalone(s%f(id)%c%omega,ibnd,ispin,inr,(ipsi == cubepsi_mlwf),caux)
       else
          ! UNK and PSINK keyword
          call s%f(id)%grid%get_qe_psink_standalone(s%f(id)%c%omega,ibnd,ik,ispin,(ipsi == cubepsi_psink),inr,caux)
       end if
       if (ipart == cubepart_real) then
          f = real(caux,8)
       elseif (ipart == cubepart_imag) then
          f = aimag(caux)
       else
          f = abs(caux)
       end if
    end if

  end subroutine rhoplot_cube_gridvalues

  !> Write the values f of a cube with origin x0 and step vectors xd
  !> (Cartesian, bohr, cell frame of crystal cr) to file, in format
  !> outform (cubefmt_*). The VASP and XSF formats take no origin or
  !> steps: they are for grids on the whole cell. If header, write
  !> only the header (f may be unallocated). ishift shifts the nodes
  !> of a text cube. Returns a non-empty errmsg if the file cannot be
  !> written.
  module subroutine rhoplot_cube_write(cr,outform,f,file,header,xd,x0,ishift,errmsg)
    use crystalmod, only: crystal
    type(crystal), intent(in) :: cr
    integer, intent(in) :: outform
    real*8, allocatable, intent(in) :: f(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: header
    real*8, intent(in) :: xd(3,3), x0(3)
    integer, intent(in) :: ishift(3)
    character(len=:), allocatable, intent(out) :: errmsg

    errmsg = ""
    if (outform == cubefmt_bincube) then
       if (header) then
          errmsg = "BINCUBE format is incompatible with HEADER"
       else
          call cr%writegrid_cube(f,file,.false.,.true.,xd,x0+cr%molx0)
       end if
    elseif (outform == cubefmt_cube) then
       if (header) then
          call cr%writegrid_cube(f,file,.true.,.false.,xd,x0+cr%molx0)
       else
          call cr%writegrid_cube(f,file,.false.,.false.,xd,x0+cr%molx0,ishift0=ishift)
       end if
    elseif (outform == cubefmt_vasp) then
       call cr%writegrid_vasp(f,file,header,.false.)
    elseif (outform == cubefmt_vaspnov) then
       call cr%writegrid_vasp(f,file,header,.true.)
    elseif (outform == cubefmt_xsf) then
       call cr%writegrid_xsf(f,file,header)
    end if

  end subroutine rhoplot_cube_write

  !> The cube format (cubefmt_*) of file, from its extension: cube,
  !> bincube, xsf, ELFCAR (VASP without the valence), and VASP
  !> otherwise.
  module function rhoplot_cube_outform(file) result(outform)
    use tools_io, only: equal
    character*(*), intent(in) :: file
    integer :: outform

    character(len=:), allocatable :: wext1

    wext1 = file(index(file,'.',.true.)+1:)
    if (equal(wext1,'cube')) then
       outform = cubefmt_cube
    elseif (equal(wext1,'bincube')) then
       outform = cubefmt_bincube
    elseif (equal(wext1,'xsf')) then
       outform = cubefmt_xsf
    elseif (equal(wext1,'ELFCAR')) then
       outform = cubefmt_vaspnov
    else
       outform = cubefmt_vasp
    end if

  end function rhoplot_cube_outform

  !> Calculate properties on a plane.
  module subroutine rhoplot_plane(line)
    use systemmod, only: sy
    use global, only: eval_next, dunit0, iunit, fileroot, iunitname0, iunit, dunit0,&
       abort_requested
    use tools_io, only: ferror, faterr, lgetword, equal, getword, &
       isexpression_or_word, fopen_write, uout, string, fclose
    use tools_math, only: plane_scale_extend, assign_ziso, niso_manual,&
       niso_lin, niso_log, niso_atan, niso_bader
    use types, only: realloc
    character*(*), intent(in) :: line

    integer :: lp2, lp, nti, id, luout, nx, ny, niso_type, niso, nn
    real*8 :: x0(3), x1(3), x2(3), xp(3), du, dv
    real*8 :: uu(3), vv(3), fmin, fmax
    real*8 :: sx0, sy0, zx0, zx1, zy0, zy1, zmin, zmax, rdum
    logical :: docontour, dorelief, docolormap, fset
    character(len=:), allocatable :: word, outfile, root0, expr
    logical :: ok
    integer :: ix, iy, cmopt
    real*8, allocatable :: ff(:,:), ziso(:)
    character(len=:), allocatable :: errmsg

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
                return
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

    ! field values on the plane
    if (.not.allocated(expr)) expr = ""
    allocate(ff(nx,ny))
    call rhoplot_plane_values(sy,x0,uu,vv,nx,ny,id,expr,nti,ff,errmsg,progress=.true.)
    if (len_trim(errmsg) > 0) &
       call ferror('rhoplot_plane',errmsg,faterr)
    if (abort_requested) then
       write (uout,'("+ The calculation was cancelled: no file was written."/)')
       return
    end if

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

  !> Calculate the values of a field on the points x0 + (ix-1)*uu +
  !> (iy-1)*vv, ix = 1..nx, iy = 1..ny, of a plane (Cartesian, bohr,
  !> cell frame of system s) in ff(nx,ny). If id >= 0, use field id
  !> and its property nti (0 = f, 1-3 = gradient components, 4 =
  !> gradient norm, 5-10 = Hessian xx, xy, xz, yy, yz, zz, 11 =
  !> Laplacian); otherwise, evaluate the expression expr. Returns a
  !> non-empty errmsg if the evaluation failed. If progress, report
  !> the progress by rows and stop when an abort is requested (the
  !> rows not calculated are zero).
  module subroutine rhoplot_plane_values(s,x0,uu,vv,nx,ny,id,expr,nti,ff,errmsg,progress)
    use systemmod, only: system
    use types, only: scalar_value, vstring, field_evaluation_avail
    use global, only: progress_start, progress_step, abort_poll
    type(system), intent(inout) :: s
    real*8, intent(in) :: x0(3), uu(3), vv(3)
    integer, intent(in) :: nx, ny, id, nti
    character*(*), intent(in) :: expr
    real*8, intent(out) :: ff(nx,ny)
    character(len=:), allocatable, intent(out) :: errmsg
    logical, intent(in), optional :: progress

    integer :: ix, iy
    logical :: progress_
    real*8 :: xp(3), rhopt
    type(scalar_value) :: res
    type(vstring) :: lerrs
    type(field_evaluation_avail) :: request
    character(len=:), allocatable :: lerrmsg

    progress_ = .false.
    if (present(progress)) progress_ = progress

    request = rhoplot_nti_request(nti)

    ! calculate the plane; the private error string is in a vstring:
    ! gfortran mishandles private deferred-length character variables
    lerrmsg = ""
    ff = 0d0
    if (progress_) call progress_start(nx,"rows")
    !$omp parallel do private(xp,res,rhopt,lerrs,iy) schedule(dynamic)
    do ix = 1, nx
       if (progress_) then
          if (abort_poll()) cycle
       end if
       do iy = 1, ny
          xp = x0 + real(ix-1,8) * uu + real(iy-1,8) * vv

          lerrs%s = ""
          if (id >= 0) then
             call s%f(id)%grd(xp,request,res,periodic=.not.s%c%ismolecule)
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
             case default
                rhopt = 0d0
             end select
          else
             rhopt = s%eval(expr,lerrs%s,xp)
          endif
          ff(ix,iy) = rhopt
          if (len_trim(lerrs%s) > 0) then
             !$omp critical (iowrite)
             lerrmsg = lerrs%s
             !$omp end critical (iowrite)
          end if
       end do
       if (progress_) call progress_step()
    end do
    !$omp end parallel do
    errmsg = lerrmsg

  end subroutine rhoplot_plane_values

  !> Contour lines of the values ff(nx,ny) on the points of a plane
  !> with sides of length r01 (the first index) and r02 (the second)
  !> forming an angle whose cosine is cosalfa, for the isovalues
  !> ziso(1:niso). Returns the nline lines, in order of isovalue:
  !> line i has the points lptr(i) to lptr(i+1)-1 of xy, which are
  !> the coordinates (bohr) in an orthogonal frame of the plane (x
  !> along the first side, y perpendicular to it), and the isovalue
  !> ziso(liso(i)). A closed line repeats its first point at the end.
  module subroutine rhoplot_contour_lines(ff,nx,ny,r01,r02,cosalfa,niso,ziso,nline,lptr,xy,liso)
    use types, only: realloc
    integer, intent(in) :: nx, ny
    real*8, intent(in) :: ff(nx,ny)
    real*8, intent(in) :: r01, r02, cosalfa
    integer, intent(in) :: niso
    real*8, intent(in) :: ziso(niso)
    integer, intent(out) :: nline
    integer, allocatable, intent(inout) :: lptr(:)
    real*8, allocatable, intent(inout) :: xy(:,:)
    integer, allocatable, intent(inout) :: liso(:)

    integer :: i, j, npt, np0
    real*8, allocatable :: x(:), y(:)
    real*8, allocatable :: xc(:), yc(:), pic(:), pjc(:)

    nline = 0
    npt = 0
    if (allocated(lptr)) deallocate(lptr)
    if (allocated(xy)) deallocate(xy)
    if (allocated(liso)) deallocate(liso)
    allocate(lptr(11),xy(2,1000),liso(10))
    lptr(1) = 1

    ! the grid in each direction
    allocate(x(nx),y(ny))
    do i = 1, nx
       x(i) = (i-1) * r01 / real(nx-1,8)
    end do
    do j = 1, ny
       y(j) = (j-1) * r02 / real(ny-1,8)
    end do

    ! the crossings for each isovalue, joined into lines
    allocate(xc(1000),yc(1000),pic(1000),pjc(1000))
    do i = 1, niso
       call hallarpuntos(ff,ziso(i),x,y,nx,ny,np0,xc,yc,pic,pjc)
       call ordenarpuntos(np0,xc,yc,pic,pjc,cosalfa,i,nline,npt,lptr,xy,liso)
    end do
    call realloc(lptr,nline+1)
    call realloc(liso,nline)
    call realloc(xy,2,npt)

  end subroutine rhoplot_contour_lines

  !> Plot of gradient paths and contours in the style of aimpac's grdvec.
  module subroutine rhoplot_grdvec()
    use systemmod, only: sy
    use global, only: fileroot, eval_next, dunit0, iunit
    use tools_io, only: uout, uin, ucopy, getline, lgetword, equal,&
       faterr, ferror, string, ioj_right, fopen_write, getword, fclose
    use tools_math, only: plane_scale_extend, assign_ziso, &
       niso_manual, niso_atan, niso_lin, niso_log, niso_bader
    use types, only: scalar_value, realloc, field_evaluation_avail
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
    integer :: n1, n2, niso
    type(scalar_value) :: res
    real*8, allocatable :: ff(:,:), ziso(:)
    type(field_evaluation_avail) :: request

    ! Header
    write (uout,'("* GRDVEC: gradient paths and contours in 2d")')

    ! Initialization
    call request%field_nder2()
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
             call sy%f(sy%iref)%grd(q0,request,res,periodic=.not.sy%c%ismolecule)
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
       request = rhoplot_nti_request(nfi)
       !$omp parallel do private(xp,res,rhopt) schedule(dynamic)
       do ix = 1, n1
          do iy = 1, n2
             xp = x0 + real(ix-1,8) * uu + real(iy-1,8) * vv
             call sy%f(sy%iref)%grd(xp,request,res,periodic=.not.sy%c%ismolecule)
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
             ff(ix,iy) = rhopt
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
    integer :: i, j, nline
    real*8 :: r012, ua, va, ub, vc
    integer, allocatable :: lptr(:), liso(:)
    real*8, allocatable :: xy(:,:), x(:), y(:)

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

    ! the contour lines: the positive isovalues in one file and the
    ! rest in the other
    isneg = minval(ff) < 0d0
    call rhoplot_contour_lines(ff,nx,ny,ub-ua,vc-va,cosalfa,niso,ziso,nline,lptr,xy,liso)
    do i = 1, nline
       j = lptr(i+1) - lptr(i)
       x = xy(1,lptr(i):lptr(i+1)-1)
       y = xy(2,lptr(i):lptr(i+1)-1)
       if (ziso(liso(i)).gt.0) then
          call linea(x,y,j,lud,ziso(liso(i)))
       else
          call linea(x,y,j,lud1,ziso(liso(i)))
       endif
    enddo
    call fclose(lud)
    call fclose(lud1)

    if (dolabels) call write_fichlabel(root0)
    if (dognu) call write_fichgnu(root0,dolabels,.true.,.false.)

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

  !> Find the points of the contour with value zc on the grid ff
  !> (coordinates x, y), using linear interpolation: np of them, in
  !> xc, yc, with the grid indices pic, pjc (half-integer = between
  !> two nodes). The arrays grow as needed.
  subroutine hallarpuntos(ff,zc,x,y,nx,ny,np,xc,yc,pic,pjc)
    use ieee_arithmetic, only: ieee_is_finite
    use param, only: zero
    use types, only: realloc
    integer, intent(in) :: nx, ny
    real*8, intent(in) :: ff(nx,ny)
    real*8, intent(in) :: zc
    real*8, intent(in) :: x(:), y(:)
    integer, intent(out) :: np
    real*8, allocatable, intent(inout) :: xc(:), yc(:), pic(:), pjc(:)

    integer :: i, j
    real*8 :: xa, ya, za, xb, yb, zb
    real*8 :: zazc, zbzc

    ! initialize
    np = 0

    ! A node on the level counts as above it (half-open rule): a
    ! crossing is an edge with one end below and the other on or above
    ! the level. Nodes exactly on the level used to be points of their
    ! own, which turned a plateau on the level (a zero field, the
    ! minimum of the data) into a lattice of spurious lines. The lines
    ! do not enter a region of non-finite values.

    ! run over columns
    do i = 1, nx
       xa = x(i)
       ya = y(1)
       za = ff(i,1)
       zazc = za-zc

       do j = 2,ny
          xb = x(i)
          yb = y(j)
          zb = ff(i,j)
          zbzc = zb-zc

          ! sign changed, interpolate and write a point
          if (ieee_is_finite(za) .and. ieee_is_finite(zb)) then
             if ((zazc < zero) .neqv. (zbzc < zero)) &
                call add_point(xb,ya-zazc*(yb-ya)/(zb-za),real(i,8),j-0.5d0)
          end if
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
          if (ieee_is_finite(za) .and. ieee_is_finite(zb)) then
             if ((zazc < zero) .neqv. (zbzc < zero)) &
                call add_point(xa-zazc*(xb-xa)/(zb-za),yb,i-0.5d0,real(j,8))
          end if
          ! reassign zazc
          xa = xb
          za = zb
          zazc = zbzc
       enddo
    enddo

  contains
    subroutine add_point(xx,yy,pi,pj)
      real*8, intent(in) :: xx, yy, pi, pj

      np = np + 1
      if (np > size(xc)) then
         call realloc(xc,2*np)
         call realloc(yc,2*np)
         call realloc(pic,2*np)
         call realloc(pjc,2*np)
      end if
      xc(np) = xx
      yc(np) = yy
      pic(np) = pi
      pjc(np) = pj

    end subroutine add_point
  end subroutine hallarpuntos

  !> Determines the connectivity of the np contour points xc, yc
  !> (grid indices pic, pjc; see hallarpuntos) and adds the lines,
  !> with isovalue index iso, to the nline lines with npt points in
  !> lptr, xy, liso (see rhoplot_contour_lines), in the orthogonal
  !> frame of the plane (calpha = cosine of the angle of its sides).
  subroutine ordenarpuntos(np,xc,yc,pic,pjc,calpha,iso,nline,npt,lptr,xy,liso)
    use param, only: one, half, zero
    use types, only: realloc
    real*8, parameter :: eps = 0.10d0

    integer, intent(in) :: np
    real*8, intent(in) :: xc(:), yc(:), pic(:), pjc(:)
    real*8, intent(in) :: calpha
    integer, intent(in) :: iso
    integer, intent(inout) :: nline, npt
    integer, allocatable, intent(inout) :: lptr(:)
    real*8, allocatable, intent(inout) :: xy(:,:)
    integer, allocatable, intent(inout) :: liso(:)

    real*8 :: salpha

    logical :: cerrada, primerabusqueda, hallado, incompleta
    integer :: puntoinicial
    integer :: puntofinal, puntoactual, salto, nptoscurva
    integer :: nptosestudiados, npuntos
    ! allocatable: np is not bounded, and these could not fit in the stack
    logical, allocatable :: lc(:)
    integer, allocatable :: punto(:)
    real*8, allocatable :: x(:), y(:)
    integer :: k, n1
    real*8  :: pi0, pj0, pi1, pj1, pi2, pj2, si1, sj1

    ! inline functions
    logical :: adyacentes
    adyacentes(pi0,pj0,pi1,pj1,si1,sj1)=(abs(abs(pi0-pi1)-si1) < eps).and.&
       (abs(abs(pj0-pj1)-sj1) < eps)

    ! initialize
    npuntos = np
    allocate(lc(np),punto(-np-1:np+1),x(np+1),y(np+1))
    lc = .false.
    nptosestudiados = 0

    ! run over all points
    x = 0d0
    y = 0d0
    n1 = 1
    do while (nptosestudiados.lt.npuntos)

       ! search for a non-used point and start a curve
       k = n1
       do while (k.le.npuntos)
          if (.not.lc(k)) exit
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

             ! (pi1,pj1) is the last found point: a crossing of a column
             ! (integer pi1) continues to the next column, a crossing of
             ! a row (integer pj1) to the next row, and either to a
             ! crossing of the other kind in the same cell
             si1 = zero
             sj1 = zero
             if (abs(pi1-nint(pi1)) <= eps) si1 = one
             if (abs(pj1-nint(pj1)) <= eps) sj1 = one
             k = n1
             do while ((.not.hallado) .and. (k.lt.npuntos))
                k = k+1
                if (.not.lc(k)) then
                   pi2 = pic(k)
                   pj2 = pjc(k)
                   if (adyacentes(pi1,pj1,pi2,pj2,half,half) .or. adyacentes(pi1,pj1,pi2,pj2,si1,sj1)) then
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
             enddo
          enddo

          ! no adyacent point found. Probe if it is a closed or open curve
          if ((puntofinal-puntoinicial+1).gt.2) then
             pi0 = pic(punto(puntoinicial))
             pj0 = pjc(punto(puntoinicial))
             pi1 = pic(punto(puntofinal))
             pj1 = pjc(punto(puntofinal))
             si1 = zero
             sj1 = zero
             if (abs(pi1-nint(pi1)) <= eps) si1 = one
             if (abs(pj1-nint(pj1)) <= eps) sj1 = one
             if (adyacentes(pi0,pj0,pi1,pj1,half,half) .or. adyacentes(pi0,pj0,pi1,pj1,si1,sj1)) &
                cerrada = .true.
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

       call add_line(nptoscurva)

       ! update number of points
       if (cerrada) nptoscurva = nptoscurva-1
       nptosestudiados = nptosestudiados+nptoscurva
    enddo

  contains
    !> Add the line in x(1:n), y(1:n) to the output.
    subroutine add_line(n)
      integer, intent(in) :: n

      integer :: kk

      if (npt + n > size(xy,2)) call realloc(xy,2,max(2*size(xy,2),npt+n))
      if (nline + 2 > size(lptr)) call realloc(lptr,2*size(lptr))
      if (nline + 1 > size(liso)) call realloc(liso,2*size(liso))
      do kk = 1, n
         xy(1,npt+kk) = x(kk)
         xy(2,npt+kk) = y(kk)
      end do
      npt = npt + n
      nline = nline + 1
      liso(nline) = iso
      lptr(nline+1) = npt + 1

    end subroutine add_line
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
    use types, only: scalar_value, gpathp, field_evaluation_avail
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
    type(field_evaluation_avail) :: request

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
    call request%field_nder2()

    write (uout,'("+ List of gradient paths traced")')
    write (uout,'("# i       xcrys        ycrys        zcrys        type    up down    pts")')
    do iorig = 1, norig
       ntotpts = 0
       grpx(:,iorig) = sy%c%x2c(grpx(:,iorig))

       if (grpatr(iorig) .eq. 1) then
          ! 3D attraction or repulsion critical points:
          do i = 1, grpup(iorig)
             phii = real(i-1,8) * 2d0 * pi / real(grpup(iorig),8)
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
             phii = real(i-1,8) * 2d0 * pi / real(grpdwn(iorig),8)
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
          call sy%f(sy%iref)%grd(xstart,request,res,periodic=.not.sy%c%ismolecule)
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
    use types, only: scalar_value, field_evaluation_avail
    integer :: i, j, k, l, ncopies
    integer :: iorde(2*indmax+1), indcell(3,(2*indmax+1)**3), iii, inum
    real*8  :: xp(3)
    real*8  :: x0(3), x1(3), uu, vv, hh, hmin
    type(scalar_value) :: res
    real*8 :: rp0(3)
    type(field_evaluation_avail) :: request
    logical :: rejected(mncritp)

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
    call request%field_nder2()
    rejected(1:newncriticp) = .false.
    do i = 1, newncriticp
       if (rejected(i)) cycle
       xp = newcriticp(:,i)

       !.discard repeated points
       do j = i+1, newncriticp
          if (rejected(j)) cycle
          if (sy%c%are_close(newcriticp(:,j),xp,epsf)) then
             write (uout,'("  CP ",A," is equivalent to ",A," -> Rejected!")') string(j), string(i)
             rejected(j) = .true.
          endif
       enddo

       !.Get the properties
       xp = sy%c%x2c(newcriticp(:,i))
       call sy%f(sy%iref)%grd(xp,request,res,periodic=.not.sy%c%ismolecule)
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
