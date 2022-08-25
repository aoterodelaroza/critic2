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

! sigmahole integration
submodule (sigmahole) proc
  implicit none

contains

  ! Driver for the Hirshfeld routines.
  module subroutine sigmahole_driver(s,line)
    use systemmod, only: system
    use fieldmod, only: type_wfn
    use global, only: iunit, dunit0, iunitname0, eval_next, fileroot
    use tools_math, only: cross
    use tools_io, only: ferror, faterr, uout, fopen_write, fclose, string, ioj_center,&
       isinteger, lgetword, equal
    use types, only: scalar_value, realloc
    use param, only: pi

    type(system), intent(inout) :: s
    character*(*), intent(in) :: line

    real*8, parameter :: brak_init = 1.0d0
    real*8, parameter :: brak_rat = 1.25d0
    real*8, parameter :: brak_max = 10d0
    real*8, parameter :: range_max = 50d0
    real*8, parameter :: bisect_eps = 1d-5

    real*8 :: minv, xp(3), uvec(3), isoval
    real*8 :: xc(3), xx(3), u, v, th, ph, xl(3)
    real*8 :: umat(3,3), rat, rini, rfin, rmid, mep, maxang
    integer :: iu, iv, i, lu, lug, nmax, iunext, iuprev
    integer :: i1st, lp, ic, ix, nu, nv
    type(scalar_value) :: res
    real*8, allocatable :: mval(:,:), rval(:,:), ms(:), ls(:), rs(:), ss(:), ds(:)
    integer, allocatable :: iassign(:,:), uvmax(:,:), ival(:)
    character*(:), allocatable :: strrs, strss, word, filename, gnuname
    logical :: doit, changed, ok

    ! checks
    if (.not.s%c%ismolecule) &
       call ferror("sigmahole_driver","sigmahole only for molecules",faterr)
    if (s%f(s%iref)%type /= type_wfn) &
       call ferror("sigmahole_driver","sigmahole only for wfn fields",faterr)
    if (.not.allocated(s%f(s%iref)%wfn%cint)) &
       call ferror('sigmahole_driver','basis set data required for sigmahole calculation',faterr)

    ! parse input
    filename = trim(fileroot) // "_sigmahole.dat"
    gnuname = trim(fileroot) // "_sigmahole.gnu"

    nu = 101
    nv = 101
    isoval = 1d-3
    maxang = 90d0
    write (uout,'("* SIGMAHOLE: charaterization of molecular sigma-hole regions")')
    lp = 1
    ok = isinteger(ic,line,lp)
    ok = ok .and. isinteger(ix,line,lp)
    if (.not.ok) then
       call ferror('sigmahole_driver','Must indicate identifiers for basal and electronegative atoms',faterr,line,syntax=.true.)
       return
    end if
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"npts")) then
          ok = isinteger(nu,line,lp)
          ok = ok .and. isinteger(nv,line,lp)
          if (.not.ok) then
             call ferror('sigmahole_driver','Error reading NPTS in SIGMAHOLE',faterr,line,syntax=.true.)
             return
          end if
          if (nu < 0 .or. nv < 0) then
             call ferror('sigmahole_driver','Invalid NPTS in SIGMAHOLE',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"isoval")) then
          ok = eval_next(isoval,line,lp)
          if (.not.ok) then
             call ferror('sigmahole_driver','Error reading ISOVAL in SIGMAHOLE',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"maxangle")) then
          ok = eval_next(maxang,line,lp)
          if (.not.ok) then
             call ferror('sigmahole_driver','Error reading MAXANGLE in SIGMAHOLE',faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror('sigmahole_driver','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! some more checks
    if (ic < 0 .or. ic > s%c%ncel) &
       call ferror('sigmahole_driver','erroneous atomic ID for basal atom',faterr)
    if (ix < 0 .or. ix > s%c%ncel) &
       call ferror('sigmahole_driver','erroneous atomic ID for hole-bearing atom',faterr)

    ! some info to the output
    write (uout,'("+ Data file: ",A)') trim(filename)
    write (uout,'("+ gnuplot file: ",A)') trim(gnuname)
    write (uout,'("  Base atom (B): ",A," (",A,")")') trim(s%c%spc(s%c%at(ic)%is)%name), string(ic)
    write (uout,'("  Electronegative atom (X): ",A," (",A,")")') trim(s%c%spc(s%c%at(ix)%is)%name), string(ix)
    write (uout,'("  Electron density isosurface (a.u.): ",A)') string(isoval,'e',decimal=4)
    write (uout,'("  Number of u-points: ",A)') string(nu)
    write (uout,'("  Number of v-points: ",A)') string(nv)
    write (uout,'("  Maximum B-X-sigma_hole angle (degrees): ",A)') string(maxang,'f',decimal=2)

    ! calculate center point and local coordinate frame
    xc = s%c%atcel(ic)%r
    xx = s%c%atcel(ix)%r
    umat(:,3) = xx - xc
    umat(:,3) = umat(:,3) / norm2(umat(:,3))
    xp = (/1d0,0d0,0d0/)
    xp = cross(xp,umat(:,3))
    if (dot_product(xp,xp) < 1d-2) then
       xp = (/0d0,1d0,0d0/)
       xp = cross(xp,umat(:,3))
    end if
    umat(:,2) = xp / norm2(xp)
    umat(:,1) = cross(umat(:,2),umat(:,3))
    umat(:,1) = umat(:,1) / norm2(umat(:,1))

    ! minimum v value
    minv = (cos(maxang * pi / 180d0) + 1d0) / 2d0
    write (uout,'("  Solid angle per uv-point (sphere = 4*pi): ",A)') string(4 * pi * (1-minv) / real(nu*nv,8),'e',decimal=4)

    allocate(mval(nu,nv),rval(nu,nv))

    ! calculate and plot
    !$omp parallel do private (u,v,ph,th,xl,uvec,rini,rat,rfin,xp,res,rmid,mep)
    do iu = 1, nu
       u = real(iu-1,8) / real(nu-1,8)
       do iv = 1, nv
          ! calculate the unitary vector
          v = (real(iv-1,8) / real(nv-1,8)) * (1-minv) + minv
          ph = acos(2*v-1)
          th = 2*pi*u
          xl(1) = cos(th)*sin(ph)
          xl(2) = sin(th)*sin(ph)
          xl(3) = cos(ph)
          uvec = xl(1) * umat(:,1) + xl(2) * umat(:,2) + xl(3) * umat(:,3)

          ! bracket the value
          rini = brak_init
          rat = brak_init
          rfin = -1d0
          do while(rat < brak_max)
             xp = xx + rat * uvec
             call s%f(s%iref)%grd(xp,0,res)
             if (res%f < isoval) then
                if (rat == brak_init) &
                   call ferror('sigmahole_driver','failed bracketing: lower than isovalue at first point',faterr)
                rfin = rat
                exit
             else
                rini = rat
             end if
             rat = rat * brak_rat
          end do
          if (rfin < 0d0) &
             call ferror('sigmahole_driver','failed bracketing: isovalue not found at rmax',faterr)

          ! bisect the value
          do while (abs(rfin-rini) > bisect_eps)
             rmid = 0.5d0 * (rini + rfin)
             xp = xx + rmid * uvec
             call s%f(s%iref)%grd(xp,0,res)
             if (res%f > isoval) then
                rini = rmid
             else
                rfin = rmid
             end if
          end do
          mep = s%f(s%iref)%wfn%mep(xp)

          !$omp critical (iowrite)
          rval(iu,iv) = rmid
          mval(iu,iv) = mep
          !$omp end critical (iowrite)
       end do
    end do

    ! write the data file
    lu = fopen_write(filename)
    write (lu,'("# u v theta phi x y z r(u,v) mep(u,v)")')
    do iu = 1, nu
       u = real(iu-1,8) / real(nu-1,8)
       do iv = 1, nv
          v = (real(iv-1,8) / real(nv-1,8)) * (1-minv) + minv
          ph = acos(2*v-1)
          th = 2*pi*u
          xl(1) = cos(th)*sin(ph)
          xl(2) = sin(th)*sin(ph)
          xl(3) = cos(ph)
          xl = xl * rval(iu,iv)

          write (lu,'(99(A," "))') string(u,'e',decimal=8), string(v,'e',decimal=8), &
             string(th,'e',decimal=8), string(ph,'e',decimal=8), &
             (string(xl(i)*dunit0(iunit),'e',decimal=8),i=1,3), &
             string(rval(iu,iv)*dunit0(iunit),'e',decimal=8), string(mval(iu,iv),'e',decimal=8)
       end do
       write (lu,*)
    end do
    call fclose(lu)

    ! write the gnuplot template file
    lug = fopen_write(gnuname)
    write (lug,'("set terminal wxt size 800,800")')
    write (lug,'("set mapping cartesian")')
    write (lug,'("")')
    write (lug,'("set style line 1 lt 1 lw 1 lc rgb ""#000000""")')
    write (lug,'("set style line 2 lt 1 lw 1 lc rgb ""#000000""")')
    write (lug,'("unset clabel")')
    write (lug,'("")')
    write (lug,'("set palette defined (-0.03 ""red"", 0 ""white"", 0.03 ""blue"")")')
    write (lug,'("")')
    write (lug,'("unset title")')
    write (lug,'("unset key")')
    write (lug,'("")')
    write (lug,'("set view equal xyz")')
    write (lug,'("set xlabel ""x (ang)""")')
    write (lug,'("set ylabel ""y (ang)""")')
    write (lug,'("set zlabel ""z (ang)""")')
    write (lug,'("set view 45")')
    write (lug,'("set pm3d depthorder # interpolate 0,0")')
    write (lug,'("set cbrange [-0.03:0.03]")')
    write (lug,'("")')
    write (lug,'("splot """,A,""" u 5:6:7:9 with pm3d")') filename
    write (lug,'("pause -1")')
    call fclose(lug)

    ! prepare for assignment
    allocate(iassign(nu,nv))
    where(mval < 0d0)
       iassign = -1
    elsewhere
       iassign = 0
    end where

    ! locate maxima
    nmax = 0
    allocate(uvmax(2,10))
    if (all(mval(1,nv) > mval(:,nv-1))) then
       ! maximum at the zenith
       nmax = 1
       uvmax(1,nmax) = 1
       uvmax(2,nmax) = nv
       if (mval(1,nv) > 0d0) iassign(:,nv) = nmax
    end if
    do iv = 2, nv-1
       do iu = 1, nu
          iunext = modulo(iu-1 +1,nu)+1
          iuprev = modulo(iu-1 -1,nu)+1

          if (mval(iu,iv) > mval(iunext,iv) .and. mval(iu,iv) > mval(iuprev,iv) .and.&
             mval(iu,iv) > mval(iu,iv-1) .and. mval(iu,iv) > mval(iu,iv+1)) then
             nmax = nmax + 1
             if (nmax > size(uvmax,2)) call realloc(uvmax,2,2*nmax)
             uvmax(1,nmax) = iu
             uvmax(2,nmax) = iv
             iassign(iu,iv) = nmax
          end if
       end do
    end do
    call realloc(uvmax,2,2*nmax)

    ! calculate maxima properties
    allocate(ms(max(nmax,1)),ls(max(nmax,1)),rs(max(nmax,1)),ss(max(nmax,1)),ds(max(nmax,1)))
    do i = 1, nmax
       iu = uvmax(1,i)
       iv = uvmax(2,i)
       u = real(iu-1,8) / real(nu-1,8)
       v = (real(iv-1,8) / real(nv-1,8)) * (1-minv) + minv
       ph = acos(2*v-1)
       th = 2*pi*u
       xl(1) = cos(th)*sin(ph)
       xl(2) = sin(th)*sin(ph)
       xl(3) = cos(ph)
       uvec = xl(1) * umat(:,1) + xl(2) * umat(:,2) + xl(3) * umat(:,3)

       ! magnitude and distance
       ms(i) = mval(iu,iv)
       ds(i) = rval(iu,iv)

       ! range, bracket and bisect
       if (ms(i) > 0d0) then
          ! bracket
          rini = rval(iu,iv)
          rat = rval(iu,iv) * brak_rat
          rfin = -1d0
          do while(rat < range_max)
             xp = xx + rat * uvec
             mep = s%f(s%iref)%wfn%mep(xp)
             if (mep < 0d0) then
                rfin = rat
                exit
             else
                rini = rat
             end if
             rat = rat * brak_rat
          end do
          if (rfin < 0d0) then
             rs(i) = -1d0
          else
             ! bisect
             do while (abs(rfin-rini) > bisect_eps)
                rmid = 0.5d0 * (rini + rfin)
                xp = xx + rmid * uvec
                mep = s%f(s%iref)%wfn%mep(xp)
                if (mep > 0d0) then
                   rini = rmid
                else
                   rfin = rmid
                end if
             end do
             rs(i) = rmid
          end if
       else
          rs(i) = -1d0
       end if

       ! linearity
       xp = xx + rval(iu,iv) * uvec - xx
       xl = xc - xx
       ls(i) = acos(min(max(dot_product(xp,xl) / norm2(xp) / norm2(xl),-1d0),1d0)) * 180d0 / pi
    end do

    ! build mapping for calculating size
    allocate(ival(4))
    changed = .true.
    do while (changed)
       changed = .false.
       do iu = 1, nu
          do iv = 1, nv-1
             if (iassign(iu,iv) == 0) then
                iunext = modulo(iu-1 +1,nu)+1
                iuprev = modulo(iu-1 -1,nu)+1
                ival(1) = iassign(iunext,iv)
                ival(2) = iassign(iuprev,iv)
                ival(3) = iassign(iu,iv+1)
                if (iv == 1) then
                   ival(4) = 0
                else
                   ival(4) = iassign(iu,iv-1)
                end if
                if (any(ival > 0)) then
                   do i = 1, 4
                      if (ival(i) > 0) then
                         i1st = ival(i)
                         exit
                      end if
                   end do

                   doit = .true.
                   do i = 1, 4
                      if (ival(i) > 0 .and. ival(i) /= i1st) then
                         where (iassign == ival(i) .or. iassign == i1st)
                            iassign = -1
                         end where
                         doit = .false.
                      end if
                   end do

                   if (doit) iassign(iu,iv) = i1st
                   changed = .true.
                end if
             end if
          end do
       end do
    end do
    deallocate(ival)
    if (iassign(1,nv) == 0) then
       allocate(ival(nu))
       ival = iassign(:,nv-1)
       i1st = 0
       do iu = 1, nu
          if (ival(iu) > 0) then
             i1st = ival(iu)
             exit
          end if
       end do

       if (i1st == 0) then
          iassign(:,nv) = -1
       else
          doit = .true.
          do iu = 1, nu
             if (ival(iu) > 0 .and. ival(iu) /= i1st) then
                where (iassign == ival(iu) .or. iassign == i1st)
                   iassign = -1
                end where
                doit = .false.
             end if
          end do
          if (doit) iassign(:,nv) = i1st
       end if
       deallocate(ival)
    end if

    ! calculate size
    ss = 0d0
    do i = 1, nmax
       if (any(iassign == i)) then
          do iu = 1, nu
             do iv = 1, nv
                if (iassign(iu,iv) == i) then
                   ss(i) = ss(i) + rval(iu,iv)**2
                end if
             end do
          end do
          ss(i) = 4 * pi * (1-minv) / real(nu*nv,8) * ss(i)
       else
          ss(i) = -1d0
       end if
    end do

    ! write results
    write (uout,'("+ List of sigma-holes found (",A,"): ")') string(nmax)
    write (uout,'("# See Kolar and Hobza, Chem. Rev. 116 (2016) 5155. doi:10.1021/acs.chemrev.5b00560 (Figure 7)")')
    write (uout,'("# dist = distance between X and the sigma-hole (",A,")")') iunitname0(iunit)
    write (uout,'("# ls (linearity) = base-X-sigma_hole angle (degree)")')
    write (uout,'("# ms (magnitude) = value of the MEP at the sigma-hole (a.u.)")')
    write (uout,'("# rs (range) = distance between X and the point at which the MEP changes sign (",A,")")') iunitname0(iunit)
    write (uout,'("# ss (size) = average aperture angle for MEP-positive region on isosurface (",A,"^2)")') iunitname0(iunit)
    write (uout,'("# ")')
    write (uout,'("#id dist     ls       ms         rs       ss")')
    do i = 1, nmax
       if (rs(i) > 0) then
          strrs = string(rs(i)*dunit0(iunit),'f',length=9,decimal=4,justify=ioj_center)
       else
          strrs = "   n/a   "
       end if
       if (rs(i) > 0) then
          strss = string(ss(i)*dunit0(iunit)**2,'f',length=9,decimal=4,justify=ioj_center)
       else
          strss = "   n/a   "
       end if

       write (uout,'(99(A," "))') string(i,2), string(ds(i)*dunit0(iunit),'f',length=8,decimal=4),&
          string(ls(i),'f',length=7,decimal=2), string(ms(i),'e',decimal=4),&
          strrs, strss
    end do
    write (uout,*)

  end subroutine sigmahole_driver

end submodule proc
