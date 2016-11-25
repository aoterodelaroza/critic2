!> Parts of code taken from pi7, by VLC, EFM, AMP, MFA, MBV, MAB.
!> (c) Victor Lua~na Cabal, Universidad de Oviedo, 1987--

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

!> Interface to aiPI (pi7) densities.
module pi_private
  private

  public :: pi_end
  public :: pi_register_struct
  public :: pi_read_ion
  public :: pi_rho2
  public :: fillinterpol

  ! parameters for PI common info
  integer, parameter :: msa = 4
  integer, parameter :: mnt = 20
  integer, parameter :: msto = 45
  integer, parameter :: maos = 45
  integer, parameter :: mstosym = 15
  integer, parameter :: mtvef = 8
  integer, parameter :: mtvex = 30
  integer, parameter :: mdens = 360

  ! private structural info
  real*8, allocatable :: renv(:,:)
  integer, allocatable :: idxenv(:), zenv(:)
  integer :: nenv
contains

  ! Terminate the pi arrays
  subroutine pi_end()

    if (allocated(renv)) deallocate(renv)
    if (allocated(idxenv)) deallocate(idxenv)
    if (allocated(zenv)) deallocate(zenv)

  end subroutine pi_end

  !> Register structural information
  subroutine pi_register_struct(nenv0,renv0,idx0,zenv0)

    integer, intent(in) :: nenv0
    real*8, intent(in) :: renv0(:,:)
    integer, intent(in) :: idx0(:), zenv0(:)

    nenv = nenv0
    if (allocated(renv)) deallocate(renv)
    allocate(renv(size(renv0,1),size(renv0,2)))
    renv = renv0
    if (allocated(idxenv)) deallocate(idxenv)
    allocate(idxenv(size(idx0)))
    idxenv = idx0
    if (allocated(zenv)) deallocate(zenv)
    allocate(zenv(size(zenv0)))
    zenv = zenv0
    nenv = nenv0
    
  end subroutine pi_register_struct

  !> From PI: read an ion description file.
  subroutine pi_read_ion(fichero,f,ni)
    use tools_io, only: fopen_read, getline_raw, ferror, faterr, fclose
    use types, only: field
    use param, only: fact, zero
    implicit none

    character*(*)  fichero
    type(field)    f
    integer        ni

    ! parameters
    integer, parameter :: mpar=5

    !.....variables locales:
    character(len=:), allocatable :: linea
    character*6 :: tition
    character(len=1024) :: chpar(mpar)

    integer, dimension(mpar) :: ipar
    integer, parameter :: very_old = -1, pi5_type = 1, pi7_type = 3
    integer :: style, isy
    logical :: openshell, ok
    integer :: lui, iao, nchpar, nipar
    integer :: i, j, k, jj, kk, ll, mm, nn
    integer :: nk2, nj, nk
    real*8 :: deg, xnum, xden, ead, zef, aef, alfaef, aexch, qclas
    real*8 :: ecorr, ekin, epot, zn, eef, eorb, dumyj, dumyk
    integer :: ileeaexch, ndensi, ntaos, ntsto, ncsh, nosh, ntef
    integer :: nref, nxef, nyef, nzef

    ! allocate
    if (.not.allocated(f%pi_used)) then
       allocate(f%pi_used(0:mnt))
       f%pi_used = .false.
    end if
    if (.not.allocated(f%piname)) allocate(f%piname(0:mnt))
    if (.not.allocated(f%nsym)) allocate(f%nsym(0:mnt))
    if (.not.allocated(f%naos)) allocate(f%naos(msa,0:mnt))
    if (.not.allocated(f%naaos)) allocate(f%naaos(msa,0:mnt))
    if (.not.allocated(f%nsto)) allocate(f%nsto(msa,0:mnt))
    if (.not.allocated(f%nasto)) allocate(f%nasto(msa,0:mnt))
    if (.not.allocated(f%nn)) allocate(f%nn(msto,0:mnt))
    if (.not.allocated(f%nelec)) allocate(f%nelec(maos,0:mnt))
    if (.not.allocated(f%z)) allocate(f%z(msto,0:mnt))
    if (.not.allocated(f%xnsto)) allocate(f%xnsto(msto,0:mnt))
    if (.not.allocated(f%c)) allocate(f%c(mstosym,maos,0:mnt))

    f%pi_used(ni) = .true.
    f%piname(ni) = fichero

    !.....abrir el fichero:
    if (ni.gt.mnt) stop 'pi(leerion): el ion por leer es >mnt !'
    if (ni.lt.0) stop 'pi(leerion): ni<0 !'
    lui=fopen_read(fichero)

    !.....Determinar la version de pi a que corresponde el fichero:
    ok = getline_raw(lui,linea,.true.)
    call buscapar (linea,chpar,nchpar,ipar,nipar)
    if (nchpar.le.0) goto 999
    linea=chpar(1)
    style = pi5_type
    if (linea(1:3).eq.'PI7') then
       style = pi7_type
    else if (linea(1:3).ne.'STO'.and.linea(1:3).ne.'GTO'.and.linea(1:4).ne.'CGTO' ) then
       style = very_old
       backspace (lui)
    else if (linea(1:3).eq.'STO') then
    else if (linea(1:3).eq.'GTO') then
       stop 'pi(leerion): sto/gto mixing!'
    else if (linea(1:4).eq.'CGTO') then
       stop 'pi(leerion): sto/gto mixing!'
    else
       call ferror('pi_read_ion','sto/gto incompatible with the global type!',faterr)
       goto 999
    endif

    !.....Leer STOs:
    ok = getline_raw(lui,linea,.true.)
    read (lui,'(a6,x,f10.0)',err=999) tition, zn
    read (lui,*,err=999) f%nsym(ni)
    if (f%nsym(ni).gt.msa) stop 'pi(leerion): nsym > msa !'
    if (f%nsym(ni).lt.0) stop 'pi(leerion): nsym < 0 !'
    read (lui,*,err=999) (f%nsto(i,ni),i=1,f%nsym(ni))
    f%nasto(1,ni)=0  
    ntsto=0
    do i=1,f%nsym(ni)
       ntsto=ntsto+f%nsto(i,ni)
       if (i.gt.1) f%nasto(i,ni)=f%nasto(i-1,ni)+f%nsto(i-1,ni)
       if (f%nsto(i,ni).gt.mstosym) stop 'pi(leerion): demasiados stos en una simetria !'
    end do
    if (ntsto.gt.msto) stop 'pi(leerion): demasiados stos !'
    read (lui,*,err=999) (f%nn(k,ni),k=1,ntsto)
    read (lui,*,err=999) (f%z(k,ni),k=1,ntsto)
    do k=1,ntsto
       nk2=f%nn(k,ni)*2
       f%xnsto(k,ni)=dsqrt((2d0*f%z(k,ni))**(nk2+1)/fact(nk2))
    end do

    !.....Leer informacion orbital:
    read (lui,*,err=999) (f%naos(i,ni),i=1,f%nsym(ni))
    f%naaos(1,ni)=0  
    ntaos=0
    do i=1,f%nsym(ni)
       ntaos=ntaos+f%naos(i,ni)
       if (i.gt.1) f%naaos(i,ni)=f%naaos(i-1,ni)+f%naos(i-1,ni)
    end do
    if (ntaos.gt.maos) stop 'pi(leerion): demasiados aos !'
    read (lui,*,err=999) (f%nelec(k,ni),k=1,ntaos)
    read (lui,*,err=999) (eorb,k=1,ntaos)
    read (lui,*,err=999) (((f%c(k,j+f%naaos(i,ni),ni),k=1,f%nsto(i,ni)),j=1,f%naos(i,ni)),i=1,f%nsym(ni))

    !.....Este fragmento depende de la version de pi:
    if (style .lt. pi7_type) then

       !........versiones previas a la 7:
       read (lui,*,err=999) (dumyj,i=1,10)
       read (lui,*,err=999) (dumyk,i=1,10)
       read (lui,*,err=999) eef,ead
       read (lui,*,err=999) zef,ntef
       if (ntef.gt.mtvef) stop 'pi(leerion): demasiados terminos vef'
       if (ntef.gt.0) then
          read (lui,*,err=999) (aef,nref,nxef,nyef,nzef,alfaef,k=1,ntef)
       endif
    else

       !........version pi7:
       openshell = .false.
       iao = 0
       do isy = 1, f%nsym(ni)
          ncsh = 0
          nosh = 0
          do i = 1, f%naos(isy,ni)
             iao = iao + 1
             deg = (2 * isy - 1) * 2
             if (f%nelec(iao,ni).eq.zero) then
                nosh = nosh + 1
             else if (f%nelec(iao,ni).eq.deg) then
                ncsh = ncsh + 1
             else
                nosh = nosh + 1
                openshell = .true.
             endif
          enddo
          if (nosh .gt. 1) stop 'leerion: more than 1 open sh. per sym!'
       enddo
       if (openshell) then
          read (lui,*,err=999) nj, nk
          do jj = 1, nj
             read (lui,*,err=999) ll,mm,nn,xnum,xden
          enddo
          do kk = 1, nk
             read (lui,*,err=999) ll,mm,nn,xnum,xden
          enddo
       endif
       read (lui,*,err=999) eef, ead, ecorr, ekin, epot
       read (lui,*,err=999) qclas, zef, ntef
       if (ntef.gt.mtvef) stop 'leerion: too many terms in Vef'
       if (ntef.gt.0) then
          do k=1,ntef
             read (lui,*,err=999) aef,nref,alfaef
          enddo
       endif
    endif

    !.....Exchange potential:
    if (style .eq. very_old) then
       ileeaexch=1
    else
       read (lui,*,err=999) ileeaexch
    endif
    ndensi=0
    do i=1,f%nsym(ni)
       ndensi=ndensi+(f%nsto(i,ni)*(f%nsto(i,ni)+1))/2
    end do
    if (ndensi.gt.mdens) stop 'pi(leerion): too many charge densities!'
    if (ndensi.gt.0.and.ileeaexch.eq.1) then
       read (lui,*,err=999) (aexch,i=1,ndensi)
    else
    endif

    !.....cerrar el fichero de datos:
    call fclose (lui)

    return
    !.....error de lectura:
999 close (lui)
    stop 'pi(leerion): read error !'

  end subroutine pi_read_ion

  !xx! PRIVATE functions and subroutines
  subroutine buscapar (line,chpar,nchpar,ipar,nipar)
    implicit none
    
    integer, parameter :: mpar=3
    integer :: ipar(mpar)
    character*(*) :: chpar(mpar)
    character*(*) :: line
    character*(1024) :: word
    character*(1) :: ch
    integer :: nchpar, nipar, nch, i, j

    nchpar=0
    nipar=0
    nch=0
    word=""
    do i=1,mpar
       chpar(i)=word
       ipar(i)=0
    end do

    !.run over all the characters:
    do i = 1, len(line)
       ch=line(i:i)
       if (ch == " ") then

          !.each word is assumed to end with a blanc:
          if (nch.gt.0) then

             !.a new word has been found only if word has any character:
             if (entero(word,j)) then

                !.an integer parameter:
                nipar=nipar+1
                if (nipar.gt.mpar) return
                ipar(nipar)=j
             else

                !.a character parameter:
                nchpar=nchpar+1
                if (nchpar.gt.mpar) return
                chpar(nchpar)=word
                if (nchpar.eq.1) then
                   if (word(1:3).eq.'tit'.or.word(1:3).eq.'rem') then
                      nchpar=nchpar+1
                      chpar(nchpar)=line(5:)
                      return
                   endif
                endif
             endif
             nch=0
             word=" "
          endif
       else
          !.store in word the non-blanc character found:
          nch=nch+1
          word(nch:nch)=ch
       endif
    end do

    ! if the line doesn't end with a blanc, word may contain some
    ! characters:
    if (nch.gt.0) then
       if (entero(word,j)) then
          ! an integer parameter:
          nipar=nipar+1
          if (nipar.gt.mpar) then 
             return
          end if
          ipar(nipar)=j
       else
          !.a character parameter:
          nchpar=nchpar+1
          if (nchpar.gt.mpar) then
             return
          end if
          chpar(nchpar)=word
          if (nchpar.eq.1) then
             if (word(1:3).eq.'tit'.or.word(1:3).eq.'rem') then
                nchpar=nchpar+1
                chpar(nchpar)=line(5:)
                return
             endif
          endif
       endif
    endif

  end subroutine buscapar

  logical function entero (palabra,ipal)
    implicit none
    
    character*(*)     palabra
    character*(1)     cero,nueve,ch
    data cero /'0'/, nueve /'9'/
    integer           ipal, i

      entero=.false.
      ipal=0
      do 10 i = 1, len(palabra)
         ch=palabra(i:i)
         if (ch.lt.cero .or. ch.gt.nueve) return
         ipal=ipal*10+(ichar(ch)-ichar(cero))
 10   continue
      entero=.true.

  end function entero

  !> Fills the interpolation grids for ion densities.
  subroutine fillinterpol(f)
    use global, only: cutrad
    use types, only: field
    implicit none

    type(field), intent(inout) :: f

    real*8, parameter :: az = exp(-6d0)
    real*8, parameter :: b = 0.002
    real*8, parameter :: pi_cutdens = 1d-12 

    integer :: i, ii, j
    real*8 :: crad, rrho, rrho1, rrho2
    logical :: done(mnt)

    allocate(f%pgrid(mnt))
    done = .false.
    do ii = 1, nenv
       i = idxenv(ii)
       if (done(i)) cycle
       done(i) = .true.

       ! determine cutoff radius (crad)
       crad = cutrad(zenv(ii))
       call rhoex1(f,i,crad,rrho,rrho1,rrho2)
       do while (rrho > pi_cutdens) 
          crad = crad * 1.05d0
          call rhoex1(f,i,crad,rrho,rrho1,rrho2)
       end do
       
       ! fill some grid info
       f%pgrid(i)%init = .true.
       f%pgrid(i)%a = az / real(zenv(ii),8)
       f%pgrid(i)%b = b
       f%pgrid(i)%ngrid = ceiling(log(crad/f%pgrid(i)%a) / f%pgrid(i)%b) + 1
       f%pgrid(i)%rmax = f%pgrid(i)%a * exp((f%pgrid(i)%ngrid - 1) * f%pgrid(i)%b)
       f%pgrid(i)%rmax2 = f%pgrid(i)%rmax * f%pgrid(i)%rmax
       allocate(f%pgrid(i)%r(f%pgrid(i)%ngrid))
       allocate(f%pgrid(i)%f(f%pgrid(i)%ngrid))
       allocate(f%pgrid(i)%fp(f%pgrid(i)%ngrid))
       allocate(f%pgrid(i)%fpp(f%pgrid(i)%ngrid))
       do j = 1, f%pgrid(i)%ngrid
          f%pgrid(i)%r(j) = f%pgrid(i)%a * exp((j-1)*f%pgrid(i)%b)
          call rhoex1(f,i,f%pgrid(i)%r(j),f%pgrid(i)%f(j),f%pgrid(i)%fp(j),f%pgrid(i)%fpp(j))
       end do
    end do

  end subroutine fillinterpol

  !> Calculates radial density and its radial derivatives for an atom.
  subroutine rhoex1(f, ni, rion0, rhoval, firstder, secondder)
    use tools_math, only: ep
    use types, only: field
    use param, only: pi, zero, two
    implicit none

    type(field), intent(inout) :: f
    integer, intent(in) :: ni
    real*8, intent(in) :: rion0
    real*8, intent(out) :: rhoval, firstder, secondder

    real*8, parameter :: pi4 = 4d0 * pi
    real*8  :: dumr, dumgr, dumgr2, gradr, grad2r
    real*8  :: or, zj, rion, rion1, rion2
    integer :: nj1, nj2, l, norb, norb1, j, j1

    real*8, parameter :: eps0 = 1d-7

    rhoval = zero
    gradr = zero
    grad2r = zero

    rion = max(rion0,eps0)
    rion1 = 1d0 / rion
    rion2 = rion1 * rion1
    do l = 1, f%nsym(ni)
       do norb = 1, f%naos(l,ni)
          norb1  = norb + f%naaos(l,ni)
          dumr   = zero
          dumgr  = zero
          dumgr2 = zero
          do j = 1, f%nsto(l,ni)
             j1  = j + f%nasto(l,ni)
             nj1 = f%nn(j1,ni)-1
             nj2 = nj1-1
             zj  = f%z(j1,ni)
             or  = ep(rion,nj1) * exp(-zj*rion) * f%c(j,norb1,ni)
             or  = or * f%xnsto(j1,ni)
             dumr = dumr + or
             dumgr = dumgr + or * (nj1 * rion1 - zj)
             dumgr2 = dumgr2 + or * (nj2*nj1*rion2-2*zj*nj1*rion1+zj*zj)
          enddo
          rhoval = rhoval + f%nelec(norb1,ni) * dumr * dumr
          gradr = gradr + f%nelec(norb1,ni) * dumgr * dumr
          grad2r = grad2r + f%nelec(norb1,ni) * (dumgr * dumgr + dumr * dumgr2)
       enddo
    enddo
    rhoval = rhoval/pi4
    firstder = gradr*two/pi4
    secondder = grad2r*two/pi4

  end subroutine rhoex1

  !> Determine the density and derivatives at a given target point
  !> (cartesian).  It is possible to use the 'approximate' method, by
  !> interpolating on a grid.  This routine is thread-safe.
  subroutine pi_rho2 (f,xpos,rho,grad,h)
    use grid1_tools, only: grid1_interp
    use tools_math, only: ep, norm
    use param, only: pi, one
    use types, only: field
    implicit none

    type(field), intent(in) :: f
    real*8, intent(in) :: xpos(3)
    real*8, intent(out) :: rho, grad(3), h(3,3)

    real*8, parameter :: pi4 = 4d0 * pi
    real*8, parameter :: eps = 1d-6
    integer :: ion
    integer :: i, j, k, j1, l
    real*8 :: xion, yion, zion, rions2, rion, rion1, rion2
    real*8 :: xr, yr, zr, x2r, y2r, z2r, xyr, xzr, yzr
    integer :: ni, llplus1, norb, norb1
    real*8 :: phi, phip, phipp, zj, or
    real*8 :: rhop, rhopp, rhofac1, xx2r(3), xxion(3), rfac, radd
    integer :: n0, nm1, nm2
    real*8 :: tmprho

    real*8, parameter :: eps0 = 1d-7

    if (f%exact) then
       ! calculate exactly the contribution of each atom
       rho = 0d0
       grad = 0d0
       h = 0d0
       !.....recorre todos los iones de la red
       !
       do ion= 1, nenv 
          rhop = 0d0
          rhopp = 0d0
          !
          xion = xpos(1) - renv(1,ion)
          yion = xpos(2) - renv(2,ion)
          zion = xpos(3) - renv(3,ion)
          xxion = (/ xion, yion, zion /)
          rions2 = xion*xion+yion*yion+zion*zion
          rion = sqrt(rions2)
          if (rion.gt.eps) then
             !.r powers
             rion1=one/rion
             rion2=rion1*rion1
             !.x powers 
             xr=xion*rion1
             yr=yion*rion1
             zr=zion*rion1
             x2r=xr*xr
             y2r=yr*yr
             z2r=zr*zr
             xx2r = (/ x2r, y2r, z2r /)
             xyr=xr*yr
             xzr=xr*zr
             yzr=yr*zr
          else
             ! trap nuclei
             rho = 0d0
             grad = (/ 0d0, 0d0, 0d0 /)
             h = 0d0
             h(1,1) = -1d0
             h(2,2) = -1d0
             h(3,3) = -1d0
             return
          endif
          !........Every atomic symmetry
          ni= idxenv(ion)
          do l = 1, f%nsym(ni)
             llplus1=l*(l-1)
             !...........every orbital
             do norb = 1, f%naos(l,ni)
                norb1 = norb + f%naaos(l,ni)
                phi = 0d0
                phip = 0d0
                phipp = 0d0

                !..............every primitive
                do j = 1, f%nsto(l,ni)
                   j1 = j + f%nasto(l,ni)
                   n0 = f%nn(j1,ni)
                   nm1 = n0-1   
                   nm2 = n0-2
                   zj = f%z(j1,ni)
                   or=ep(rion,nm1)*exp(-zj*rion)*f%c(j,norb1,ni)
                   or=or*f%xnsto(j1,ni)
                   phi = phi + or
                   phip = phip + or * (nm1*rion1-zj)
                   phipp = phipp + or * (nm2*nm1*rion2-2*zj*nm1*rion1+zj*zj)
                enddo
                rho = rho + f%nelec(norb1,ni) * phi * phi
                rhop = rhop + f%nelec(norb1,ni) * phi * phip
                rhopp = rhopp + f%nelec(norb1,ni) * (phip*phip + phi*phipp)
             enddo
          enddo
          grad = grad + rion1 * rhop * (/ xion, yion, zion /)
          rhofac1 = (rhopp - rhop * rion1)
          do i = 1, 3
             h(i,i) = h(i,i) + rhop * rion1 + xx2r(i) * rhofac1
             do j = i+1,3
                h(i,j) = h(i,j) + xxion(i) * xxion(j) * rion2 * rhofac1
             end do
          end do
       enddo
       rho = rho / pi4
       grad = grad * 2d0 / pi4
       h = h * 2d0 / pi4
       h(2,1) = h(1,2)
       h(3,1) = h(1,3)
       h(3,2) = h(2,3)
    else
       ! use the density grids
       rho = 0d0
       grad = 0d0
       h = 0d0
       do i = 1, nenv
          xxion = xpos - renv(:,i)
          rion = max(norm(xxion),eps0)
          rion1 = 1d0 / rion
          rion2 = rion1 * rion1
          ni = idxenv(i)
          call grid1_interp(f%pgrid(ni),rion,tmprho,rhop,rhopp)
          rho = rho + tmprho
          grad = grad + rhop * xxion * rion1
          rfac = (rhopp - rhop * rion1)
          do j = 1, 3
             h(j,j) = h(j,j) + rhop * rion1 + rfac * rion2 * xxion(j) * xxion(j)
             do k = 1, j-1
                radd = rfac * rion2 * xxion(j) * xxion(k)
                h(j,k) = h(j,k) + radd
                h(k,j) = h(k,j) + radd
             end do
          end do
       end do
       h(2,1) = h(1,2)
       h(3,1) = h(1,3)
       h(3,2) = h(2,3)
    endif

  end subroutine pi_rho2

end module pi_private

