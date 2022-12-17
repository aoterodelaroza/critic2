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

!> Interface to aiPI (pi7) densities and structures.
submodule (pi_private) proc
  implicit none

  !xx! private procedures
  ! subroutine read_ion(f,fichero,ni,ti)
  ! subroutine buscapar(line,chpar,nchpar,ipar,nipar)
  ! function entero (palabra,ipal)
  ! subroutine rhoex1(f,ni,rion0,rhoval,firstder,secondder)

  ! parameters for PI common info
  integer, parameter :: msa = 4
  integer, parameter :: msto = 45
  integer, parameter :: maos = 45
  integer, parameter :: mstosym = 15
  integer, parameter :: mtvef = 8
  integer, parameter :: mtvex = 30
  integer, parameter :: mdens = 360

contains

  ! Terminate the pi arrays
  module subroutine pi_end(f)
    class(piwfn), intent(inout) :: f

    if (allocated(f%bas)) deallocate(f%bas)
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    if (f%isealloc) then
       if (associated(f%e)) deallocate(f%e)
    end if
    nullify(f%e)
    f%isealloc = .false.

  end subroutine pi_end

  !> Build a pi field from external file
  module subroutine pi_read(f,nfile,piat,file,env,errmsg,ti)
    use global, only: cutrad
    use tools_io, only: isinteger, equali
    use param, only: mlen
    class(piwfn), intent(inout) :: f
    integer, intent(in) :: nfile
    character*10, intent(in) :: piat(:)
    character(len=mlen), intent(in) :: file(:)
    type(environ), intent(in), target :: env
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, ithis
    logical :: iok, found
    real*8 :: crad, rrho, rrho1, rrho2

    real*8, parameter :: az = exp(-6d0)
    real*8, parameter :: b = 0.002
    real*8, parameter :: pi_cutdens = 1d-12

    errmsg = ""
    ! restart and allocate all fields
    if (allocated(f%bas)) deallocate(f%bas)
    allocate(f%bas(env%nspc))

    do i = 1, nfile
       iok = isinteger(ithis,piat(i))
       found = .false.
       do j = 1, env%nspc
          if (equali(piat(i),env%spc(j)%name)) then
             call read_ion(f,file(i),j,errmsg,ti=ti)
             if (len_trim(errmsg) > 0) goto 999
             found = .true.
          else if (iok) then
             if (ithis == j) then
                call read_ion(f,file(i),j,errmsg,ti=ti)
                if (len_trim(errmsg) > 0) goto 999
                found = .true.
             end if
          end if
       end do
       if (.not.found) then
          errmsg = "unknown species for pi ion file: " // trim(file(i))
          goto 999
       end if
    end do

    ! find the individual species cutoffs and maximum cutoff
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    allocate(f%spcutoff(env%nspc,2))
    f%spcutoff = 0d0
    f%globalcutoff = -1d0
    do i = 1, env%nspc
       f%spcutoff(i,2) = cutrad(env%spc(i)%z)
       f%globalcutoff = max(f%globalcutoff,f%spcutoff(i,2))
    end do

    ! pointer to the environment
    if (f%isealloc) then
       if (associated(f%e)) deallocate(f%e)
    end if
    nullify(f%e)
    if (f%globalcutoff >= env%dmax0 .and..not.env%ismolecule) then
       ! Create a new environment to satisfy all searches.
       ! The environment contains all the atoms in molecules anyway.
       f%isealloc = .true.
       nullify(f%e)
       allocate(f%e)
       call f%e%extend(env,f%globalcutoff)
    else
       ! keep a pointer to the environment
       f%isealloc = .false.
       f%e => env
    end if

    ! fill the interpolation tables
    do i = 1, env%nspc
       if (.not.f%bas(i)%pi_used) cycle
       ! determine cutoff radius (crad)
       crad = cutrad(env%spc(i)%z)
       call rhoex1(f,i,crad,rrho,rrho1,rrho2)
       do while (rrho > pi_cutdens)
          crad = crad * 1.05d0
          call rhoex1(f,i,crad,rrho,rrho1,rrho2)
       end do

       ! fill some grid info
       f%bas(i)%pgrid%z = env%spc(i)%z
       f%bas(i)%pgrid%qat = 0
       f%bas(i)%pgrid%isinit = .true.
       f%bas(i)%pgrid%a = az / real(env%spc(i)%z,8)
       f%bas(i)%pgrid%b = b
       f%bas(i)%pgrid%ngrid = ceiling(log(crad/f%bas(i)%pgrid%a) / f%bas(i)%pgrid%b) + 1
       f%bas(i)%pgrid%rmax = f%bas(i)%pgrid%a * exp((f%bas(i)%pgrid%ngrid - 1) * f%bas(i)%pgrid%b)
       f%bas(i)%pgrid%rmax2 = f%bas(i)%pgrid%rmax * f%bas(i)%pgrid%rmax
       allocate(f%bas(i)%pgrid%r(f%bas(i)%pgrid%ngrid))
       allocate(f%bas(i)%pgrid%f(f%bas(i)%pgrid%ngrid))
       allocate(f%bas(i)%pgrid%fp(f%bas(i)%pgrid%ngrid))
       allocate(f%bas(i)%pgrid%fpp(f%bas(i)%pgrid%ngrid))
       do j = 1, f%bas(i)%pgrid%ngrid
          f%bas(i)%pgrid%r(j) = f%bas(i)%pgrid%a * exp((j-1)*f%bas(i)%pgrid%b)
          call rhoex1(f,i,f%bas(i)%pgrid%r(j),f%bas(i)%pgrid%f(j),f%bas(i)%pgrid%fp(j),f%bas(i)%pgrid%fpp(j))
       end do
    end do

    return
999 continue ! error condition
    call f%end()

  end subroutine pi_read

  !> Determine the density and derivatives at a given target point
  !> (cartesian).  It is possible to use the 'approximate' method, by
  !> interpolating on a grid.  This routine is thread-safe.
  module subroutine rho2(f,xpos,exact,rho,grad,h)
    use tools_math, only: ep
    use param, only: pi, one, icrd_cart
    class(piwfn), intent(in) :: f
    real*8, intent(in) :: xpos(3)
    logical, intent(in) :: exact
    real*8, intent(out) :: rho, grad(3), h(3,3)

    real*8, parameter :: pi4 = 4d0 * pi
    real*8, parameter :: eps = 1d-6
    real*8, parameter :: eps0 = 1d-7

    integer :: ion
    integer :: i, j, k, j1, l
    real*8 :: xion, yion, zion, rions2, rion, rion1, rion2
    real*8 :: xr, yr, zr, x2r, y2r, z2r, xyr, xzr, yzr
    integer :: ni, llplus1, norb, norb1
    real*8 :: phi, phip, phipp, zj, or
    real*8 :: rhop, rhopp, rhofac1, xx2r(3), xxion(3), rfac, radd
    integer :: n0, nm1, nm2, nenv, ierr
    real*8 :: tmprho
    integer, allocatable :: eid(:)

    ! calculate the environment of the input point
    rho = 0d0
    grad = 0d0
    h = 0d0
    call f%e%list_near_atoms(xpos,icrd_cart,.false.,nenv,ierr,eid,up2dsp=f%spcutoff)
    if (ierr > 0) return ! could happen if in a molecule and very far -> zero

    if (exact) then
       ! calculate exactly the contribution of each atom
       !.....recorre todos los iones de la red
       do ion= 1, nenv
          ni = f%e%at(eid(ion))%is
          if (.not.f%bas(ni)%pi_used) cycle
          rhop = 0d0
          rhopp = 0d0
          !
          xion = xpos(1) - f%e%at(eid(ion))%r(1)
          yion = xpos(2) - f%e%at(eid(ion))%r(2)
          zion = xpos(3) - f%e%at(eid(ion))%r(3)
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
          do l = 1, f%bas(ni)%nsym
             llplus1=l*(l-1)
             !...........every orbital
             do norb = 1, f%bas(ni)%naos(l)
                norb1 = norb + f%bas(ni)%naaos(l)
                phi = 0d0
                phip = 0d0
                phipp = 0d0

                !..............every primitive
                do j = 1, f%bas(ni)%nsto(l)
                   j1 = j + f%bas(ni)%nasto(l)
                   n0 = f%bas(ni)%nn(j1)
                   nm1 = n0-1
                   nm2 = n0-2
                   zj = f%bas(ni)%z(j1)
                   or=ep(rion,nm1)*exp(-zj*rion)*f%bas(ni)%c(j,norb1)
                   or=or*f%bas(ni)%xnsto(j1)
                   phi = phi + or
                   phip = phip + or * (nm1*rion1-zj)
                   phipp = phipp + or * (nm2*nm1*rion2-2*zj*nm1*rion1+zj*zj)
                enddo
                rho = rho + f%bas(ni)%nelec(norb1) * phi * phi
                rhop = rhop + f%bas(ni)%nelec(norb1) * phi * phip
                rhopp = rhopp + f%bas(ni)%nelec(norb1) * (phip*phip + phi*phipp)
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
       do i = 1, nenv
          ni = f%e%at(eid(i))%is
          if (.not.f%bas(ni)%pi_used) cycle

          xxion = xpos - f%e%at(eid(i))%r
          rion = max(norm2(xxion),eps0)
          rion1 = 1d0 / rion
          rion2 = rion1 * rion1
          call f%bas(ni)%pgrid%interp(rion,tmprho,rhop,rhopp)
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

  end subroutine rho2

  !xx! private procedures

  !> Read a PI ion description file from file fichero and species ni.
  subroutine read_ion(f,fichero,ni,errmsg,ti)
    use tools_io, only: fopen_read, getline_raw, fclose
    use param, only: fact, zero
    type(piwfn), intent(inout) :: f
    character*(*) :: fichero
    integer :: ni
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

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
    errmsg = "error reading files"
    f%bas(ni)%pi_used = .true.
    f%bas(ni)%piname = trim(fichero)
    if (allocated(f%bas(ni)%naos)) deallocate(f%bas(ni)%naos)
    if (allocated(f%bas(ni)%naaos)) deallocate(f%bas(ni)%naaos)
    if (allocated(f%bas(ni)%nsto)) deallocate(f%bas(ni)%nsto)
    if (allocated(f%bas(ni)%nasto)) deallocate(f%bas(ni)%nasto)
    if (allocated(f%bas(ni)%nn)) deallocate(f%bas(ni)%nn)
    if (allocated(f%bas(ni)%z)) deallocate(f%bas(ni)%z)
    if (allocated(f%bas(ni)%xnsto)) deallocate(f%bas(ni)%xnsto)
    if (allocated(f%bas(ni)%c)) deallocate(f%bas(ni)%c)
    if (allocated(f%bas(ni)%nelec)) deallocate(f%bas(ni)%nelec)
    allocate(f%bas(ni)%naos(msa))
    allocate(f%bas(ni)%naaos(msa))
    allocate(f%bas(ni)%nsto(msa))
    allocate(f%bas(ni)%nasto(msa))
    allocate(f%bas(ni)%nn(msto))
    allocate(f%bas(ni)%z(msto))
    allocate(f%bas(ni)%xnsto(msto))
    allocate(f%bas(ni)%c(mstosym,maos))
    allocate(f%bas(ni)%nelec(maos))

    !.....abrir el fichero:
    lui=fopen_read(fichero,ti=ti)
    if (lui < 0) goto 999

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
       errmsg = 'sto/gto incompatible with the global type!'
       goto 999
    endif

    !.....Leer STOs:
    ok = getline_raw(lui,linea,.true.)
    read (lui,'(a6,x,f10.0)',err=999,end=999) tition, zn
    read (lui,*,err=999,end=999) f%bas(ni)%nsym
    if (f%bas(ni)%nsym.gt.msa) stop 'pi(leerion): nsym > msa !'
    if (f%bas(ni)%nsym.lt.0) stop 'pi(leerion): nsym < 0 !'
    read (lui,*,err=999,end=999) (f%bas(ni)%nsto(i),i=1,f%bas(ni)%nsym)
    f%bas(ni)%nasto(1)=0
    ntsto=0
    do i=1,f%bas(ni)%nsym
       ntsto=ntsto+f%bas(ni)%nsto(i)
       if (i.gt.1) f%bas(ni)%nasto(i)=f%bas(ni)%nasto(i-1)+f%bas(ni)%nsto(i-1)
       if (f%bas(ni)%nsto(i).gt.mstosym) stop 'pi(leerion): demasiados stos en una simetria !'
    end do
    if (ntsto.gt.msto) stop 'pi(leerion): demasiados stos !'
    read (lui,*,err=999,end=999) (f%bas(ni)%nn(k),k=1,ntsto)
    read (lui,*,err=999,end=999) (f%bas(ni)%z(k),k=1,ntsto)
    do k=1,ntsto
       nk2=f%bas(ni)%nn(k)*2
       f%bas(ni)%xnsto(k)=dsqrt((2d0*f%bas(ni)%z(k))**(nk2+1)/fact(nk2))
    end do

    !.....Leer informacion orbital:
    read (lui,*,err=999,end=999) (f%bas(ni)%naos(i),i=1,f%bas(ni)%nsym)
    f%bas(ni)%naaos(1)=0
    ntaos=0
    do i=1,f%bas(ni)%nsym
       ntaos=ntaos+f%bas(ni)%naos(i)
       if (i.gt.1) f%bas(ni)%naaos(i)=f%bas(ni)%naaos(i-1)+f%bas(ni)%naos(i-1)
    end do
    if (ntaos.gt.maos) stop 'pi(leerion): demasiados aos !'
    read (lui,*,err=999,end=999) (f%bas(ni)%nelec(k),k=1,ntaos)
    read (lui,*,err=999,end=999) (eorb,k=1,ntaos)
    read (lui,*,err=999,end=999) (((f%bas(ni)%c(k,j+f%bas(ni)%naaos(i)),k=1,f%bas(ni)%nsto(i)),&
       j=1,f%bas(ni)%naos(i)),i=1,f%bas(ni)%nsym)

    !.....Este fragmento depende de la version de pi:
    if (style .lt. pi7_type) then

       !........versiones previas a la 7:
       read (lui,*,err=999,end=999) (dumyj,i=1,10)
       read (lui,*,err=999,end=999) (dumyk,i=1,10)
       read (lui,*,err=999,end=999) eef,ead
       read (lui,*,err=999,end=999) zef,ntef
       if (ntef.gt.mtvef) stop 'pi(leerion): demasiados terminos vef'
       if (ntef.gt.0) then
          read (lui,*,err=999,end=999) (aef,nref,nxef,nyef,nzef,alfaef,k=1,ntef)
       endif
    else

       !........version pi7:
       openshell = .false.
       iao = 0
       do isy = 1, f%bas(ni)%nsym
          ncsh = 0
          nosh = 0
          do i = 1, f%bas(ni)%naos(isy)
             iao = iao + 1
             deg = (2 * isy - 1) * 2
             if (f%bas(ni)%nelec(iao).eq.zero) then
                nosh = nosh + 1
             else if (f%bas(ni)%nelec(iao).eq.deg) then
                ncsh = ncsh + 1
             else
                nosh = nosh + 1
                openshell = .true.
             endif
          enddo
          if (nosh .gt. 1) stop 'leerion: more than 1 open sh. per sym!'
       enddo
       if (openshell) then
          read (lui,*,err=999,end=999) nj, nk
          do jj = 1, nj
             read (lui,*,err=999,end=999) ll,mm,nn,xnum,xden
          enddo
          do kk = 1, nk
             read (lui,*,err=999,end=999) ll,mm,nn,xnum,xden
          enddo
       endif
       read (lui,*,err=999,end=999) eef, ead, ecorr, ekin, epot
       read (lui,*,err=999,end=999) qclas, zef, ntef
       if (ntef.gt.mtvef) stop 'leerion: too many terms in Vef'
       if (ntef.gt.0) then
          do k=1,ntef
             read (lui,*,err=999,end=999) aef,nref,alfaef
          enddo
       endif
    endif

    !.....Exchange potential:
    if (style .eq. very_old) then
       ileeaexch=1
    else
       read (lui,*,err=999,end=999) ileeaexch
    endif
    ndensi=0
    do i=1,f%bas(ni)%nsym
       ndensi=ndensi+(f%bas(ni)%nsto(i)*(f%bas(ni)%nsto(i)+1))/2
    end do
    if (ndensi.gt.mdens) stop 'pi(leerion): too many charge densities!'
    if (ndensi.gt.0.and.ileeaexch.eq.1) then
       read (lui,*,err=999,end=999) (aexch,i=1,ndensi)
    else
    endif

    !.....cerrar el fichero de datos:
    call fclose (lui)
    errmsg = ""

    return

    !.....error de lectura:
999 continue
    if (lui > 0) close (lui)

  end subroutine read_ion

  subroutine buscapar(line,chpar,nchpar,ipar,nipar)
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

  function entero(palabra,ipal)
    logical :: entero
    character*(*) :: palabra
    integer :: ipal

    character*(1)     cero,nueve,ch
    data cero /'0'/, nueve /'9'/
    integer           i

      entero=.false.
      ipal=0
      do 10 i = 1, len(palabra)
         ch=palabra(i:i)
         if (ch.lt.cero .or. ch.gt.nueve) return
         ipal=ipal*10+(ichar(ch)-ichar(cero))
 10   continue
      entero=.true.

  end function entero

  !> Calculates radial density and its radial derivatives for an atom.
  subroutine rhoex1(f,ni,rion0,rhoval,firstder,secondder)
    use tools_math, only: ep
    use param, only: pi, zero, two
    class(piwfn), intent(in) :: f
    integer, intent(in) :: ni
    real*8, intent(in) :: rion0
    real*8, intent(out) :: rhoval, firstder, secondder

    real*8, parameter :: pi4 = 4d0 * pi
    real*8, parameter :: eps0 = 1d-7

    real*8  :: dumr, dumgr, dumgr2, gradr, grad2r
    real*8  :: or, zj, rion, rion1, rion2
    integer :: nj1, nj2, l, norb, norb1, j, j1

    rhoval = zero
    gradr = zero
    grad2r = zero

    rion = max(rion0,eps0)
    rion1 = 1d0 / rion
    rion2 = rion1 * rion1
    do l = 1, f%bas(ni)%nsym
       do norb = 1, f%bas(ni)%naos(l)
          norb1  = norb + f%bas(ni)%naaos(l)
          dumr   = zero
          dumgr  = zero
          dumgr2 = zero
          do j = 1, f%bas(ni)%nsto(l)
             j1  = j + f%bas(ni)%nasto(l)
             nj1 = f%bas(ni)%nn(j1)-1
             nj2 = nj1-1
             zj  = f%bas(ni)%z(j1)
             or  = ep(rion,nj1) * exp(-zj*rion) * f%bas(ni)%c(j,norb1)
             or  = or * f%bas(ni)%xnsto(j1)
             dumr = dumr + or
             dumgr = dumgr + or * (nj1 * rion1 - zj)
             dumgr2 = dumgr2 + or * (nj2*nj1*rion2-2*zj*nj1*rion1+zj*zj)
          enddo
          rhoval = rhoval + f%bas(ni)%nelec(norb1) * dumr * dumr
          gradr = gradr + f%bas(ni)%nelec(norb1) * dumgr * dumr
          grad2r = grad2r + f%bas(ni)%nelec(norb1) * (dumgr * dumgr + dumr * dumgr2)
       enddo
    enddo
    rhoval = rhoval/pi4
    firstder = gradr*two/pi4
    secondder = grad2r*two/pi4

  end subroutine rhoex1

end submodule proc
