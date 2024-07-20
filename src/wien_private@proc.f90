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

!> Interface to WIEN2k densities and structures
submodule (wien_private) proc
  implicit none

  !xx! private procedures
  ! subroutine wien_read_struct(f,file,errmsg,ti)
  ! subroutine readslm(f,lu,errmsg)
  ! subroutine readk(f,lu,errmsg)
  ! subroutine gbass(rbas,gbas)
  ! subroutine rotdef(f)
  ! subroutine gener(f,errmsg)
  ! subroutine sternb(f)
  ! subroutine rotator(vt,iz,tau,a)
  ! subroutine rotato(vt,iz,tau)
  ! subroutine rotat(vt,rotloc)
  ! subroutine reduc(v,atp,npos,pos,iat,rmt)
  ! function vnorm(v)
  ! subroutine charge(f,nder,chg,grad,hess,ir,r,jatom,v,natnr)
  ! subroutine radial(f,rlm,rho,rho1,rho2,r,ir,jatom)
  ! subroutine rhoout(f,v,nder,chg,grad,hess)
  ! subroutine ylm(v,lmax,y)

  integer, parameter  :: ncom = 122 !< maximum number of LM pairs
  integer, parameter  :: nrad = 781 !< maximum value of npt
  integer, parameter  :: nsym = 48 !< max. number of symmetry operations
  integer, parameter  :: lmax2 = 14 !< maximum l for LM lattice harmonics expansion of charte inside muffins.

  integer, parameter :: nsa = 3, nsb = 3, nsc = 3
  integer, parameter :: nnpos = (2*nsa+1)*(2*nsb+1)*(2*nsc+1) !< maximum number of cell origins

  ! pw cutoff
  real*8, parameter :: pwcutoff = 1d-30

  ! radial grid derivation formulas
  real*8, parameter :: coef1(6,3) = reshape((/&
     -274,  600, -600,  400,  -150,  24,&
        6,  -60,  -40,  120,   -30,   4,&
      -24, 150, -400,  600,  -600, 274 /),shape(coef1))
  integer, parameter :: noef(6,3) = reshape((/&
      0,  1,  2,  3,  4,  5,&
     -2, -1,  0,  1,  2,  3,&
     -5, -4, -3, -2, -1,  0/),shape(noef))
  real*8, parameter :: coef2(6,3) = reshape((/&
     225, -770,  1070,  -780,   305,   -50,&
      -5,   80,  -150,    80,    -5,     0,&
     -50,  305,  -780 ,  1070,  -770,   225/),shape(coef2))
  real*8, parameter :: coef3(6,3) = reshape((/&
     -85,  355,  -590,   490,  -205,    35,&
     -5 ,  -5 ,  50  ,  -70 ,   35 ,    -5,&
     -35,  205,  -490,  590 ,  -355,    85/),shape(coef3))
  real*8, parameter :: coef4(6,3) = reshape((/&
     15,  -70,   130,  -120,    55,   -10,&
      5,  -20,    30,   -20,     5,     0,&
    -10,   55,  -120,   130,   -70,    15/),shape(coef4))
  real*8, parameter :: fac1=1d0/120d0
  real*8, parameter :: fac2=2d0/120d0
  real*8, parameter :: fac3=6d0/120d0
  real*8, parameter :: fac4=24d0/120d0

contains

  module subroutine wien_end(f)
    class(wienwfn), intent(inout) :: f

    if (allocated(f%lm)) deallocate(f%lm)
    if (allocated(f%lmmax)) deallocate(f%lmmax)
    if (allocated(f%slm)) deallocate(f%slm)
    if (allocated(f%sk)) deallocate(f%sk)
    if (allocated(f%ski)) deallocate(f%ski)
    if (allocated(f%tauk)) deallocate(f%tauk)
    if (allocated(f%tauki)) deallocate(f%tauki)
    if (allocated(f%rotloc)) deallocate(f%rotloc)
    if (allocated(f%rnot)) deallocate(f%rnot)
    if (allocated(f%rmt)) deallocate(f%rmt)
    if (allocated(f%dx)) deallocate(f%dx)
    if (allocated(f%jri)) deallocate(f%jri)
    if (allocated(f%multw)) deallocate(f%multw)
    if (allocated(f%iatnr)) deallocate(f%iatnr)
    if (allocated(f%pos)) deallocate(f%pos)
    if (allocated(f%iop)) deallocate(f%iop)
    if (allocated(f%iz)) deallocate(f%iz)
    if (allocated(f%tau)) deallocate(f%tau)
    if (allocated(f%krec)) deallocate(f%krec)
    if (allocated(f%inst)) deallocate(f%inst)

  end subroutine wien_end

  ! Given the field f, returns the Rmt of the atom at crystallographic
  ! position x.
  module function rmt_atom(f,x)
    class(wienwfn), intent(in) :: f
    real*8, intent(in) :: x(3)
    real*8 :: rmt_atom

    integer :: ipos, iat, ii, jatom
    real*8 :: vt(3), r, v(3)

    jatom = 0
    do ii = 1, 3
       v(ii) = f%br1(1,ii)*x(1)+f%br1(2,ii)*x(2)+f%br1(3,ii)*x(3)
    end do
    do ipos=1,f%npos
       do iat=1,f%ndat
          do ii=1,3
             vt(ii)=v(ii)-(f%atp(ii,ipos)+f%pos(ii,iat))
          enddo
          r=sqrt(vt(1)**2+vt(2)**2+vt(3)**2)
          jatom=iabs(f%iatnr(iat))
          if(r .lt. f%rmt(jatom)) then
             goto 11
          endif
       enddo
    enddo
11  continue

    rmt_atom = f%rmt(jatom)

  end function rmt_atom

  ! Read a clmsum (file) and struct (file2) file and returns the
  ! wien2k field f.
  module subroutine read_clmsum(f,file,file2,errmsg,ti)
    use tools_io, only: fopen_read, fclose
    class(wienwfn), intent(inout) :: f
    character*(*), intent(in) :: file, file2
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, ntot, i, j

    errmsg = ""

    ! load field data from the struct file and save it in f
    call wien_read_struct(f,file2,errmsg,ti)
    if (len_trim(errmsg) > 0) return

    ! read the clmsum
    lu = fopen_read(file,ti=ti)

    call readslm(f,lu,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    call readk(f,lu,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    call fclose(lu)

    ntot = 0
    do i = 1, f%nwav
       do j = 1, f%inst(i)
          ntot = ntot + 1
       end do
    end do

    return
999 continue ! error condition
    if (lu > 0) call fclose(lu)

  end subroutine read_clmsum

  !> Calculate the density and its derivatives at point v0
  !> (cryst. coord.) up to the nder derivative. Returns
  !> the density (rho), gradient (grad), and Hessian (h).
  !> This routine is thread-safe.
  module subroutine rho2(f,v0,nder,rho,grad,h)
    class(wienwfn), intent(in) :: f
    real*8, dimension(3), intent(in) :: v0
    integer, intent(in) :: nder
    real*8, intent(out) :: rho, grad(3), h(3,3)

    logical                :: insphere
    real*8, dimension(3)   :: vt,vt1,v(3)
    real*8, dimension(3,3) :: mat, temp, mhess
    real*8                 :: gx,gy,gz
    real*8                 :: r
    integer                :: ir
    real*8                 :: chg, gradd(3), hess(6)
    integer                :: ipos,iat,jatom,jatom1,ii,i,j,j1,k,m,n,l
    integer                :: index

    !.transform v to cartesian coordinates
    do i = 1, 3
       v(i) = f%br1(1,i)*v0(1)+f%br1(2,i)*v0(2)+f%br1(3,i)*v0(3)
    end do

    !.test if v falls in a muffin tin sphere
    insphere=.false.
    ! loop over cell origins
    do ipos=1,f%npos
       !    loop over different atoms
       do iat=1,f%ndat
          !       check distance to atom iat in cell with origin ipos
          !       mark insphere as true if v is inside muffin tin
          do  ii=1,3
             vt(ii)=v(ii)-(f%atp(ii,ipos)+f%pos(ii,iat))
          enddo
          r=sqrt(vt(1)**2+vt(2)**2+vt(3)**2)
          jatom=iabs(f%iatnr(iat))
          if(r .lt. f%rmt(jatom)) then
             insphere=.true.
             goto 11
          endif
       enddo
    enddo
11  continue

    rho = 0d0
    grad = 0d0
    hess = 0d0
    gradd = 0d0
    mhess = 0d0
    !.calculate properties according to the position of v
    ! point in interstitial
    if (.not.insphere) then
       call rhoout(f,v,nder,chg,gradd,hess)
       if (nder >= 2) then
          mhess(1,1)=hess(1)
          mhess(2,1)=hess(2)
          mhess(3,1)=hess(3)
          mhess(2,2)=hess(4)
          mhess(3,2)=hess(5)
          mhess(3,3)=hess(6)
          mhess(1,2)=mhess(2,1)
          mhess(1,3)=mhess(3,1)
          mhess(2,3)=mhess(3,2)
       end if
       ! point in muffin tin of iatnr(iat)
    else
       !    identify type of atom
       jatom=iabs(f%iatnr(iat))
       vt = v
       !    i am going to transfer coordinates from sphere iat to jatom
       !    rotate vt so as to make it fall in the equivalent point
       !    inside first of equivalent's atom muffin tin, instead of iat.

       !    case of orthogonal unit cell axes
       if(f%ortho) then

          !       This can be done only for ortho, because unit cell axes coincide
          !       with orthogonal axes, i. e. br1 = id*a, br3 = id*1/a.
          !       beware! rotator accepts a 3x3 array as iz and a 3x1 array as tau
          call rotator(vt,f%iz(1,1,f%iop(iat)),f%tau(1,f%iop(iat)),f%a)

          !       store the transpose of the matrix used in mat
          do m=1,3
             do n=1,3
                mat(m,n)=f%iz(n,m,f%iop(iat))
             enddo
          enddo
          !    case of non-orthogonal unit cell axes
       else
          !       convert v to crystallographic coordinates (tocrsytal())
          do j1=1,3
             vt1(j1)=f%br3(j1,1)*vt(1)+f%br3(j1,2)*vt(2)+f%br3(j1,3)*vt(3)
          enddo
          !       mat = br3
          do m=1,3
             do n=1,3
                mat(m,n)=f%br3(m,n)
             enddo
          enddo
          !       rotate v so as to fall in the equivalent point inside the
          !       muffin of an image of the first equivalent atom of the list (jatom = iatnr(iat)).

          call rotato(vt1,f%iz(:,:,f%iop(iat)),f%tau(:,f%iop(iat)))

          !       mat = (iz)^t * br3
          do m=1,3
             do n=1,3
                temp(m,n)=0d0
                do l=1,3
                   temp(m,n)=temp(m,n)+f%iz(l,m,f%iop(iat))*mat(l,n)
                enddo
             enddo
          enddo
          do m=1,3
             do n=1,3
                mat(m,n)=temp(m,n)
             enddo
          enddo
          !       convert v back to cartesian coordinates

          do  j1=1,3
             vt(j1)=f%br1(1,j1)*vt1(1)+f%br1(2,j1)*vt1(2)+f%br1(3,J1)*vt1(3)
          enddo

          !       mat = (br1)^t * (iz)^t * br3
          do m=1,3
             do n=1,3
                temp(m,n)=0d0
                do l=1,3
                   temp(m,n)=temp(m,n)+f%br1(l,m)*mat(l,n)
                enddo
             enddo
          enddo
          do m=1,3
             do n=1,3
                mat(m,n)=temp(m,n)
             enddo
          enddo
       end if

       !....find the unit cell's atom list index for the first atom in the equivalent
       !    list (jatom), and assign its value to index.
       index=1
       jatom1=1
       do while (jatom1.ne.jatom)
          index=index+f%multw(jatom1)
          jatom1=jatom1+1
       enddo

       !....v is inside the muffin tin of atom index, convert v (cartesian coords) to
       !    a vector referred to the system of reference of the closest atom.
       !          v = v - pos(closest atom)
       !    this doesn't affect the gradient

       call reduc(vt,f%atp,f%npos,f%pos,index,f%rmt(jatom))

       !....transfer coordinates into local coordinate system
       !    using local rotation matrix

       call rotat(vt,f%rotloc(:,:,jatom))

       !... mat = (rotloc) * (br1)^t * (iz)^t * (br3)
       !
       do m=1,3
          do n=1,3
             temp(m,n)=0d0
             do l=1,3
                temp(m,n)=temp(m,n)+f%rotloc(m,l,jatom)*mat(l,n)
             enddo
          enddo
       enddo
       do m=1,3
          do n=1,3
             mat(m,n)=temp(m,n)
          enddo
       enddo

       !....calculate charge density properties inside muffin tin.
       !    if r is very close, r = r0 and ir=1. Else, use the
       !    formula in the above comments
       if(r.lt.f%rnot(jatom)) then
          ir=1
          r=f%rnot(jatom)
       else
          ir=1+int(log(r/f%rnot(jatom))/f%dx(jatom))
          if(ir.lt.1) ir=1
       endif

       call charge(f,nder,chg,gradd,hess,ir,r,jatom,vt,f%iatnr(iat))

       if (nder >= 1) then
          !....reconstruct gradient and hessian in original setting.
          !    x' = mat * x  -->  gr = (mat)^t * gr'
          if(r < f%rnot(jatom)) then
             gradd = 0d0
          else
             gx=mat(1,1)*gradd(1)+mat(2,1)*gradd(2)+mat(3,1)*gradd(3)
             gy=mat(1,2)*gradd(1)+mat(2,2)*gradd(2)+mat(3,2)*gradd(3)
             gz=mat(1,3)*gradd(1)+mat(2,3)*gradd(2)+mat(3,3)*gradd(3)
             gradd(1)=gx
             gradd(2)=gy
             gradd(3)=gz
          end if

          if (nder >= 2) then
             !    store the values for the hessian
             mhess(1,1)=hess(1)
             mhess(2,1)=hess(2)
             mhess(3,1)=hess(3)
             mhess(2,2)=hess(4)
             mhess(3,2)=hess(5)
             mhess(3,3)=hess(6)
             mhess(1,2)=mhess(2,1)
             mhess(1,3)=mhess(3,1)
             mhess(2,3)=mhess(3,2)

             ! (amp: it can be optimized further)
             ! x' = mat * x --> h = (mat)^t * h' * mat
             do i=1,3
                do j=1,i
                   temp(i,j)=0d0
                   do k=1,3
                      do l=1,3
                         temp(i,j)=temp(i,j)+mhess(k,l)*mat(k,i)*mat(l,j)
                      enddo
                   enddo
                   if (j.ne.i) temp(j,i)=temp(i,j)
                enddo
             enddo
             do i=1,3
                do j=1,3
                   mhess(i,j)=temp(i,j)
                enddo
             enddo
          end if
       end if
    endif

    rho = chg
    if (nder >= 1) grad = matmul(f%br1,gradd)
    if (nder >= 2) h = matmul(matmul(f%br1,mhess),transpose(f%br1))

  end subroutine rho2

  !> Convert a given wien2k scalar field into its Laplacian.
  module subroutine tolap(f)
    use param, only: tpi2
    class(wienwfn), intent(inout) :: f

    integer :: jatom, lmmx, ilm, ir
    real*8 :: rho, rho1, rho2
    integer :: l, lp1, i, kpp, niz
    real*8 :: r, r2, tpik2
    real*8, allocatable :: rgrid(:)
    real*8 :: factor(3)

    ! atomic spheres
    allocate(rgrid(size(f%slm,1)))
    do jatom = 1, f%nat
       lmmx = f%lmmax(jatom)

       do ilm = 1, lmmx
          l=abs(f%lm(1,ilm,jatom))
          lp1 = l + 1

          ! save the old radial grid
          rgrid = f%slm(:,ilm,jatom)

          ! run over all radial points
          do ir = 1, f%jri(jatom)
             r = f%rnot(jatom)*exp((ir-1)*f%dx(jatom))
             r2 = r * r
             call radial(f,rgrid,rho,rho1,rho2,r,ir,jatom)
             f%slm(ir,ilm,jatom) = -l*lp1*rho + 2d0*r*rho1 + rho2 * r2
          end do
       end do
    end do
    deallocate(rgrid)

    ! interstitial
    if (f%ortho) then
       factor(1)=1d0/f%a(1)
       factor(2)=1d0/f%a(2)
       factor(3)=1d0/f%a(3)
    else
       factor(1)=1d0
       factor(2)=1d0
       factor(3)=1d0
    endif
    kpp = 0
    niz=f%lastniz
    do i = f%lastind,1,-1
       kpp = kpp + 1
       tpik2 = -dot_product(f%krec(:,i),f%krec(:,i)) * tpi2
       if(kpp > f%inst(niz)) then
          kpp=1
          niz=niz-1
       end if
       if (kpp /= 1) cycle
       tpik2 = -dot_product(f%krec(:,i)*factor,f%krec(:,i)*factor) * tpi2
       f%sk(niz) = f%sk(niz) * tpik2
       if (f%cmpl) then
          f%ski(niz) = f%ski(niz) * tpik2
       end if
    end do

  end subroutine tolap

  !xx! private procedures

  !> Read the crystal structure from a STRUCT file, with logical unit lut.
  !> Store WIEN2k variables (rmt, etc.) in the module.
  subroutine wien_read_struct(f,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose
    use param, only: pi
    class(wienwfn), intent(inout) :: f
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lut
    character*80 :: titel
    character*4 :: cform
    real*8 :: alphatemp,betatemp,gammatemp
    real*8 :: localcalph, localcbeta, localcgamm
    real*8 :: localsalph, localsbeta, localsgamm
    integer :: index
    integer :: jatom, mu
    character*10 :: aname
    real*8 :: znuc, cosg1, gamma0, ar2, car2
    integer :: i, j, i1, i2, j1, nato, ndif
    real*8 :: aix, aiy, aiz
    character*4 :: lattic

    integer, dimension(:), allocatable :: temp_iatnr
    real*8, dimension(:,:), allocatable :: temp_pos

    errmsg = "Error reading file: " // trim(file)
    lut = fopen_read(file,ti=ti)
    if (lut < 0) goto 999
    READ(lut,102,err=999,end=999) TITEL
    READ(lut,103,err=999,end=999) LATTIC,f%NAT,cform
    f%ishlat = (lattic.eq.'H   ')

102 FORMAT(A80)
103 FORMAT(A4,23X,I3,1x,a4,/,4X,4X) ! new

    nato = f%nat
    ndif = 48*nato
    if (allocated(f%rotloc)) deallocate(f%rotloc)
    if (allocated(f%rnot)) deallocate(f%rnot)
    if (allocated(f%rmt)) deallocate(f%rmt)
    if (allocated(f%dx)) deallocate(f%dx)
    if (allocated(f%jri)) deallocate(f%jri)
    if (allocated(f%multw)) deallocate(f%multw)
    allocate (f%rotloc(3,3,nato))
    allocate (f%rnot(nato),f%rmt(nato),f%dx(nato))
    allocate (f%jri(nato),f%multw(nato))
    allocate (temp_iatnr(ndif))
    allocate (temp_pos(3,ndif))

    READ(lut,100,err=999,end=999) f%a,alphatemp,betatemp,gammatemp
100 FORMAT(6F10.5)
    if(gammatemp == 0.d0) gammatemp=90.d0

    ! write cell parameters
    if(gammatemp .eq. 0.d0) then
       gammatemp=90.d0
    end if
    localcalph = cos(alphatemp*pi/180d0)
    localcbeta = cos(betatemp*pi/180d0)
    localcgamm = cos(gammatemp*pi/180d0)
    localsalph = sin(alphatemp*pi/180d0)
    localsbeta = sin(betatemp*pi/180d0)
    localsgamm = sin(gammatemp*pi/180d0)

    !.calculate br1, br2, ortho, centering vectors
    f%br1 = 0d0
    f%br2 = 0d0
    f%br3 = 0d0
    IF(LATTIC(1:1).EQ.'S'.OR.LATTIC(1:1).EQ.'P') THEN
       cosg1=(localcgamm-localcalph*localcbeta)/(localsalph*localsbeta)
       gamma0=acos(cosg1)
       F%BR2(1,1)=f%A(1)*sin(gamma0)*localsbeta
       f%br2(1,2)=f%a(1)*cos(gamma0)*localsbeta
       f%br2(1,3)=f%a(1)*localcbeta
       F%BR2(2,2)=f%A(2)*localsalph
       F%BR2(2,3)=f%A(2)*localcalph
       F%BR2(3,3)=f%A(3)
       F%BR1(1,1)=f%A(1)*sin(gamma0)*localsbeta
       F%BR1(1,2)=f%A(1)*cos(gamma0)*localsbeta
       f%br1(1,3)=f%a(1)*localcbeta
       F%BR1(2,2)=f%A(2)*localsalph
       F%BR1(2,3)=f%A(2)*localcalph
       F%BR1(3,3)=f%A(3)
       f%ortho=.true.
       if(gammatemp.ne.90.d0) f%ortho=.false.
       if(alphatemp.ne.90.d0) f%ortho=.false.
       if(betatemp .ne.90.d0) f%ortho=.false.
    ELSE IF(LATTIC(1:1).EQ.'F') THEN
       F%BR2(1,1)=0.5*f%A(1)
       F%BR2(2,1)=0.5*f%A(1)
       F%BR2(2,2)=0.5*f%A(2)
       F%BR2(3,2)=0.5*f%A(2)
       F%BR2(1,3)=0.5*f%A(3)
       F%BR2(3,3)=0.5*f%A(3)
       F%BR1(1,1)=f%A(1)
       F%BR1(2,2)=f%A(2)
       F%BR1(3,3)=f%A(3)
       f%ortho=.true.
    ELSE IF(LATTIC(1:1).EQ.'B') THEN
       F%BR2(1,1)=-0.5d0*F%A(1)
       F%BR2(2,1)=0.5d0*F%A(1)
       F%BR2(3,1)=0.5d0*F%A(1)
       F%BR2(1,2)=0.5d0*F%A(2)
       F%BR2(2,2)=-0.5d0*F%A(2)
       F%BR2(3,2)=0.5d0*F%A(2)
       F%BR2(1,3)=0.5d0*F%A(3)
       F%BR2(2,3)=0.5d0*F%A(3)
       F%BR2(3,3)=-0.5d0*F%A(3)
       F%BR1(1,1)=F%A(1)
       F%BR1(2,2)=F%A(2)
       F%BR1(3,3)=F%A(3)
       f%ortho=.true.
    ELSE IF(LATTIC(1:1).EQ.'H') THEN
       F%BR1(1,1)=SQRT(3.d0)/2.d0*F%A(1)
       F%BR1(1,2)=-0.5d0*F%A(2)
       F%BR1(2,2)=F%A(2)
       F%BR1(3,3)=F%A(3)
       F%BR2(1,1)=SQRT(3.d0)/2.d0*F%A(1)
       F%BR2(1,2)=-0.5d0*F%A(2)
       F%BR2(2,2)=F%A(2)
       F%BR2(3,3)=F%A(3)
       f%ortho=.false.
    ELSE IF(LATTIC(1:1).EQ.'R') THEN
       F%BR1(1,1)=1.d0/SQRT(3.d0)/2.d0*F%A(1)
       F%BR1(1,2)=-0.5d0*F%A(2)
       F%BR1(1,3)=1.d0/3.d0*F%A(3)
       F%BR1(2,1)=1.d0/SQRT(3.d0)/2.d0*F%A(1)
       F%BR1(2,2)=0.5d0*F%A(2)
       F%BR1(2,3)=1.d0/3.d0*F%A(3)
       F%BR1(3,1)=-1.d0/SQRT(3.d0)*F%A(1)
       F%BR1(3,2)=0.d0
       F%BR1(3,3)=1.d0/3.d0*F%A(3)
       F%BR2(1,1)=1.d0/SQRT(3.d0)/2.d0*F%A(1)
       F%BR2(1,2)=-0.5d0*F%A(2)
       F%BR2(1,3)=1.d0/3.d0*F%A(3)
       F%BR2(2,1)=1.d0/SQRT(3.d0)/2.d0*F%A(1)
       F%BR2(2,2)=0.5d0*F%A(2)
       F%BR2(2,3)=1.d0/3.d0*F%A(3)
       F%BR2(3,1)=-1.d0/SQRT(3.d0)*F%A(1)
       F%BR2(3,2)=0.d0
       F%BR2(3,3)=1.d0/3.d0*F%A(3)
       f%ortho=.false.
       ar2 = (3d0 * F%A(1)*F%A(1) + F%A(3)*F%A(3)) / 9d0
       car2 = (-1.5d0 * F%A(1)*F%A(1) + F%A(3)*F%A(3)) / (9d0 * ar2)
    ELSE IF(LATTIC(1:3).EQ.'CXY') THEN
       F%BR2(1,1)=0.5d0*F%A(1)
       F%BR2(2,1)=0.5d0*F%A(1)
       F%BR2(3,1)=0.d0
       F%BR2(1,2)=0.5d0*F%A(2)
       F%BR2(2,2)=-0.5d0*F%A(2)
       F%BR2(3,2)=0.d0
       F%BR2(1,3)=0.d0
       F%BR2(2,3)=0.d0
       F%BR2(3,3)=F%A(3)
       F%BR1(1,1)=F%A(1)
       F%BR1(2,2)=F%A(2)
       F%BR1(3,3)=F%A(3)
       f%ortho=.true.
    ELSE IF(LATTIC(1:3).EQ.'CYZ') THEN
       F%BR2(1,1)=F%A(1)
       F%BR2(2,1)=0.d0
       F%BR2(3,1)=0.d0
       F%BR2(1,2)=0.d0
       F%BR2(2,2)=-0.5d0*F%A(2)
       F%BR2(3,2)=0.5d0*F%A(2)
       F%BR2(1,3)=0.d0
       F%BR2(2,3)=0.5d0*F%A(3)
       F%BR2(3,3)=0.5d0*F%A(3)
       F%BR1(1,1)=F%A(1)
       F%BR1(2,2)=F%A(2)
       F%BR1(3,3)=F%A(3)
       f%ortho=.true.
    ELSE IF(LATTIC(1:3).EQ.'CXZ') THEN
       F%BR2(1,1)=0.5d0*F%A(1)*localsgamm
       F%BR2(1,2)=0.5d0*F%A(1)*localcgamm
       F%BR2(1,3)=-0.5d0*F%A(3)
       F%BR2(2,2)=F%A(2)
       F%BR2(3,1)=0.5d0*F%A(1)*localsgamm
       F%BR2(3,2)=0.5d0*F%A(1)*localcgamm
       F%BR2(3,3)=0.5d0*F%A(3)
       F%BR1(1,1)=F%A(1)*localsgamm
       f%br1(1,2)=f%a(1)*localcgamm
       F%BR1(2,2)=F%A(2)
       F%BR1(3,3)=F%A(3)
       f%ortho=.false.
    ELSE
       STOP 'LATTIC NOT DEFINED'
    END IF

    !.calculate br3
    call gbass(f%br1,f%br3)

    !.read positions and set nneq, pos, iatnr, multw, equivalent variables in
    ! critic, etc.
    INDEX=0
    DO JATOM=1,f%NAT
       INDEX=INDEX+1
       READ(lut,1012,err=999,end=999) temp_IATNR(INDEX),&
          TEMP_POS(1,INDEX),TEMP_POS(2,INDEX),TEMP_POS(3,INDEX),f%MULTW(JATOM)
       DO MU=1,f%MULTW(JATOM)-1
          INDEX=INDEX+1
          READ(lut,1013,err=999,end=999) temp_IATNR(INDEX),&
             TEMP_POS(1,INDEX),TEMP_POS(2,INDEX),TEMP_POS(3,INDEX)
       end DO

       READ(lut,113,err=999,end=999) ANAME,f%JRI(JATOM),f%RNOT(JATOM),f%RMT(JATOM),Znuc
       if (f%jri(jatom) > nrad) then
          errmsg = "Max. radial points (nrad) exceeded"
          goto 999
       end if
       f%dx(jatom) = log(f%rmt(jatom)/f%rnot(jatom)) / (f%jri(jatom)-1)
       READ(lut,1051,err=999,end=999) ((f%ROTLOC(I1,J1,JATOM),I1=1,3),J1=1,3)
    end DO
113 FORMAT(A10,5X,I5,5X,F10.5,5X,F10.5,5X,F5.2)
1012 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/15X,I2) ! new
1013 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7) ! new
1051 FORMAT(20X,3F10.8)

    !.total number of atoms in unit cell
    f%NDAT=INDEX
    if (allocated(f%iatnr)) deallocate(f%iatnr)
    if (allocated(f%iop)) deallocate(f%iop)
    if (allocated(f%pos)) deallocate(f%pos)
    allocate (f%iatnr(f%ndat),f%iop(f%ndat))
    allocate (f%pos(3,f%ndat))
    f%iatnr(:) = temp_iatnr(1:f%ndat)
    f%pos(:,:) = temp_pos(:,1:f%ndat)
    deallocate(temp_iatnr,temp_pos)

    !.read number of symmetry operations, sym. operations
    READ(lut,114,err=999,end=999) f%nIORD
114 FORMAT(I4)

    if (allocated(f%iz)) deallocate(f%iz)
    if (allocated(f%tau)) deallocate(f%tau)
    allocate(f%iz(3,3,nsym))
    allocate(f%tau(3,nsym))
    DO I=1,f%nIORD
       READ(lut,115,err=999,end=999) ((f%IZ(I1,I2,I),I1=1,3),f%TAU(I2,I),I2=1,3)
    end DO
115 FORMAT(3(3I2,F10.5,/))

    !.find iop matrices
    CALL ROTDEF(f,lattic)

    !.transform pos to cartesian coordinates
    DO INDEX=1,f%NDAT
       AIX=f%POS(1,INDEX)
       AIY=f%POS(2,INDEX)
       AIZ=f%POS(3,INDEX)
       DO J=1,3
          f%POS(J,INDEX)=F%BR1(1,J)*AIX+F%BR1(2,J)*AIY+F%BR1(3,J)*AIZ
       end DO
    end DO

    !.generate origins of next unit cells (atp)
    CALL GENER(f,errmsg)

    call fclose(lut)

    errmsg = ""
    return
999 continue
    if (lut > 0) call fclose(lut)

  end subroutine wien_read_struct

  !> Read the atomic sphere part of a clmsum style file
  subroutine readslm(f,lu,errmsg)
    use param, only: sqfp
    class(wienwfn), intent(inout) :: f
    integer, intent(in) :: lu !< Input logical unit
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: jatom, jrj, l, l1, l2, ll, i
    integer :: mlmmax, mjrj
    integer, allocatable :: alm(:,:,:)
    integer, allocatable :: almmax(:)
    real*8, allocatable :: aslm(:,:,:)

    errmsg = ""
    if (allocated(f%lm)) deallocate(f%lm)
    if (allocated(f%lmmax)) deallocate(f%lmmax)
    if (allocated(f%slm)) deallocate(f%slm)
    allocate(f%lm(2,ncom,f%nat),f%lmmax(f%nat),f%slm(nrad,ncom,f%nat))

    READ(lu,2032,err=999,end=999)

    mjrj = 0
    mlmmax = 0
    f%slm = 0d0
    DO JATOM=1,f%NAT
       JRJ=f%JRI(JATOM)
       mjrj = max(mjrj,jrj)
       READ(lu,118,err=999,end=999) LL
       mlmmax = max(mlmmax,LL)
       f%lmmax(jatom)=LL
       if (f%lmmax(jatom) > ncom) then
          errmsg = "Max. LM pairs (ncom) exceeded"
          return
       end if
       do l=1,ll
          READ(lu,2010,err=999,end=999) l1, l2
          f%lm(1,l,jatom) = l1
          f%lm(2,l,jatom) = l2
          if (f%lm(1,l,jatom) > lmax2) then
             errmsg = "Max. L in LM list (lmax2) exceeded"
             return
          end if
          READ(lu,2021,err=999,end=999) (f%slm(I,l,jatom),I=1,JRJ)
          READ(lu,2031,err=999,end=999)
          if (l == 1) then
             f%slm(1:jrj,l,jatom)=f%slm(1:jrj,l,jatom)/sqfp
          end if
       end do
       READ(lu,2033,err=999,end=999)
    end DO
    f%cnorm = .true.  ! we have density normalization (/sqfp above)

    ! scale back
    allocate(alm(2,mlmmax,f%nat))
    alm = f%lm(1:2,1:mlmmax,1:f%nat)
    call move_alloc(alm,f%lm)
    allocate(almmax(f%nat))
    almmax = f%lmmax(1:f%nat)
    call move_alloc(almmax,f%lmmax)
    allocate(aslm(mjrj,mlmmax,f%nat))
    aslm = f%slm(1:mjrj,1:mlmmax,1:f%nat)
    call move_alloc(aslm,f%slm)

118 FORMAT(/,15X,I3,//) ! new
2032 FORMAT(//)
2010 FORMAT(15X,I3,5X,I2,/) ! new
2021 format(3x,4e19.12)
2031 FORMAT(/)
2033 FORMAT(///)

999 continue

  end subroutine readslm

  !> Read the plane wave part of a clmsum-style file
  subroutine readk(f,lu,errmsg)
    class(wienwfn), intent(inout) :: f
    integer, intent(in) :: lu !< Input logical unit
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: ind, i
    real*8 :: ttmpk, ttmpki
    real*8 :: vt1(3), vt(3)
    integer :: j1, ii, jj, j
    integer, allocatable :: k2(:,:)
    real*8 :: rkmod
    integer :: indmax, k1(3), nst, istg(3,nsym)
    real*8, allocatable :: taup(:)
    real*8, allocatable :: taupi(:)

    ! read number of pw and allocate space for kvectors
    errmsg = ""
    IND=0
    READ(lu,117,err=999,end=999) f%nwav
    if (allocated(f%krec)) deallocate(f%krec)
    if (allocated(f%inst)) deallocate(f%inst)
    if (allocated(k2)) deallocate(k2)
    if (allocated(f%sk)) deallocate(f%sk)
    if (allocated(f%ski)) deallocate(f%ski)
    allocate(f%krec(3,f%nwav*nsym),f%inst(f%nwav))
    allocate(k2(3,f%nwav))
    allocate(f%sk(f%nwav))
    f%sk = 0d0

    if (allocated(f%ski)) deallocate(f%ski)
    f%cmpl = .false.
    do i = 1, f%nwav
       read(lu,2071,END=36,err=999) (K2(J,i),J=1,3), ttmpk, ttmpki
       if (abs(ttmpki) > pwcutoff .and..not.f%cmpl) then
          f%cmpl = .true.
          allocate(f%ski(f%nwav))
          f%ski = 0d0
       end if
       f%sk(i) = ttmpk
       if (f%cmpl) f%ski(i) = ttmpki
    end do
36  continue

    ! continue
    if (allocated(f%tauk)) deallocate(f%tauk)
    if (allocated(f%tauki)) deallocate(f%tauki)
    allocate (taup(nsym))
    allocate (f%tauk(f%nwav*nsym))
    allocate (taupi(nsym))
    allocate (f%tauki(f%nwav*nsym))

    f%lastind = 0
    do i = 1, f%nwav
       do j=1,3
          k1(j)=k2(j,i)
       enddo
       CALL STERNB(f,k1,nst,istg,taup,taupi)
       rkmod = 0d0
       f%sk(i) = f%sk(i) / nst
       rkmod = max(rkmod,abs(f%sk(i)))
       if (f%cmpl) then
          f%ski(i) = f%ski(i) / nst
          rkmod = max(rkmod,abs(f%ski(i)))
       end if
       ! copy taup to tauk, istg to krec, nst to inst. Index is converted to int,
       ! which assumes values 1...nst(1),nst(1)+1,...,nst(1)+nst(2),...,nst*nk
       DO J=1,NST
          IND=IND+1
          f%tauk(ind) = taup(j)
          if (f%cmpl) f%tauki(IND)=taupi(J)
          IF(IND.GT.f%NWAV*NSYM) then
             errmsg = " ind > nwav*nsym"
             return
          end IF
          DO JJ=1,3
             f%KREC(JJ,IND)=ISTG(JJ,J)
          end DO
       end DO
       f%INST(I)=NST
       if (rkmod > pwcutoff) then
          f%lastind = ind
          f%lastniz = i
       end if
    end do

    INDMAX=IND
    !.for non-orthogonal lattice convert to orthogonal coordinates
    IF(.not.f%ortho) THEN
       DO II=1,INDMAX
          VT(1)=f%krec(1,ii)
          VT(2)=f%krec(2,ii)
          VT(3)=f%krec(3,ii)
          IF(f%ishlat) THEN
             F%KREC(1,II)=F%KREC(1,II)*2.d0/SQRT(3.d0)+F%KREC(2,II)/SQRT(3.d0)
          ELSE
             F%KREC(1,II)=vt(1)*1.d0/SQRT(3.d0)+vt(2)*1./SQRT(3.d0)-2.d0/SQRT(3.d0)*vt(3)
             F%KREC(2,II)=-vt(1)+vt(2)
             F%KREC(3,II)=vt(1)+vt(2)+vt(3)
          END IF
          do J1=1,3
             Vt1(J1)=f%BR3(1,j1)*vt(1)+f%BR3(2,j1)*vt(2)+f%BR3(3,J1)*vt(3)
             f%krec(j1,ii)=vt1(j1)
          end do
       end DO
    END IF

    deallocate(taup)
    deallocate(taupi)
    deallocate(k2)

117 FORMAT(//,13X,I6)
2071 FORMAT(3X,3I5,2E19.12)

    return
999 continue ! error condition
    errmsg = "Error reading file (k-points, clmsum)"

  end subroutine readk

  ! PRIVATE from wien2k
  ! The following routines have been adapted from the WIEN2k source by
  ! P. Blaha, K. Schwarz et al.
  subroutine gbass(rbas,gbas)
    real*8, dimension(3,3) :: rbas, gbas
    real*8 :: det
    integer :: i, j

    GBAS(1,1)=RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)
    GBAS(2,1)=RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)
    GBAS(3,1)=RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)
    GBAS(1,2)=RBAS(2,3)*RBAS(3,1)-RBAS(3,3)*RBAS(2,1)
    GBAS(2,2)=RBAS(3,3)*RBAS(1,1)-RBAS(1,3)*RBAS(3,1)
    GBAS(3,2)=RBAS(1,3)*RBAS(2,1)-RBAS(2,3)*RBAS(1,1)
    GBAS(1,3)=RBAS(2,1)*RBAS(3,2)-RBAS(3,1)*RBAS(2,2)
    GBAS(2,3)=RBAS(3,1)*RBAS(1,2)-RBAS(1,1)*RBAS(3,2)
    GBAS(3,3)=RBAS(1,1)*RBAS(2,2)-RBAS(2,1)*RBAS(1,2)
    DET=0.D0
    DO I=1,3
       DET=DET+GBAS(I,1)*RBAS(I,1)
    end DO
    DO I=1,3
       DO J=1,3
          GBAS(I,J)=GBAS(I,J)/DET
       end DO
    end DO
  end subroutine gbass

  subroutine rotdef(f,lattic)
    ! new version of rotdef, adapted from rotdef1.f wien2k_08.3
    ! replaces old rotdef (rotdef.f)
    class(wienwfn), intent(inout) :: f
    character*4, intent(in) :: lattic

    real*8 :: toler, one, toler2
    integer :: index, ncount, jatom, index1, m, i, j, i1
    real*8 :: x, y, z, x1, y1, z1

    DATA TOLER/ 1.d-4/,ONE/1.d0/
    toler2=toler/2.d0
    INDEX=0
    NCOUNT=0
    DO 20 JATOM=1,f%NAT
         INDEX1=INDEX+1
         DO 30 M=1,f%MULTW(JATOM)
            INDEX=INDEX+1
            DO 25 I=1,f%nIORD
            X=0.
            DO J=1,3
               X=X+f%IZ(J,1,I)*f%POS(J,INDEX)
            END DO
            X=X+f%TAU(1,I) + 5.d0
            X= MOD (X+toler2,ONE)-toler2
            Y=0.
            DO J=1,3
               Y=Y+f%IZ(J,2,I)*f%POS(J,INDEX)
            END DO
            Y=Y+f%TAU(2,I) + 5.d0
            Y= MOD (Y+toler2,ONE)-toler2
            Z=0.
            DO J=1,3
               Z=Z+f%IZ(J,3,I)*f%POS(J,INDEX)
            END DO
            Z=Z+f%TAU(3,I) + 5.d0
            Z= MOD (Z+toler2,ONE)-toler2
            X1=ABS(X-f%POS(1,INDEX1))
            Y1=ABS(Y-f%POS(2,INDEX1))
            Z1=ABS(Z-f%POS(3,INDEX1))
            x1=min(x1,abs(x1-one))
            y1=min(y1,abs(y1-one))
            z1=min(z1,abs(z1-one))
!         WRITE(*,*) 'JATOM,INDEX1,INDEX,I',JATOM,INDEX1,INDEX,I
!         WRITE(*,*) X1,Y1,Z1
!         WRITE(*,*) X,Y,Z
            IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN
            NCOUNT=NCOUNT+1
            f%IOP(INDEX)=I
            GOTO 30
            END IF
            if(lattic(1:1).eq.'B') then
              x1=mod(x1+0.500000001d0,one)
              y1=mod(y1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN
              NCOUNT=NCOUNT+1
              f%IOP(INDEX)=I
              GOTO 30
              END IF
            end if
           if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
              x1=mod(x1+0.500000001d0,one)
              y1=mod(y1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN
              NCOUNT=NCOUNT+1
              f%IOP(INDEX)=I
              GOTO 30
              END IF
              x1=mod(x1+0.5d0,one)
              y1=mod(y1+0.5d0,one)
            endif
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
              x1=mod(x1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN
              NCOUNT=NCOUNT+1
              f%IOP(INDEX)=I
              GOTO 30
              END IF
              x1=mod(x1+0.5d0,one)
              z1=mod(z1+0.5d0,one)
            endif
            if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
              y1=mod(y1+0.500000001d0,one)
              z1=mod(z1+0.500000001d0,one)
              IF(X1.LT.TOLER.AND.Y1.LT.TOLER.AND.Z1.LT.TOLER) THEN
              NCOUNT=NCOUNT+1
              f%IOP(INDEX)=I
              GOTO 30
              END IF
            end if
 25         CONTINUE
            WRITE(*,*) 'NO SYMMETRY OPERATION FOUND IN ROTDEF'
            WRITE(*,*) JATOM,INDEX
            WRITE(*,*) (f%POS(I1,JATOM),I1=1,3)
            WRITE(*,*) (f%POS(I1,INDEX),I1=1,3)
            STOP 'ROTDEF'

!
 30      CONTINUE
 20   CONTINUE
      IF(NCOUNT.NE.INDEX) THEN
         WRITE(6,1000) NCOUNT
         STOP ' ROTDEF: NCOUNT NE INDEX'
      ENDIF
      RETURN
 1000 FORMAT(///,3X,'ERROR IN ROTDEF: ROTIJ NOT DEFINED FOR ALL ',       &
       'ATOMS OF BASIS',/,20X,'NCOUNT=',I2)
   end subroutine rotdef

   SUBROUTINE GENER(f,errmsg)
    class(wienwfn), intent(inout) :: f
    character(len=:), allocatable, intent(out) :: errmsg

    integer             :: ja, jb, jc
    integer             :: i

    errmsg = ""
    f%NPOS=0
    JA=-NSA-1
 10 JA=JA+1
    IF(JA.GT.NSA) GOTO 1
    JB=-NSB-1
 11 JB=JB+1
    IF(JB.GT.NSB) GOTO 10
    JC=-NSC-1
 12 JC=JC+1
    IF(JC.GT.NSC) GOTO 11
    f%NPOS=f%NPOS+1
    if (f%npos > nnpos) then
       errmsg = "Max. cell origins (nnpos) exceeded"
       return
    end if
    DO I=1,3
       f%ATP(I,f%NPOS)=f%BR2(1,I)*JA+f%BR2(2,I)*JB+f%BR2(3,I)*JC
    END DO
    GOTO 12
  1 CONTINUE
    RETURN
  END SUBROUTINE GENER

  SUBROUTINE STERNB(f,k1,nst,istg,taup,taupi)
    use param, only: tpi
    class(wienwfn), intent(in) :: f
    integer, intent(in) :: k1(3)
    integer, intent(out) :: nst
    integer, intent(out) :: istg(3,nsym)
    real*8, intent(inout) :: taup(*)
    real*8, intent(inout) :: taupi(*)

    ! generate star of g
    ! taken from wien2k_08.3
    ! replaces old stern.
    integer :: ind(nsym)
    integer :: i, j, k, l, m
    real*8 :: tk

    NST=0
    DO I=1,f%nIORD
       TK=0.
       DO J=1,3
          TK=TK+f%TAU(J,I)*k1(J)*TPI
          K=0
          DO L=1,3
             K=f%IZ(J,L,I)*k1(L)+K
          end do
          ISTG(J,I)=K
       end DO
       IF(NST.EQ.0) GOTO 7
       DO M=1,NST
          DO J=1,3
             IF(ISTG(J,M).NE.ISTG(J,I)) GOTO 4
          end DO
          IND(M)=IND(M)+1
          taup(m) = taup(m) + cos(tk)
          if (f%cmpl) taupi(m) = taupi(M) + sin(tk)
          GOTO 1
4         continue
       end DO
7      NST=NST+1
       DO J=1,3
          ISTG(J,NST)=ISTG(J,I)
       end DO
       IND(NST)=1
       taup(nst) = cos(tk)
       if (f%cmpl) taupi(NST)=sin(tk)
1      continue
    end do

    DO I=1,NST
       taup(i) = taup(i) / ind(i)
       if (f%cmpl) taupi(i)=taupi(i)/ind(i)
    end DO

  end subroutine sternb

  SUBROUTINE ROTATOR(VT,IZ,TAU,A)
    integer :: jc, j
    real*8 :: vt(3), tau(3), vn(3), a(3), dotpro
    integer :: iz(3,3)
    ! SYNOPSIS
    !     Rotates the real vector VT by integer rotation matrix iz.
    !     Then, translates according to tau, a translation vector in
    !     crystallographic coords, scaled to atomic units through a, which
    !     usually contains cell parameters
    !        (vt(1) vt(2) vt(3)) = IZ * (VT(1) VT(2) VT(3))^t + <tau * a>
    !        vt(i) = (IZ*VT)(i) + tau(i) * a(i)
    !
    ! INPUT
    !    vt                       vector to rotate, in crystallographic coords.
    !    iz                       integer rotation matrix
    !    tau                      translation vector in crystallographic coords.
    !    a                        dimension vector
    !
    ! OUTPUT
    !    vt                       rotated vector
    !
    ! LOCAL VARIABLES
    !    vn                       temp vt
    !    dotpro                   dot product
    !    j,jc                     indexes
    !
      DO JC=1,3
         DOTPRO=0.0
         DO J=1,3
            DOTPRO=DOTPRO+VT(J)*IZ(J,Jc)
         END DO
         VN(JC)=DOTPRO+TAU(JC)*A(JC)
      END DO
      DO J=1,3
         VT(J)= VN(J)
      END DO
      RETURN
  END SUBROUTINE ROTATOR

  SUBROUTINE ROTATO(VT,IZ,TAU)
    real*8 :: vt(3), tau(3),vn(3), dotpro
    integer :: iz(3,3), jc, j

    ! SYNOPSIS
    !     Rotates the real vector VT
    !     vt = IZ * VT + tau
    !
    ! INPUT
    !     vt                      vector to rotate, in crystallograhpic coords.
    !     iz                      integer rotation matrix
    !     tau                     additional vector, in crystallograhpic coords.
    !
    ! OUTPUT
    !     vt                      rotated vector, in crystallograhpic coords.
    !
    ! LOCAL VARIABLES
    !     vn                      temp vt
    !     dotpro                  dot product
    !     j,jc                    indexes
    !
    DO JC=1,3
       DOTPRO=0.0
       DO J=1,3
          DOTPRO=DOTPRO+VT(J)*IZ(J,Jc)
       END DO
       VN(JC)=DOTPRO+TAU(JC)
    END DO
    DO J=1,3
       VT(J)= VN(J)
    END DO
    RETURN
  END SUBROUTINE ROTATO

  SUBROUTINE ROTAT(VT,ROTLOC)
    ! SYNOPSIS
    !     Rotates real vector vt with local rotation matrix
    !     uses transposed matrix, because r(-1) is read in struct
    !     vt = rotloc * VT
    !
    ! INPUT
    !     vt                      vector to rotate
    !     rotloc                  local rotation matrix
    ! OUTPUT
    !     vt                      rotated vector
    !
    ! LOCAL VARIABLES
    !     vn                      temp vt
    !     dotpro                  dot product
    !     j,jc                    indexes
    !
    real*8 :: VT(3),ROTLOC(3,3),VN(3), dotpro
    integer :: jc, j
      DO JC=1,3
         DOTPRO=0.0
         DO J=1,3
            DOTPRO=DOTPRO+VT(J)*ROTLOC(JC,J)
         END DO
         VN(JC)=DOTPRO
      END DO
      DO J=1,3
         VT(J)= VN(J)
      END DO
      RETURN
  END SUBROUTINE ROTAT

  SUBROUTINE REDUC(V,ATP,NPOS,POS,IAT,rmt)
    use tools_io, only: uout
    ! INPUT
    !     v                       vector to reduce (cartesian coordinates). Must belong to
    !     some atom's muffin tin.
    !     atp                     see wien.inc
    !     npos                    see wien.inc
    !     pos                     see wien.inc
    !     iat                     index (ndif) of atom in the unit cell
    !     rmt                     see wien.inc
    !
    ! OUTPUT
    !     v                       reduced vector
    !
    ! LOCAL VARIABLES
    !     vhelp                   temp variable for v
    !     vtest                   temp variable for testing v.
    !
    !
    integer :: j, i, npos, iat
    real*8 :: V(3),ATP(3,*),VHELP(3),POS(3,*),VTEST(3), rmt, r, rtest

      DO J=1,3
         VHELP(J)=V(J)
      END DO
!      R=VNORM(V)
      R=99999.
      DO 10 I=1,NPOS
      DO J=1,3
         VTEST(J)=V(J)-ATP(J,I)-POS(J,IAT)
      END DO
      RTEST=VNORM(VTEST)
      IF(RTEST.LT.R) THEN
      R=RTEST
      DO J=1,3
         VHELP(J)=VTEST(J)
      END DO
      END IF
  10  CONTINUE
      DO J=1,3
         V(J)=VHELP(J)
      END DO
      if(r.gt.rmt+1d-4) then
         write(uout,*) ' Reduction failed, r-reduc gt.rmt',r,rmt, abs(r-rmt)
      endif
      RETURN
  END SUBROUTINE REDUC

  FUNCTION VNORM(V)
    REAL*8 :: vnorm
! INPUT
!     v                       input vector
!
! OUTPUT
!     vnorm                   norm of v
!
! LOCAL VARIABLES
!     i                       iteration counter
!
      real*8 :: V(3)
      integer :: i

      VNORM=0.0
      DO I=1,3
         VNORM=VNORM+V(I)**2
      END DO
      VNORM=SQRT(VNORM)
      RETURN
  END FUNCTION VNORM

  ! Calculate the density and derivatives in the muffins, at
  ! point v (cryst. coords.). Calculate derivatives up to nder.
  ! Return density (chg), gradient (grad) and Hessian (hess) or
  ! zero if nder is not high enough.
  subroutine charge(f,nder,chg,grad,hess,ir,r,jatom,v,natnr)
    use tools_math, only: ylmderiv
    use tools_io, only: uout
    use param, only: c_kub
    class(wienwfn), intent(in) :: f
    integer, intent(in) :: nder
    real*8, intent(out) :: chg
    real*8, dimension(3), intent(out) :: grad
    real*8, dimension(6), intent(out) :: hess
    integer, intent(in) :: ir
    real*8,intent(in) :: r
    integer, intent(in) :: jatom
    real*8, dimension(3), intent(in) :: v
    integer, intent(in) :: natnr

    real*8 :: sqin2
    integer :: i, j
    integer :: l, m, m1, minu
    integer :: lmmx
    integer :: ilm
    integer :: idx
    real*8 :: c1, c2, c3
    real*8 :: frho, frho1, frho2

    real*8, dimension(ncom) :: rho, ang
    real*8, dimension(ncom) :: rho1, rho2

    complex*16 :: imag0
    complex*16, dimension(3) :: ygrad, ygrad1, ygrad2, ygradtemp
    complex*16, dimension(6) :: yhess, yhess1, yhess2, yhesstemp

    complex*16, dimension((lmax2+1)**2) :: yl !< spherical harmonics

    sqin2=1d0/sqrt(2d0)

    !.obtain all spherical harmonics and derivatives needed.
    call ylm(v,lmax2,yl)
    lmmx=f%lmmax(jatom)

    !.obtain radial clm's and derivatives, for every lm in lm list.
    do ilm=1,lmmx
       call radial(f,f%slm(:,ilm,jatom),rho(ilm),rho1(ilm),rho2(ilm),r,ir,jatom)
    end do

    !.initiailze output variables
    chg = 0d0
    grad = 0d0
    hess = 0d0

    !.noncubic lattice
    if (natnr.le.0) then
       !    run over lm pairs
       do ilm=1,lmmx

          ! orig: it will be welcome to adapt this routine to more simple forms
          l=f%lm(1,ilm,jatom)
          m=f%lm(2,ilm,jatom)
          m1=(-1)**m
          if (l.ge.0) then
             imag0=(1d0,0d0)*sqin2*m1       !(-1)**m 1/sq2(y+(-1)**m Y*)
             minu=m1                        !(-1)**m 1/isq2(y-(-1)**m Y*)
          else
             imag0=(0d0,-1d0)*sqin2*m1
             minu=-m1
          endif
          l=abs(l)
          idx=l*(l+1)+1

          ! for m=0, ylm = slm. for m!=0, calculate real harmonics.
          if (m.eq.0) then
             ang(ilm)=dble(yl(idx))       !yl0
             call ylmderiv(yl,r,l,m,rho(ilm),rho1(ilm),rho2(ilm),nder,ygrad,yhess)
          else
             ! real harmonics, slm and derivatives
             ang(ilm)=dble(imag0*(yl(idx+m)+minu*yl(idx-m)))  !slm
             call ylmderiv(yl,r,l, m,rho(ilm),rho1(ilm),rho2(ilm),nder,ygrad,yhess)
             call ylmderiv(yl,r,l,-m,rho(ilm),rho1(ilm),rho2(ilm),nder,ygrad1,yhess1)
             !real harmonics derivatives transformation
             do i=1,6
                if (i.le.3) ygrad(i)=imag0*(ygrad(i)+minu*ygrad1(i))
                yhess(i)=            imag0*(yhess(i)+minu*yhess1(i))
             enddo
          endif
          ! add to density, gradient and hessian
          chg=chg+rho(ilm)*ang(ilm)
          do i=1,6
             if (i.le.3) grad(i)=grad(i)+dble(ygrad(i))
             hess(i)=hess(i)+dble(yhess(i))
          enddo
       enddo
       !.cubic lattice. wien treats cubic lattices in a special way.
       ! see kara and kurki-suonio, acta cryst. a 37 (1981) 201.
    else

       ! calculate real harmonics from spherical harmonics
       do i=1,lmmx
          l=f%lm(1,i,jatom)
          m=f%lm(2,i,jatom)
          m1=(-1)**m
          if (l.ge.0) then
             imag0=(1d0,0d0)*sqin2*m1       !(-1)**m 1/sq2(y+(-1)**m Y*)
             minu=m1                        !(-1)**m 1/isq2(y-(-1)**m Y*)
          else
             imag0=(0d0,-1d0)*sqin2*m1
             minu=-m1
          endif
          l=abs(l)
          idx=l*(l+1)+1
          if (m == 0) then
             ang(i)=dble(yl(idx))
          else
             ang(i)=dble(imag0*(yl(idx+m)+minu*yl(idx-m)))
          end if
       end do

       ! build up rho (adapted from lapw5, wien2k_06), grad and hess
       i = 1
       do while (i <= lmmx)
          l = f%lm(1,i,jatom)
          m = f%lm(2,i,jatom)
          m1=(-1)**m
          if (l.ge.0) then
             imag0=(1d0,0d0)*sqin2*m1       !(-1)**m 1/sq2(y+(-1)**m Y*)
             minu=m1                        !(-1)**m 1/isq2(y-(-1)**m Y*)
          else
             imag0=(0d0,-1d0)*sqin2*m1
             minu=-m1
          endif
          if(l.eq.0 .and. m.eq.0) then
             call ylmderiv(yl,r,0,0,rho(i),rho1(i),rho2(i),nder,ygrad,yhess)
             chg=chg + (rho(i) * ang(i))
             do j=1,6
                if (j.le.3) grad(j)=grad(j)+dble(ygrad(j))
                hess(j)=hess(j)+dble(yhess(j))
             enddo
             i=i+1
          elseif (l.eq.-3 .and. m.eq.2) then
             call ylmderiv(yl,r,3, 2,rho(i),rho1(i),rho2(i),nder,ygrad,yhess)
             call ylmderiv(yl,r,3,-2,rho(i),rho1(i),rho2(i),nder,ygradtemp,yhesstemp)
             ygrad = imag0*(ygrad+minu*ygradtemp)
             yhess = imag0*(yhess+minu*yhesstemp)

             chg=chg + (rho(i) * ang(i))
             grad = grad + dble(ygrad)
             hess = hess + dble(yhess)

             i=i+1
          elseif (l.eq.4 .or. l.eq.6 .or. l.eq.-7 .or. l.eq.-9) then
             c1=c_kub(abs(l),m)
             c2=c_kub(abs(l),m+4)

             frho  = c1*rho(i) + c2*rho(i+1)
             frho1 = c1*rho1(i) + c2*rho1(i+1)
             frho2 = c1*rho2(i) + c2*rho2(i+1)

             call ylmderiv(yl,r,abs(l),m,frho,frho1,frho2,nder,ygrad,yhess)
             if (m /= 0) then
                call ylmderiv(yl,r,abs(l),-m,frho,frho1,frho2,nder,ygradtemp,yhesstemp)
                ygrad = imag0*(ygrad+minu*ygradtemp)
                yhess = imag0*(yhess+minu*yhesstemp)
             end if
             call ylmderiv(yl,r,abs(l),m+4,frho,frho1,frho2,nder,ygrad1,yhess1)
             if (m+4 /= 0) then
                call ylmderiv(yl,r,abs(l),-(m+4),frho,frho1,frho2,nder,ygradtemp,yhesstemp)
                ygrad1 = imag0*(ygrad1+minu*ygradtemp)
                yhess1 = imag0*(yhess1+minu*yhesstemp)
             end if
             chg=chg + (c1*ang(i) + c2*ang(i+1)) * frho
             grad = grad + dble(c1*ygrad + c2*ygrad1)
             hess = hess + dble(c1*yhess + c2*yhess1)
             i=i+2
          elseif (l.eq.8 .or. l.eq.10) then
             c1=c_kub(abs(l),m)
             c2=c_kub(abs(l),m+4)
             c3=c_kub(abs(l),m+8)

             frho  = c1*rho(i) + c2*rho(i+1) + c3*rho(i+2)
             frho1 = c1*rho1(i) + c2*rho1(i+1) + c3*rho1(i+2)
             frho2 = c1*rho2(i) + c2*rho2(i+1) + c3*rho2(i+2)

             call ylmderiv(yl,r,abs(l),m,frho,frho1,frho2,nder,ygrad,yhess)
             if (m /= 0) then
                call ylmderiv(yl,r,abs(l),-m,frho,frho1,frho2,nder,ygradtemp,yhesstemp)
                ygrad = imag0*(ygrad+minu*ygradtemp)
                yhess = imag0*(yhess+minu*yhesstemp)
             end if
             call ylmderiv(yl,r,abs(l),m+4,frho,frho1,frho2,nder,ygrad1,yhess1)
             if (m+4 /= 0) then
                call ylmderiv(yl,r,abs(l),-(m+4),frho,frho1,frho2,nder,ygradtemp,yhesstemp)
                ygrad1 = imag0*(ygrad1+minu*ygradtemp)
                yhess1 = imag0*(yhess1+minu*yhesstemp)
             end if
             call ylmderiv(yl,r,abs(l),m+8,frho,frho1,frho2,nder,ygrad2,yhess2)
             if (m+8 /= 0) then
                call ylmderiv(yl,r,abs(l),-(m+8),frho,frho1,frho2,nder,ygradtemp,yhesstemp)
                ygrad2 = imag0*(ygrad2+minu*ygradtemp)
                yhess2 = imag0*(yhess2+minu*yhesstemp)
             end if
             chg=chg + (c1*ang(i) + c2*ang(i+1) + c3*ang(i+2)) * (c1*rho(i) + c2*rho(i+1) + c3*rho(i+2))
             grad = grad + dble(c1*ygrad + c2*ygrad1 + c3*ygrad2)
             hess = hess + dble(c1*yhess + c2*yhess1 + c3*yhess2)
             i=i+3
          else
             write(uout,*) 'incorrect lm list for cubic structure'
             write(uout,'(a2,i4,a3,i4)') 'l=',l,' m=',m
             stop
          endif
       end do
    endif

    !.capture nuclei
    if (r.lt.f%rnot(jatom)+1.d-10) then
       grad = 0d0
       hess = 0d0
       hess(1)=-1d15
       hess(4)=-1d15
       hess(6)=-1d15
    endif

  end subroutine charge

  ! Evaluate the radial density and derivatives
  subroutine radial(f,rlm,rho,rho1,rho2,r,ir,jatom)
    ! note: rlm = clm, but name is changed to avoid collisions with wien.inc
    !.obtain by linear (oops) interpolation radial coefficients
    ! rho, and first and second derivatives rho1,rho2,at requested
    ! points from radial mesh
    ! rho=rlm/r^2
    ! rho1=d(rlm/r2)/dr
    ! rho2=d2(rlm/r2)/dr2
    class(wienwfn), intent(in) :: f
    real*8, dimension(nrad), intent(in) :: rlm
    real*8, intent(out) :: rho, rho1, rho2
    real*8, intent(in) :: r
    integer, intent(in) :: ir, jatom

    integer :: i, j
    real*8, dimension(4) :: c, cp, cpp
    real*8, dimension(4) :: r1, dr1
    real*8, dimension(4,4) :: x1dr12
    real*8 :: delta, delta2
    real*8, dimension(4) :: x1r1, x1r2, x1r3, x1r4
    integer :: ii, jj, ic
    real*8 :: prod
    integer :: temp_ir
    real*8 :: rrlm(4,0:2)

    ! careful with grid limits.
    if (ir <= 2) then
       temp_ir = 2
    else if (ir >= (f%jri(jatom) - 2)) then
       temp_ir = f%jri(jatom) - 2
    else
       temp_ir = ir
    end if

    delta=1.d0/f%dx(jatom)
    delta2=delta*delta

    rrlm(:,0) = rlm(temp_ir-1:temp_ir+2)
    x1dr12 = 0d0
    do i = 1, 4
       ii = temp_ir - 2 + i
       if (ii <= 2) then
          ic = 1
       else if ( ii >= (f%jri(jatom)-2)) then
          ic = 3
       else
          ic = 2
       end if

       rrlm(i,1:2) = 0d0
       do j = 1, 6
          jj = ii + noef(j,ic)
          rrlm(i,1) = rrlm(i,1) + coef1(j,ic) * rlm(jj)
          rrlm(i,2) = rrlm(i,2) + coef2(j,ic) * rlm(jj)
       end do
       rrlm(i,1) = rrlm(i,1) * fac1
       rrlm(i,2) = rrlm(i,2) * fac2

       ! calculate factors and distances
       r1(i) = f%rnot(jatom)*exp((ii-1)*f%dx(jatom))
       dr1(i) = r - r1(i)
       do j = 1, i-1
          x1dr12(i,j) = 1.d0 / (r1(i) - r1(j))
          x1dr12(j,i) = -x1dr12(i,j)
       end do
       x1r1(i) = 1.d0 / r1(i)
       x1r2(i) = x1r1(i) * x1r1(i)
       x1r3(i) = x1r1(i) * x1r2(i)
       x1r4(i) = x1r1(i) * x1r3(i)

       ! calculate radial derivatives
       c(i) = rrlm(i,0) * x1r2(i)
       cp(i) = (rrlm(i,1) * delta - 2.d0 * rrlm(i,0)) * x1r3(i)
       cpp(i) = (rrlm(i,2) * delta2 - 5.d0 * rrlm(i,1) * delta + 6.d0 * rrlm(i,0)) * x1r4(i)
    end do

    ! interpolate, lagrange 3rd order, 4 nodes
    rho = 0.d0
    rho1 = 0.d0
    rho2 = 0.d0
    do i = 1, 4
       prod = 1.d0
       do j = 1 ,4
          if (i == j) then
             cycle
          end if
          prod = prod * dr1(j) * x1dr12(i,j)
       end do
       rho = rho + c(i) * prod
       rho1 = rho1 + cp(i) * prod
       rho2 = rho2 + cpp(i) * prod
    end do

  end subroutine radial

  ! Calculate the density and derivatives in the interstitial, at
  ! point v (cryst. coords.). Calculate derivatives up to nder.
  ! Return density (chg), gradient (grad) and Hessian (hess) or
  ! zero if nder is not high enough.
  subroutine rhoout(f,v,nder,chg,grad,hess)
    use param, only: tpi, tpi2
    class(wienwfn), intent(in) :: f
    real*8, dimension(3), intent(in) :: v
    integer, intent(in) :: nder
    real*8, intent(out) :: chg
    real*8, dimension(3), intent(out) :: grad
    real*8, dimension(6), intent(out) :: hess

    real*8, dimension(3) :: factor
    real*8, dimension(3) :: vt
    integer :: kpp, niz, ntot
    integer :: i
    real*8 :: tpiarg, expo, dexpo, ddexpo
    real*8 :: ro, rhoexpo, rhodexpo, rhoddexpo
    real*8 :: aro, aroi
    complex*16 :: roc

    !.initialize charge, grad, hessian
    chg = 0d0
    grad = 0d0
    hess = 0d0

    !.to crystallographic if cell is orthogonal
    if (f%ortho) then
       factor = 1d0 / f%a
    else
       factor = 1d0
    endif
    vt = v * factor

    !.run over planewaves
    kpp=0
    niz=f%lastniz
    if (f%cmpl) then
       do i = f%lastind,1,-1
          ! kpp is number of equivalent plane wave
          ! niz is number of non-equivalent plane wave
          kpp=kpp+1
          if(kpp.gt.f%inst(niz)) then
             kpp=1
             niz=niz-1
          end if
          if ((abs(f%sk(niz))+abs(f%ski(niz))) < pwcutoff) then
             cycle
          endif

          tpiarg=tpi*dot_product(vt,f%krec(:,i))

          aro = f%sk(niz)
          aroi = f%ski(niz)
          roc = cmplx(aro,aroi,8)
          roc = roc * cmplx(f%tauk(i),f%tauki(i),8)
          roc = roc * cmplx(cos(tpiarg),sin(tpiarg),8)

          rhoexpo=real(roc,8)
          chg=chg+rhoexpo
          if (nder <= 0) cycle

          rhodexpo=-tpi*aimag(roc)
          grad = grad + rhodexpo * f%krec(:,i)
          if (nder <= 1) cycle

          rhoddexpo=-tpi2*rhoexpo
          hess(1)=hess(1)+rhoddexpo*f%krec(1,i)*f%krec(1,i)
          hess(2)=hess(2)+rhoddexpo*f%krec(1,i)*f%krec(2,i)
          hess(3)=hess(3)+rhoddexpo*f%krec(1,i)*f%krec(3,i)
          hess(4)=hess(4)+rhoddexpo*f%krec(2,i)*f%krec(2,i)
          hess(5)=hess(5)+rhoddexpo*f%krec(2,i)*f%krec(3,i)
          hess(6)=hess(6)+rhoddexpo*f%krec(3,i)*f%krec(3,i)
       enddo
    else
       ntot = 0
       do  i=f%lastind,1,-1
          ! kpp is number of equivalent plane wave
          ! niz is number of non-equivalent plane wave
          kpp=kpp+1
          if(kpp.gt.f%inst(niz)) then
             kpp=1
             niz=niz-1
          end if
          if (abs(f%sk(niz)) < pwcutoff) cycle
          ntot = ntot + 1

          tpiarg=tpi*dot_product(vt,f%krec(:,i))

          aro = f%sk(niz)
          ro = aro * f%tauk(i)
          expo=cos(tpiarg)
          rhoexpo=ro*expo
          chg = chg + rhoexpo
          if (nder <= 0) cycle

          dexpo=-sin(tpiarg)*tpi
          rhodexpo=ro*dexpo
          grad = grad + rhodexpo * f%krec(:,i)
          if (nder <= 1) cycle

          ddexpo=-expo*tpi2
          rhoddexpo=ro*ddexpo
          hess(1)=hess(1)+rhoddexpo*f%krec(1,i)*f%krec(1,i)
          hess(2)=hess(2)+rhoddexpo*f%krec(1,i)*f%krec(2,i)
          hess(3)=hess(3)+rhoddexpo*f%krec(1,i)*f%krec(3,i)
          hess(4)=hess(4)+rhoddexpo*f%krec(2,i)*f%krec(2,i)
          hess(5)=hess(5)+rhoddexpo*f%krec(2,i)*f%krec(3,i)
          hess(6)=hess(6)+rhoddexpo*f%krec(3,i)*f%krec(3,i)
       enddo
    end if

    !.change back to cartesian coordinates for ortho structures
    if (nder >= 1) grad = grad * factor
    if (nder >= 2) then
       hess(1)=hess(1)*factor(1)*factor(1)
       hess(2)=hess(2)*factor(1)*factor(2)
       hess(3)=hess(3)*factor(1)*factor(3)
       hess(4)=hess(4)*factor(2)*factor(2)
       hess(5)=hess(5)*factor(2)*factor(3)
       hess(6)=hess(6)*factor(3)*factor(3)
    end if

  end subroutine rhoout

  !> Spherical harmonics up to lmax
  SUBROUTINE YLM(V,LMAX,Y)
    use param, only: pi
! -----------------------------------------------------------------
! This subroutine is taken from SRC_LAPW2 and subsitutes the old
! version of ylm.f which was found to be numerically unstable for
! arguments close to but not on the z-axis (P. Blaha, priv. comm.)
! -----------------------------------------------------------------

      INTEGER            LMAX
      double precision   V(3)
      COMPLEX*16         Y(*)
!
!     ..................................................................
! 1.     PROGRAM UNIT 'YLM'
!           Calculates spherical harmonics
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           The spherical harmonics (Condon and Shortley convention)
!             Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!           for vector V (given in Cartesian coordinates)
!           are calculated. In the Condon Shortley convention the
!           spherical harmonics are defined as
!                             +------+
!                        m    |   1     m              im(Phi)
!           Y(l,m) = (-1)  -+ | -----  P (cos(Theta)) e
!                            \| 2(Pi)   l
!                  m
!           where P (cos(Theta)) is the normalized Associated Legendre
!                  l
!           function. Thus,
!                                          m      *
!                            Y(l,-m) = (-1) Y(l,m)
!
!
! 3.     USAGE
!           DOUBLE PRECISION V(3), Y(5*5)
!           V(1) = ...
!           V(2) = ...
!           V(3) = ...
!           CALL YLM(V,4,Y)
!
!        ARGUMENT-DESCRIPTION
!           V      - DOUBLE PRECISION vector, dimension 3        (input)
!                    Must be given in Cartesian coordinates.
!                    Conversion of V to polar coordinates gives the
!                    angles Theta and Phi necessary for the calculation
!                    of the spherical harmonics.
!           LMAX   - INTEGER value                               (input)
!                    upper bound of L for which spherical harmonics
!                    will be calculated
!                    constraint:
!                       LMAX .GE. 0 (not checked)
!           Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
!                    contains the calculated spherical harmonics
!                    Y(1)                   for L .EQ. 0 (M = 0)
!                    Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!                    ...
!                    Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!                                           for L .EQ. LMAX
!                                              (M = -L,...,L)
!                    constraint:
!                       Dimension of Y .GE. (LMAX+1)**2 (not checked)
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           Type COMPLEX*16 is used which does not conform to the
!           FORTRAN 77 standard.
!           Also the non-standard type conversion function DCMPLX()
!           is used which combines two double precision values into
!           one double complex value.
!
! 4.     REMARKS
!           none
!
! 5.     METHOD
!           The basic algorithm used to calculate the spherical
!           harmonics for vector V is as follows:
!
!           Y(0,0)
!           Y(1,0)
!           Y(1,1)
!           Y(1,-1) = -Y(1,1)
!           DO L = 2, LMAX
!              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
!              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!              DO M = L-2, 0, -1
!                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
!                 Y(L,-M)= (-1)**M*Y(L,M)
!              ENDDO
!           ENDDO
!
!           In the following the necessary recursion formulas and
!           starting values are given:
!
!        Start:
!                        +------+
!                        |   1
!           Y(0,0) =  -+ | -----
!                       \| 4(Pi)
!
!                                   +------+
!                                   |   3
!           Y(1,0) =  cos(Theta) -+ | -----
!                                  \| 4(Pi)
!
!                                     +------+
!                                     |   3    i(Phi)
!           Y(1,1) =  - sin(Theta) -+ | ----- e
!                                    \| 8(Pi)
!
!        Formula 1:
!
!           Y(l,l) =
!                           +--------+
!                           | (2l+1)   i(Phi)
!            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!                          \|   2l
!
!        Formula 2:
!                                  +---------------+
!                                  |  (2l-1)(2l+1)
!           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!                                 \|   (l-m)(l+m)
!
!                                    +--------------------+
!                                    |(l-1+m)(l-1-m)(2l+1)
!                              -  -+ |-------------------- Y(l-2,m)
!                                   \|  (2l-3)(l-m)(l+m)
!
!        Formula 3: (not used in the algorithm because of the division
!                    by sin(Theta) which may be zero)
!
!                                    +--------------+
!                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
!           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
!                      sin(Theta)   \| (l+m+1)(l-m)
!
!                                    +--------------+
!                                    |(l-m-1)(l+m+2)  -2i(Phi)
!                              -  -+ |-------------- e        Y(l,m+2)
!                                   \| (l-m)(l+m+1)
!
! 6.     DATE
!           26. April 1994                                   Version 1.2
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
!
      INTEGER            I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      DOUBLE PRECISION   A, B, C, AB, ABC, ABMAX, ABCMAX
      DOUBLE PRECISION   D4LL1C, D2L13
      DOUBLE PRECISION   COSTH, SINTH, COSPH, SINPH
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION   YLLR, YLL1R, YL1L1R, YLMR
      DOUBLE PRECISION   YLLI, YLL1I, YL1L1I, YLMI

!
!        Y(0,0)
!
      YLLR = 1.0D+0/SQRT(4.0D+0*PI)
      YLLI = 0.0D+0
      Y(1) = DCMPLX(YLLR,YLLI)
!
!        continue only if spherical harmonics for (L .GT. 0) are desired
!
      IF (LMAX .LE. 0) then
         GOTO 999
      end IF
!
!        calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
!        Theta, Phi ... polar angles of vector V
!
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. 0.0D+0) THEN
         A = V(1)/ABMAX
         B = V(2)/ABMAX
         AB = SQRT(A*A+B*B)
         COSPH = A/AB
         SINPH = B/AB
      ELSE
         COSPH = 1.0D+0
         SINPH = 0.0D+0
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. 0.0D+0) THEN
         A = V(1)/ABCMAX
         B = V(2)/ABCMAX
         C = V(3)/ABCMAX
         AB = A*A + B*B
         ABC = SQRT(AB + C*C)
         COSTH = C/ABC
         SINTH = SQRT(AB)/ABC
      ELSE
         COSTH = 1.0D+0
         SINTH = 0.0D+0
      ENDIF
!
!        Y(1,0)
!
      Y(3) = DCMPLX(SQRT(3.0D+0)*YLLR*COSTH,0.0D+0)
!
!        Y(1,1) ( = -DCONJG(Y(1,-1)))
!
      TEMP1 = -SQRT(1.5D+0)*YLLR*SINTH
      Y(4) = DCMPLX(TEMP1*COSPH,TEMP1*SINPH)
      Y(2) = -DCONJG(Y(4))
!
      DO 20 L = 2, LMAX
         INDEX  = L*L+1
         INDEX2 = INDEX + 2*L
         MSIGN  = 1 - 2*MOD(L,2)
!
!        YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         YL1L1R = DBLE(Y(INDEX-1))
         YL1L1I = DIMAG(Y(INDEX-1))
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTH
         YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
         YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
         Y(INDEX2) = DCMPLX(YLLR,YLLI)
         Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
!        YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!               (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP2 = SQRT(DBLE(2*L+1))*COSTH
         YLL1R = TEMP2*YL1L1R
         YLL1I = TEMP2*YL1L1I
         Y(INDEX2) = DCMPLX(YLL1R,YLL1I)
         Y(INDEX)  = -MSIGN*DCONJG(Y(INDEX2))
         INDEX2 = INDEX2 - 1
         INDEX  = INDEX  + 1
!
         I4L2 = INDEX2 - 4*L + 2
         I2L  = INDEX2 - 2*L
         D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
         D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!
         DO 10 M = L - 2, 0, -1
!
!        YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!
            TEMP1 = 1.0D+0/SQRT(DBLE((L+M)*(L-M)))
            TEMP2 = D4LL1C*TEMP1
            TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
            YLMR = TEMP2*DBLE(Y(I2L))  + TEMP3*DBLE(Y(I4L2))
            YLMI = TEMP2*DIMAG(Y(I2L)) + TEMP3*DIMAG(Y(I4L2))
            Y(INDEX2) = DCMPLX(YLMR,YLMI)
            Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
!
            MSIGN  = -MSIGN
            INDEX2 = INDEX2 - 1
            INDEX  = INDEX  + 1
            I4L2   = I4L2   - 1
            I2L    = I2L    - 1
   10    CONTINUE
   20 CONTINUE
!
  999 RETURN
!
!        End of 'YLM'
!
  END SUBROUTINE YLM

end submodule proc
