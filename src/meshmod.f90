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

! meshmod: becke-style meshes for molecular integrals.
! The meshmod module has been adapted from postg: 
! Copyright (c) 2013 Alberto Otero de la Roza
! <aoterodelaroza@ucmerced.edu>, Felix Kannemann
! <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider
! <hs7@post.queensu.ca>, and Axel D. Becke <axel.becke@dal.ca>
module meshmod
  use param
  implicit none

  private
  public :: genmesh
  public :: fillmesh
  private :: rmesh
  private :: z2nr
  private :: z2nang
  private :: bhole
  private :: xfuncs

contains

  !> Generate a Becke-style molecular 
  function genmesh(c) result(mesh)
    use struct_basic
    use tools_math
    use tools_io
    use types

    type(crystal), intent(in) :: c
    type(molmesh) :: mesh

    real*8 :: rr(c%ncel,c%ncel), rmid, r, r1, r2, hypr, vp0, vpsum, vpi
    integer :: i, j, k, kk
    real*8, allocatable :: rads(:), wrads(:), xang(:), yang(:), zang(:), wang(:)
    integer :: nr, nang, ir, il, istat, mang, mr, iz, iz2
    real*8 :: cutoff(c%ncel,c%ncel), x(3)
    real*8, allocatable :: meshrl(:,:,:), meshx(:,:,:,:)

    if (.not.c%ismolecule) &
       call ferror("genmesh","genmesh can only work with molecules",faterr)
    
    ! interatomic distances
    rr = 0d0
    do i = 1, c%ncel
       do j = i+1, c%ncel
          rr(i,j) = sqrt((c%atcel(i)%r(1)-c%atcel(j)%r(1))**2+(c%atcel(i)%r(2)-c%atcel(j)%r(2))**2+(c%atcel(i)%r(3)-c%atcel(j)%r(3))**2)
          rr(j,i) = rr(i,j)
       enddo
    enddo

    ! allocate space for the mesh
    mesh%n = 0
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1) cycle
       mesh%n = mesh%n + z2nr(iz) * z2nang(iz)
    enddo
    allocate(mesh%w(mesh%n),mesh%x(3,mesh%n),stat=istat)

    ! allocate work arrays
    mr = -1
    mang = -1
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1) cycle
       mang = max(mang,z2nang(iz))
       mr = max(mr,z2nr(iz))
    end do
    allocate(meshrl(mang,mr,c%ncel),meshx(3,mang,mr,c%ncel))
    allocate(rads(mr),wrads(mr),stat=istat)
    if (istat /= 0) call ferror('genmesh','could not allocate memory for radial meshes',faterr)
    allocate(xang(mang),yang(mang),zang(mang),wang(mang),stat=istat)
    if (istat /= 0) call ferror('readwfn','could not allocate memory for angular meshes',faterr)
    
    ! Precompute the mesh weights with multiple threads. The job has to be
    ! split in two because the nodes have to be positioned in the array in 
    ! the correct order 
    !$omp parallel do private(nr,nang,rmid,ir,r,il,x,j,k,r1,r2,hypr,&
    !$omp cutoff,vp0,vpsum,vpi) firstprivate(rads,wrads,xang,yang,zang,wang)
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1) cycle

       ! radial mesh
       nr = z2nr(iz)
       nang = z2nang(iz)
       rmid = 1d0/real(iz,8)**third
       call rmesh(nr,rmid,rads,wrads)

       ! angular mesh
       call good_lebedev(nang)
       call select_lebedev(nang,xang,yang,zang,wang)
       wang = wang / fourpi

       ! 3d mesh, do not parallelize to get the nodes in order
       do ir = 1, nr
          r = rads(ir)
          do il = 1, nang
             x = c%atcel(i)%r + r * (/xang(il),yang(il),zang(il)/)
             do j = 2, c%ncel
                iz = c%at(c%atcel(j)%idx)%z 
                if (iz < 1) cycle
                do k = 1, j-1
                   iz2 = c%at(c%atcel(k)%idx)%z 
                   if (iz2 < 1) cycle
                   r1 = sqrt((x(1)-c%atcel(j)%r(1))**2+(x(2)-c%atcel(j)%r(2))**2+(x(3)-c%atcel(j)%r(3))**2)
                   r2 = sqrt((x(1)-c%atcel(k)%r(1))**2+(x(2)-c%atcel(k)%r(2))**2+(x(3)-c%atcel(k)%r(3))**2)
                   hypr = (r1-r2) / rr(j,k)
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   hypr = 1.5d0*hypr-0.5d0*hypr**3
                   cutoff(j,k) = (1d0-hypr) / 2d0
                   cutoff(k,j) = (1d0+hypr) / 2d0
                enddo
                cutoff(j,j) = 1d0
             enddo
             cutoff(1,1) = 1d0
             vp0 = 1d0
             vpsum = 0d0
             do j = 1, c%ncel
                iz = c%at(c%atcel(j)%idx)%z 
                if (iz < 1) cycle
                vp0=vp0*cutoff(i,j)
                vpi=1d0
                do k = 1, c%ncel
                   iz2 = c%at(c%atcel(k)%idx)%z 
                   if (iz2 < 1) cycle
                   vpi = vpi * cutoff(j,k)
                enddo
                vpsum = vpsum + vpi
             enddo
             !$omp critical (mmesh)
             meshrl(il,ir,i) = vp0/vpsum * wrads(ir) * wang(il)
             meshx(:,il,ir,i) = x
             !$omp end critical (mmesh)
          enddo
       enddo
    end do
    !$omp end parallel do

    ! clean up
    if (allocated(rads)) deallocate(rads)
    if (allocated(wrads)) deallocate(wrads)
    if (allocated(xang)) deallocate(xang)
    if (allocated(yang)) deallocate(yang)
    if (allocated(zang)) deallocate(zang)
    if (allocated(wang)) deallocate(wang)

    ! fill the 3d mesh
    kk = 0
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1) cycle
       nr = z2nr(iz)
       nang = z2nang(iz)
       do ir = 1, nr
          do il = 1, nang
             kk = kk + 1
             mesh%w(kk) = meshrl(il,ir,i)
             mesh%x(:,kk) = meshx(:,il,ir,i)
          enddo
       enddo
    enddo

  end function genmesh

  !> Calculate one or more scalar fields on the molecular mesh (m)
  !> using field f. id and prop are one-dimensional arrays of the same
  !> size. id contains the index in the mesh scalar field array and
  !> prop is the integer label for the property (param module).
  subroutine fillmesh(m,ff,id,prop)
    use fields
    use tools_io
    use types

    type(molmesh), intent(inout) :: m
    type(field), intent(inout) :: ff
    integer, intent(in) :: id(:)
    integer, intent(in) :: prop(:)
    
    real*8, parameter :: bsmall = 1d-10

    type(scalar_value) :: res
    integer :: i, j, n
    real*8 :: x(3), rhos, drho2, d2rho, taup, dsigs, quads

    if (size(id) /= size(prop)) &
       call ferror("fillmesh","incongruent id and prop arrays",faterr)
    if (.not.ff%init) &
       call ferror("fillmesh","field not initialized",faterr)

    n = size(id)
    if (allocated(m%f)) then
       if (size(m%f,2) /= n) call realloc(m%f,m%n,n)
    else
       allocate(m%f(m%n,n))
    end if

    do i = 1, m%n
       call grd(ff,m%x(:,i),2,res,.false.)
       do j = 1, n
          select case(prop(j))
          case(im_rho)
             m%f(i,j) = res%f
          case(im_gradrho)
             m%f(i,j) = res%gfmod
          case(im_gkin)
             m%f(i,j) = res%gkin
          case(im_b)
             if (res%f > bsmall) then
                rhos = 0.5d0 * res%f
                drho2 = 0.5d0 * res%gfmod * res%gfmod
                d2rho = 0.5d0 * res%del2f
                taup = 0.5d0 * res%gkin
                dsigs = taup - 0.25d0 * drho2 / max(rhos,1d-30)
                quads = (d2rho - 2d0 * dsigs) / 6d0
                call bhole(rhos,quads,1d0,m%f(i,j))
             endif
          case default
             m%f(i,j) = 0d0
          end select
       end do
    end do

  end subroutine fillmesh

  !> Radial mesh for integration of exponential-like functions.
  subroutine rmesh(n,rmid,r,wintr)
    ! Radial mesh and integration weights,
    ! derivatives of variable r with respect to variable q.
    !
    ! The q-mesh is uniform on the interval (0,+1). Transformation is
    !
    !                    Q
    !     R  =  RMID ---------
    !                ( 1 - Q )
    !
    ! Also, 3 and 7-point finite difference matrices for d2(wrt)r on q-mesh.

    integer, intent(in) :: n
    real*8, intent(in) :: rmid
    real*8, intent(out) :: r(n), wintr(n)

    real*8 :: h, q
    integer :: i

    h = 1d0/real(n+1,8)
    do i=1,n
       q=h*i
       r(i) = rmid * q / (1.d0-q)
       wintr(i) = fourpi * h * r(i)**2 * rmid / (1.d0-q)**2
    enddo

  end subroutine rmesh

  function z2nr(z)
    integer, intent(in) :: z
    integer :: z2nr

    z2nr = 40
    if (z > 2) z2nr = 60
    if (z > 10) z2nr = 80
    if (z > 18) z2nr = 100
    if (z > 36) z2nr = 120
    if (z > 54) z2nr = 140
    if (z > 86) z2nr = 160

  endfunction z2nr

  function z2nang(z)
    integer, intent(in) :: z
    integer :: z2nang

    z2nang = 194

  endfunction z2nang

  subroutine bhole(rho,quad,hnorm,b)
    use param

    real*8, intent(in) :: rho, quad, hnorm
    real*8, intent(out) :: b 

    real*8 :: rhs, x0, shift, x, x1, expo, prefac, alf, f, df
    integer :: i

    rhs=twothird*(pi*rho/hnorm)**twothird*rho/quad
    x0=2.d0
    shift=1.d0
    if(rhs.lt.0.d0)go to 10
    if(rhs.gt.0.d0)go to 20
10  do i=1,16
       x=x0-shift
       call xfuncs(x,rhs,f,df)
       if(f.lt.0.d0)go to 88
       shift=0.1d0*shift
    enddo
    write(*,1002)
    stop
20  do i=1,16
       x=x0+shift
       call xfuncs(x,rhs,f,df)
       if(f.gt.0.d0)go to 88
       shift=0.1d0*shift
    enddo
    write(*,1002)
    stop
88  continue
    do i=1,100
       call xfuncs(x,rhs,f,df)
       x1=x-f/df
       if(dabs(x1-x).lt.1.d-10)go to 111
       x=x1
    enddo
    write(*,1001)
    stop
111 x=x1
    expo=dexp(-x)
    prefac=rho/expo
    alf=(8.d0*pi*prefac/hnorm)**third
    b=x/alf
    return
1001 format(' ','bhole: newton algorithm fails to converge!')
1002 format(' ','bhole: newton algorithm fails to initialize!')
  end subroutine bhole
  
  subroutine xfuncs(x,rhs,f,df)
    real*8, intent(in) :: x, rhs
    real*8, intent(out) :: f, df

    real*8 :: expo23

    expo23=dexp(-2.d0/3.d0*x)
    f = x*expo23/(x-2.d0) - rhs
    df=2.d0/3.d0*(2.d0*x-x**2-3.d0)/(x-2.d0)**2*expo23
  end subroutine xfuncs

end module meshmod
