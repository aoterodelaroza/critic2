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

! meshmod: becke-style meshes for molecular integrals.
! The meshmod module has been adapted from postg: 
! Copyright (c) 2013 Alberto Otero de la Roza
! <aoterodelaroza@ucmerced.edu>, Felix Kannemann
! <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
! Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider
! <hs7@post.queensu.ca>, and Axel D. Becke <axel.becke@dal.ca>
module meshmod
  implicit none

  private
  public :: genmesh
  private :: genmesh_becke
  private :: genmesh_franchini
  public :: fillmesh
  private :: rmesh_postg
  private :: z2nr_postg
  private :: z2nr_franchini
  private :: z2nang_postg
  private :: z2nang_franchini
  private :: bhole
  private :: xfuncs

contains

  !> Driver for the generation of a molecular mesh. Uses the global
  !> MESH_type to decide the type and quality of the mesh.
  function genmesh(c) result(mesh)
    use crystalmod, only: crystal
    use global, only: mesh_type
    use types, only: molmesh
    type(crystal), intent(inout) :: c
    type(molmesh) :: mesh

    integer :: tmesh

    tmesh = MESH_type
    if (.not.c%ismolecule .and. MESH_type == 0) then
       tmesh = 3
    end if

    if (tmesh == 0) then
       mesh = genmesh_becke(c)
    else
       mesh = genmesh_franchini(c,tmesh)
    end if
    
  end function genmesh

  !> Generate a Becke-style molecular mesh. Only for molecules.
  function genmesh_becke(c) result(mesh)
    use crystalmod, only: crystal
    use tools_math, only: good_lebedev, select_lebedev
    use tools_io, only: ferror, faterr
    use types, only: molmesh
    use param, only: fourpi, maxzat
    type(crystal), intent(in) :: c
    type(molmesh) :: mesh

    real*8 :: rr(c%ncel,c%ncel), r, r1, r2, hypr, vp0, vpsum, vpi
    integer :: i, j, k, kk
    real*8, allocatable :: rads(:), wrads(:), xang(:), yang(:), zang(:), wang(:)
    integer :: nr, nang, ir, il, istat, mang, mr, iz, iz2
    real*8 :: cutoff(c%ncel,c%ncel), x(3)
    real*8, allocatable :: meshrl(:,:,:), meshx(:,:,:,:)

    if (.not.c%ismolecule) &
       call ferror("genmesh_becke","genmesh_becke can only work with molecules",faterr)
    
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
       if (iz < 1 .or. iz > maxzat) cycle
       mesh%n = mesh%n + z2nr_postg(iz) * z2nang_postg(iz)
    enddo
    allocate(mesh%w(mesh%n),mesh%x(3,mesh%n),stat=istat)

    ! allocate work arrays
    mr = -1
    mang = -1
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       mang = max(mang,z2nang_postg(iz))
       mr = max(mr,z2nr_postg(iz))
    end do
    allocate(meshrl(mang,mr,c%ncel),meshx(3,mang,mr,c%ncel))
    allocate(rads(mr),wrads(mr),stat=istat)
    if (istat /= 0) call ferror('genmesh_becke','could not allocate memory for radial meshes',faterr)
    allocate(xang(mang),yang(mang),zang(mang),wang(mang),stat=istat)
    if (istat /= 0) call ferror('genmesh_becke','could not allocate memory for angular meshes',faterr)
    
    ! Precompute the mesh weights with multiple threads. The job has to be
    ! split in two because the nodes have to be positioned in the array in 
    ! the correct order 
    !$omp parallel do private(nr,nang,ir,r,il,x,j,k,r1,r2,hypr,&
    !$omp cutoff,vp0,vpsum,vpi,iz,iz2) firstprivate(rads,wrads,xang,yang,zang,wang)
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1 .or. iz > maxzat) cycle

       ! radial mesh
       nr = z2nr_postg(iz)
       nang = z2nang_postg(iz)
       call rmesh_postg(nr,iz,rads,wrads)

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
                if (iz < 1 .or. iz > maxzat) cycle
                do k = 1, j-1
                   iz2 = c%at(c%atcel(k)%idx)%z 
                   if (iz2 < 1 .or. iz2 > maxzat) cycle
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
                if (iz < 1 .or. iz > maxzat) cycle
                vp0=vp0*cutoff(i,j)
                vpi=1d0
                do k = 1, c%ncel
                   iz2 = c%at(c%atcel(k)%idx)%z 
                   if (iz2 < 1 .or. iz2 > maxzat) cycle
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
       if (iz < 1 .or. iz > maxzat) cycle
       nr = z2nr_postg(iz)
       nang = z2nang_postg(iz)
       do ir = 1, nr
          do il = 1, nang
             kk = kk + 1
             mesh%w(kk) = meshrl(il,ir,i)
             mesh%x(:,kk) = meshx(:,il,ir,i)
          enddo
       enddo
    enddo

  end function genmesh_becke

  !> Generate a Becke-style molecular mesh, Franchini weights
  !> J. Comput. Chem. 34 (2013) 1819.
  !> The lvl parameter controls the quality:
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  !> This mesh is good for periodic systems because the calculation of the 
  !> weights does not involve a double sum over atoms.
  function genmesh_franchini(c,lvl) result(mesh)
    use crystalmod, only: crystal
    use tools_math, only: good_lebedev, select_lebedev
    use tools_io, only: faterr, ferror
    use types, only: molmesh
    use param, only: maxzat, fourpi
    type(crystal), intent(in) :: c
    integer, intent(in) :: lvl
    type(molmesh) :: mesh

    real*8 :: r, r1, vp0, vpsum
    integer :: i, j, kk
    real*8, allocatable :: rads(:), wrads(:), xang(:), yang(:), zang(:), wang(:)
    integer :: nr, nang, ir, il, istat, mang, mr, iz, iz2
    real*8 :: x(3), fscal, fscal2
    real*8, allocatable :: meshrl(:,:,:), meshx(:,:,:,:)

    ! allocate space for the mesh
    mesh%n = 0
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       mesh%n = mesh%n + z2nr_franchini(iz,lvl) * z2nang_franchini(iz,lvl)
    enddo
    allocate(mesh%w(mesh%n),mesh%x(3,mesh%n),stat=istat)

    ! allocate work arrays
    mr = -1
    mang = -1
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       mang = max(mang,z2nang_franchini(iz,lvl))
       mr = max(mr,z2nr_franchini(iz,lvl))
    end do
    allocate(meshrl(mang,mr,c%ncel),meshx(3,mang,mr,c%ncel))
    allocate(rads(mr),wrads(mr),stat=istat)
    if (istat /= 0) call ferror('genmesh_franchini','could not allocate memory for radial meshes',faterr)
    allocate(xang(mang),yang(mang),zang(mang),wang(mang),stat=istat)
    if (istat /= 0) call ferror('genmesh_franchini','could not allocate memory for angular meshes',faterr)
    
    ! Precompute the mesh weights with multiple threads. The job has to be
    ! split in two because the nodes have to be positioned in the array in 
    ! the correct order 
    !$omp parallel do private(iz,fscal,nr,nang,r,vp0,x,vpsum,iz2,fscal2,r1) &
    !$omp firstprivate(rads,wrads,xang,yang,zang,wang)
    do i = 1, c%ncel
       iz = c%at(c%atcel(i)%idx)%z 
       if (iz < 1 .or. iz > maxzat) then
          cycle
       elseif (iz == 1) then
          fscal = 0.3d0
       else
          fscal = 1d0
       end if

       ! radial mesh
       nr = z2nr_franchini(iz,lvl)
       nang = z2nang_franchini(iz,lvl)
       call rmesh_franchini(nr,iz,rads,wrads)

       ! angular mesh
       call good_lebedev(nang)
       call select_lebedev(nang,xang,yang,zang,wang)
       wang = wang / fourpi

       ! 3d mesh, do not parallelize to get the nodes in order
       do ir = 1, nr
          r = rads(ir)
          vp0 = fscal * exp(-2d0 * r) / max(r,1d-10)**3

          do il = 1, nang
             x = c%atcel(i)%r + r * (/xang(il),yang(il),zang(il)/)

             vpsum = 0d0
             do j = 1, c%nenv
                iz2 = c%at(c%atenv(j)%idx)%z 
                if (iz2 < 1 .or. iz2 > maxzat) then
                   cycle
                elseif (iz2 == 1) then
                   fscal2 = 0.3d0
                else
                   fscal2 = 1d0
                end if
                r1 = sqrt((x(1)-c%atenv(j)%r(1))**2+(x(2)-c%atenv(j)%r(2))**2+(x(3)-c%atenv(j)%r(3))**2)
                vpsum = vpsum + fscal2 * exp(-2d0 * r1) / max(r1,1d-10)**3
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
       if (iz < 1 .or. iz > maxzat) cycle
       nr = z2nr_franchini(iz,lvl)
       nang = z2nang_franchini(iz,lvl)
       do ir = 1, nr
          do il = 1, nang
             kk = kk + 1
             mesh%w(kk) = meshrl(il,ir,i)
             mesh%x(:,kk) = meshx(:,il,ir,i)
          enddo
       enddo
    enddo

  end function genmesh_franchini

  !> Calculate one or more scalar fields on the molecular mesh (m)
  !> using field f. id and prop are one-dimensional arrays of the same
  !> size. id contains the index in the mesh scalar field array and
  !> prop is the integer label for the property (param module).
  !> If periodic, assume the mesh is for a crystal (calculate the 
  !> properties by moving the points back to the main cell.
  subroutine fillmesh(m,ff,id,prop,periodic)
    use fields, only: grd
    use tools_io, only: faterr, ferror
    use types, only: molmesh, field, scalar_value_noalloc, realloc
    use param, only: im_rho, im_gradrho, im_gkin, im_b
    type(molmesh), intent(inout) :: m
    type(field), intent(inout) :: ff
    integer, intent(in) :: id(:)
    integer, intent(in) :: prop(:)
    logical, intent(in) :: periodic
    
    real*8, parameter :: bsmall = 1d-10

    type(scalar_value_noalloc) :: res
    integer :: i, j, n
    real*8 :: rhos, drho2, d2rho, taup, dsigs, quads
    real*8 :: fval

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

    !$omp parallel do private(fval,res,rhos,drho2,d2rho,taup,dsigs,quads)
    do i = 1, m%n
       call grd(ff,m%x(:,i),2,periodic,res0_noalloc=res)
       do j = 1, n
          select case(prop(j))
          case(im_rho)
             fval = res%f
          case(im_gradrho)
             fval = res%gfmod
          case(im_gkin)
             fval = res%gkin
          case(im_b)
             if (res%f > bsmall) then
                rhos = 0.5d0 * res%fval
                drho2 = 0.25d0 * res%gfmodval * res%gfmodval
                d2rho = 0.5d0 * res%del2fval
                taup = res%gkin
                dsigs = taup - 0.25d0 * drho2 / max(rhos,1d-30)
                quads = (d2rho - 2d0 * dsigs) / 6d0
                call bhole(rhos,quads,1d0,fval)
             endif
          case default
             fval = 0d0
          end select
          !$omp critical (save)
          m%f(i,j) = fval
          !$omp end critical (save)
       end do
    end do
    !$omp end parallel do

  end subroutine fillmesh

  !> Radial mesh for integration of exponential-like functions.
  !> From postg.
  subroutine rmesh_postg(n,iz,r,wintr)
    ! Radial mesh and integration weights,
    ! derivatives of variable r with respect to variable q.
    ! The q-mesh is uniform on the interval (0,+1). Transformation is
    ! r = rmid * q / (1-q)
    use param, only: third, fourpi
    integer, intent(in) :: n
    integer, intent(in) :: iz
    real*8, intent(out) :: r(n), wintr(n)

    real*8 :: h, q, rmid
    integer :: i

    rmid = 1d0/real(iz,8)**third
    h = 1d0/real(n+1,8)
    do i=1,n
       q=h*i
       r(i) = rmid * q / (1.d0-q)
       wintr(i) = fourpi * h * r(i)**2 * rmid / (1.d0-q)**2
    enddo

  end subroutine rmesh_postg

  !> Radial mesh for integration of exponential-like functions.
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> The transformation is:
  !>   r = zeta/ln(2) * (1+q) * ln(2/(1-q))
  !> where zeta is a function of the atomic number.
  subroutine rmesh_franchini(n,iz,r,wintr)
    use tools_math, only: gauleg
    use param, only: log2, fourpi
    integer, intent(in) :: n
    integer, intent(in) :: iz
    real*8, intent(out) :: r(n), wintr(n)

    !> From the SI in Franchini's.
    real*8, parameter :: zeta(103) = (/&
       0.8, 0.9, 1.8, 1.4, 1.3, 1.1, 0.9, 0.9, 0.9, 0.9,& !  1-10
       1.4, 1.3, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0, 1.5, 1.4,& ! 11-20
       1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1, 1.1, 1.1,& ! 21-30
       1.1, 1.0, 0.9, 0.9, 0.9, 0.9, 1.4, 1.4, 1.1, 1.3,& ! 31-40
       1.0, 1.2, 0.9, 0.9, 0.9, 1.0, 0.9, 1.0, 1.0, 1.3,& ! 41-50
       1.2, 1.2, 0.9, 1.0, 1.7, 1.5, 1.5, 1.3, 1.3, 1.4,& ! 51-60
       1.8, 1.4, 1.2, 1.3, 1.3, 1.4, 1.1, 1.1, 1.2, 1.6,& ! 61-70
       1.4, 1.3, 1.2, 1.0, 1.0, 0.9, 1.3, 1.2, 1.2, 1.0,& ! 71-80
       1.2, 1.2, 1.1, 1.2, 1.1, 2.1, 2.2, 1.8, 1.7, 1.3,& ! 81-90
       1.4, 1.2, 1.2, 1.3, 1.4, 1.4, 1.7, 1.9, 1.9, 2.0,& ! 91-100
       2.0, 1.6, 2.0/) ! 101-103

    real*8 :: h, q
    integer :: i

    call gauleg(-1d0,1d0,r,wintr,n)
    h = 1d0/real(n+1,8)
    do i = 1, n
       q = r(i)
       r(i) = zeta(iz)/log2 * (1+q) * log(2/(1-q))
       wintr(i) = fourpi * r(i)**2 * wintr(i) * zeta(iz)/log2 * (log(2/(1-q)) + (1+q)/(1-q))
    enddo

  end subroutine rmesh_franchini

  !> Atomic number to radial point number. From postg - really good
  !> accuracy but slow.
  function z2nr_postg(z) result(nr)
    integer, intent(in) :: z
    integer :: nr

    nr = 40
    if (z > 2) nr = 60
    if (z > 10) nr = 80
    if (z > 18) nr = 100
    if (z > 36) nr = 120
    if (z > 54) nr = 140
    if (z > 86) nr = 160

  endfunction z2nr_postg

  !> Atomic number to radial point number. 
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  !> Some changes to the choice of ni0
  function z2nr_franchini(z,lvl) result(nr)
    use tools_io, only: faterr, ferror
    integer, intent(in) :: z
    integer, intent(in) :: lvl
    integer :: nr

    nr = 15
    if (z > 2) nr = 20
    if (z > 10) nr = 25
    if (z > 18) nr = 35
    if (z > 36) nr = 60
    if (z > 54) nr = 85
    if (z > 86) nr = 110

    if (lvl == 1) then
       nr = ceiling(nr * 2.37)
    elseif (lvl == 2) then
       nr = ceiling(nr * 3.08)
    elseif (lvl == 3) then
       nr = ceiling(nr * 3.42)
    elseif (lvl == 4) then
       nr = ceiling(nr * 4.27)
    elseif (lvl == 5) then
       nr = ceiling(nr * 6.72)
    else
       call ferror("z2nr_franchini","unknown accuracy level",faterr)
    end if

  endfunction z2nr_franchini

  !> Select the number of points for the angular (Lebedev) quadrature.
  !> From postg.
  function z2nang_postg(z) result(nang)
    integer, intent(in) :: z
    integer :: nang

    nang = 194

  endfunction z2nang_postg

  !> Select the number of points for the angular (Lebedev) quadrature.
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  function z2nang_franchini(z,lvl) result(nang)
    use tools_io, only: ferror, faterr
    integer, intent(in) :: z
    integer, intent(in) :: lvl
    integer :: nang

    if (lvl == 1) then
       nang = 110
    elseif (lvl == 2) then
       nang = 194
    elseif (lvl == 3) then
       nang = 302
    elseif (lvl == 4) then
       nang = 590
    elseif (lvl == 5) then
       nang = 770
    else
       call ferror("z2nang_franchini","unknown accuracy level",faterr)
    end if
       
  endfunction z2nang_franchini

  subroutine bhole(rho,quad,hnorm,b)
    use param, only: twothird, pi, third

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
