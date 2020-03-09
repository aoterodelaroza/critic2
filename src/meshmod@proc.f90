! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (meshmod) proc
  implicit none

  !xx! private procedures
  ! subroutine rmesh_postg(n,iz,r,wintr)
  ! subroutine rmesh_franchini(n,iz,r,wintr)
  ! function z2nr(z,lvl) result(nr)
  ! function z2nang(z,lvl) result(nang)
  ! subroutine bhole(rho,quad,hnorm,b)
  ! subroutine xfuncs(x,rhs,f,df)

contains

  !> Deallocate and uninitalize
  module subroutine endmesh(m)
    class(mesh), intent(inout) :: m

    m%n = 0
    if (allocated(m%w)) deallocate(m%w)
    if (allocated(m%x)) deallocate(m%x)
    if (allocated(m%f)) deallocate(m%f)

  end subroutine endmesh

  !> Driver for the generation of a molecular mesh. Uses the global
  !> MESH_type to decide the type and quality of the mesh.
  module subroutine genmesh(m,c,type,lvl)
    use crystalmod, only: crystal
    class(mesh), intent(inout) :: m
    type(crystal), intent(inout) :: c
    integer, intent(in), optional :: type
    integer, intent(in), optional :: lvl

    integer :: tmesh, lmesh

    if (.not.c%ismolecule) then
       tmesh = mesh_type_franchini
    else if (present(type)) then
       tmesh = type
    else
       tmesh = mesh_type_becke
    end if

    if (present(lvl)) then
       lmesh = lvl
    else
       lmesh = mesh_level_good
    end if

    if (tmesh == mesh_type_becke) then
       call m%gen_becke(c,lmesh)
    else
       call m%gen_franchini(c,lmesh)
    end if
    m%type = tmesh
    m%lvl = lmesh
    
  end subroutine genmesh

  !> Generate a Becke-style molecular mesh. Only for molecules.
  module subroutine genmesh_becke(m,c,lvl)
    use crystalmod, only: crystal
    use tools_math, only: good_lebedev, select_lebedev
    use tools_io, only: ferror, faterr
    use param, only: fourpi, maxzat
    class(mesh), intent(inout) :: m
    type(crystal), intent(in) :: c
    integer, intent(in) :: lvl

    real*8 :: rr(c%ncel,c%ncel), r, r1, r2, hypr, vp0, vpsum, vpi
    integer :: i, j, k, kk
    real*8, allocatable :: rads(:), wrads(:), xang(:), yang(:), zang(:), wang(:)
    integer :: nr, nang, ir, il, istat, mang, mr, iz, iz2
    real*8 :: cutoff(c%ncel,c%ncel), x(3)
    real*8, allocatable :: meshrl(:,:,:), meshx(:,:,:,:)

    if (.not.c%ismolecule) &
       call ferror("genmesh_becke","Becke mesh only for molecules",faterr)
    
    ! reset the arrays
    call m%end()

    ! interatomic distances
    rr = 0d0
    do i = 1, c%ncel
       do j = i+1, c%ncel
          rr(i,j) = sqrt((c%atcel(i)%r(1)-c%atcel(j)%r(1))**2+(c%atcel(i)%r(2)-c%atcel(j)%r(2))**2&
             + (c%atcel(i)%r(3)-c%atcel(j)%r(3))**2)
          rr(j,i) = rr(i,j)
       enddo
    enddo

    ! allocate space for the mesh
    m%n = 0
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       m%n = m%n + z2nr(iz,lvl) * z2nang(iz,lvl)
    enddo
    allocate(m%w(m%n),m%x(3,m%n),stat=istat)

    ! allocate work arrays
    mr = -1
    mang = -1
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       mang = max(mang,z2nang(iz,lvl))
       mr = max(mr,z2nr(iz,lvl))
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
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) cycle

       ! radial mesh
       nr = z2nr(iz,lvl)
       nang = z2nang(iz,lvl)
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
                iz = c%spc(c%atcel(j)%is)%z 
                if (iz < 1 .or. iz > maxzat) cycle
                do k = 1, j-1
                   iz2 = c%spc(c%atcel(k)%is)%z 
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
                iz = c%spc(c%atcel(j)%is)%z 
                if (iz < 1 .or. iz > maxzat) cycle
                vp0=vp0*cutoff(i,j)
                vpi=1d0
                do k = 1, c%ncel
                   iz2 = c%spc(c%atcel(k)%is)%z 
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
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       nr = z2nr(iz,lvl)
       nang = z2nang(iz,lvl)
       do ir = 1, nr
          do il = 1, nang
             kk = kk + 1
             m%w(kk) = meshrl(il,ir,i)
             m%x(:,kk) = meshx(:,il,ir,i)
          enddo
       enddo
    enddo

  end subroutine genmesh_becke

  !> Generate a Becke-style molecular mesh, Franchini weights
  !> J. Comput. Chem. 34 (2013) 1819.
  !> The lvl parameter controls the quality:
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  !> This mesh is good for periodic systems because the calculation of the 
  !> weights does not involve a double sum over atoms.
  module subroutine genmesh_franchini(m,c,lvl)
    use crystalmod, only: crystal
    use environmod, only: environ
    use tools_math, only: good_lebedev, select_lebedev
    use tools_io, only: faterr, ferror
    use param, only: maxzat, fourpi, icrd_cart
    class(mesh), intent(inout) :: m
    type(crystal), intent(in) :: c
    integer, intent(in) :: lvl

    real*8 :: r, vp0, vpsum
    integer :: i, j, kk
    real*8, allocatable :: rads(:), wrads(:), xang(:), yang(:), zang(:), wang(:)
    integer, allocatable :: eid(:)
    integer :: nr, nang, ir, il, istat, mang, mr, iz, izmr, iz2, nat, lvec(3), ierr
    real*8 :: x(3), fscal, fscal2, xnuc(3), rmax
    real*8, allocatable :: meshrl(:,:,:), meshx(:,:,:,:), dist(:)
    logical :: isealloc
    type(environ), allocatable :: env

    real*8, parameter :: rthres = 12d0 ! contribution to weight: 2e-14

    ! reset the arrays
    call m%end()

    ! allocate space for the mesh
    m%n = 0
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       m%n = m%n + z2nr(iz,lvl) * z2nang(iz,lvl)
    enddo
    allocate(m%w(m%n),m%x(3,m%n),stat=istat)

    ! allocate work arrays
    rmax = 0d0
    izmr = -1
    mr = -1
    mang = -1
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       mang = max(mang,z2nang(iz,lvl))
       nr = z2nr(iz,lvl)
       if (nr > mr .or. nr == mr .and. iz > izmr) then
          mr = nr
          izmr = iz
       end if
    end do
    allocate(meshrl(mang,mr,c%ncel),meshx(3,mang,mr,c%ncel))
    allocate(rads(mr),wrads(mr),stat=istat)
    if (istat /= 0) call ferror('genmesh_franchini','could not allocate memory for radial meshes',faterr)
    allocate(xang(mang),yang(mang),zang(mang),wang(mang),stat=istat)
    if (istat /= 0) call ferror('genmesh_franchini','could not allocate memory for angular meshes',faterr)
    
    ! calculate the maximum r
    call rmesh_franchini(mr,izmr,rads,wrads)
    rmax = max(rads(mr),rthres)

    ! pointer to the environment
    if (rmax >= c%env%dmax0.and..not.c%ismolecule) then
       isealloc = .true.
       allocate(env)
       call env%extend(c%env,rmax)
    else
       ! keep a pointer to the environment
       isealloc = .false.
    end if
    
    ! Precompute the mesh weights with multiple threads. The job has to be
    ! split in two because the nodes have to be positioned in the array in 
    ! the correct order 
    !$omp parallel do private(iz,fscal,nr,nang,r,vp0,x,vpsum,iz2,fscal2,xnuc,nat,lvec,ierr) &
    !$omp firstprivate(rads,wrads,xang,yang,zang,wang,eid,dist)
    do i = 1, c%ncel
       xnuc = c%x2xr(c%atcel(i)%x)
       xnuc = xnuc - floor(xnuc)
       xnuc = c%xr2c(xnuc)
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) then
          cycle
       elseif (iz == 1) then
          fscal = 0.3d0
       else
          fscal = 1d0
       end if

       ! radial mesh
       nr = z2nr(iz,lvl)
       nang = z2nang(iz,lvl)
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
             x = xnuc + r * (/xang(il),yang(il),zang(il)/)

             ! find all atoms within a distance = rthres from the mesh point
             if (isealloc) then
                call env%list_near_atoms(x,icrd_cart,.false.,nat,eid,dist,lvec,ierr,up2d=rmax)
             else
                call c%env%list_near_atoms(x,icrd_cart,.false.,nat,eid,dist,lvec,ierr,up2d=rmax)
             end if
             if (ierr > 0) then
                call ferror('genmesh_franchini','could not find environment of a mesh point',faterr)
             end if

             vpsum = 0d0
             do j = 1, nat
                if (isealloc) then
                   iz2 = c%spc(env%at(eid(j))%is)%z 
                else
                   iz2 = c%spc(c%env%at(eid(j))%is)%z 
                end if
                if (iz2 < 1 .or. iz2 > maxzat) then
                   cycle
                elseif (iz2 == 1) then
                   fscal2 = 0.3d0
                else
                   fscal2 = 1d0
                end if
                vpsum = vpsum + fscal2 * exp(-2d0 * dist(j)) / max(dist(j),1d-10)**3
             enddo
             vpsum = max(vp0,vpsum)
                
             !$omp critical (mmesh)
             meshrl(il,ir,i) = vp0/max(vpsum,1d-40) * wrads(ir) * wang(il)
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
    if (isealloc) deallocate(env)

    ! fill the 3d mesh
    kk = 0
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z 
       if (iz < 1 .or. iz > maxzat) cycle
       nr = z2nr(iz,lvl)
       nang = z2nang(iz,lvl)
       do ir = 1, nr
          do il = 1, nang
             kk = kk + 1
             m%w(kk) = meshrl(il,ir,i)
             m%x(:,kk) = meshx(:,il,ir,i)
          enddo
       enddo
    enddo

  end subroutine genmesh_franchini

  !> Calculate one or more scalar fields on the molecular mesh (m)
  !> using field f. prop is a one-dimensional array containing an
  !> integer label for the property (param module).  If prop is 100+n,
  !> then calculate molecular orbital number n (only for if ff is a
  !> molecular wavefunction).  If periodic, assume the mesh is for a
  !> crystal (calculate the properties by moving the points back to
  !> the main cell.
  module subroutine fillmesh(m,ff,prop,periodic)
    use fieldmod, only: field
    use tools_io, only: faterr, ferror
    use tools_math, only: bhole
    use types, only: scalar_value, realloc
    use param, only: im_volume, im_rho, im_gradrho, im_gkin, im_b, im_null
    class(mesh), intent(inout) :: m
    type(field), intent(inout) :: ff
    integer, intent(in) :: prop(:)
    logical, intent(in) :: periodic
    
    type(scalar_value) :: res
    integer :: i, j, n, nder
    real*8 :: fval, rhos, laps, tau, drhos2, dsigs, quads, br_b, br_alf, br_a
    character*10 :: fder

    if (.not.ff%isinit) &
       call ferror("fillmesh","field not initialized",faterr)

    n = size(prop)
    if (allocated(m%f)) then
       if (size(m%f,2) /= n) call realloc(m%f,m%n,n)
    else
       allocate(m%f(m%n,n))
    end if

    ! calcualte the maximum derivative
    nder = -1
    do j = 1, n
       select case(prop(j))
       case(im_volume,im_null)
          nder = -1
       case(im_rho)
          nder = max(nder,0)
       case(im_gradrho,im_gkin)
          nder = max(nder,1)
       case(im_b)
          nder = max(nder,2)
       end select
    end do

    !$omp parallel do private(fval,res,fder,rhos,laps,tau,drhos2,dsigs,quads,br_b,br_alf,br_a)
    do i = 1, m%n
       if (nder >= 0) then
          call ff%grd(m%x(:,i),nder,res,periodic=periodic)
       end if
       do j = 1, n
          if (prop(j) == im_volume) then
             fval = 1d0
          else if (prop(j) == im_rho) then
             fval = res%f
          else if (prop(j) == im_gradrho) then
             fval = res%gfmod
          else if (prop(j) == im_gkin) then
             fval = res%gkin
          else if (prop(j) == im_b) then
             rhos = 0.5d0 * res%f
             laps = 0.5d0 * res%del2f
             tau = res%gkin
             drhos2 = (0.5d0 * res%gfmod)
             drhos2 = drhos2 * drhos2
             dsigs = tau - 0.25d0 * drhos2 / max(rhos,1d-30)
             quads = (laps - 2d0 * dsigs) / 6d0
             call bhole(rhos,quads,1d0,br_b,br_alf,br_a)
             fval = br_b
          else if (prop(j) > 100) then
             write (fder,'(I10)') prop(j) - 100
             call ff%grd(m%x(:,i),-1,res,fder=fder,periodic=periodic)
             fval = res%fspc
          else if (prop(j) == im_null) then
             fval = 0d0
          else
             call ferror("fillmesh","unknown quantity to calculate in fillmesh",faterr)
          end if
          !$omp critical (save)
          m%f(i,j) = fval
          !$omp end critical (save)
       end do
    end do
    !$omp end parallel do

  end subroutine fillmesh

  !> Write information about the mesh to the standard output.
  module subroutine report(m)
    use tools_io, only: uout, string
    class(mesh), intent(inout) :: m
    
    write (uout,'(2X,"Mesh size      ",A)') string(m%n)
    if (m%type == mesh_type_becke) then
       write (uout,'(2X,"Mesh type      Becke")')
    elseif (m%type == mesh_type_franchini) then
       write (uout,'(2X,"Mesh type      Franchini")')
    end if
    if (m%lvl == mesh_level_small) then
       write (uout,'(2X,"Mesh level     small")')
    elseif (m%lvl == mesh_level_normal) then
       write (uout,'(2X,"Mesh level     normal")')
    elseif (m%lvl == mesh_level_good) then
       write (uout,'(2X,"Mesh level     good")')
    elseif (m%lvl == mesh_level_vgood) then
       write (uout,'(2X,"Mesh level     very good")')
    elseif (m%lvl == mesh_level_amazing) then
       write (uout,'(2X,"Mesh level     amazing")')
    end if

  end subroutine report

  !xx! private procedures

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

  !> Atomic number to radial point number. 
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  !> Some changes to the choice of ni0
  function z2nr(z,lvl) result(nr)
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

    if (lvl == mesh_level_small) then
       nr = ceiling(nr * 2.37)
    elseif (lvl == mesh_level_normal) then
       nr = ceiling(nr * 3.08)
    elseif (lvl == mesh_level_good) then
       nr = ceiling(nr * 3.42)
    elseif (lvl == mesh_level_vgood) then
       nr = ceiling(nr * 4.27)
    elseif (lvl == mesh_level_amazing) then
       nr = ceiling(nr * 6.72)
    else
       call ferror("z2nr","unknown accuracy level",faterr)
    end if

  endfunction z2nr

  !> Select the number of points for the angular (Lebedev) quadrature.
  !> After Franchini et al., J. Comput. Chem. 34 (2013) 1819.
  !> lvl = 1 (small), 2 (normal), 3 (good), 4(very good), 5 (excellent)
  function z2nang(z,lvl) result(nang)
    use tools_io, only: ferror, faterr
    integer, intent(in) :: z
    integer, intent(in) :: lvl
    integer :: nang

    if (lvl == mesh_level_small) then
       nang = 110
    elseif (lvl == mesh_level_normal) then
       nang = 194
    elseif (lvl == mesh_level_good) then
       nang = 302
    elseif (lvl == mesh_level_vgood) then
       nang = 590
    elseif (lvl == mesh_level_amazing) then
       nang = 770
    else
       call ferror("z2nang","unknown accuracy level",faterr)
    end if
       
  endfunction z2nang

end submodule proc
