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

!> Surface and mini-surface user-defined types and tools to work with them. 
module surface
  implicit none

  private

  public :: minisurf_init, minisurf_close, minisurf_clean
  public :: minisurf_spheresphere, minisurf_spheretriang, minisurf_spherecub
  public :: minisurf_transform
  public :: minisurf_write3dmodel
  public :: minisurf_writebasin
  public :: minisurf_writedbasin
  public :: minisurf_writeint, minisurf_readint

  integer, parameter :: surf_fpv = 8 !< Maximum faces per vertex
  integer, parameter :: surf_vpf = 4 !< Maximum vertex per face

contains

  !> Initialize a minisurface
  subroutine minisurf_init(s,m,f)
    use types, only: minisurf
    type(minisurf), intent(inout) :: s
    integer, intent(in) :: m, f
    
    allocate(s%r(m),s%th(m),s%ph(m))
    if (f /= 0) then
       allocate(s%f(f))
       s%init = 2
    else
       s%init = 1
    end if
    s%mv = m
    s%mf = f
    s%rgb = (/128,128,128/)

  end subroutine minisurf_init

  !> Destroy a minisurface
  subroutine minisurf_close(s)
    use types, only: minisurf
    type(minisurf), intent(inout) :: s
    
    deallocate(s%r,s%th,s%ph)
    if (s%init > 1) deallocate(s%f)
    
    s%init = 0
    
  end subroutine minisurf_close

  !> Clean the information in the minisurface.
  subroutine minisurf_clean(s)
    use types, only: minisurf
    type(minisurf), intent(inout) :: s

    s%n = (/ 0d0, 0d0, 0d0 /)
    s%r = -1d0
    s%th = 0d0
    s%ph = 0d0
    if (s%init > 1) then
       s%f(:)%v(1) = 0
       s%f(:)%v(2) = 0
       s%f(:)%v(3) = 0
       s%f(:)%v(4) = 0
    end if
    s%nv = 0
    s%nf = 0
    s%rgb = (/128,128,128/)

  end subroutine minisurf_clean

  !> Triangulate the unit sphere using an incremental pole->equator
  !> method to place the points. Analogous to surf_spheresphere,
  !> minisurface version.
  subroutine minisurf_spheresphere(s,xnuc,ntheta,nphi)
    use tools_io, only: ferror, faterr
    use param, only: pi, two, zero
    use types, only: minisurf
    type(minisurf), intent(inout) :: s
    real*8, dimension(3), intent(in) :: xnuc
    integer, intent(in) :: ntheta, nphi
    
    real*8 :: delta_theta, delta_phi, theta, phi
    integer :: nvertex, nphiact
    integer :: i, j, nn, npoints, first, second

    call minisurf_clean(s)
    
    ! Fill origin
    s%n = xnuc

    s%nv = 2*nphi*(2**ntheta-1)+2
    if (s%nv > s%mv) then
       call ferror('minisurf_spheresphere','max. number of vertex exceeded',faterr)
    end if

    s%nf = 6*nphi*(2**(ntheta-1)-1)+nphi+nphi
    if (s%nf > s%mf) then
       call ferror('minisurf_spheresphere','maximum number of faces exceeded',faterr)
    end if
    
    nvertex = 1
    delta_theta=pi/two/ntheta

    ! upper semisphere
    s%th(1) = 0d0
    s%ph(1) = 0d0
    s%r(1) = 1d0

    theta = delta_theta
    nphiact = nphi
    do i = 1, ntheta
       delta_phi = two*pi/nphiact
       phi = zero
       do j = 1, nphiact
          nvertex = nvertex + 1
          s%th(nvertex) = theta 
          s%ph(nvertex) = phi
          s%r(nvertex) = 1d0
          
          phi=phi+delta_phi
       end do
       theta = theta + delta_theta
       nphiact = nphiact + nphiact
    end do
    
    ! lower semisphere
    nn = s%nv / 2
    do i = 1,nn
       s%th(nn+i) = pi - s%th(i)
       s%ph(nn+i) = s%ph(i)
       s%r(nn+i) = 1d0
    end do

    ! Save face information
    s%nf = 0
    nn = s%nv / 2
    if (s%init > 1) then
       ! upper and lower cap
       do i = 1, nphi-1
          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = 1
          s%f(s%nf)%v(2) = i+1
          s%f(s%nf)%v(3) = i+2

          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = nn+1
          s%f(s%nf)%v(2) = i+nn+1
          s%f(s%nf)%v(3) = i+nn+2
       end do
       s%nf = s%nf + 1
       s%f(s%nf)%nv = 3
       s%f(s%nf)%v(1) = 1
       s%f(s%nf)%v(2) = nphi+1
       s%f(s%nf)%v(3) = 2

       s%nf = s%nf + 1
       s%f(s%nf)%nv = 3
       s%f(s%nf)%v(1) = nn+1
       s%f(s%nf)%v(2) = nphi+nn+1
       s%f(s%nf)%v(3) = nn+2

       ! body
       npoints = nphi
       first = 1
       second = first + npoints
       do i = 1, ntheta -1
          j = 1
          do while (j<npoints)
             s%nf = s%nf + 1
             s%f(s%nf)%nv = 3
             s%f(s%nf)%v(1) = second+1
             s%f(s%nf)%v(2) = first+1
             s%f(s%nf)%v(3) = second+2

             s%nf = s%nf + 1
             s%f(s%nf)%nv = 3
             s%f(s%nf)%v(1) = second+2
             s%f(s%nf)%v(2) = first+1
             s%f(s%nf)%v(3) = first+2

             s%nf = s%nf + 1
             s%f(s%nf)%nv = 3
             s%f(s%nf)%v(1) = second+2
             s%f(s%nf)%v(2) = first+2
             s%f(s%nf)%v(3) = second+3

             s%nf = s%nf + 1
             s%f(s%nf)%nv = 3
             s%f(s%nf)%v(1) = second+nn+1
             s%f(s%nf)%v(2) = first+nn+1
             s%f(s%nf)%v(3) = second+nn+2

             s%nf = s%nf + 1
             s%f(s%nf)%nv = 3
             s%f(s%nf)%v(1) = second+nn+2
             s%f(s%nf)%v(2) = first+nn+1
             s%f(s%nf)%v(3) = first+nn+2

             s%nf = s%nf + 1
             s%f(s%nf)%nv = 3
             s%f(s%nf)%v(1) = second+nn+2
             s%f(s%nf)%v(2) = first+nn+2
             s%f(s%nf)%v(3) = second+nn+3

             second = second + 2
             first = first + 1
             j = j + 1
          end do
          ! last three
          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = second+1
          s%f(s%nf)%v(2) = first+1
          s%f(s%nf)%v(3) = second+2

          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = second+2
          s%f(s%nf)%v(2) = first+1
          s%f(s%nf)%v(3) = first-npoints+2

          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = second+2
          s%f(s%nf)%v(2) = first-npoints+2
          s%f(s%nf)%v(3) = first+2

          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = second+nn+1
          s%f(s%nf)%v(2) = first+nn+1
          s%f(s%nf)%v(3) = second+nn+2

          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = second+nn+2
          s%f(s%nf)%v(2) = first+nn+1
          s%f(s%nf)%v(3) = first-npoints+nn+2

          s%nf = s%nf + 1
          s%f(s%nf)%nv = 3
          s%f(s%nf)%v(1) = second+nn+2
          s%f(s%nf)%v(2) = first-npoints+nn+2
          s%f(s%nf)%v(3) = first+nn+2
          npoints = npoints + npoints
          first = first + 1
          second = first + npoints
       end do
    end if

  end subroutine minisurf_spheresphere

  !> Triangulate the unit sphere by recursively subdividing an
  !> octahedron level times. Minisurface version.
  subroutine minisurf_spheretriang(s,xnuc,level)
    use tools_io, only: faterr, ferror, warning, uout
    use types, only: minisurf
    type(minisurf), intent(inout) :: s
    real*8, dimension(3), intent(in) :: xnuc
    integer, intent(in) :: level
    
    integer, parameter :: mlevel = 7

    ! nvert                       number of vertex
    ! nface                       number of faces in each workspace (0,1)
    ! iface(i,j,k)                vertex number i of face j in workspace k
    ! isw0, isw1                  switchers between workspaces
    ! i{xyz}(i)                   integer coordinate for vertex i
    ! i{123456}                   temp. for new vertex
    ! ilev                        level iteration counter
    ! i{xyz}{456}                 temp. for new vertex coordinates
    ! new                         flag for new vertex

    integer ::   nvert               
    integer ::   nface(0:1)          
    integer ::   iface(3,s%mf,0:1)
    integer ::   ix(s%mv), iy(s%mv), iz(s%mv) 
    integer ::   i1, i2, i3, i4, i5, i6, ilev
    integer ::   ix4, iy4, iz4, ix5, iy5, iz5, ix6, iy6, iz6
    integer ::   isw0, isw1          
    integer ::   nf1
    logical ::   new
    integer ::   i
    
    call minisurf_clean(s)

    ! Fill origin
    s%n = xnuc

    !.....First, let us compute a convex polyhedra that approximately
    !     triangulates the sphere. We will start with a octahedron and
    !     succesively divide each triangular facet into four new ones.
    !     The coordinates for the vertices will be maintained as integer
    !     values to get best precision and efficiency.
    nvert = 6
    if (nvert > s%mv) &
       call ferror('minisurf_spheretriang','max. number of vertex exceeded',faterr)
    ix(1) =  1024
    iy(1) =     0
    iz(1) =     0
    ix(2) =     0
    iy(2) =  1024
    iz(2) =     0
    ix(3) =     0
    iy(3) =     0
    iz(3) =  1024
    ix(4) = -1024
    iy(4) =     0
    iz(4) =     0
    ix(5) =     0
    iy(5) = -1024
    iz(5) =     0
    ix(6) =     0
    iy(6) =     0
    iz(6) = -1024
    !.....store initial faces of the octahedron:
    isw0 = 1
    isw1 = 0
    nface(0) = 8
    iface(1,1,0) = 1
    iface(2,1,0) = 2
    iface(3,1,0) = 3
    iface(1,2,0) = 1
    iface(2,2,0) = 3
    iface(3,2,0) = 5
    iface(1,3,0) = 5
    iface(2,3,0) = 3
    iface(3,3,0) = 4
    iface(1,4,0) = 2
    iface(2,4,0) = 4
    iface(3,4,0) = 3
    !
    iface(1,5,0) = 6
    iface(2,5,0) = 1
    iface(3,5,0) = 5
    iface(1,6,0) = 6
    iface(2,6,0) = 2
    iface(3,6,0) = 1
    iface(1,7,0) = 6
    iface(2,7,0) = 4
    iface(3,7,0) = 2
    iface(1,8,0) = 6
    iface(2,8,0) = 5
    iface(3,8,0) = 4

    if (level > mlevel) then
       call ferror('surf_spheretriang','subdivision level too large',warning)
       write (uout,'(A,I2)') "  -> The new level is : ", mlevel
    end if
    do ilev = 1, min(level,mlevel)

       !........exchange bins for old and new levels:
       !
       isw0 = mod(isw0+1,2)
       isw1 = mod(isw0+1,2)
       nface(isw1) = 0
       !
       !........iterate over the facets of the old level:
       !
       do i = 1, nface(isw0)
          i1 = iface(1,i,isw0)
          i2 = iface(2,i,isw0)
          i3 = iface(3,i,isw0)
          !
          !...........get new points to divide this triangle:
          !
          ix4 = (ix(i1)+ix(i2))/2
          iy4 = (iy(i1)+iy(i2))/2
          iz4 = (iz(i1)+iz(i2))/2
          ix5 = (ix(i2)+ix(i3))/2
          iy5 = (iy(i2)+iy(i3))/2
          iz5 = (iz(i2)+iz(i3))/2
          ix6 = (ix(i3)+ix(i1))/2
          iy6 = (iy(i3)+iy(i1))/2
          iz6 = (iz(i3)+iz(i1))/2
          !
          !...........store only new points:
          !
          new = .true.
          i4 = 0
          do while (new .and. i4.lt.nvert)
             i4 = i4 + 1
             if (ix(i4) .eq. ix4) then
                if (iy(i4) .eq. iy4) then
                   if (iz(i4) .eq. iz4) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spheretriang','max. number of vertex exceeded',faterr)
             i4 = nvert
             ix(i4) = ix4
             iy(i4) = iy4
             iz(i4) = iz4
          endif
          new = .true.
          i5 = 0
          do while (new .and. i5.lt.nvert)
             i5 = i5 + 1
             if (ix(i5) .eq. ix5) then
                if (iy(i5) .eq. iy5) then
                   if (iz(i5) .eq. iz5) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spheretriang','max. number of vertex exceeded',faterr)
             i5 = nvert
             ix(i5) = ix5
             iy(i5) = iy5
             iz(i5) = iz5
          endif
          new = .true.
          i6 = 0
          do while (new .and. i6.lt.nvert)
             i6 = i6 + 1
             if (ix(i6) .eq. ix6) then
                if (iy(i6) .eq. iy6) then
                   if (iz(i6) .eq. iz6) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spheretriang','max. number of vertex exceeded',faterr)
             i6 = nvert
             ix(i6) = ix6
             iy(i6) = iy6
             iz(i6) = iz6
          endif
          !
          !...........store the four new faces that result from the division of
          !           the old face:
          !
          nf1 = nface(isw1)
          iface(1,nf1+1,isw1) = i1
          iface(2,nf1+1,isw1) = i6
          iface(3,nf1+1,isw1) = i4
          iface(1,nf1+2,isw1) = i6
          iface(2,nf1+2,isw1) = i5
          iface(3,nf1+2,isw1) = i4
          iface(1,nf1+3,isw1) = i6
          iface(2,nf1+3,isw1) = i3
          iface(3,nf1+3,isw1) = i5
          iface(1,nf1+4,isw1) = i4
          iface(2,nf1+4,isw1) = i5
          iface(3,nf1+4,isw1) = i2
          nface(isw1) = nf1 + 4
       enddo
    enddo
    
    ! Save vertex information
    s%nv = nvert

    !.....Determine the spherical coordinates of the vertices, and store
    !     sin(theta)sin(phi), sin(theta)cos(phi), and cos(theta):
    do i = 1, s%nv
       s%r(i) = sqrt(real(ix(i)*ix(i)+iy(i)*iy(i)+iz(i)*iz(i),8))
       s%th(i) = acos(real(iz(i),8) /s%r(i)) 
       s%ph(i) = atan2(real(iy(i),8), real(ix(i),8))
       s%r(i) = 1d0
    enddo

    ! Save face information
    s%nf = nface(isw1)
    if (s%nf > s%mf) then
       call ferror('minisurf_spheretriang','maximum number of faces exceeded',faterr)
    end if
    if (s%init > 1) then
       do i=1,s%nf
          !! Vertices, id
          s%f(i)%nv   = 3
          s%f(i)%v(1) = iface(1,i,isw1)
          s%f(i)%v(2) = iface(2,i,isw1)
          s%f(i)%v(3) = iface(3,i,isw1)
       end do
    end if

  end subroutine minisurf_spheretriang

  !> Triangulate the unit sphere by recursively subdividing a cube
  !> level times. Minisurface version.
  subroutine minisurf_spherecub(s,xnuc,level)
    use tools_io, only: ferror, warning, uout, faterr
    use types, only: minisurf
    type(minisurf), intent(inout) :: s
    real*8, dimension(3), intent(in) :: xnuc
    integer, intent(in) :: level
    
    integer, parameter :: mlevel = 7

    ! nvert                       number of vertex
    ! nface                       number of faces in each workspace (0,1)
    ! iface(i,j,k)                vertex number i of face j in workspace k
    ! isw0, isw1                  switchers between workspaces
    ! i{xyz}(i)                   integer coordinate for vertex i
    ! i{123456}                   temp. for new vertex
    ! ilev                        level iteration counter
    ! i{xyz}{456}                 temp. for new vertex coordinates
    ! new                         flag for new vertex

    integer ::   nvert               
    integer ::   nface(0:1)          
    integer ::   iface(4,s%mf,0:1)
    integer ::   ix(s%mv), iy(s%mv), iz(s%mv) 
    integer ::   ix5, iy5, iz5, ix6, iy6, iz6
    integer ::   ix7, iy7, iz7, ix8, iy8, iz8, ix9, iy9, iz9
    integer ::   i, i1, i2, i3, i4, i5, i6, i7, i8, i9, ilev
    integer ::   isw0, isw1          
    integer ::   nf1
    logical ::   new
    
    call minisurf_clean(s)

    ! Fill origin
    s%n = xnuc

    !.First, let us compute a convex polyhedra that approximately
    ! cubicate the sphere. We will start with a cube and
    ! succesively divide each square facet into four new ones.
    ! The coordinates for the vertices will be maintained as integer
    ! values to get best precision and efficiency.
    !
    nvert = 8
    ix(1) =  1024
    iy(1) =  1024
    iz(1) =  1024
    !
    ix(2) = -1024
    iy(2) =  1024
    iz(2) =  1024
    !
    ix(3) = -1024
    iy(3) = -1024
    iz(3) =  1024
    !
    ix(4) =  1024
    iy(4) = -1024
    iz(4) =  1024
    !
    ix(5) =  1024
    iy(5) =  1024
    iz(5) = -1024
    !
    ix(6) = -1024
    iy(6) =  1024
    iz(6) = -1024
    !
    ix(7) = -1024
    iy(7) = -1024
    iz(7) = -1024
    !
    ix(8) =  1024
    iy(8) = -1024
    iz(8) = -1024
    !
    !.....store initial faces of the cube:
    !
    isw0 = 1
    isw1 = 0
    nface(0) = 6
    iface(1,1,0) = 1
    iface(2,1,0) = 5
    iface(3,1,0) = 6
    iface(4,1,0) = 2
    !
    iface(1,2,0) = 2
    iface(2,2,0) = 6
    iface(3,2,0) = 7
    iface(4,2,0) = 3
    !
    iface(1,3,0) = 3
    iface(2,3,0) = 7
    iface(3,3,0) = 8
    iface(4,3,0) = 4
    !
    iface(1,4,0) = 4
    iface(2,4,0) = 8
    iface(3,4,0) = 5
    iface(4,4,0) = 1
    !
    iface(1,5,0) = 3
    iface(2,5,0) = 4
    iface(3,5,0) = 1
    iface(4,5,0) = 2
    !
    iface(1,6,0) = 8
    iface(2,6,0) = 7
    iface(3,6,0) = 6
    iface(4,6,0) = 5
    !
    if (level > mlevel) then
       call ferror('surf_spherecub','subdivision level too large',warning)
       write (uout,'(A,I2)') "  -> The new level is : ", mlevel
    end if
    do ilev = 1, min(level,mlevel)
       !
       !........exchange bins for old and new levels:
       !
       isw0 = mod(isw0+1,2)
       isw1 = mod(isw0+1,2)
       nface(isw1) = 0
       !
       !........iterate over the facets of the old level:
       !
       do i = 1, nface(isw0)
          i1 = iface(1,i,isw0)
          i2 = iface(2,i,isw0)
          i3 = iface(3,i,isw0)
          i4 = iface(4,i,isw0)
          !
          !...........get new points to divide this square:
          !
          ix5 = (ix(i1)+ix(i2))/2
          iy5 = (iy(i1)+iy(i2))/2
          iz5 = (iz(i1)+iz(i2))/2
          !
          ix6 = (ix(i2)+ix(i3))/2
          iy6 = (iy(i2)+iy(i3))/2
          iz6 = (iz(i2)+iz(i3))/2
          !
          ix7 = (ix(i3)+ix(i4))/2
          iy7 = (iy(i3)+iy(i4))/2
          iz7 = (iz(i3)+iz(i4))/2
          !
          ix8 = (ix(i1)+ix(i4))/2
          iy8 = (iy(i1)+iy(i4))/2
          iz8 = (iz(i1)+iz(i4))/2
          !
          ix9 = (ix(i2)+ix(i4))/2
          iy9 = (iy(i2)+iy(i4))/2
          iz9 = (iz(i2)+iz(i4))/2
          !
          !...........store only new points:
          !
          new = .true.
          i5 = 0
          do while (new .and. i5.lt.nvert)
             i5 = i5 + 1
             if (ix(i5) .eq. ix5) then
                if (iy(i5) .eq. iy5) then
                   if (iz(i5) .eq. iz5) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spherecub','max. number of vertex exceeded',faterr)
             i5 = nvert
             ix(i5) = ix5
             iy(i5) = iy5
             iz(i5) = iz5
          endif
          new = .true.
          i6 = 0
          do while (new .and. i6.lt.nvert)
             i6 = i6 + 1
             if (ix(i6) .eq. ix6) then
                if (iy(i6) .eq. iy6) then
                   if (iz(i6) .eq. iz6) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spherecub','max. number of vertex exceeded',faterr)
             i6 = nvert
             ix(i6) = ix6
             iy(i6) = iy6
             iz(i6) = iz6
          endif
          new = .true.
          i7 = 0
          do while (new .and. i7.lt.nvert)
             i7 = i7 + 1
             if (ix(i7) .eq. ix7) then
                if (iy(i7) .eq. iy7) then
                   if (iz(i7) .eq. iz7) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spherecub','max. number of vertex exceeded',faterr)
             i7 = nvert
             ix(i7) = ix7
             iy(i7) = iy7
             iz(i7) = iz7
          endif
          new = .true.
          i8 = 0
          do while (new .and. i8.lt.nvert)
             i8 = i8 + 1
             if (ix(i8) .eq. ix8) then
                if (iy(i8) .eq. iy8) then
                   if (iz(i8) .eq. iz8) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spherecub','max. number of vertex exceeded',faterr)
             i8 = nvert
             ix(i8) = ix8
             iy(i8) = iy8
             iz(i8) = iz8
          endif
          new = .true.
          i9 = 0
          do while (new .and. i9.lt.nvert)
             i9 = i9 + 1
             if (ix(i9) .eq. ix9) then
                if (iy(i9) .eq. iy9) then
                   if (iz(i9) .eq. iz9) then
                      new = .false.
                   endif
                endif
             endif
          enddo
          if (new) then
             nvert = nvert + 1
             if (nvert > s%mv) &
                call ferror('surf_spherecub','max. number of vertex exceeded',faterr)
             i9 = nvert
             ix(i9) = ix9
             iy(i9) = iy9
             iz(i9) = iz9
          endif
          !
          !...........store the four new faces that result from the division of
          !           the old face:
          !
          nf1 = nface(isw1)
          iface(1,nf1+1,isw1) = i1
          iface(2,nf1+1,isw1) = i5
          iface(3,nf1+1,isw1) = i9
          iface(4,nf1+1,isw1) = i8
          !
          iface(1,nf1+2,isw1) = i8
          iface(2,nf1+2,isw1) = i9
          iface(3,nf1+2,isw1) = i7
          iface(4,nf1+2,isw1) = i4
          !
          iface(1,nf1+3,isw1) = i5
          iface(2,nf1+3,isw1) = i2
          iface(3,nf1+3,isw1) = i6
          iface(4,nf1+3,isw1) = i9
          !
          iface(1,nf1+4,isw1) = i9
          iface(2,nf1+4,isw1) = i6
          iface(3,nf1+4,isw1) = i3
          iface(4,nf1+4,isw1) = i7
          nface(isw1) = nf1 + 4
       enddo
    enddo
    
    ! Save vertex information
    s%nv = nvert

    !.....Determine the spherical coordinates of the vertices, and store
    !     sin(theta)sin(phi), sin(theta)cos(phi), and cos(theta):
    do i = 1, s%nv
       s%r(i) = sqrt(real(ix(i)*ix(i)+iy(i)*iy(i)+iz(i)*iz(i),8))
       s%th(i) = acos(real(iz(i),8) /s%r(i)) 
       s%ph(i) = atan2(real(iy(i),8), real(ix(i),8)) 
       s%r(i) = 1d0
    enddo

    ! Save face information
    s%nf = nface(isw1)
    if (s%nf > s%mf) then
       call ferror('minisurf_spherecub','maximum number of faces exceeded',faterr)
    end if
    if (s%init > 1) then
       do i=1,s%nf
          !! Vertices, id
          s%f(i)%nv   = 4
          s%f(i)%v(1) = iface(1,i,isw1)
          s%f(i)%v(2) = iface(2,i,isw1)
          s%f(i)%v(3) = iface(3,i,isw1)
          s%f(i)%v(4) = iface(4,i,isw1)
       end do
    end if

  end subroutine minisurf_spherecub

  !> Transform the points in the minisurface s using the symmetry operation op 
  !> and the translation vector tvec.
  subroutine minisurf_transform(s,op,tvec)
    use struct_basic, only: cr
    use types, only: minisurf
    type(minisurf) :: s
    integer, intent(in) :: op
    real*8, intent(in) :: tvec(3)

    integer :: i
    real*8 :: x(3), n(3), r

    n = cr%x2c(matmul(cr%rotm(:,1:3,op),cr%c2x(s%n)) + cr%rotm(:,4,op) + tvec)

    do i = 1, s%nv
       x = s%n + (/ s%r(i) * sin(s%th(i)) * cos(s%ph(i)),&
                    s%r(i) * sin(s%th(i)) * sin(s%ph(i)),&
                    s%r(i) * cos(s%th(i)) /)
       x = cr%c2x(x)
       x = cr%x2c(matmul(cr%rotm(:,1:3,op),x) + cr%rotm(:,4,op) + tvec)
       x = x - n
       r = sqrt(dot_product(x,x))
       s%th(i) = acos(x(3)/r) 
       s%ph(i) = atan2(x(2),x(1))
    end do
    s%n = n

  end subroutine minisurf_transform

  !> Write the minisurface s to the OFF file offfile.
  subroutine minisurf_write3dmodel(s,fmt,file,expr)
    use fields, only: f, grd, fields_fcheck, fields_feval
    use global, only: refden
    use graphics, only: graphics_open, graphics_surf, graphics_close
    use arithmetic, only: eval
    use tools_io, only: faterr, ferror
    use types, only: minisurf, scalar_value
    type(minisurf), intent(inout) :: s
    character*3, intent(in) :: fmt
    character*(*), intent(in) :: file
    character*(*), intent(in), optional :: expr
    
    type(scalar_value) :: res
    integer :: lu1, lu2
    integer :: i
    real*8 :: x(3)
    real*8, allocatable :: fsurf(:)
    logical :: iok

    if (s%init <= 1) &
       call ferror ('minisurf_write3dmodel','No face information in minisurf',faterr)

    ! calculate the field value on the basin
    if (present(expr)) then
       allocate(fsurf(s%nv))
       do i = 1, s%nv
          x(1) = s%n(1) + s%r(i) * sin(s%th(i)) * cos(s%ph(i))
          x(2) = s%n(2) + s%r(i) * sin(s%th(i)) * sin(s%ph(i))
          x(3) = s%n(3) + s%r(i) * cos(s%th(i))             
          call grd(f(refden),x,0,res)
          fsurf(i) = eval(expr,.true.,iok,x,fields_fcheck,fields_feval)
       end do
    else
       allocate(fsurf(3))
       fsurf = real(s%rgb,8) / 255d0
    end if

    call graphics_open(fmt,file,lu1,lu2)
    call graphics_surf(fmt,lu1,s,fsurf)
    call graphics_close(fmt,lu1,lu2)
    deallocate(fsurf)

  end subroutine minisurf_write3dmodel

  !> Write the surface s to the BASIN file offfile. Minisurface version.
  subroutine minisurf_writebasin (s,offfile,doprops)
    use global, only: refden
    use fields, only: integ_prop, f, nprops, grd, grdall
    use struct_basic, only: cr
    use tools_io, only: faterr, ferror, fopen_write, fclose
    use types, only: minisurf, scalar_value
    type(minisurf), intent(inout) :: s
    character*(*), intent(in) :: offfile
    logical, intent(in) :: doprops
    
    integer :: lud
    integer :: i, j
    real*8 :: lprop(Nprops), x(3), rr
    type(scalar_value) :: res

    if (s%init <= 1) then
       call ferror ('minisurf_writebasin','No face information in minisurf',faterr)
    end if

    ! open BASIN output file:
    lud = fopen_write(offfile)

    x = cr%c2x(s%n)
    write (lud,'("# POS(cryst) ",3(E22.14,X))') x
    write (lud,'("# CRYS2CART ")') 
    rr = 0d0
    do i = 1, 3
       write (lud,'("# ",3(E22.14,X),E10.2)') cr%crys2car(i,1:3), rr
    end do
    write (lud,'("# ",3(E22.14,X),E10.2)') 0d0, 0d0, 0d0, 0d0
    write (lud,'("# CART2CRYS ")') 
    rr = 0d0
    do i = 1, 3
       write (lud,'("# ",3(E22.14,X),E10.2)') cr%car2crys(i,1:3), rr
    end do
    write (lud,'("# ",3(E22.14,X),E10.2)') 0d0, 0d0, 0d0, 0d0

    write (lud,105) s%nv, s%nf, s%nv + s%nf - 2
    if (doprops) then
       write (lud,105) Nprops + 5
       write (lud,106) 'f','fval','|gradf|','lapf','lapfval',(integ_prop(i)%prop_name, i=1,Nprops)
    else
       write (lud,105) 1
       write (lud,106) 'f'
    end if

    ! Use the density and plot the properties
    do i = 1, s%nv
       x = s%n + (/ s%r(i) * sin(s%th(i)) * cos(s%ph(i)),&
                    s%r(i) * sin(s%th(i)) * sin(s%ph(i)),&
                    s%r(i) * cos(s%th(i)) /)
       call grd(f(refden),x,2,res)
       if (doprops) then
          call grdall(x,lprop)
          write (lud,110) x, res%f, res%fval, res%gfmod, res%del2f,&
             res%del2fval, lprop
       else
          write (lud,110) x, res%f
       end if
    end do
    
    do i = 1, s%nf
       write (lud,115) s%f(i)%nv,(s%f(i)%v(j)-1,j=1,s%f(i)%nv)
    end do

    call fclose(lud)

105 format (3i7)
106 format (999(a22,1x))
110 format (1p,999(e22.15,1x))
115 format (999(i7))

  end subroutine minisurf_writebasin

  !> Write the minisurface s and scalar field information on the
  !> basin rays to the DBASIN file offfile.
  subroutine minisurf_writedbasin (s,npoint,offfile)
    use global, only: refden
    use fields, only: f, grd
    use tools_io, only: ferror, faterr, fopen_write, fclose
    use types, only: minisurf, scalar_value
    type(minisurf), intent(inout) :: s
    integer, intent(in) :: npoint
    character*(*), intent(in) :: offfile
    
    integer :: lud
    integer :: i, j
    real*8  :: xxx(3)
    real*8 :: rdelta
    real*8, allocatable :: fpoint(:)
    type(scalar_value) :: res

    if (s%init <= 1) then
       call ferror ('minisurf_writedbasin','No face information in minisurf',faterr)
    end if

    allocate(fpoint(0:npoint))
    
    lud = fopen_write(offfile)

    write (lud,305) s%nv, s%nf, s%nv + s%nf - 2

    xxx = s%n
    call grd(f(refden),xxx,0,res)
    fpoint(0) = res%f
    write (lud,306) npoint, s%n(1), s%n(2), s%n(3), res%f

    do i = 1, s%nv
       rdelta = s%r(i) / npoint
       do j = 1, npoint
          xxx(1) = s%n(1) + rdelta * j * sin(s%th(i)) * cos(s%ph(i))
          xxx(2) = s%n(2) + rdelta * j * sin(s%th(i)) * sin(s%ph(i))
          xxx(3) = s%n(3) + rdelta * j * cos(s%th(i))             
          call grd(f(refden),xxx,0,res)
          fpoint(j) = res%f
       end do
       write (lud,310) xxx, (fpoint(j), j = 1, npoint)
    end do

    do i = 1, s%nf
       write (lud,315) s%f(i)%nv, (s%f(i)%v(j)-1,j=1,s%f(i)%nv)
    end do

    call fclose(lud)

    deallocate(fpoint)

305 format (3i7)
306 format (i7, 3f12.6, 1p, e14.6)
310 format (3f12.6, 1p, 10e14.6, (/ (10e14.6)))
315 format (999(i7))

  end subroutine minisurf_writedbasin

  !> Write the radii of the rays of the minisurface s to an INT file.
  !> These files are used to write and read the IAS information
  !> for a particular quadrature method (meth) with n1 / n2 nodes.
  subroutine minisurf_writeint (s,n1,n2,meth,offfile)
    use tools_io, only: fopen_write, fclose
    use types, only: minisurf
    type(minisurf), intent(inout) :: s
    integer, intent(in) :: n1, n2, meth
    character*(*), intent(in) :: offfile
    
    integer :: lud
    integer :: i

    lud = fopen_write(offfile)

    write (lud,'(I10,x,I10,x,I2)') n1, n2, meth
    write (lud,'(3(E23.15,x))') s%n(1), s%n(2), s%n(3)
    write (lud,'(3(E23.15,x))') (s%r(i),i=1,s%nv)

    call fclose(lud)

  end subroutine minisurf_writeint

  !> Read the IAS radii of the rays in the INT file into the
  !> minisurface s. The method fingerprint (meth, n1, n2) must match
  !> the one contained in the file.
  subroutine minisurf_readint (s,n1,n2,meth,offfile,ierr)
    use global, only: int_gauleg
    use tools_io, only: fopen_read, ferror, warning, string, fclose
    use types, only: minisurf
    type(minisurf), intent(inout) :: s
    integer, intent(in) :: n1, n2, meth
    character*(*), intent(in) :: offfile
    integer, intent(out) :: ierr
    
    integer :: lud
    integer :: n1_, n2_, method_
    real*8 :: xn(3)
    integer :: i

    ierr = 0
    lud = fopen_read(offfile)

    read (lud,'(I10,x,I10,x,I2)') n1_, n2_, method_
    if (n1 /= n1_ .or. meth /= method_ .or. (meth == INT_gauleg .and. n2 /= n2_)) then
       call ferror ('minisurf_readint', 'Wrong ntheta, nphi or method',warning,string(offfile)) 
       ierr = 1
       return
    end if
    read (lud,'(3(E23.15,x))') xn
    if (any(abs(xn - s%n) > 1d-13)) then
       call ferror ('minisurf_readint', 'Wrong CP position in integrals file',warning,string(offfile)) 
       ierr = 1
       return
    end if
    read (lud,'(3(E23.15,x))') (s%r(i),i=1,s%nv)

    call fclose(lud)

  end subroutine minisurf_readint

!  !> Write an OFF file from an xyz and face list.
!  subroutine xyz_writeoff (xyz,nface,face,offfile)
!    
!    real*8 :: xyz(3,:)
!    integer :: nface(:), face(:,:)
!    character*40, intent(in) :: offfile
!    
!    integer :: lud
!    integer :: nios, i, j, nv, nf
!    integer :: ipid
!
!    ! open OFF output file:
!    lud=fopen_write(offfile)
!
!    nv = size(xyz,2)
!    nf = size(nface,1)
!    write (lud,200) 
!    write (lud,205) nv, nf, nv+nf-2
!    
!    do i = 1, nv
!       write (lud,210) xyz(1:3,i)
!    end do
!
!    do i = 1, nf
!       write (lud,215) nface(i),(face(j,i)-1,j=1,nface(i))
!    end do
!
!    call fclose(lud)
!
!200 format ('OFF')
!205 format (3(i10,x))
!210 format (3(e22.14,x))
!215 format (999(i10,x))
!
!  end subroutine xyz_writeoff

end module surface

