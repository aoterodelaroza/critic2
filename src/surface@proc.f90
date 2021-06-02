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

submodule (surface) proc
  implicit none

  integer, parameter :: surf_fpv = 8 !< Maximum faces per vertex
  integer, parameter :: surf_vpf = 4 !< Maximum vertex per face

contains

  !> Initialize a minisurface
  module subroutine minisurf_begin(s,m,f)
    class(minisurf), intent(inout) :: s
    integer, intent(in) :: m, f

    if (allocated(s%r)) deallocate(s%r)
    if (allocated(s%th)) deallocate(s%th)
    if (allocated(s%ph)) deallocate(s%ph)
    if (allocated(s%f)) deallocate(s%f)
    allocate(s%r(m),s%th(m),s%ph(m))
    if (f /= 0) then
       allocate(s%f(f))
       s%isinit = 2
    else
       s%isinit = 1
    end if
    s%mv = m
    s%mf = f
    s%rgb = (/128,128,128/)

  end subroutine minisurf_begin

  !> Destroy a minisurface
  module subroutine minisurf_close(s)
    class(minisurf), intent(inout) :: s

    if (allocated(s%r)) deallocate(s%r)
    if (allocated(s%th)) deallocate(s%th)
    if (allocated(s%ph)) deallocate(s%ph)
    if (allocated(s%f)) deallocate(s%f)
    s%isinit = 0

  end subroutine minisurf_close

  !> Clean the information in the minisurface.
  module subroutine minisurf_clean(s)
    class(minisurf), intent(inout) :: s

    s%n = (/ 0d0, 0d0, 0d0 /)
    s%r = -1d0
    s%th = 0d0
    s%ph = 0d0
    if (s%isinit > 1) then
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
  module subroutine spheresphere(s,xnuc,ntheta,nphi)
    use tools_io, only: ferror, faterr
    use param, only: pi, two, zero
    class(minisurf), intent(inout) :: s
    real*8, dimension(3), intent(in) :: xnuc
    integer, intent(in) :: ntheta, nphi

    real*8 :: delta_theta, delta_phi, theta, phi
    integer :: nvertex, nphiact
    integer :: i, j, nn, npoints, first, second

    call s%clean()

    ! Fill origin
    s%n = xnuc

    s%nv = 2*nphi*(2**ntheta-1)+2
    if (s%nv > s%mv) then
       call ferror('spheresphere','max. number of vertex exceeded',faterr)
    end if

    s%nf = 6*nphi*(2**(ntheta-1)-1)+nphi+nphi
    if (s%nf > s%mf) then
       call ferror('spheresphere','maximum number of faces exceeded',faterr)
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
    if (s%isinit > 1) then
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

  end subroutine spheresphere

  !> Triangulate the unit sphere by recursively subdividing an
  !> octahedron level times. Minisurface version.
  module subroutine spheretriang(s,xnuc,level)
    use tools_io, only: faterr, ferror, warning, uout
    class(minisurf), intent(inout) :: s
    real*8, dimension(3), intent(in) :: xnuc
    integer, intent(in) :: level

    integer, parameter :: mlevel = 7

    integer :: nvert
    integer :: nface(0:1)
    integer :: iface(3,s%mf,0:1)
    integer :: ix(s%mv), iy(s%mv), iz(s%mv)
    integer :: i1, i2, i3, i4, i5, i6, ilev
    integer :: ix4, iy4, iz4, ix5, iy5, iz5, ix6, iy6, iz6
    integer :: isw0, isw1
    integer :: nf1
    logical :: new
    integer :: i

    call s%clean()

    ! Fill origin
    s%n = xnuc

    !.....First, let us compute a convex polyhedra that approximately
    !     triangulates the sphere. We will start with a octahedron and
    !     succesively divide each triangular facet into four new ones.
    !     The coordinates for the vertices will be maintained as integer
    !     values to get best precision and efficiency.
    nvert = 6
    if (nvert > s%mv) &
       call ferror('spheretriang','max. number of vertex exceeded',faterr)
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
       call ferror('spheretriang','maximum number of faces exceeded',faterr)
    end if
    if (s%isinit > 1) then
       do i=1,s%nf
          !! Vertices, id
          s%f(i)%nv   = 3
          s%f(i)%v(1) = iface(1,i,isw1)
          s%f(i)%v(2) = iface(2,i,isw1)
          s%f(i)%v(3) = iface(3,i,isw1)
       end do
    end if

  end subroutine spheretriang

  !> Triangulate the unit sphere by recursively subdividing a cube
  !> level times. Minisurface version.
  module subroutine spherecub(s,xnuc,level)
    use tools_io, only: ferror, warning, uout, faterr
    class(minisurf), intent(inout) :: s
    real*8, dimension(3), intent(in) :: xnuc
    integer, intent(in) :: level

    integer, parameter :: mlevel = 7

    integer :: nvert
    integer :: nface(0:1)
    integer :: iface(4,s%mf,0:1)
    integer :: ix(s%mv), iy(s%mv), iz(s%mv)
    integer :: ix5, iy5, iz5, ix6, iy6, iz6
    integer :: ix7, iy7, iz7, ix8, iy8, iz8, ix9, iy9, iz9
    integer :: i, i1, i2, i3, i4, i5, i6, i7, i8, i9, ilev
    integer :: isw0, isw1
    integer :: nf1
    logical :: new

    call s%clean()

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
       call ferror('spherecub','maximum number of faces exceeded',faterr)
    end if
    if (s%isinit > 1) then
       do i=1,s%nf
          !! Vertices, id
          s%f(i)%nv   = 4
          s%f(i)%v(1) = iface(1,i,isw1)
          s%f(i)%v(2) = iface(2,i,isw1)
          s%f(i)%v(3) = iface(3,i,isw1)
          s%f(i)%v(4) = iface(4,i,isw1)
       end do
    end if

  end subroutine spherecub

  !> Write the radii of the rays of the minisurface s to an INT file.
  !> These files are used to write and read the IAS information
  !> for a particular quadrature method (meth) with n1 / n2 nodes.
  module subroutine writeint (s,n1,n2,meth,offfile)
    use tools_io, only: fopen_write, fclose
    class(minisurf), intent(inout) :: s
    integer, intent(in) :: n1, n2, meth
    character*(*), intent(in) :: offfile

    integer :: lud
    integer :: i

    lud = fopen_write(offfile)

    write (lud,'(I10,x,I10,x,I2)') n1, n2, meth
    write (lud,'(3(E23.15,x))') s%n(1), s%n(2), s%n(3)
    write (lud,'(3(E23.15,x))') (s%r(i),i=1,s%nv)

    call fclose(lud)

  end subroutine writeint

  !> Read the IAS radii of the rays in the INT file into the
  !> minisurface s. The method fingerprint (meth, n1, n2) must match
  !> the one contained in the file.
  module subroutine readint (s,n1,n2,meth,offfile,ierr)
    use global, only: int_gauleg
    use tools_io, only: fopen_read, ferror, warning, string, fclose
    class(minisurf), intent(inout) :: s
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
       call ferror ('readint', 'Wrong ntheta, nphi or method',warning,string(offfile))
       ierr = 1
       return
    end if
    read (lud,'(3(E23.15,x))') xn
    if (any(abs(xn - s%n) > 1d-13)) then
       call ferror ('readint', 'Wrong CP position in integrals file',warning,string(offfile))
       ierr = 1
       return
    end if
    read (lud,'(3(E23.15,x))') (s%r(i),i=1,s%nv)

    call fclose(lud)

  end subroutine readint

  !> Sets the angular nodes for the quadrature of a sphere by treating
  !> ntheta and nphi as 1d Gauss-Legendre quadratures. For a given
  !> theta (polar angle), the number of phi (azimuthal angle) points
  !> is int(nphi*abs(sin(theta)))+1, where nphi is provided in the input.
  !> minisurface version
  module subroutine gauleg_nodes(srf,ntheta,nphi)
    use param, only: pi
    use tools_math, only: gauleg
    class(minisurf), intent(inout) :: srf !< Surface where the node info is placed
    integer, intent(in) :: ntheta !< Number of phi points
    integer, intent(in) :: nphi   !< Number of azimuthal points

    integer :: j, k
    integer :: realnphi
    ! Nodes and weights
    real*8, allocatable, dimension(:) :: tpoints
    real*8, allocatable, dimension(:) :: tweights
    real*8, allocatable, dimension(:) :: ppoints
    real*8, allocatable, dimension(:) :: pweights

    allocate(tpoints(ntheta),tweights(ntheta))
    allocate(ppoints(nphi+1),pweights(nphi+1))
    ! place the gauss-legendre nodes on the surface
    call gauleg(0d0,pi,tpoints,tweights,ntheta)
    srf%nv = 0
    do j = 1, ntheta
       realnphi = int(nphi*abs(sin(tpoints(j))))+1
       call gauleg(0d0,2d0*pi,ppoints,pweights,realnphi)
       do k = 1, realnphi
          srf%nv = srf%nv + 1
          srf%r(srf%nv)  = -1d0
          srf%th(srf%nv) = tpoints(j)
          srf%ph(srf%nv) = ppoints(k)
       end do
    end do

    deallocate(tpoints,tweights,ppoints,pweights)

  end subroutine gauleg_nodes

  !> Sets the angular nodes for the quadrature of a sphere using
  !> Lebedev method (see above for references), where npts is the
  !> number of points, that must be one for which a cubature
  !> exists (6, 14, 26, 38, 50, 74,
  !> 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974,
  !> 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,
  !> 5294, 5810). The (theta,phi) information is written to the surface.
  !> minisurface version.
  module subroutine lebedev_nodes(srf,npts)
    use tools_math, only: select_lebedev
    class(minisurf), intent(inout) :: srf !< Surface where the node info is placed
    integer, intent(in) :: npts !< Number of quadrature points

    integer :: j
    ! Nodes and weights
    real*8 :: xleb(5810)
    real*8 :: yleb(5810)
    real*8 :: zleb(5810)
    real*8 :: wleb(5810)

    call select_lebedev(npts,xleb,yleb,zleb,wleb)

    srf%nv = 0
    do j = 1, npts
       srf%nv = srf%nv + 1
       srf%r(srf%nv)  = -1d0
       srf%th(srf%nv) = acos(zleb(j))
       srf%ph(srf%nv) = atan2(yleb(j),xleb(j))
    end do

  end subroutine lebedev_nodes

end submodule proc
