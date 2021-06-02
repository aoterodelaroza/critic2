! This module contains Bader integration-on-a-grid as proposed by
! Henkelman and collaborators. The code in this module has been
! adapted from the 'bader' program (version 0.28a, 07/12/12), that
! can be located at the Henkelman's group webpage:
!    http://theory.cm.utexas.edu/
!    http://theory.cm.utexas.edu/vasp/bader/
! The authors of bader are Wenjie Tang, Andri Arnaldsson, Samuel T.
! Chill, and Graeme Henkelman.
! Also, see:
!   Comput. Mater. Sci. 36, 254-360 (2006).
!   J. Comput. Chem. 28, 899-908 (2007).
!   J. Phys.: Condens. Matter 21, 084204 (2009)

! Copyright 2009
! Wenjie Tang, Andri Arnaldsson, Samuel T. Chill, and Graeme Henkelman
!
! Bader is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

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

submodule (bader) proc
  implicit none

  !xx! private procedure
  ! subroutine refine_edge(f,irefine_edge,ref_itrs)
  ! subroutine max_neargrid(f,p)
  ! subroutine step_neargrid(f,p)
  ! subroutine step_ongrid(f,p)
  ! function rho_grad_dir(f,p) result(res)
  ! function is_max(f,p)
  ! subroutine pbc(p)
  ! function rho_val(ff,p1,p2,p3)
  ! function volnum_val(p1,p2,p3)
  ! subroutine assign_surrounding_pts(p)
  ! subroutine known_volnum_ongrid(p)
  ! function is_vol_edge(p)
  ! subroutine reassign_volnum_ongrid2(p)

  ! private to the module, initialized at the beginning of bader_integrate
  integer, allocatable :: volnum(:,:,:) !< Bader volume identifier
  integer, allocatable :: known(:,:,:) !< Is the point known?
  integer, allocatable :: path(:,:) !< A path through the grid
  integer :: pnum !< number of points in the path
  integer :: nbasin !< number of bader maxima
  real*8, dimension(-1:1,-1:1,-1:1) :: lat_dist !< distance between neighbor grid points
  real*8, dimension(-1:1,-1:1,-1:1) :: lat_i_dist !< inverse of that
  real*8 :: lat2car(3,3) !< from integer to cartesian
  real*8 :: car2lat(3,3) !< from cartesian to integer
  integer :: n(3) !< grid dimensions

contains

  !> Do a grid integration using the BADER method. s is the system and
  !> id is the field id. Return the number of basins (nbasin0), their
  !> coordinates (crystallographic corods, xcoord). volnum0 gives the
  !> id of the basin (from 1 to nbasin0) on the lattice. If the
  !> arithmetic expression discexpr is not empty, then apply that
  !> expression to the basin attractors. If the expression is
  !> non-zero, discard the attractor. If atexist is true, then the
  !> code is aware of the presence of atoms, which are added as
  !> attractors at the beginning of the run. Two attractors are
  !> considered equal if they are within a ditsance of ratom (bohr).
  module subroutine bader_integrate(s,bas,iref)
    use systemmod, only: system
    use tools_io, only: faterr, ferror
    use tools_math, only: matinv
    use arithmetic, only: eval
    use param, only: vsmall, icrd_crys
    use types, only: realloc, basindat
    type(system), intent(inout) :: s
    type(basindat), intent(inout) :: bas
    integer, intent(in) :: iref

    integer :: i, j, k, l, path_volnum, p(3)
    integer :: ptemp(3), ref_itrs, irefine_edge, nid
    real*8 :: dlat(3), dcar(3), dv(3), x(3), fval
    integer :: bat(s%c%ncel)
    logical :: isassigned, ok

    if (.not.s%isinit) &
       call ferror("bader_integrate","system not initialized",faterr)
    if (.not.associated(s%c)) &
       call ferror("bader_integrate","system does not have crystal",faterr)

    ! deallocate the arguments and private globals
    if (allocated(bas%idg)) deallocate(bas%idg)
    if (allocated(volnum)) deallocate(volnum)
    if (allocated(known)) deallocate(known)
    if (allocated(path)) deallocate(path)

    ! Pre-allocate atoms and nnm as maxima
    allocate(bas%xattr(3,s%f(iref)%ncpcel))
    bas%xattr = 0d0
    nbasin = 0
    if (bas%atexist) then
       nbasin = s%f(iref)%ncpcel
       do i = 1, s%f(iref)%ncpcel
          bas%xattr(:,i) = s%f(iref)%cpcel(i)%x
       end do
    end if

    ! initialize
    bat = 0

    ! metrics
    do i = 1, 3
       n(i) = size(bas%f,i)
       lat2car(:,i) = s%c%m_x2c(:,i) / n(i)
    end do
    car2lat = lat2car
    call matinv(car2lat,3)

    ! distance between neighboring points
    lat_i_dist = 0d0
    do i= -1, 1
       dlat(1)=real(i,8)
       do j= -1, 1
          dlat(2)=real(j,8)
          do k= -1, 1
             dlat(3)=real(k,8)
             dcar = matmul(lat2car,dlat)
             lat_dist(i,j,k) = sqrt(sum(dcar*dcar))
             if ((i == 0).and.(j == 0).and.(k == 0)) cycle
             lat_i_dist(i,j,k) = 1d0 / lat_dist(i,j,k)
          end do
       end do
    end do

    allocate(volnum(n(1),n(2),n(3)),known(n(1),n(2),n(3)),path(3,10))
    volnum = 0
    known = 0

    do i = 1, n(1)
       do j = 1, n(2)
          do k = 1, n(3)
             p = (/i, j, k/)
             if (volnum(i,j,k) == 0) then
                call max_neargrid(bas%f,p)
                path_volnum = volnum(p(1),p(2),p(3))

                ! maximum
                if (path_volnum == 0) then
                   dv = real(p-1,8) / n

                   ! check if it is an atom (use ratom)
                   isassigned = .false.
                   if (bas%atexist) then
                      nid = s%c%identify_atom(dv,icrd_crys,distmax=bas%ratom)
                      if (nid > 0) then
                         path_volnum = nid
                         isassigned = .true.
                      end if
                   end if
                   ! check if it is a known nnm
                   if (.not.isassigned .and. bas%ratom > vsmall) then
                      do l = 1, nbasin
                         if (s%c%are_lclose(dv,bas%xattr(:,l),bas%ratom)) then
                            path_volnum = l
                            isassigned = .true.
                            exit
                         end if
                      end do
                   end if
                   ! well, it must be a new attractor then
                   if (.not.isassigned) then
                      ok = .true.
                      if (len_trim(bas%expr) > 0) then
                         x = s%c%x2c(dv)
                         fval = s%eval(bas%expr,.false.,ok,x)
                         if (.not.ok) &
                            call ferror("yt","invalid DISCARD expression",faterr)
                         ok = (abs(fval) < 1d-30)
                      end if
                      if (ok) then
                         nbasin = nbasin + 1
                         if (nbasin > size(bas%xattr,2)) call realloc(bas%xattr,3,2*nbasin)
                         path_volnum = nbasin
                         bas%xattr(:,nbasin) = dv
                      end if
                   end if
                end if

                ! assign all points along the trajectory
                do l = 1, pnum
                   ptemp = path(:,l)
                   if (volnum(ptemp(1),ptemp(2),ptemp(3)) /= -1) then
                      volnum(ptemp(1),ptemp(2),ptemp(3)) = path_volnum
                   end if
                   if (known(ptemp(1),ptemp(2),ptemp(3)) /= 2) then
                      known(ptemp(1),ptemp(2),ptemp(3)) = 0
                   end if
                   call assign_surrounding_pts(ptemp)
                end do
             end if
          end do
       end do
    end do

    ! refine
    ref_itrs = 1
    irefine_edge = -1
    do while (.true.)
       call refine_edge(bas%f,irefine_edge,ref_itrs)
       if (irefine_edge == 0) exit
       ref_itrs = ref_itrs + 1
    end do

    ! wrap up
    deallocate(known,path)
    call realloc(bas%xattr,3,nbasin)
    call move_alloc(volnum,bas%idg)
    bas%nattr = nbasin
    if (allocated(known)) deallocate(known)
    if (allocated(path)) deallocate(path)

  end subroutine bader_integrate

  !> Remap the attractors from a bader calculation
  module subroutine bader_remap(s,bas,nattn,idg1,ilvec,iatt)
    use types, only: realloc, basindat
    type(system), intent(in) :: s
    type(basindat), intent(in) :: bas
    integer, intent(out) :: nattn
    integer, allocatable, intent(inout) :: iatt(:), ilvec(:,:), idg1(:,:,:)

    integer :: i, m1, m2, m3, p(3)
    real*8 :: x(3), xs(3), d2
    logical :: found

    if (allocated(iatt)) deallocate(iatt)
    allocate(iatt(bas%nattr))
    nattn = bas%nattr
    do i = 1, bas%nattr
       iatt(i) = i
    enddo
    if (allocated(ilvec)) deallocate(ilvec)
    allocate(ilvec(3,bas%nattr))
    ilvec = 0

    if (allocated(idg1)) deallocate(idg1)
    allocate(idg1(bas%n(1),bas%n(2),bas%n(3)))
    do m3 = 1, bas%n(3)
       do m2 = 1, bas%n(2)
          do m1 = 1, bas%n(1)
             idg1(m1,m2,m3) = bas%idg(m1,m2,m3)
             p = (/m1,m2,m3/)
             x = real(p-1,8) / bas%n - bas%xattr(:,bas%idg(m1,m2,m3))
             xs = x
             call s%c%shortest(xs,d2)
             p = nint(x - s%c%c2x(xs))
             if (any(p /= 0)) then
                found = .false.
                do i = bas%nattr+1, nattn
                   if (iatt(i) == bas%idg(m1,m2,m3) .and. all(p == ilvec(:,i))) then
                      found = .true.
                      idg1(m1,m2,m3) = i
                      exit
                   end if
                end do
                if (.not.found) then
                   nattn = nattn + 1
                   if (nattn > size(ilvec,2)) then
                      call realloc(ilvec,3,2*nattn)
                      call realloc(iatt,2*nattn)
                   end if
                   ilvec(:,nattn) = p
                   idg1(m1,m2,m3) = nattn
                   iatt(nattn) = bas%idg(m1,m2,m3)
                end if
             end if
          end do
       end do
    end do
    call realloc(ilvec,3,nattn)
    call realloc(iatt,nattn)

  endsubroutine bader_remap

  !xx! private procedure

  subroutine refine_edge(f,irefine_edge,ref_itrs)
    use tools_io, only: faterr, ferror
    real*8, intent(in) :: f(:,:,:)
    integer, intent(inout) :: irefine_edge
    integer, intent(inout) :: ref_itrs

    integer :: p(3), pt(3)
    integer :: n1, n2, n3, path_volnum, bvolnum, i
    integer :: num_edge, num_reassign, num_check
    integer :: d1, d2, d3

    if(ref_itrs == 1) then
       num_edge = 0
       do n1 = 1, n(1)
          do n2 = 1, n(2)
             do n3 = 1, n(3)
                p = (/n1,n2,n3/)
                ! change for calculating the vacuum volume
                if (volnum(n1,n2,n3) == nbasin + 1) cycle
                if (is_vol_edge(p) .and. (.not.is_max(f,p))) then
                   num_edge = num_edge + 1
                   volnum(p(1),p(2),p(3)) = -volnum(p(1),p(2),p(3))
                   known(p(1),p(2),p(3)) = 0
                   call reassign_volnum_ongrid2(p)
                end if
             end do
          end do
       end do
    else
       num_check=0
       do n1 = 1, n(1)
          do n2 = 1, n(2)
             do n3 = 1, n(3)
                p = (/n1,n2,n3/)
                ! change for calculating the vacuum volume
                if (volnum(n1,n2,n3) == nbasin+1) cycle

                if(volnum(n1,n2,n3) < 0 .and. known(n1,n2,n3) /=-1) then
                   do d1 = -1,1
                      do d2 = -1,1
                         do d3 = -1,1
                            pt = p + (/d1,d2,d3/)
                            call pbc(pt)
                            ! change for calculating the vacuum volume
                            if (volnum(pt(1),pt(2),pt(3)) == nbasin+1) cycle
                            if(.not.is_max(f,pt)) then
                               if(volnum(pt(1),pt(2),pt(3)) > 0) then
                                  volnum(pt(1),pt(2),pt(3)) = -volnum(pt(1),pt(2),pt(3))
                                  known(pt(1),pt(2),pt(3)) = -1
                                  num_check=num_check+1
                               else if(volnum(pt(1),pt(2),pt(3))<0 .and. known(pt(1),pt(2),pt(3)) == 0) then
                                  known(pt(1),pt(2),pt(3)) = -2
                                  num_check = num_check + 1
                               end if
                            end if
                         end do
                      end do
                   end do
                   num_check = num_check - 1
                   if (known(pt(1),pt(2),pt(3)) /= -2) then
                      volnum(p(1),p(2),p(3)) = abs(volnum(p(1),p(2),p(3)))
                   end if
                   ! end of mark
                end if
             end do
          end do
       end do

       ! make the surrounding points unknown
       do n1 = 1, n(1)
          do n2 = 1, n(2)
             do n3 = 1, n(3)
                p = (/n1,n2,n3/)
                bvolnum = volnum(n1,n2,n3)

                if (bvolnum < 0) then
                   do d1 = -1,1
                      do d2 = -1,1
                         do d3 = -1,1
                            pt = p + (/d1,d2,d3/)
                            call pbc(pt)
                            if(known(pt(1),pt(2),pt(3)) == 2) known(pt(1),pt(2),pt(3)) = 0
                         end do
                      end do
                   end do
                end if
             end do
          end do
       end do
    end if

    num_reassign = 0
    do n1 = 1, n(1)
       do n2 = 1, n(2)
          do n3 = 1, n(3)
             p = (/n1,n2,n3/)
             bvolnum = volnum(n1,n2,n3)
             if (bvolnum < 0) then
                call max_neargrid(f,p)
                path_volnum = volnum(p(1),p(2),p(3))
                if (path_volnum < 0 .or. path_volnum > nbasin) then
                   call ferror('refine_edge','should be no new maxima in edge refinement',faterr)
                end if
                volnum(n1,n2,n3) = path_volnum
                if (abs(bvolnum) /= path_volnum) then
                   num_reassign = num_reassign + 1
                   volnum(n1,n2,n3) = -path_volnum
                end if
                do i = 1,pnum
                   pt = path(:,i)
                   if (known(pt(1),pt(2),pt(3)) /= 2) then
                      known(pt(1),pt(2),pt(3)) = 0
                   end if
                end do
             end if
          end do
       end do
    end do

    ! flag to indicate that we are done refining
    if (num_reassign==0) irefine_edge = 0

  end subroutine refine_edge

  ! max_neargrid
  ! From the point p do a maximization on the charge density grid and
  ! assign the maximum found to the volnum array.
  subroutine max_neargrid(f,p)
    use types, only: realloc
    real*8, intent(in) :: f(:,:,:)
    integer, dimension(3), intent(inout) :: p

    pnum = 1
    path(:,pnum) = p

    do while (.true.)
       call step_neargrid(f,p)

       ! if we didn't move, we're at a maximum
       if (all(p == path(:,pnum))) exit

       ! otherwise, add point to path
       pnum = pnum + 1
       if (pnum > size(path,2)) call realloc(path,3,2*pnum)
       path(:,pnum) = p

       ! quit at a known point
       if (known(p(1),p(2),p(3)) == 2) exit
    end do

  end subroutine max_neargrid

  ! step_neargrid
  ! Do a single iteration of a maximization on the charge density
  ! grid from the point (px,py,pz).
  subroutine step_neargrid(f,p)
    real*8, intent(in) :: f(:,:,:)
    integer, intent(inout) :: p(3)

    real*8, save :: dr(3)
    real*8 :: gradrl(3), coeff
    integer :: pm(3)

    if (pnum == 1) then
      dr = 0d0
    end if
    gradrl = rho_grad_dir(f,p)

    if (maxval(abs(gradrl)) < 1d-30) then
       dr = 0d0
       if (is_max(f,p)) then
          return
       else
          pm = p
          call step_ongrid(f,pm)
       end if
    else
       coeff = 1d0/maxval(abs(gradrl))
       gradrl = coeff*gradrl
       pm = p + nint(gradrl)
       dr = dr + gradrl - nint(gradrl)
       pm = pm + nint(dr)
       dr = dr - nint(dr)
    end if
    known(p(1),p(2),p(3)) = 1

    call pbc(pm)
    if (known(pm(1),pm(2),pm(3)) == 1) then
       pm=p
       call step_ongrid(f,pm)
       dr = 0d0
    end if
    p = pm

  end subroutine step_neargrid

  ! step_ongrid
  ! Do a single iteration of a maximization on the charge density
  ! grid from the point (px,py,pz).  Return a logical indicating
  ! if the current  point is a charge density maximum.
  subroutine step_ongrid(f,p)
    real*8, intent(in) :: f(:,:,:)
    integer, intent(inout) :: p(3)

    integer :: pm(3), pt(3)
    real*8 :: rho_max, rho_tmp, rho_ctr
    integer :: d1, d2, d3

    pm = p
    rho_ctr = rho_val(f,p(1),p(2),p(3))
    rho_max = rho_ctr
    do d1 = -1, 1
       do d2 = -1, 1
          do d3 = -1, 1
             pt = p+(/d1,d2,d3/)
             rho_tmp = rho_val(f,pt(1),pt(2),pt(3))
             rho_tmp = rho_ctr+(rho_tmp-rho_ctr)*lat_i_dist(d1,d2,d3)
             if (rho_tmp > rho_max) then
                rho_max = rho_tmp
                pm = pt
             end if
          end do
       end do
    end do
    call pbc(pm)
    p = pm

  end subroutine step_ongrid

  !rho_grad_dir
  ! Return the direction of the gradient in lattice vectors
  ! at the grid position p
  function rho_grad_dir(f,p) result(res)
    real*8, intent(in) :: f(:,:,:)
    integer, intent(in) :: p(3)
    real*8 :: res(3)

    integer :: p1, p2, p3
    real*8 :: rho000, rho001, rho010, rho100, rho00_1, rho_100, rho0_10
    real*8 :: rho_grad_lat(3), rho_grad_car(3)

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)

    rho000 = rho_val(f,p1,p2,p3)
    rho001 = rho_val(f,p1,p2,p3+1)
    rho010 = rho_val(f,p1,p2+1,p3)
    rho100 = rho_val(f,p1+1,p2,p3)
    rho00_1 = rho_val(f,p1,p2,p3-1)
    rho_100 = rho_val(f,p1-1,p2,p3)
    rho0_10 = rho_val(f,p1,p2-1,p3)

    rho_grad_lat(1) = (rho100-rho_100)/2d0
    rho_grad_lat(2) = (rho010-rho0_10)/2d0
    rho_grad_lat(3) = (rho001-rho00_1)/2d0

    if (rho100 < rho000.and.rho_100 < rho000) rho_grad_lat(1) = 0d0
    if (rho010 < rho000.and.rho0_10 < rho000) rho_grad_lat(2) = 0d0
    if (rho001 < rho000.and.rho00_1 < rho000) rho_grad_lat(3) = 0d0

    ! convert to cartesian coordinates
    rho_grad_car = matmul(rho_grad_lat,car2lat)

    ! express this vector in direct coordinates
    res = matmul(car2lat,rho_grad_car)

  end function rho_grad_dir

  ! is_max
  ! return .true. if the grid point is a maximum of charge density
  function is_max(f,p)
    real*8, intent(in) :: f(:,:,:)
    integer, intent(in) :: p(3)
    logical :: is_max

    real*8 :: rho
    integer :: d1, d2, d3, p1, p2, p3

    is_max=.true.
    p1 = p(1)
    p2 = p(2)
    p3 = p(3)
    rho=rho_val(f,p1,p2,p3)
    do d1=-1,1
      p1=p(1)+d1
      do d2=-1,1
        p2=p(2)+d2
        do d3=-1,1
          p3=p(3)+d3
          if(rho_val(f,p1,p2,p3) > rho) then
            is_max = .false.
          end if
        end do
      end do
    end do

  end function is_max

  ! pbc
  ! Wrap the point (p(1),p(2),p(3)) to the boundary conditions [0,pmax].
  subroutine pbc(p)
    integer, intent(inout) :: p(3)

    integer :: i

    do i = 1, 3
      do
        if(p(i) > 0) exit
        p(i) = p(i) + n(i)
      end do
      do
        if(p(i) <= n(i)) exit
        p(i) = p(i) - n(i)
      end do
    end do

  end subroutine pbc

  function rho_val(ff,p1,p2,p3)
    real*8, intent(in) :: ff(:,:,:)
    integer, intent(in) :: p1, p2, p3
    real*8 :: rho_val

    integer :: p(3), i

    p=(/p1,p2,p3/)
    do i = 1, 3
       do while (.true.)
          if(p(i) >= 1) exit
          p(i) = p(i) + n(i)
       end do
       do while (.true.)
          if(p(i) <= n(i)) exit
          p(i) = p(i) - n(i)
       end do
    end do
    rho_val = ff(p(1),p(2),p(3))

  end function rho_val

  function volnum_val(p1,p2,p3)
    integer, intent(in) :: p1, p2, p3
    integer :: volnum_val

    integer :: p(3), i

    p=(/p1,p2,p3/)
    do i = 1, 3
       do while (.true.)
          if(p(i) >= 1) exit
          p(i) = p(i) + n(i)
       end do
       do while (.true.)
          if(p(i) <= n(i)) exit
          p(i) = p(i) - n(i)
       end do
    end do
    volnum_val = volnum(p(1),p(2),p(3))

  end function volnum_val

  ! assign_surrounding_pts
  ! check the surrounding points of p to see if their volnum
  ! is known
  subroutine assign_surrounding_pts(p)
    integer, intent(in) :: p(3)
    integer :: pt(3)

    pt = p + (/1,0,0/)
    call pbc(pt)
    if (known(pt(1),pt(2),pt(3)) /= 2) then
       call known_volnum_ongrid(pt)
    end if
    pt = p + (/-1,0,0/)
    call pbc(pt)
    if(known(pt(1),pt(2),pt(3)) /= 2) then
       call known_volnum_ongrid(pt)
    end if
    pt = p + (/0,1,0/)
    call pbc(pt)
    if(known(pt(1),pt(2),pt(3)) /= 2) then
       call known_volnum_ongrid(pt)
    end if
    pt = p + (/0,-1,0/)
    call pbc(pt)
    if(known(pt(1),pt(2),pt(3)) /= 2) then
       call known_volnum_ongrid(pt)
    end if
    pt = p + (/0,0,1/)
    call pbc(pt)
    if(known(pt(1),pt(2),pt(3)) /= 2) then
       call known_volnum_ongrid(pt)
    end if
    pt = p + (/0,0,-1/)
    call pbc(pt)
    if(known(pt(1),pt(2),pt(3)) /= 2) then
       call known_volnum_ongrid(pt)
    end if

  end subroutine assign_surrounding_pts

  ! known_volnum_ongrid
  ! return number of the associated bader volnum if nearest
  ! grid points are known to be associated with the same bader volnum
  subroutine known_volnum_ongrid(p)
    integer, intent(in) :: p(3)

    integer :: volnum_, p1, p2, p3

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)

    volnum_ = volnum_val(p1,p2,p3)
    if(volnum_ <= 0) return

    if (volnum_val(p1,p2,p3+1) /= volnum_) return
    if (volnum_val(p1,p2,p3-1) /= volnum_) return
    if (volnum_val(p1,p2+1,p3) /= volnum_) return
    if (volnum_val(p1,p2-1,p3) /= volnum_) return
    if (volnum_val(p1+1,p2,p3) /= volnum_) return
    if (volnum_val(p1-1,p2,p3) /= volnum_) return

    known(p1,p2,p3) = 2

  end subroutine known_volnum_ongrid

  ! is_vol_edge
  ! return .true. if the grid point is on the edge of a Bader volume.
  function is_vol_edge(p)
    logical :: is_vol_edge
    integer, intent(in) :: p(3)

    integer :: d1, d2, d3, volnum_, volnbr, pt(3)

    volnum_ = volnum(p(1),p(2),p(3))
    is_vol_edge = .false.
    neighborloop: do d1 = -1,1
       do d2 = -1,1
          do d3 = -1,1
             pt = p + (/d1,d2,d3/)
             call pbc(pt)
             volnbr = volnum(pt(1),pt(2),pt(3))
             if (abs(volnbr) /= abs(volnum_)) then
                is_vol_edge = .true.
                exit neighborloop
             end if
          end do
       end do
    end do neighborloop

  end function is_vol_edge

  ! reassign_volnum_ongrid
  ! reassign the surrounding points of a edge point as unknown points
  subroutine reassign_volnum_ongrid2(p)
    integer, intent(in) :: p(3)

    integer :: d1, d2, d3, pt(3)

    do d1 = -1,1
       do d2 = -1,1
          do d3 = -1,1
             pt = p + (/d1,d2,d3/)
             call pbc(pt)
             known(pt(1),pt(2),pt(3)) = 0
          end do
       end do
    end do

  end subroutine reassign_volnum_ongrid2

end submodule proc
