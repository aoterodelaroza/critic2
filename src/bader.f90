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
module bader
  implicit none

  private

  public :: bader_integrate

  ! private to the module, initialized at the beginning of bader_integrate
  integer, allocatable :: volnum(:,:,:) !< Bader volume identifier
  integer, allocatable :: known(:,:,:) !< Is the point known?
  integer, allocatable :: path(:,:) !< A path through the grid
  integer :: pnum !< number of points in the path
  integer :: bnum !< number of bader maxima
  real*8, dimension(-1:1,-1:1,-1:1) :: lat_dist !< distance between neighbor grid points
  real*8, dimension(-1:1,-1:1,-1:1) :: lat_i_dist !< inverse of that
  real*8 :: lat2car(3,3) !< from integer to cartesian
  real*8 :: car2lat(3,3) !< from cartesian to integer
  integer :: n(3) !< grid dimensions

contains

  !> Do a grid integration using the BADER method. cr is the crystal
  !> and f is the field. Return the number of basins (bnum0), their
  !> coordinates (crystallographic corods, volpos_lat) and the basin
  !> id of every point in the grid (idatt, positive corresponds to an
  !> atom and negative to an nnm). volnum0 gives the id of the basin
  !> (from 1 to bnum0) on the lattice.  ratom is the minimum distance
  !> for atom recognition.
  subroutine bader_integrate(c,f,bnum0,volpos_lat,idatt,volnum0,ratom)
    use global
    use struct_basic
    use tools_math
    use types

    type(crystal), intent(in) :: c
    type(field), intent(in) :: f
    integer, intent(out) :: bnum0
    real*8, allocatable, intent(inout) :: volpos_lat(:,:)
    integer, allocatable, intent(inout) :: idatt(:)
    integer, allocatable, intent(inout) :: volnum0(:,:,:)
    real*8, intent(in) :: ratom

    integer :: i, j, k, l, path_volnum, p(3)
    integer :: ptemp(3), ref_itrs, irefine_edge, nid, lvec(3)
    real*8 :: dlat(3), dcar(3), dist, dv(3)
    integer :: bat(c%ncel), nnnm

    ! initialize
    bat = 0

    ! metrics
    n = f%n
    do i = 1, 3
       lat2car(:,i) = c%crys2car(:,i) / n(i)
    end do
    car2lat = matinv(lat2car)

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

    allocate(volnum(n(1),n(2),n(3)),volpos_lat(3,10),idatt(10),known(n(1),n(2),n(3)),path(3,10))
    volnum = 0
    known = 0
    bnum = 0
    nnnm = 0

    do i = 1, n(1)
       do j = 1, n(2)
          do k = 1, n(3)
             p = (/i, j, k/)
             if (volnum(i,j,k) == 0) then
                call max_neargrid(f,p)
                path_volnum = volnum(p(1),p(2),p(3))

                ! maximum
                if (path_volnum == 0) then
                   dv = real(p-1,8) / n
                   nid = 0
                   call c%nearest_atom(dv,nid,dist,lvec)
                   if (dist <= ratom) then
                      if (bat(nid) > 0) then
                         path_volnum = bat(nid)
                      else
                         bnum = bnum + 1
                         if (bnum > size(volpos_lat,2)) call realloc(volpos_lat,3,2*bnum)
                         if (bnum > size(idatt)) call realloc(idatt,2*bnum)
                         bat(nid) = bnum
                         path_volnum = bnum
                         volpos_lat(:,bnum) = c%atcel(nid)%x * n
                         idatt(bnum) = nid
                      endif
                   else
                      bnum = bnum + 1
                      nnnm = nnnm + 1
                      if (bnum > size(volpos_lat,2)) call realloc(volpos_lat,3,2*bnum)
                      if (bnum > size(idatt)) call realloc(idatt,2*bnum)
                      path_volnum = bnum
                      volpos_lat(:,bnum) = real(p-1,8)
                      idatt(bnum) = -nnnm
                   endif
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
       call refine_edge(f,irefine_edge,ref_itrs)
       if (irefine_edge == 0) exit
       ref_itrs = ref_itrs + 1
    end do

    ! wrap up
    deallocate(known,path)
    call realloc(volpos_lat,3,bnum)
    do i = 1, bnum
       volpos_lat(:,i) = volpos_lat(:,i) / n
    end do
    call realloc(idatt,bnum)
    call move_alloc(volnum,volnum0)
    bnum0 = bnum

  end subroutine bader_integrate

  subroutine refine_edge(f,irefine_edge,ref_itrs)
    use tools_io
    use types

    type(field), intent(in) :: f
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
                if (volnum(n1,n2,n3) == bnum + 1) cycle
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
                if (volnum(n1,n2,n3) == bnum+1) cycle
                
                if(volnum(n1,n2,n3) < 0 .and. known(n1,n2,n3) /=-1) then
                   do d1 = -1,1
                      do d2 = -1,1
                         do d3 = -1,1
                            pt = p + (/d1,d2,d3/)
                            call pbc(pt)
                            ! change for calculating the vacuum volume
                            if (volnum(pt(1),pt(2),pt(3)) == bnum+1) cycle
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

       ! make the surrounding points unkown
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
                if (path_volnum < 0 .or. path_volnum > bnum) then
                   call ferror('refine_edge','should be no new maxima in edge refinement',2)
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
    use types

    type(field), intent(in) :: f
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
    use types
    
    type(field), intent(in) :: f
    integer,dimension(3),intent(inout) :: p

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
  SUBROUTINE step_ongrid(f,p)
    use types

    type(field), intent(in) :: f
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

  END SUBROUTINE step_ongrid

  !rho_grad_dir
  ! Return the direction of the gradient in lattice vectors
  ! at the grid position p
  function rho_grad_dir(f,p)
    use types

    type(field), intent(in) :: f
    integer, intent(in) :: p(3)
    real*8 :: rho_grad_dir(3)

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
    rho_grad_dir = matmul(car2lat,rho_grad_car)

  end function rho_grad_dir

  ! is_max
  ! return .true. if the grid point is a maximum of charge density
  function is_max(f,p)
    use types

    type(field), intent(in) :: f
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

  function rho_val(f,p1,p2,p3)
    use types

    type(field), intent(in) :: f
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
    rho_val = f%f(p(1),p(2),p(3))

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

end module bader
