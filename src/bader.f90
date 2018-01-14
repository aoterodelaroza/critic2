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
module bader
  implicit none

  private

  public :: bader_integrate
  
  interface
     module subroutine bader_integrate(s,ff,discexpr,atexist,ratom,nbasin0,xcoord,volnum0)
       use systemmod, only: system
       type(system), intent(inout) :: s
       real*8, intent(in) :: ff(:,:,:)
       character*(*), intent(in) :: discexpr
       logical, intent(in) :: atexist
       real*8, intent(in) :: ratom
       integer, intent(out) :: nbasin0
       real*8, allocatable, intent(inout) :: xcoord(:,:)
       integer, allocatable, intent(inout) :: volnum0(:,:,:)
     end subroutine bader_integrate
     module subroutine refine_edge(f,irefine_edge,ref_itrs)
       use tools_io, only: faterr, ferror
       real*8, intent(in) :: f(:,:,:)
       integer, intent(inout) :: irefine_edge
       integer, intent(inout) :: ref_itrs
     end subroutine refine_edge
     module subroutine max_neargrid(f,p)
       use types, only: realloc
       real*8, intent(in) :: f(:,:,:)
       integer, dimension(3), intent(inout) :: p
     end subroutine max_neargrid
     module subroutine step_neargrid(f,p)
       real*8, intent(in) :: f(:,:,:)
       integer, intent(inout) :: p(3)
     end subroutine step_neargrid
     module subroutine step_ongrid(f,p)
       real*8, intent(in) :: f(:,:,:)
       integer, intent(inout) :: p(3)
     end subroutine step_ongrid
     module function rho_grad_dir(f,p) result(res)
       real*8, intent(in) :: f(:,:,:)
       integer, intent(in) :: p(3)
       real*8 :: res(3)
     end function rho_grad_dir
     module function is_max(f,p)
       real*8, intent(in) :: f(:,:,:)
       integer, intent(in) :: p(3)
       logical :: is_max
     end function is_max
     module subroutine pbc(p)
       integer, intent(inout) :: p(3)
     end subroutine pbc
     module function rho_val(ff,p1,p2,p3)
       real*8, intent(in) :: ff(:,:,:)
       integer, intent(in) :: p1, p2, p3
       real*8 :: rho_val
     end function rho_val
     module function volnum_val(p1,p2,p3)
       integer, intent(in) :: p1, p2, p3
       integer :: volnum_val
     end function volnum_val
     module subroutine assign_surrounding_pts(p)
       integer, intent(in) :: p(3)
       integer :: pt(3)
     end subroutine assign_surrounding_pts
     module subroutine known_volnum_ongrid(p)
       integer, intent(in) :: p(3)
     end subroutine known_volnum_ongrid
     module function is_vol_edge(p)
       logical :: is_vol_edge
       integer, intent(in) :: p(3)
     end function is_vol_edge
     module subroutine reassign_volnum_ongrid2(p)
       integer, intent(in) :: p(3)
     end subroutine reassign_volnum_ongrid2
  end interface

end module bader
