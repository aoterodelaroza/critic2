!> Abinit structure and scalar field reader.
!> Most routines in this module were adapted from abinit.
! COPYRIGHT
! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MKV, MT)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .
! For the initials of contributors, see ~abinit/doc/developers/contributor

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

! abinit tools and readers
module abinit_private
  implicit none

  private
  public :: hdr_type, hdr_io ! abinit

  ! The hdr_type to be exported
  type hdr_type
     real*8 :: rprimd(3,3)
     integer :: natom
     integer :: ntypat
     real*8, allocatable :: xred(:,:)
     real*8, allocatable :: znucltypat(:)
     integer, allocatable :: typat(:)
     integer :: ngfft(3)
  end type hdr_type

  interface
     module subroutine hdr_io(fform,hdr,rdwr,unitfi,errmsg)
       integer, intent(inout) :: fform
       integer, intent(in) :: rdwr,unitfi
       type(hdr_type), intent(inout) :: hdr
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine hdr_io
  end interface

end module abinit_private
