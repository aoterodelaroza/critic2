! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! Driver routines for structure operations
module struct_drivers
  use systemmod, only: system
  use crystalmod, only: crystal
  use crystalseedmod, only: crystalseed
  implicit none

  private
  public :: struct_crystal_input
  public :: struct_clearsym
  public :: struct_charges
  public :: struct_write
  public :: struct_atomlabel
  public :: struct_powder
  public :: struct_rdf
  public :: struct_compare
  public :: struct_environ
  public :: struct_coord
  public :: struct_packing
  public :: struct_newcell
  public :: struct_molcell
  public :: struct_identify
  
  interface
     module subroutine struct_crystal_input(line,mol0,allownofile,verbose,s0,cr0,seed0)
       character*(*), intent(in) :: line
       integer, intent(in) :: mol0
       logical, intent(in) :: allownofile
       logical, intent(in) :: verbose
       type(system), intent(inout), optional :: s0
       type(crystal), intent(inout), optional :: cr0
       type(crystalseed), intent(inout), optional :: seed0
     end subroutine struct_crystal_input
     module subroutine struct_clearsym(s)
       type(system), intent(inout) :: s
     end subroutine struct_clearsym
     module subroutine struct_charges(s,line,oksyn)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
       logical, intent(out) :: oksyn
     end subroutine struct_charges
     module subroutine struct_write(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_write
     module subroutine struct_atomlabel(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_atomlabel
     module subroutine struct_powder(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_powder
     module subroutine struct_rdf(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_rdf
     module subroutine struct_compare(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_compare
     module subroutine struct_environ(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_environ
     module subroutine struct_coord(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_coord
     module subroutine struct_packing(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_packing
     module subroutine struct_newcell(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_newcell
     module subroutine struct_molcell(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_molcell
     module subroutine struct_identify(s,line0,lp)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line0
       integer, intent(inout) :: lp
     end subroutine struct_identify
  end interface

end module struct_drivers
