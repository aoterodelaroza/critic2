! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Driver routines for structural calculations
module struct_drivers
  use systemmod, only: system
  use crystalmod, only: crystal
  use crystalseedmod, only: crystalseed
  implicit none

  private
  public :: struct_crystal_input
  public :: struct_sym
  public :: struct_charges
  public :: struct_write
  public :: struct_atomlabel
  public :: struct_powder
  public :: struct_rdf
  public :: struct_amd
  public :: struct_compare
  public :: struct_comparevc
  public :: struct_environ
  public :: struct_econ
  public :: struct_coord
  public :: struct_polyhedra
  public :: struct_packing
  public :: struct_vdw
  public :: struct_edit
  public :: struct_newcell
  public :: struct_vibrations
  public :: struct_molcell
  public :: struct_identify
  public :: struct_makemols_neighcrys
  public :: struct_molreorder
  public :: struct_molmove
  public :: struct_kpoints
  public :: struct_bz

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
     module subroutine struct_sym(s,line,verbose)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
       logical, intent(in) :: verbose
     end subroutine struct_sym
     module subroutine struct_charges(s,line,oksyn)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
       logical, intent(out) :: oksyn
     end subroutine struct_charges
     module subroutine struct_write(s,line,usexyznames)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
       logical, intent(in) :: usexyznames
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
       type(system), intent(inout), target :: s
       character*(*), intent(in) :: line
     end subroutine struct_rdf
     module subroutine struct_amd(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_amd
     module subroutine struct_compare(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_compare
     module subroutine struct_comparevc(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_comparevc
     module subroutine struct_environ(s,line)
       type(system), intent(inout), target :: s
       character*(*), intent(in) :: line
     end subroutine struct_environ
     module subroutine struct_econ(s)
       type(system), intent(inout) :: s
     end subroutine struct_econ
     module subroutine struct_coord(s,line)
       type(system), intent(inout), target :: s
       character*(*), intent(in) :: line
     end subroutine struct_coord
     module subroutine struct_polyhedra(s,line)
       type(system), intent(inout), target :: s
       character*(*), intent(in) :: line
     end subroutine struct_polyhedra
     module subroutine struct_packing(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_packing
     module subroutine struct_vdw(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_vdw
     module subroutine struct_edit(s,verbose)
       type(system), intent(inout) :: s
       logical, intent(in) :: verbose
     end subroutine struct_edit
     module subroutine struct_newcell(s,line,verbose)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
       logical, intent(in) :: verbose
     end subroutine struct_newcell
     module subroutine struct_vibrations(s,line,verbose)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
       logical, intent(in) :: verbose
     end subroutine struct_vibrations
     module subroutine struct_molcell(s,line)
       type(system), intent(inout) :: s
       character*(*), intent(in) :: line
     end subroutine struct_molcell
     module subroutine struct_identify(s,line0,lp)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line0
       integer, intent(inout) :: lp
     end subroutine struct_identify
     module subroutine struct_makemols_neighcrys(line0,lp)
       character*(*), intent(in) :: line0
       integer, intent(inout) :: lp
     end subroutine struct_makemols_neighcrys
     module subroutine struct_molreorder(line,lp)
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp
     end subroutine struct_molreorder
     module subroutine struct_molmove(line,lp)
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp
     end subroutine struct_molmove
     module subroutine struct_kpoints(s,line)
       type(system), intent(in) :: s
       character*(*), intent(in) :: line
     end subroutine struct_kpoints
     module subroutine struct_bz(s)
       type(system), intent(in) :: s
     end subroutine struct_bz
  end interface

end module struct_drivers
