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

! sysmod: system class and associated routines
module systemmod
  use hashmod, only: hash
  use types, only: integrable, pointpropable
  use fieldmod, only: field
  use crystalmod, only: crystal
  implicit none

  private

  public :: systemmod_init
  public :: systemmod_end
  private :: field_fcheck
  private :: field_feval
  private :: field_cube

  type system
     logical :: isinit = .false. !< Is the system initialized?
     type(crystal), pointer :: c => null() !< Crystal structure (always allocated)
     integer :: nf = -1 !< Number of fields
     type(field), allocatable :: f(:) !< Fields for this system
     logical :: refset = .false. !< Has the reference been set?
     integer :: iref = 0 !< Reference field
     integer :: npropi = 0 !< Number of integrable properties
     type(integrable), allocatable :: propi(:) !< Integrable properties
     integer :: npropp = 0 !< Number of properties at points
     type(pointpropable), allocatable :: propp(:) !< Properties at points
     type(hash) :: fh !< Hash of function aliases
     procedure(field_fcheck), nopass, pointer :: fcheck => null() !< Pointer to functions for arithmetics
     procedure(field_feval), nopass, pointer :: feval => null() !< Pointer to functions for arithmetics
     procedure(field_cube), nopass, pointer :: cube => null() !< Pointer to functions for arithmetics
   contains
     procedure :: end => system_end !< Terminate a system object
     procedure :: init => system_init !< Allocate space for crystal structure
     procedure :: reset_fields !< Reset fields, properties, and aliases to promolecular
     procedure :: set_reference !< Set a given field as reference
     procedure :: set_default_integprop !< Reset to default integrable properties
     procedure :: set_default_pointprop !< Reset to default point properties
     procedure :: report !< Write information about the system to the stdout
     procedure :: aliasstring !< A string containing the aliases of a given field
     procedure :: new_from_seed !< Build a system from a crystal seed
     procedure :: load_field_string !< Load a field using a command string
     procedure :: goodfield !< Returns true if the field is initialized
     procedure :: fieldname_to_idx !< Find the field ID from the alias
     procedure :: getfieldnum !< Find an open slot for a new field
     procedure :: field_copy !< Copy a field from one slot to another
     procedure :: unload_field !< Unload a field
     procedure :: new_integrable_string !< Define a field as integrable from a command
     procedure :: new_pointprop_string !< Define a field as point prop from a command
     procedure :: eval => system_eval !< Evaluate an arithmetic expression using the system's fields
     procedure :: propty !< Calculate the properties of a field or all fields at a point
     procedure :: grdall !< Calculate all integrable properties at a point
     procedure :: addcp !< Add a critical point to a field's CP list, maybe with discarding expr
  end type system
  public :: system

  ! Text-mode operation. Only one crystal and one system at a time.
  type(system), allocatable, target :: sy_(:)
  type(system), pointer :: sy => null()
  public :: sy

  ! integrable properties enumerate
  integer, parameter, public :: itype_v = 1
  integer, parameter, public :: itype_f = 2
  integer, parameter, public :: itype_fval = 3
  integer, parameter, public :: itype_gmod = 4
  integer, parameter, public :: itype_lap = 5
  integer, parameter, public :: itype_lapval = 6
  integer, parameter, public :: itype_expr = 7
  integer, parameter, public :: itype_mpoles = 8
  integer, parameter, public :: itype_deloc = 9
  character*10, parameter, public :: itype_names(9) = (/&
     "Volume    ","Field     ","Field (v) ","Gradnt mod","Laplacian ",&
     "Laplcn (v)","Expression","Multipoles","Deloc indx"/)
  
  interface
     module subroutine systemmod_init(isy)
       integer, intent(in) :: isy
     end subroutine systemmod_init
     module subroutine systemmod_end()
     end subroutine systemmod_end
     module subroutine system_end(s)
       class(system), intent(inout) :: s
     end subroutine system_end
     module subroutine system_init(s)
       class(system), intent(inout) :: s
     end subroutine system_init
     module subroutine reset_fields(s)
       class(system), intent(inout) :: s
     end subroutine reset_fields
     module subroutine set_reference(s,id,maybe)
       class(system), intent(inout) :: s
       integer, intent(in) :: id
       logical, intent(in) :: maybe
     end subroutine set_reference
     module subroutine set_default_integprop(s)
       class(system), intent(inout) :: s
     end subroutine set_default_integprop
     module subroutine set_default_pointprop(s)
       class(system), intent(inout) :: s
     end subroutine set_default_pointprop
     module subroutine report(s,lcrys,lfield,lpropi,lpropp,lalias,lzpsp,lcp)
       class(system), intent(inout) :: s
       logical, intent(in) :: lcrys
       logical, intent(in) :: lfield
       logical, intent(in) :: lpropi
       logical, intent(in) :: lpropp
       logical, intent(in) :: lalias
       logical, intent(in) :: lzpsp
       logical, intent(in) :: lcp
     end subroutine report
     module subroutine aliasstring(s,id,nal,str)
       class(system), intent(in) :: s
       integer, intent(in) :: id
       integer, intent(out) :: nal
       character(len=:), allocatable, intent(out) :: str
     end subroutine aliasstring
     module subroutine new_from_seed(s,seed)
       use crystalseedmod, only: crystalseed
       class(system), intent(inout) :: s
       type(crystalseed), intent(in) :: seed
     end subroutine new_from_seed
     module subroutine load_field_string(s,line,id,errmsg)
       class(system), intent(inout), target :: s
       character*(*), intent(in) :: line
       integer, intent(out) :: id
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine load_field_string
     module function goodfield(s,id,key) result(ok)
       class(system), intent(in) :: s
       integer, intent(in), optional :: id
       character*(*), intent(in), optional :: key
       logical :: ok
     end function goodfield
     module function fieldname_to_idx(s,id) result(fid)
       class(system), intent(in) :: s
       character*(*), intent(in) :: id
       integer :: fid
     end function fieldname_to_idx
     module function getfieldnum(s) result(id)
       use fieldmod, only: realloc_field
       use tools_io, only: string
       class(system), intent(inout) :: s
       integer :: id
     end function getfieldnum
     module subroutine field_copy(s,id0,id1)
       use fieldmod, only: realloc_field
       use tools_io, only: string
       class(system), intent(inout) :: s
       integer, intent(in) :: id0
       integer, intent(in) :: id1
     end subroutine field_copy
     module function field_fcheck(sptr,id,iout)
       use iso_c_binding, only: c_ptr
       type(c_ptr), intent(in) :: sptr
       character*(*), intent(in) :: id
       integer, intent(out), optional :: iout
       logical :: field_fcheck
     end function field_fcheck
     recursive module function field_feval(sptr,id,nder,fder,x0,periodic)
       use iso_c_binding, only: c_ptr
       use types, only: scalar_value
       type(scalar_value) :: field_feval
       type(c_ptr), intent(in) :: sptr
       character*(*), intent(in) :: id
       integer, intent(in) :: nder
       character*(*), intent(in) :: fder
       real*8, intent(in) :: x0(3)
       logical, intent(in), optional :: periodic
     end function field_feval
     module function field_cube(sptr,n,id,fder,dry,ifail) result(q)
       use iso_c_binding, only: c_ptr
       type(c_ptr), intent(in) :: sptr
       character*(*), intent(in) :: id
       integer, intent(in) :: n(3)
       character*(*), intent(in) :: fder
       logical, intent(in) :: dry
       logical, intent(out) :: ifail
       real*8 :: q(n(1),n(2),n(3))
     end function field_cube
     module subroutine unload_field(s,id)
       class(system), intent(inout) :: s
       integer, intent(in) :: id
     end subroutine unload_field
     module subroutine new_integrable_string(s,line,errmsg)
       class(system), intent(inout) :: s
       character*(*), intent(in) :: line
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine new_integrable_string
     module subroutine new_pointprop_string(s,line0,errmsg)
       class(system), intent(inout) :: s
       character*(*), intent(in) :: line0
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine new_pointprop_string
     module function system_eval(s,expr,hardfail,iok,x0) 
       class(system), intent(inout), target :: s
       character(*), intent(in) :: expr
       logical, intent(in) :: hardfail
       logical, intent(out) :: iok
       real*8, intent(in), optional :: x0(3)
       real*8 :: system_eval
     end function system_eval
     module subroutine propty(s,id,x0,res,verbose,allfields)
       use types, only: scalar_value
       class(system), intent(inout) :: s
       integer, intent(in) :: id
       real*8, dimension(:), intent(in) :: x0
       type(scalar_value), intent(out) :: res
       logical, intent(in) :: verbose
       logical, intent(in) :: allfields
     end subroutine propty
     module subroutine grdall(s,xpos,lprop,pmask)
       class(system), intent(inout) :: s
       real*8, intent(in) :: xpos(3)
       real*8, intent(out) :: lprop(s%npropi)
       logical, intent(in), optional :: pmask(s%npropi)
     end subroutine grdall
     module subroutine addcp(s,id,x0,discexpr,cpeps,nuceps,nucepsh,itype)
       class(system), intent(inout) :: s
       integer, intent(in) :: id
       real*8, intent(in) :: x0(3)
       character*(*), intent(in) :: discexpr
       real*8, intent(in) :: cpeps
       real*8, intent(in) :: nuceps
       real*8, intent(in) :: nucepsh
       integer, intent(in), optional :: itype
     end subroutine addcp
  end interface

end module systemmod
