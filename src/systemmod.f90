! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! system class and associated routines
module systemmod
  use iso_c_binding, only: c_ptr
  use grid1mod, only: grid1
  use hashmod, only: hash
  use types, only: integrable, pointpropable, discard_cp_expr
  use fieldmod, only: field
  use crystalmod, only: crystal
  use types, only: thread_info
  implicit none

  private

  public :: systemmod_init
  public :: systemmod_end

  ! The keywords of POINTPROP that define a point property of the
  ! reference field by name (besides STRESS, a tensor), and what each
  ! one is
  integer, parameter, public :: npointprop_keywords = 15
  character(len=7), parameter, public :: pointprop_keywords(npointprop_keywords) = (/&
     character(len=7) :: "gtf","vtf","htf","gtf_kir","vtf_kir","htf_kir","gkin","kkin","lag",&
     "elf","vir","he","lol","lol_kir","rdg"/)
  character(len=96), parameter, public :: pointprop_keyword_desc(npointprop_keywords) = (/&
     character(len=96) :: &
     "Thomas-Fermi kinetic energy density",&
     "Potential energy density from the Thomas-Fermi kinetic energy density and the virial theorem",&
     "Total energy density from the Thomas-Fermi kinetic energy density and the virial theorem",&
     "Thomas-Fermi kinetic energy density with the Kirzhnits gradient correction",&
     "Potential energy density from gtf_kir and the local virial theorem",&
     "Total energy density from gtf_kir and the local virial theorem",&
     "Kinetic energy density, G (positive definite) version",&
     "Kinetic energy density, K (Schrodinger) version",&
     "Lagrangian density (-1/4 of the Laplacian)",&
     "Electron localization function (ELF)",&
     "Electronic potential energy density (virial field)",&
     "Energy density (G + V)",&
     "Localized-orbital locator (LOL)",&
     "Localized-orbital locator with the Kirzhnits kinetic energy density",&
     "Reduced density gradient"/)

  ! The system class. A system contains:
  ! - One crystal structure (%c)
  ! - One or more fields (%nf fields in %f(:))
  ! - A reference field with id %iref
  ! - Several integrable properties (%npropi and %propi(:))
  ! - Several point properties (%npropp and %propp(:))
  ! - A set of field aliases (%fh)
  ! Some of the methods in this object require pass a C pointer down
  ! the hierarchy over to the arithmetic module. If one of those methods needs
  ! to be used, the system itself must be declared as TARGET.
  type system
     logical :: isinit = .false. !< Is the system initialized?
     type(crystal), allocatable :: c !< Crystal structure (always allocated)
     integer :: nf = -1 !< Number of fields
     type(field), allocatable :: f(:) !< Fields for this system
     integer :: fieldgen = 0 !< Generation counter, bumped whenever the field set changes (loaded, copied, unloaded, reset)
     logical :: refset = .false. !< Has the reference been set?
     integer :: iref = 0 !< Reference field
     integer :: npropi = 0 !< Number of integrable properties
     type(integrable), allocatable :: propi(:) !< Integrable properties
     integer :: npropp = 0 !< Number of properties at points
     type(pointpropable), allocatable :: propp(:) !< Properties at points
     type(hash) :: fh !< Hash of function aliases
   contains
     procedure :: end => system_end !< Terminate a system object
     procedure :: init => system_init !< Allocate space for crystal structure
     procedure :: clearsym !< Clear symmetry in the system's structure and the CP list
     procedure :: reset_fields !< Reset fields, properties, and aliases to promolecular
     procedure :: set_reference !< Set a given field as reference
     procedure :: set_default_integprop !< Reset to default integrable properties
     procedure :: set_default_pointprop !< Reset to default point properties
     procedure :: report !< Write information about the system to the stdout
     procedure :: aliasstring !< A string containing the aliases of a given field
     procedure :: add_field !< Add a field to the system, given as an argument
     procedure :: new_from_seed !< Build a system from a crystal seed
     procedure :: load_field_string !< Load a field using a command string
     procedure :: goodfield !< Returns true if the field is initialized
     procedure :: fieldname_to_idx !< Find the field ID from the alias
     procedure :: propi_field !< Slot of the field of an integrable property (-1 if gone)
     procedure :: getfieldnum !< Find an open slot for a new field
     procedure :: field_copy !< Copy a field from one slot to another
     procedure :: unload_field !< Unload a field
     procedure :: new_integrable_string !< Define a field as integrable from a command
     procedure :: new_pointprop_string !< Define a field as point prop from a command
     procedure :: delete_pointprop !< Remove a point property from the list
     procedure :: eval => system_eval_expression !< Evaluate an arithmetic expression using the system's fields
     procedure :: check_expression => system_check_expression !< Validate an expression by evaluating it at a probe point
     procedure :: propty !< Calculate the properties of a field or all fields at a point
     procedure :: grdall !< Calculate all integrable properties at a point
     procedure :: addcp !< Add a critical point to a field's CP list, maybe with discarding expr
  end type system
  public :: system

  !> A reference to a field of a system: its slot, and the unique
  !> identifier of the field that was in the slot when the reference
  !> was set. The reference is good only while that same field is
  !> there, not after it is unloaded or another field takes the slot.
  type field_ref
     integer :: id = -1 !< slot of the field (s%f)
     integer*8 :: uid = 0 !< unique identifier of that field (0 = the slot was empty)
   contains
     procedure :: set => field_ref_set !< point to the field in a slot
     procedure :: ok => field_ref_ok !< the field pointed to is still there
  end type field_ref
  public :: field_ref

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
  integer, parameter, public :: itype_deloc_wnr = 9
  integer, parameter, public :: itype_deloc_psink = 10
  integer, parameter, public :: itype_deloc_sijchk = 11
  integer, parameter, public :: itype_deloc_fachk = 12
  integer, parameter, public :: itype_hirshfeld_ovpop = 13
  character*10, parameter, public :: itype_names(13) = (/&
     "Volume    ","Field     ","Field (v) ","Gradnt mod","Laplacian ",&
     "Laplcn (v)","Expression","Multipoles","Deloc indx","Deloc indx",&
     "Deloc indx","Deloc indx","Overlp pop"/)

  interface
     module subroutine systemmod_init(isy)
       integer, intent(in) :: isy
     end subroutine systemmod_init
     module subroutine systemmod_end()
     end subroutine systemmod_end
     module subroutine clearsym(s)
       class(system), intent(inout) :: s
     end subroutine clearsym
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
     module subroutine add_field(s,sptr,f,verbose,id,errmsg)
       class(system), intent(inout) :: s
       type(c_ptr), intent(in) :: sptr
       type(field), intent(in) :: f
       logical, intent(in) :: verbose
       integer, intent(out) :: id
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine add_field
     module subroutine new_from_seed(s,seed,errmsg,ti)
       use crystalseedmod, only: crystalseed
       class(system), intent(inout) :: s
       type(crystalseed), intent(in) :: seed
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
     end subroutine new_from_seed
     module subroutine load_field_string(s,line,verbose,id,errmsg,ti,readchk,autointerp)
       class(system), intent(inout), target :: s
       character*(*), intent(in) :: line
       logical, intent(in) :: verbose
       integer, intent(out) :: id
       character(len=:), allocatable, intent(out) :: errmsg
       type(thread_info), intent(in), optional :: ti
       logical, intent(in), optional :: readchk
       logical, intent(in), optional :: autointerp
     end subroutine load_field_string
     module function goodfield(s,id,key,type,n,idout,uid) result(ok)
       use fieldmod, only: type_grid
       use tools_io, only: ferror, faterr
       class(system), intent(in) :: s
       integer, intent(in), optional :: id
       character*(*), intent(in), optional :: key
       integer, intent(in), optional :: type
       integer, intent(in), optional :: n(3)
       integer, intent(out), optional :: idout
       integer*8, intent(in), optional :: uid
       logical :: ok
     end function goodfield
     module function fieldname_to_idx(s,id) result(fid)
       class(system), intent(in) :: s
       character*(*), intent(in) :: id
       integer :: fid
     end function fieldname_to_idx
     module function propi_field(s,i) result(fid)
       class(system), intent(in) :: s
       integer, intent(in) :: i
       integer :: fid
     end function propi_field
     module function getfieldnum(s) result(id)
       use fieldmod, only: realloc_field
       use tools_io, only: string
       class(system), intent(inout) :: s
       integer :: id
     end function getfieldnum
     module subroutine field_copy(s,id0,id1,keepuid)
       use fieldmod, only: realloc_field
       use tools_io, only: string
       class(system), intent(inout) :: s
       integer, intent(in) :: id0
       integer, intent(in) :: id1
       logical, intent(in), optional :: keepuid
     end subroutine field_copy
     module subroutine field_ref_set(fr,s,id)
       class(field_ref), intent(inout) :: fr
       type(system), intent(in) :: s
       integer, intent(in) :: id
     end subroutine field_ref_set
     module function field_ref_ok(fr,s) result(ok)
       class(field_ref), intent(in) :: fr
       type(system), intent(in) :: s
       logical :: ok
     end function field_ref_ok
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
       class(system), intent(inout), target :: s
       character*(*), intent(in) :: line0
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine new_pointprop_string
     module subroutine delete_pointprop(s,i)
       class(system), intent(inout) :: s
       integer, intent(in) :: i
     end subroutine delete_pointprop
     module function system_eval_expression(s,expr,errmsg,x0,toklist)
       use arithmetic, only: token
       class(system), intent(inout), target :: s
       character(*), intent(in) :: expr
       character(len=:), allocatable, intent(inout) :: errmsg
       real*8, intent(in), optional :: x0(3)
       type(token), intent(in), optional :: toklist(:)
       real*8 :: system_eval_expression
     end function system_eval_expression
     module subroutine system_check_expression(s,expr,errmsg)
       class(system), intent(inout), target :: s
       character(*), intent(in) :: expr
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine system_check_expression
     module subroutine propty(s,id,x0,res,resinput,verbose,allfields)
       use types, only: scalar_value
       class(system), intent(inout) :: s
       integer, intent(in) :: id
       real*8, dimension(:), intent(in) :: x0
       type(scalar_value), intent(inout) :: res
       logical, intent(in) :: resinput
       logical, intent(in) :: verbose
       logical, intent(in) :: allfields
     end subroutine propty
     module subroutine grdall(s,xpos,lprop,pmask)
       class(system), intent(inout) :: s
       real*8, intent(in) :: xpos(3)
       real*8, intent(out) :: lprop(s%npropi)
       logical, intent(in), optional :: pmask(s%npropi)
     end subroutine grdall
     module subroutine addcp(s,id,x0,discard,cpeps,nuceps,nucepsh,gfnormeps,itype,typeok)
       class(system), intent(inout) :: s
       integer, intent(in) :: id
       real*8, intent(in) :: x0(3)
       type(discard_cp_expr), allocatable, intent(in) :: discard(:)
       real*8, intent(in) :: cpeps
       real*8, intent(in) :: nuceps
       real*8, intent(in) :: nucepsh
       real*8, intent(in) :: gfnormeps
       integer, intent(in), optional :: itype
       logical, intent(in), optional :: typeok(4)
     end subroutine addcp
  end interface

end module systemmod
