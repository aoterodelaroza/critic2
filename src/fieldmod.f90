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

!> Field class
module fieldmod
  use fieldseedmod, only: fieldseed
  use crystalmod, only: crystal
  use fragmentmod, only: fragment
  use elk_private, only: elkwfn
  use wien_private, only: wienwfn
  use pi_private, only: piwfn
  use grid3mod, only: grid3
  use wfn_private, only: molwfn
  use dftb_private, only: dftbwfn
  use param, only: maxzat0, mlen, mmlen
  use types, only: cp_type, scalar_value, gpathp
  use hashmod, only: hash
  use iso_c_binding, only: c_ptr, c_null_ptr
  implicit none
  
  private

  public :: realloc_field

  ! pointers for the arithmetic module
  interface
     !> Check that the id is a grid and is a sane field
     function fcheck(sptr,id,iout)
       import c_ptr
       logical :: fcheck
       type(c_ptr), intent(in) :: sptr
       character*(*), intent(in) :: id
       integer, intent(out), optional :: iout
     end function fcheck
     !> Evaluate the field at a point
     function feval(sptr,id,nder,fder,x0,periodic)
       import c_ptr, scalar_value
       type(scalar_value) :: feval
       type(c_ptr), intent(in) :: sptr
       character*(*), intent(in) :: id
       integer, intent(in) :: nder
       character*(*), intent(in) :: fder
       real*8, intent(in) :: x0(3)
       logical, intent(in), optional :: periodic
     end function feval
  end interface

  !> Scalar field types
  integer, parameter, public :: type_uninit = -1 !< uninitialized
  integer, parameter, public :: type_promol = 0 !< promolecular density
  integer, parameter, public :: type_grid = 1 !< grid format
  integer, parameter, public :: type_wien = 2 !< wien2k format
  integer, parameter, public :: type_elk  = 3 !< elk format
  integer, parameter, public :: type_pi   = 4 !< pi format
  integer, parameter, public :: type_wfn  = 6 !< molecular wavefunction format
  integer, parameter, public :: type_dftb = 7 !< DFTB+ wavefunction
  integer, parameter, public :: type_promol_frag = 8 !< promolecular density from a fragment
  integer, parameter, public :: type_ghost = 9 !< a ghost field

  !> Definition of the field class
  type field
     ! parent structure information
     type(crystal), pointer :: c => null() !< crsytal
     integer :: id !< field ID
     ! general information
     logical :: isinit = .false. !< is this field initialized?
     integer :: type = type_uninit !< field type
     logical :: usecore = .false. !< augment with core densities
     logical :: numerical = .false. !< numerical derivatives
     logical :: exact = .false. !< exact or approximate calc
     integer :: typnuc = -3 !< type of nuclei
     character(len=mlen) :: name = "" !< field name
     character(len=mlen) :: file = "" !< file name
     ! scalar field types
     type(elkwfn) :: elk
     type(wienwfn) :: wien
     type(piwfn) :: pi
     type(grid3) :: grid
     type(molwfn) :: wfn
     type(dftbwfn) :: dftb
     ! promolecular and core densities
     type(fragment) :: fr
     integer :: zpsp(maxzat0)
     ! ghost field
     character(len=mmlen) :: expr
     type(hash), pointer :: fh => null()
     type(c_ptr) :: sptr = c_null_ptr
     procedure(fcheck), pointer, nopass :: fcheck => null()
     procedure(feval), pointer, nopass :: feval => null()
     ! critical point list
     integer :: ncp = 0
     type(cp_type), allocatable :: cp(:)
     integer :: ncpcel = 0
     type(cp_type), allocatable :: cpcel(:)
   contains
     procedure :: end => field_end !< Deallocate data and uninitialize
     procedure :: set_default_options => field_set_default_options !< Sets field default options
     procedure :: set_options => field_set_options !< Set field options from a command string
     procedure :: field_new !< Creates a new field from a field seed.
     procedure :: load_promolecular !< Loads a promolecular density field
     procedure :: load_as_fftgrid !< Loads as a transformation of a 3d grid
     procedure :: load_ghost !< Loads a ghost field
     procedure :: grd !< Calculate field value and its derivatives at a point
     procedure :: grd0 !< Calculate only the field value at a given point
     procedure :: der1i !< Numerical first derivatives of the field
     procedure :: der2ii !< Numerical second derivatives (diagonal)
     procedure :: der2ij !< Numerical second derivatives (mixed)
     procedure :: typestring !< Return a string identifying the field type
     procedure :: printinfo !< Print field information to stdout
     procedure :: init_cplist !< Initialize the CP list
     procedure :: nearest_cp !< Given a point, find the nearest CP of a certain type
     procedure :: identify_cp !< Identify the CP given the position
     procedure :: testrmt !< Test for MT discontinuities
     procedure :: benchmark !< Test the speed of field evaluation
     procedure :: newton !< Newton-Raphson search for a CP
     procedure :: addcp !< Add a new CP to the CP list
     procedure :: sortcps !< Sort the CP list by field value
     procedure :: gradient !< Calculate a gradient path
  end type field 
  public :: field

  integer, parameter :: ndif_jmax = 10

  interface
     module subroutine realloc_field(a,nnew)
       type(field), intent(inout), allocatable :: a(:)
       integer, intent(in) :: nnew
     end subroutine realloc_field
     module subroutine field_end(f)
       class(field), intent(inout) :: f
     end subroutine field_end
     module subroutine field_set_default_options(ff)
       class(field), intent(inout) :: ff
     end subroutine field_set_default_options
     module subroutine field_set_options(ff,line,errmsg)
       class(field), intent(inout) :: ff
       character*(*), intent(in) :: line
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine field_set_options
     module subroutine field_new(f,seed,c,id,fh,sptr,fcheck,feval,cube,errmsg)
       class(field), intent(inout) :: f
       type(fieldseed), intent(in) :: seed 
       type(crystal), intent(in), target :: c
       integer, intent(in) :: id
       type(hash), intent(in) :: fh
       type(c_ptr), intent(in) :: sptr
       character(len=:), allocatable, intent(out) :: errmsg
       interface
          function fcheck(sptr,id,iout)
            import c_ptr
            logical :: fcheck
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(out), optional :: iout
          end function fcheck
          function feval(sptr,id,nder,fder,x0,periodic)
            import c_ptr, scalar_value
            type(scalar_value) :: feval
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(in) :: nder
            character*(*), intent(in) :: fder
            real*8, intent(in) :: x0(3)
            logical, intent(in), optional :: periodic
          end function feval
          subroutine cube(sptr,n,id,fder,dry,ifail,q)
            use iso_c_binding, only: c_ptr
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(in) :: n(3)
            character*(*), intent(in) :: fder
            logical, intent(in) :: dry
            logical, intent(out) :: ifail
            real*8, intent(out) :: q(n(1),n(2),n(3))
          end subroutine cube
       end interface
     end subroutine field_new
     module subroutine load_ghost(f,c,id,name,expr,sptr,fh,fcheck,feval)
       class(field), intent(inout) :: f
       type(crystal), intent(in), target :: c
       integer, intent(in) :: id
       character*(*), intent(in) :: name
       character*(*), intent(in) :: expr
       type(c_ptr), intent(in) :: sptr
       type(hash), intent(in), target :: fh 
       interface
          function fcheck(sptr,id,iout)
            import c_ptr
            logical :: fcheck
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(out), optional :: iout
          end function fcheck
          function feval(sptr,id,nder,fder,x0,periodic)
            import c_ptr, scalar_value
            type(scalar_value) :: feval
            type(c_ptr), intent(in) :: sptr
            character*(*), intent(in) :: id
            integer, intent(in) :: nder
            character*(*), intent(in) :: fder
            real*8, intent(in) :: x0(3)
            logical, intent(in), optional :: periodic
          end function feval
       end interface
     end subroutine load_ghost
     module subroutine load_promolecular(f,c,id,name,fr)
       class(field), intent(inout) :: f
       type(crystal), intent(in), target :: c
       integer, intent(in) :: id
       character*(*), intent(in) :: name
       type(fragment), intent(in), optional :: fr
     end subroutine load_promolecular
     module subroutine load_as_fftgrid(f,c,id,name,g,ityp,isry_)
       class(field), intent(inout) :: f
       type(crystal), intent(in), target :: c
       integer, intent(in) :: id
       character*(*), intent(in) :: name
       type(grid3), intent(in) :: g
       integer, intent(in) :: ityp
       logical, intent(in), optional :: isry_
     end subroutine load_as_fftgrid
     recursive module subroutine grd(f,v,nder,res,fder,periodic)
       class(field), intent(inout) :: f
       real*8, intent(in) :: v(3)
       integer, intent(in) :: nder
       type(scalar_value), intent(out) :: res
       character*(*), intent(in), optional :: fder
       logical, intent(in), optional :: periodic
     end subroutine grd
     recursive module function grd0(f,v,periodic)
       class(field), intent(inout) :: f
       real*8, dimension(3), intent(in) :: v
       real*8 :: grd0
       logical, intent(in), optional :: periodic
     end function grd0
     recursive module function der1i(f,dir,x,h,errcnv,pool,periodic)
       class(field), intent(inout) :: f
       real*8, intent(in) :: dir(3)
       real*8, intent(in) :: x(3), h, errcnv
       real*8, intent(inout) :: pool(-ndif_jmax:ndif_jmax)
       logical, intent(in), optional :: periodic
       real*8 :: der1i
     end function der1i
     recursive module function der2ii(f,dir,x,h,errcnv,pool,periodic)
       class(field), intent(inout) :: f
       real*8 :: der2ii
       real*8, intent(in) :: dir(3)
       real*8, intent(in) :: x(3), h, errcnv
       real*8, intent(inout) :: pool(-ndif_jmax:ndif_jmax)
       logical, intent(in), optional :: periodic
     end function der2ii
     recursive module function der2ij(f,dir1,dir2,x,h1,h2,errcnv,periodic)
       class(field), intent(inout) :: f
       real*8 :: der2ij
       real*8, intent(in) :: dir1(3), dir2(3)
       real*8, intent(in) :: x(3), h1, h2, errcnv
       logical, intent(in), optional :: periodic
     end function der2ij
     module function typestring(f,short) result(s)
       class(field), intent(in) :: f
       character(len=:), allocatable :: s
       logical, intent(in) :: short
     endfunction typestring
     module subroutine printinfo(f,isload,isset)
       class(field), intent(in) :: f
       logical, intent(in) :: isload
       logical, intent(in) :: isset
     end subroutine printinfo
     module subroutine init_cplist(f)
       class(field), intent(inout) :: f
     end subroutine init_cplist
     module subroutine nearest_cp(f,xp,nid,dist,lvec,type,nid0,id0,nozero)
       class(field), intent(in) :: f
       real*8, intent(in) :: xp(:)
       integer, intent(out) :: nid
       real*8, intent(out) :: dist
       integer, intent(out), optional :: lvec(3)
       integer, intent(in), optional :: type
       integer, intent(in), optional :: nid0
       integer, intent(in), optional :: id0
       logical, intent(in), optional :: nozero
     end subroutine nearest_cp
     module function identify_cp(f,x0,eps)
       class(field), intent(in) :: f
       integer :: identify_cp
       real*8, intent(in) :: x0(3)
       real*8, intent(in) :: eps
     end function identify_cp
     module subroutine testrmt(f,ilvl,errmsg)
       class(field), intent(inout) :: f
       integer, intent(in) :: ilvl
       character(len=:), allocatable, intent(out) :: errmsg
     end subroutine testrmt
     module subroutine benchmark(f,npts)
       class(field), intent(inout) :: f
       integer, intent(in) :: npts
     end subroutine benchmark
     module subroutine newton(f,r,gfnormeps,ier)
       class(field), intent(inout) :: f
       real*8, dimension(3), intent(inout) :: r
       integer, intent(out) :: ier
       real*8, intent(in) :: gfnormeps
     end subroutine newton
     module subroutine addcp(f,x0,cpeps,nuceps,nucepsh,itype)
       class(field), intent(inout) :: f
       real*8, intent(in) :: x0(3)
       real*8, intent(in) :: cpeps
       real*8, intent(in) :: nuceps
       real*8, intent(in) :: nucepsh
       integer, intent(in), optional :: itype
     end subroutine addcp
     module subroutine sortcps(f,cpeps)
       class(field), intent(inout) :: f
       real*8, intent(in) :: cpeps
     end subroutine sortcps
     module subroutine gradient(fid,xpoint,iup,nstep,ier,up2beta,plen,path,prune,pathini)
       class(field), intent(inout) :: fid
       real*8, dimension(3), intent(inout) :: xpoint
       integer, intent(in) :: iup
       integer, intent(out) :: nstep
       integer, intent(out) :: ier
       logical, intent(in) :: up2beta
       real*8, intent(out) :: plen
       type(gpathp), intent(inout), allocatable, optional :: path(:)
       real*8, intent(in), optional :: prune
       real*8, intent(in), optional :: pathini(3)
     end subroutine gradient
  end interface
  
end module fieldmod
