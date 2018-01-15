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
  use crystalmod, only: crystal
  use fragmentmod, only: fragment
  use elk_private, only: elkwfn
  use wien_private, only: wienwfn
  use pi_private, only: piwfn
  use grid3mod, only: grid3
  use wfn_private, only: molwfn
  use dftb_private, only: dftbwfn
  use param, only: maxzat0
  use types, only: cp_type, scalar_value, gpathp
  use hashmod, only: hash
  use iso_c_binding, only: c_ptr, c_null_ptr
  implicit none
  
  private

  public :: realloc_field
  private :: adaptive_stepper
  private :: stepper_euler1
  private :: stepper_heun
  private :: stepper_bs
  private :: stepper_rkck
  private :: stepper_dp

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
     character*(255) :: name = "" !< field name
     character*(255) :: file = "" !< file name
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
     character*(2048) :: expr
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

  ! eps to move to the main cell
  real*8, parameter :: flooreps = 1d-4 ! border around unit cell

  ! numerical differentiation parameters
  real*8, parameter :: derw = 1.4d0, derw2 = derw*derw, big = 1d30, safe = 2d0
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
     module subroutine nearest_cp(f,xp,nid,dist,type,idx,nozero)
       class(field), intent(in) :: f
       real*8, intent(in) :: xp(:)
       integer, intent(out) :: nid
       real*8, intent(out) :: dist
       integer, intent(in), optional :: type
       integer, intent(in), optional :: idx
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
     module function adaptive_stepper(fid,xpoint,h0,maxstep,eps,res)
       logical :: adaptive_stepper
       type(field), intent(inout) :: fid
       real*8, intent(inout) :: xpoint(3)
       real*8, intent(inout) :: h0
       real*8, intent(in) :: maxstep, eps
       type(scalar_value), intent(inout) :: res
     end function adaptive_stepper
     module subroutine stepper_euler1(xpoint,grdt,h0,xout)
       real*8, intent(in) :: xpoint(3), h0, grdt(3)
       real*8, intent(out) :: xout(3)
     end subroutine stepper_euler1
     module subroutine stepper_heun(fid,xpoint,grdt,h0,xout,res)
       type(field), intent(inout) :: fid
       real*8, intent(in) :: xpoint(3), h0, grdt(3)
       real*8, intent(out) :: xout(3)
       type(scalar_value), intent(inout) :: res
     end subroutine stepper_heun
     module subroutine stepper_bs(fid,xpoint,grdt,h0,xout,xerr,res)
       type(field), intent(inout) :: fid
       real*8, intent(in) :: xpoint(3), h0, grdt(3)
       real*8, intent(out) :: xout(3), xerr(3)
       type(scalar_value), intent(inout) :: res
     end subroutine stepper_bs
     module subroutine stepper_rkck(fid,xpoint,grdt,h0,xout,xerr,res)
       type(field), intent(inout) :: fid
       real*8, intent(in) :: xpoint(3), grdt(3), h0
       real*8, intent(out) :: xout(3), xerr(3)
       type(scalar_value), intent(inout) :: res
     end subroutine stepper_rkck
     module subroutine stepper_dp(fid,xpoint,grdt,h0,xout,xerr,res)
       type(field), intent(inout) :: fid
       real*8, intent(in) :: xpoint(3), grdt(3), h0
       real*8, intent(out) :: xout(3), xerr(3)
       type(scalar_value), intent(inout) :: res
     end subroutine stepper_dp
  end interface
  
contains

  !> Load a new field using the given field seed and the crystal
  !> structure pointer. The ID of the field in the system is also
  !> required.
  subroutine field_new(f,seed,c,id,fh,sptr,fcheck,feval,cube,errmsg)
    use types, only: realloc
    use fieldseedmod, only: fieldseed
    use arithmetic, only: eval
    use tools_io, only: equal, isinteger
    use param, only: ifformat_unknown, ifformat_wien, ifformat_elk, ifformat_pi,&
       ifformat_cube, ifformat_abinit, ifformat_vasp, ifformat_vaspchg, ifformat_qub,&
       ifformat_xsf, ifformat_elkgrid, ifformat_siestagrid, ifformat_dftb, ifformat_chk,&
       ifformat_wfn, ifformat_wfx, ifformat_fchk, ifformat_molden, ifformat_as,&
       ifformat_as_promolecular, ifformat_as_core, ifformat_as_lap, ifformat_as_grad,&
       ifformat_as_pot, ifformat_as_clm, ifformat_as_clm_sub, ifformat_as_ghost, &
       ifformat_copy, ifformat_promolecular, ifformat_promolecular_fragment
    use hashmod, only: hash
    use iso_c_binding, only: c_ptr
    class(field), intent(inout) :: f !< Input field
    type(fieldseed), intent(in) :: seed 
    type(crystal), intent(in), target :: c
    integer, intent(in) :: id
    type(hash), intent(in) :: fh
    type(c_ptr), intent(in) :: sptr
    character(len=:), allocatable, intent(out) :: errmsg

    interface
       !> Check that the id is a grid and is a sane field
       function fcheck(sptr,id,iout)
         use iso_c_binding, only: c_ptr
         logical :: fcheck
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(out), optional :: iout
       end function fcheck
       !> Evaluate the field at a point
       function feval(sptr,id,nder,fder,x0,periodic)
         use types, only: scalar_value
         use iso_c_binding, only: c_ptr
         type(scalar_value) :: feval
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(in) :: nder
         character*(*), intent(in) :: fder
         real*8, intent(in) :: x0(3)
         logical, intent(in), optional :: periodic
       end function feval
       !> Return the grid from field id
       function cube(sptr,n,id,fder,dry,ifail) result(q)
         use iso_c_binding, only: c_ptr
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(in) :: n(3)
         character*(*), intent(in) :: fder
         logical, intent(in) :: dry
         logical, intent(out) :: ifail
         real*8 :: q(n(1),n(2),n(3))
       end function cube
    end interface

    character(len=:), allocatable :: ofile
    integer :: i, j, k, iz, n(3), ithis
    type(fragment) :: fr
    real*8 :: xdelta(3,3), x(3), rho
    logical :: iok, found

    errmsg = ""
    if (.not.c%isinit) then
       errmsg = "crystal not initialized"
       return
    end if
    call f%end()
    f%c => c
    f%id = id
    f%name = adjustl(trim(seed%fid))

    ! set the default field flags
    call f%set_default_options()

    ! inherit the pseudopotential charges from the crystal
    f%zpsp = c%zpsp

    ! interpret the seed and load the field
    if (seed%iff == ifformat_unknown) then
       errmsg = "unknown seed format"
       call f%end()
       return

    elseif (seed%iff == ifformat_wien) then
       call f%wien%end()
       call f%wien%read_clmsum(seed%file(1),seed%file(2))
       f%type = type_wien
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_elk) then
       if (seed%nfile == 1) then
          call f%grid%end()
          call f%grid%read_elk(seed%file(1))
          f%type = type_grid
          f%file = seed%file(1)
       elseif (seed%nfile == 2) then
          call f%elk%end()
          call f%elk%read_out(seed%file(1),seed%file(2))
          f%type = type_elk
          f%file = seed%file(1)
       else
          call f%elk%end()
          call f%elk%read_out(seed%file(1),seed%file(2),seed%file(3))
          f%type = type_elk
          f%file = seed%file(3)
       endif

    elseif (seed%iff == ifformat_pi) then
       call f%pi%end()
       do i = 1, seed%nfile
          iok = isinteger(ithis,seed%piat(i))
          found = .false.
          do j = 1, c%nneq
             if (equal(seed%piat(i),c%spc(c%at(j)%is)%name)) then
                call f%pi%read_ion(seed%file(i),j)
                found = .true.
             else if (iok) then
                if (ithis == j) then
                   call f%pi%read_ion(seed%file(i),j)
                   found = .true.
                end if
             end if
          end do
          if (.not.found) then
             errmsg = "unknown atom for pi ion file: " // trim(seed%file(i))
             call f%end()
             return
          end if
       end do

       call f%pi%register_struct(f%c%nenv,f%c%spc,f%c%atenv(1:f%c%nenv))
       call f%pi%fillinterpol()
       f%type = type_pi
       f%file = "<pi ion files>"

    elseif (seed%iff == ifformat_cube) then
       call f%grid%end()
       call f%grid%read_cube(seed%file(1))
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_abinit) then
       call f%grid%end()
       call f%grid%read_abinit(seed%file(1))
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_vasp) then
       call f%grid%end()
       call f%grid%read_vasp(seed%file(1),f%c%omega)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_vaspchg) then
       call f%grid%end()
       call f%grid%read_vasp(seed%file(1),1d0)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_qub) then
       call f%grid%end()
       call f%grid%read_qub(seed%file(1))
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_xsf) then
       call f%grid%end()
       call f%grid%read_xsf(seed%file(1))
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_elkgrid) then
       call f%grid%end()
       call f%grid%read_elk(seed%file(1))
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_siestagrid) then
       call f%grid%end()
       call f%grid%read_siesta(seed%file(1))
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_dftb) then
       call f%dftb%end()
       call f%dftb%read(seed%file(1),seed%file(2),seed%file(3),f%c%atcel(1:f%c%ncel),f%c%spc(1:f%c%nspc))
       call f%dftb%register_struct(f%c%crys2car,f%c%atenv(1:f%c%nenv),f%c%spc(1:f%c%nspc))
       f%type = type_dftb
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_chk) then
       call f%grid%end()

       if (seed%nfile == 1) then
          ofile = ""
       else
          ofile = seed%file(2)
       end if
       f%grid%wan%useu = .not.seed%nou
       f%grid%wan%sijchk = seed%sijchk
       f%grid%wan%fachk = seed%fachk
       f%grid%wan%haschk = .false.
       f%grid%wan%cutoff = seed%wancut
       f%type = type_grid
       f%file = seed%file(1)

       if (len_trim(seed%unkgen) > 0 .and. len_trim(seed%evc) > 0) then
          call f%grid%read_unkgen(seed%file(1),ofile,seed%unkgen,seed%evc,&
             f%c%omega,seed%sijchk)
       else
          call f%grid%read_unk(seed%file(1),ofile,f%c%omega,seed%nou,&
             seed%sijchk)
       end if

    elseif (seed%iff == ifformat_wfn) then
       call f%wfn%end()
       call f%wfn%read_wfn(seed%file(1))
       call f%wfn%register_struct(f%c%ncel,f%c%atcel)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_wfx) then
       call f%wfn%end()
       call f%wfn%read_wfx(seed%file(1))
       call f%wfn%register_struct(f%c%ncel,f%c%atcel)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_fchk) then
       call f%wfn%end()
       call f%wfn%read_fchk(seed%file(1),seed%readvirtual)
       call f%wfn%register_struct(f%c%ncel,f%c%atcel)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_molden) then
       call f%wfn%end()
       call f%wfn%read_molden(seed%file(1),seed%readvirtual)
       call f%wfn%register_struct(f%c%ncel,f%c%atcel)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_promolecular) then
       call f%load_promolecular(f%c,id,"<promolecular>")

    elseif (seed%iff == ifformat_promolecular_fragment) then
       fr = f%c%identify_fragment_from_xyz(seed%file(1))
       if (fr%nat == 0) then
          errmsg = "fragment contains unknown atoms"
          call f%end()
          return
       end if
       call f%load_promolecular(f%c,id,trim(seed%file(1)),fr)

    elseif (seed%iff == ifformat_as_promolecular.or.seed%iff == ifformat_as_core) then
       if (seed%iff == ifformat_as_promolecular) then
          if (seed%nfile > 0) then
             fr = c%identify_fragment_from_xyz(seed%file(1))
             if (fr%nat == 0) then
                errmsg = "zero atoms in the fragment"
                call f%end()
                return
             end if
             call c%promolecular_grid(f%grid,seed%n,fr=fr)
          else
             call c%promolecular_grid(f%grid,seed%n)
          end if
       else
          call c%promolecular_grid(f%grid,seed%n,zpsp=c%zpsp)
       end if
       f%type = type_grid
       f%file = ""

    elseif (seed%iff == ifformat_as_ghost) then
       call f%load_ghost(c,id,"<ghost>",seed%expr,sptr,fh,fcheck,feval)

    elseif (seed%iff == ifformat_as) then
       f%type = type_grid
       f%file = ""
       n = seed%n
       call f%grid%new_eval(sptr,n,seed%expr,fh,cube)
       if (.not.f%grid%isinit) then
          call f%grid%end()
          f%grid%n = n
          allocate(f%grid%f(n(1),n(2),n(3)))

          do i = 1, 3
             xdelta(:,i) = 0d0
             xdelta(i,i) = 1d0 / real(n(i),8)
          end do

          !$omp parallel do private(x,rho) schedule(dynamic)
          do k = 1, n(3)
             do j = 1, n(2)
                do i = 1, n(1)
                   x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
                   x = c%x2c(x)
                   rho = eval(seed%expr,.true.,iok,x,sptr,fh,fcheck,feval,.true.)
                   !$omp critical(write)
                   f%grid%f(i,j,k) = rho
                   !$omp end critical(write)
                end do
             end do
          end do
          !$omp end parallel do
          f%grid%isinit = .true.
       end if

    elseif (seed%iff == ifformat_copy .or. seed%iff == ifformat_as_lap .or.&
       seed%iff == ifformat_as_pot .or. seed%iff == ifformat_as_grad .or. &
       seed%iff == ifformat_as_clm .or. seed%iff == ifformat_as_clm_sub) then
       errmsg = "error in file format for field_new"
       call f%end()
       return
    else
       errmsg = "unknown seed format"
       call f%end()
       return
    end if

    ! set the rest of the variables passed with the field
    call f%set_options(seed%elseopt,errmsg)
    if (len_trim(errmsg) > 0) then
       call f%end()
       return
    end if

    f%isinit = .true.
    call f%init_cplist()

  end subroutine field_new

  !> Load a ghost field.
  subroutine load_ghost(f,c,id,name,expr,sptr,fh,fcheck,feval)
    use grid3mod, only: grid3
    use fragmentmod, only: fragment
    use hashmod, only: hash
    use iso_c_binding, only: c_ptr
    class(field), intent(inout) :: f !< Input/output field
    type(crystal), intent(in), target :: c
    integer, intent(in) :: id
    character*(*), intent(in) :: name
    character*(*), intent(in) :: expr
    type(c_ptr), intent(in) :: sptr
    type(hash), intent(in), target :: fh 
    interface
       !> Check that the id is a grid and is a sane field
       function fcheck(sptr,id,iout)
         use iso_c_binding, only: c_ptr
         logical :: fcheck
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(out), optional :: iout
       end function fcheck
       !> Evaluate the field at a point
       function feval(sptr,id,nder,fder,x0,periodic)
         use types, only: scalar_value
         use iso_c_binding, only: c_ptr
         type(scalar_value) :: feval
         type(c_ptr), intent(in) :: sptr
         character*(*), intent(in) :: id
         integer, intent(in) :: nder
         character*(*), intent(in) :: fder
         real*8, intent(in) :: x0(3)
         logical, intent(in), optional :: periodic
       end function feval
    end interface

    if (.not.c%isinit) return
    call f%end()
    f%c => c
    f%id = id
    f%isinit = .true.
    f%type = type_ghost
    f%usecore = .false. 
    f%numerical = .true. 
    f%exact = .false. 
    f%name = adjustl(name)
    f%file = ""
    f%zpsp = c%zpsp
    f%expr = expr
    f%fh => fh
    f%sptr = sptr
    f%fcheck => fcheck
    f%feval => feval
    call f%init_cplist()

  end subroutine load_ghost

end module fieldmod
