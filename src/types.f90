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

!> User-defined types and overloaded reallocation procedures.
module types
  implicit none

  private
  public :: grid1, atom, celatom, anyatom, fragment
  public :: cp_type
  public :: scalar_value, field
  public :: integrable, pointpropable
  public :: minisurf, miniface
  public :: neighstar
  public :: realloc

  ! overloaded functions
  interface realloc
     module procedure realloc_field
     module procedure realloc_atom
     module procedure realloc_celatom
     module procedure realloc_anyatom
     module procedure realloc_fragment
     module procedure realloc_cp
     module procedure realloc1l
     module procedure realloc1r
     module procedure realloc2r
     module procedure realloc3r
     module procedure realloc4r
     module procedure realloc5r
     module procedure realloc1i
     module procedure realloc2i
     module procedure realloc1c
     module procedure realloc1cmplx4
     module procedure realloc2cmplx4
     module procedure realloc4cmplx4
     module procedure realloc1cmplx8
  end interface

  !> Radial grid type.
  type grid1
     logical :: init !< Is initialized?
     real*8 :: a !< Logarithmic grid parameter ri = a * exp(b * (i-1))
     real*8 :: b !< Logarithmic grid parameter ri = a * exp(b * (i-1))
     real*8 :: rmax !< Max. grid distance
     real*8 :: rmax2 !< Squared max. grid distance
     integer :: ngrid !< Number of nodes
     real*8, allocatable :: r(:) !< Node positions
     real*8, allocatable :: f(:) !< Grid values, f = 4*pi*r^2*rho
     real*8, allocatable :: fp(:) !< First derivative of f
     real*8, allocatable :: fpp(:) !< Second derivative of f 
     integer :: z 
     integer :: qat
  end type grid1

  !> Non-equivalent atom list type (nneq)
  type atom
     real*8 :: x(3)   !< coordinates (crystallographic)
     real*8 :: r(3)   !< coordinates (cartesian)
     character*(10) :: name = "" !< name
     integer :: z = 0 !< atomic number
     integer :: zpsp = -1 !< pseudopotential ionic charge (Z-core)
     integer :: qat = 0 !< ionic charge for promolecular densities and Ewald
     integer :: mult  !< multiplicity
     real*8 :: rnn2   !< half the nearest neighbor distance
  end type atom
  
  !> Equivalent atom list (type)
  type celatom
     real*8 :: x(3) !< coordinates (crystallographic)
     real*8 :: r(3) !< coordinates (cartesian)
     integer :: idx !< corresponding atom from the non-equivalent atom list
     integer :: cidx !< corresponding equivalent atom from the complete atom list
     integer :: ir  !< rotation matrix to the representative equivalent atom
     integer :: ic  !< translation vector to the representative equivalent atom
     integer :: lvec(3) !< lattice vector to the representative equivalent atom
  end type celatom

  !> Any atom in the crystal (not necessarily in the main cell), type
  type anyatom
     real*8 :: x(3) !< coordinates (crystallographic)
     real*8 :: r(3) !< coordinates (cartesian)
     integer :: idx !< corresponding atom from the non-equivalent atom list
     integer :: cidx !< corresponding atom from the complete atom list
     integer :: lvec(3) !< lattice vector to the atom in the complete atom list
     integer :: z !< atomic number
  end type anyatom

  !> Type for a fragment of the crystal
  type fragment
     integer :: nat !< Number of atoms in the fragment
     type(anyatom), allocatable :: at(:) !< Atoms in the fragment
  end type fragment

  !> Critical point type
  type cp_type
     ! Position and type
     real*8 :: x(3) !< Position (cryst. coords.)
     real*8 :: r(3) !< Position (cart. coords.)
     integer :: typ  !< Type of CP (-3,-1,1,3)
     integer :: typind !< Same as type, but mapped to (0,1,2,3)
     integer :: mult !< Multiplicity
     character*3 :: pg !< Point group symbol
     character*(10) :: name !< Name
     logical :: isdeg !< Is it a degenerate CP?
     logical :: isnuc !< Is it a nucleus?
     logical :: isnnm !< Is it a non-nuclear attractor/repulsor?

     ! beta-sphere 
     real*8 :: rbeta !< beta-sphere radius

     ! Properties at the CP
     real*8 :: rho  !< Density
     real*8 :: gmod !< Norm of the gradient
     real*8 :: lap  !< Laplacian
     logical :: deferred = .true. !< Calculate the rho and lap values later

     ! BCP and RCP ias properties
     integer :: ipath(2) !< Associated attractor (bcp) or repulsor (rcp), complete list
     integer :: ilvec(3,2) !< Lattice vector to shift the cp_(ipath) position of the actual attractor
     real*8 :: brdist(2) !< If b or r, distance to attractor/repulsor
     real*8 :: brang !< If b or r, angle wrt attractors/repulsors

     ! In shell-structured scalar fields,
     integer :: innuc !< associated nucleus 
     integer :: inshell !< associated shell 
     real*8 :: dist2nuc !< distance to nucleus

     ! Complete -> reduced CP list index and conversion
     integer :: idx !< Complete to non-equivalent list index
     integer :: ir !< Rotation matrix to the neq cp list
     integer :: ic !< Translation vector to the neq cp list 
     integer :: lvec(3) !< Lattice vector to the neq cp list 
  end type cp_type

  !> Information about a field
  type field
     ! all types/more than one type
     logical :: init = .false. !< is this field initialized?
     integer :: type !< field type
     logical :: usecore = .false. !< augment with core densities
     logical :: numerical = .false. !< numerical derivatives
     integer :: n(3) !< grid points for elk and grids
     integer :: typnuc !< type of nuclei
     character*(255) :: name = "" !< field name
     character*(255) :: file = "" !< file name
     ! grids
     integer :: mode
     real*8, allocatable :: f(:,:,:)
     real*8, allocatable, dimension(:,:,:,:) :: c2
     real*8 :: c2x(3,3), x2c(3,3)
     integer :: nwan(3) ! number of wannier vectors
     real*8, allocatable :: fwan(:,:,:,:,:) ! from wannier xsf
     ! wien2k 
     logical :: cnorm
     integer, allocatable :: lm(:,:,:)
     integer, allocatable :: lmmax(:)
     real*8, allocatable :: slm(:,:,:)
     real*8, allocatable :: sk(:)
     real*8, allocatable :: ski(:)
     real*8, allocatable :: tauk(:)
     real*8, allocatable :: tauki(:)
     integer :: lastind, lastniz
     ! wien2k, structural
     logical :: ishlat
     integer :: nat
     real*8, dimension(:,:,:), allocatable  :: rotloc
     real*8, dimension(:), allocatable :: rnot 
     real*8, dimension(:), allocatable :: rmt 
     real*8, dimension(:), allocatable :: dx
     integer, dimension(:), allocatable :: jri
     integer, dimension(:), allocatable :: multw
     real*8, dimension(3,3) :: br1
     real*8, dimension(3,3) :: br2
     real*8, dimension(3,3) :: br3
     logical :: ortho
     integer :: ndat
     integer, dimension(:), allocatable :: iatnr
     real*8, dimension(:,:), allocatable :: pos 
     integer, dimension(:), allocatable :: iop 
     integer :: niord
     integer, dimension(:,:,:), allocatable :: iz 
     real*8, dimension(:,:), allocatable :: tau
     integer :: npos 
     real*8  :: atp(3,343)
     integer :: nwav
     real*8, dimension(:,:), allocatable :: krec
     integer, dimension(:), allocatable :: inst
     logical :: cmpl
     real*8 :: a(3)
     ! elk
     real*8, allocatable :: rhomt(:,:,:)
     complex*16, allocatable :: rhok(:)
     ! elk, structural
     ! --> note rmt(:) in wien_structural
     integer :: lmaxvr
     real*8, allocatable :: spr(:,:)
     real*8, allocatable :: spr_a(:)
     real*8, allocatable :: spr_b(:)
     integer, allocatable :: nrmt(:)
     integer :: ngvec
     real*8, allocatable :: vgc(:,:) 
     integer, allocatable :: igfft(:)
     integer :: ncel0
     real*8, allocatable :: xcel(:,:)
     integer, allocatable :: iesp(:)
     ! pi
     logical :: exact
     character*(20), allocatable :: piname(:)
     logical, allocatable :: pi_used(:)
     integer, allocatable :: nsym(:), naos(:,:), naaos(:,:), nsto(:,:)
     integer, allocatable :: nasto(:,:), nn(:,:)
     real*8, allocatable :: z(:,:), xnsto(:,:), c(:,:,:), nelec(:,:)
     type(grid1), allocatable :: pgrid(:)
     ! wfn
     logical :: useecp
     integer :: nmo, nalpha
     integer :: npri
     integer :: wfntyp
     integer, allocatable :: icenter(:)
     integer, allocatable :: itype(:)
     real*8, allocatable :: e(:)
     real*8, allocatable :: occ(:)
     real*8, allocatable :: cmo(:,:)
     ! promolecular from fragment
     type(fragment) :: fr
     ! ghost field
     character*(2048) :: expr
     integer :: nf 
     integer, allocatable :: fused(:)
  end type field

  !> Result of the evaluation of a scalar field
  type scalar_value
     ! basic
     real*8 :: f, fval, gf(3), hf(3,3), gfmod, del2f, del2fval
     ! kinetic energy density
     real*8 :: gkin
     ! schrodinger stress tensor
     real*8 :: stress(3,3)
     ! electronic potential energy density, virial field
     real*8 :: vir
     ! additional local properties
     real*8 :: hfevec(3,3), hfeval(3)
     integer :: r, s
     ! is it a nuclear position?
     logical :: isnuc
  end type scalar_value

  !> Information about an integrable field
  type integrable
     logical :: used
     integer :: itype
     integer :: fid
     character*(10) :: prop_name
     character*(2048) :: expr
     integer :: lmax
  end type integrable

  !> Information about an point-property field
  type pointpropable
     character*(10) :: name
     character*(2048) :: expr
     integer :: ispecial ! 0 = expression, 1 = stress
     integer :: nf 
     integer, allocatable :: fused(:)
  end type pointpropable

  !> A face on the surface
  type miniface
     integer, dimension(4) :: v !< Pointers to vertices
     integer :: nv !< Number of vertices
  end type miniface

  !> Surface type
  type minisurf 
     integer :: init !< was initialized? 0 no, 1 only v, 2 all
     integer :: rgb(3) !< optionally used for global color
     real*8, dimension(3) :: n  !< center on which the polyhedron is based
     real*8, dimension(:), allocatable :: r, th, ph !< vertices
     type(miniface), allocatable :: f(:)
     integer :: nv !< nv
     integer :: nf !< number of vertex and faces
     integer :: mv !< mv
     integer :: mf !< maximum number of vertex and faces
  end type minisurf

  !> Neighbor star -> to analyze the connectivity in a molecular crystal
  type neighstar
     integer :: ncon = 0 !< number of neighbor for this atom
     integer, allocatable :: idcon(:) !< id (atcel) of the connected atom
     integer, allocatable :: lcon(:,:) !< lattice vector of the connected atom
  end type neighstar

contains
  
  !> Adapt the size of an allocatable 1D type(field) array
  subroutine realloc_field(a,nnew)
    use tools_io

    type(field), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(field), allocatable :: temp(:)
    integer :: nold, i, l1, u1

    if (.not.allocated(a)) &
       call ferror('realloc_field','array not allocated',faterr)
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    allocate(temp(l1:nnew))

    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc_field

  !> Adapt the size of an allocatable 1D type(atom) array
  subroutine realloc_atom(a,nnew)
    use tools_io

    type(atom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(atom), allocatable :: temp(:)
    integer :: nold, i

    if (.not.allocated(a)) &
       call ferror('realloc_atom','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)
    do i = nold+1,nnew
       a(i)%name = ""
       a(i)%z = 0
       a(i)%zpsp = -1
       a(i)%qat = 0
       a(i)%rnn2 = 0d0
    end do

  end subroutine realloc_atom

  !> Adapt the size of an allocatable 1D type(celatom) array
  subroutine realloc_celatom(a,nnew)
    use tools_io

    type(celatom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(celatom), allocatable :: temp(:)
    integer :: nold

    if (.not.allocated(a)) &
       call ferror('realloc_celatom','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_celatom

  !> Adapt the size of an allocatable 1D type(celatom) array
  subroutine realloc_anyatom(a,nnew)
    use tools_io

    type(anyatom), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(anyatom), allocatable :: temp(:)
    integer :: nold

    if (.not.allocated(a)) &
       call ferror('realloc_anyatom','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_anyatom

  !> Adapt the size of an allocatable 1D type(fragment) array
  subroutine realloc_fragment(a,nnew)
    use tools_io

    type(fragment), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(fragment), allocatable :: temp(:)
    integer :: nold

    if (.not.allocated(a)) &
       call ferror('realloc_fragment','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_fragment

  !> Adapt the size of an allocatable 1D type(atom) array
  subroutine realloc_cp(a,nnew)
    use tools_io

    type(cp_type), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew
    
    type(cp_type), allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc_atom','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_cp

  !> Adapt the size of an allocatable 1D logical array
  subroutine realloc1l(a,nnew)
    use tools_io, only: ferror, faterr

    logical, intent(inout), allocatable :: a(:) !< Input array, logical, 1D
    integer, intent(in) :: nnew !< new dimension
    
    logical, allocatable :: temp(:)
    integer :: l1, u1
    
    if (.not.allocated(a)) &
       call ferror('realloc1l','array not allocated',faterr)
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    allocate(temp(l1:nnew))
    
    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc1l

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc1r(a,nnew)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    real*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1r','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1r

  !> Adapt the size of an allocatable 2D real*8 array
  subroutine realloc2r(a,n1,n2)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    real*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call ferror('realloc2r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    allocate(temp(n1,n2))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc3r(a,n1,n2,n3)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3 !< new dimension
    
    real*8, allocatable :: temp(:,:,:)
    integer :: nold(3)
    
    if (.not.allocated(a)) &
       call ferror('realloc3r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    allocate(temp(n1,n2,n3))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)))
    call move_alloc(temp,a)

  end subroutine realloc3r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc4r(a,n1,n2,n3,n4)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3, n4 !< new dimension
    
    real*8, allocatable :: temp(:,:,:,:)
    integer :: nold(4)
    
    if (.not.allocated(a)) &
       call ferror('realloc4r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    nold(4) = size(a,4)
    allocate(temp(n1,n2,n3,n4))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)))
    call move_alloc(temp,a)

  end subroutine realloc4r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc5r(a,n1,n2,n3,n4,n5)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:,:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3, n4, n5 !< new dimension
    
    real*8, allocatable :: temp(:,:,:,:,:)
    integer :: nold(5)
    
    if (.not.allocated(a)) &
       call ferror('realloc5r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    nold(4) = size(a,4)
    nold(5) = size(a,5)
    allocate(temp(n1,n2,n3,n4,n5))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5)))
    call move_alloc(temp,a)

  end subroutine realloc5r

  !> Adapt the size of an allocatable 1D integer array
  subroutine realloc1i(a,nnew)
    use tools_io, only: ferror, faterr

    integer, intent(inout), allocatable :: a(:) !< Input array, integer, 1D
    integer, intent(in) :: nnew !< New dimension
    
    integer, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1i','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1i

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc2i(a,n1,n2)
    use tools_io, only: ferror, faterr

    integer, intent(inout), allocatable :: a(:,:) !< Input array, integer, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    integer, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call ferror('realloc2i','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    allocate(temp(n1,n2))
    
    temp = 0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2i

  !> Adapt the size of an allocatable 1D character array
  subroutine realloc1c(a,nnew)
    use tools_io, only: ferror, faterr

    character*(*), intent(inout), allocatable :: a(:) !< Input array, character, 1D
    integer, intent(in) :: nnew !< New dimension
    
    character*(len(a)), allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1c','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1c

  !> Adapt the size of an allocatable 1D complex*8 array
  subroutine realloc1cmplx4(a,nnew)
    use tools_io, only: ferror, faterr

    complex*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    complex*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1cmplx4','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1cmplx4

  !> Adapt the size of an allocatable 2D complex*8 array
  subroutine realloc2cmplx4(a,n1,n2)
    use tools_io, only: ferror, faterr

    complex*8, intent(inout), allocatable :: a(:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    complex*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call ferror('realloc2cmplx4','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    allocate(temp(n1,n2))
    
    temp = cmplx(0,8)
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2cmplx4

  !> Adapt the size of an allocatable 2D complex*8 array
  subroutine realloc4cmplx4(a,n1,n2,n3,n4)
    use tools_io, only: ferror, faterr

    complex*8, intent(inout), allocatable :: a(:,:,:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2, n3, n4 !< new dimension
    
    complex*8, allocatable :: temp(:,:,:,:)
    integer :: nold(4), i
    
    if (.not.allocated(a)) &
       call ferror('realloc4cmplx4','array not allocated',faterr)
    do i = 1, 4
       nold(i) = size(a,i)
    end do
    allocate(temp(n1,n2,n3,n4))
    
    temp = cmplx(0,8)
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)))
    call move_alloc(temp,a)

  end subroutine realloc4cmplx4

  !> Adapt the size of an allocatable 1D complex*16 array
  subroutine realloc1cmplx8(a,nnew)
    use tools_io, only: ferror, faterr

    complex*16, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    complex*16, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1cmplx8','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1cmplx8

end module types
