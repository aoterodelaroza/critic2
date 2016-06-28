!> Part of the following code has been adapted from postg
!> Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
!> Felix Kannemann <felix.kannemann@dal.ca>, Erin R. Johnson <ejohnson29@ucmerced.edu>,
!> Ross M. Dickson <ross.dickson@dal.ca>, Hartmut Schmider <hs7@post.queensu.ca>,
!> and Axel D. Becke <axel.becke@dal.ca>

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

! molecular wavefunction readers and tools
module wfn_private
  implicit none
  
  private

  public :: wfn_end
  public :: wfn_read_xyz_geometry
  public :: wfn_read_wfn_geometry
  public :: wfn_read_wfx_geometry
  public :: wfn_read_fchk_geometry
  public :: wfn_read_molden_geometry
  public :: wfn_read_wfn
  public :: wfn_read_wfx
  public :: wfn_register_struct
  public :: wfn_rho2
  public :: wfx_read_integers
  public :: wfx_read_reals1

  ! wfn type identifier
  integer, parameter, public :: wfn_rhf = 0
  integer, parameter, public :: wfn_uhf = 1
  integer, parameter, public :: wfn_frac = 2

  ! atomic positions, internal copy
  integer :: nat
  real*8, allocatable :: xat(:,:)

contains

  !> Terminate the wfn arrays
  subroutine wfn_end()

    ! if (allocated(xat)) deallocate(xat)

  end subroutine wfn_end

  !> Read the molecular geometry from an xyz file
  subroutine wfn_read_xyz_geometry(file,n,at)
    use types
    use tools_io
    use param

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    type(atom), allocatable, intent(out) :: at(:) !< Atoms

    character*4 :: atsym
    integer :: lu
    integer :: i, lp
    logical :: ok

    ! read the number of atoms
    lu = fopen_read(file)
    read (lu,*) n
    read (lu,*)
    if (allocated(at)) deallocate(at)
    allocate(at(n))

    do i = 1, n
       read (lu,*) atsym, at(i)%x
       at(i)%z = zatguess(atsym)
       if (at(i)%z <= 0) then
          ! maybe it's a number
          lp = 1
          ok = isinteger(at(i)%z,atsym,lp)
          if (.not.ok) &
             call ferror('wfn_read_xyz_geometry','could not determine atomic number',faterr)
       end if
       at(i)%name = trim(adjustl(atsym))
       at(i)%x = at(i)%x / bohrtoa
    end do
    call fclose(lu)

  end subroutine wfn_read_xyz_geometry

  !> Read the molecular geometry from a wfn file
  subroutine wfn_read_wfn_geometry(file,n,at)
    use types
    use tools_io

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    type(atom), allocatable, intent(out) :: at(:) !< Atoms

    character*4 :: atsym, orbtyp
    integer :: lu
    integer :: i, i1, i2
    real*8 :: zreal

    lu = fopen_read(file)

    ! read the number of atoms
    read (lu,*)
    read (lu,101) orbtyp, i1, i2, n
    if (n <= 0) &
       call ferror('wfn_read_wfn_geometry','wrong number of atoms',faterr)
    if (allocated(at)) deallocate(at)
    allocate(at(n))

    ! read the geometry
    do i = 1, n
       read(lu,106) atsym, at(i)%x, zreal
       at(i)%z = zatguess(atsym)
       at(i)%name = trim(adjustl(atsym))
    end do
101 format (4X,A4,10X,3(I5,15X))
106 format(2X,A2,20X,3F12.8,10X,F5.1)

    call fclose(lu)

  end subroutine wfn_read_wfn_geometry

  !> Read the molecular geometry from a wfx file
  subroutine wfn_read_wfx_geometry(file,n,at)
    use types
    use tools_io

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    type(atom), allocatable, intent(out) :: at(:) !< Atoms

    integer :: lu
    character(len=:), allocatable :: line, line2
    integer, allocatable :: iz(:)
    real*8, allocatable :: x(:,:)
    integer :: i

    lu = fopen_read(file)

    ! read the number of atoms
    do while (getline_raw(lu,line,.true.))
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Number of Nuclei>") then
             read (lu,*) n
             exit
          endif
       endif
    enddo
    if (n == 0) &
       call ferror("wfn_read_wfx_geometry","Number of Nuclei tag not found",faterr)
    if (allocated(at)) deallocate(at)
    allocate(at(n),iz(n),x(3,n))

    ! read the geometry
    rewind(lu)
    do while (.true.)
       read(lu,'(A)',end=20) line
       line2 = adjustl(line)
       line = line2
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Atomic Numbers>") then
             iz = wfx_read_integers(lu,n)
             do i = 1, n
                at(i)%z = iz(i)
                at(i)%name = nameguess(iz(i))
             end do
          elseif (trim(line) == "<Nuclear Cartesian Coordinates>") then
             x = reshape(wfx_read_reals1(lu,3*n),shape(x))
             do i = 1, n
                at(i)%x = x(:,i)
             end do
          endif
       endif
    enddo
20  continue
    deallocate(iz,x)

    call fclose(lu)

  end subroutine wfn_read_wfx_geometry

  !> Read the molecular geometry from a fchk file
  subroutine wfn_read_fchk_geometry(file,n,at)
    use types
    use tools_io

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    type(atom), allocatable, intent(out) :: at(:) !< Atoms

    integer :: lu, lp, i, j
    character(len=:), allocatable :: line
    logical :: ok
    real*8, allocatable :: xat(:)
    integer, allocatable :: zat(:)

    lu = fopen_read(file)

    ! read the number of atoms
    n = 0
    do while (getline_raw(lu,line,.true.))
       lp = 45
       if (line(1:15) == "Number of atoms") then
          ok = isinteger(n,line,lp)
          if (.not.ok) call ferror('wfn_read_fchk_geometry','could not read number of atoms',faterr)
          exit
       endif
    enddo
    if (n == 0) &
       call ferror("wfn_read_fchk_geometry","error reading number of atoms",faterr)
    if (allocated(at)) deallocate(at)
    allocate(at(n))

    ! read the geometry
    allocate(zat(n),xat(3*n))
    rewind(lu)
    do while (getline_raw(lu,line))
       lp = 45
       if (line(1:29) == "Current cartesian coordinates") then
          do i = 0, (3*n-1)/5
             read(lu,'(5E16.8)') (xat(5*i+j),j=1,min(5,3*n-5*i))
          enddo
       elseif (line(1:14) == "Atomic numbers") then
          do i = 0, (n-1)/6
             read(lu,'(6I12)') (zat(6*i+j),j=1,min(6,n-6*i))
          enddo
       endif
    enddo
    do i = 1, n
       at(i)%x = xat((i-1)*3+1:(i-1)*3+3)
       at(i)%z = zat(i)
       at(i)%name = nameguess(zat(i))
    end do
    
    ! clean up
    deallocate(xat,zat)
    call fclose(lu)

  end subroutine wfn_read_fchk_geometry

  !> Read the molecular geometry from a molden file (only tested with
  !> new psi4 molden files).
  subroutine wfn_read_molden_geometry(file,n,at)
    use types
    use tools_io
    use param

    character*(*), intent(in) :: file !< Input file name
    integer, intent(out) :: n !< Number of atoms
    type(atom), allocatable, intent(out) :: at(:) !< Atoms

    integer :: lu, lp, idum, i
    character(len=:), allocatable :: line, keyword, word1, word2
    character*(1024) :: fixword
    logical :: ok, isang

    lu = fopen_read(file)

    ! read the number of atoms
    n = 0
    do while(next_keyword())
       if (trim(lower(keyword)) == "atoms") then
          n = 0
          ok = getline_raw(lu,line,.true.)
          do while(index(lower(line),"[") == 0 .and. len(trim(line)) > 0)
             n = n + 1
             ok = getline_raw(lu,line,.true.)
          end do
          exit
       end if
    end do
    if (n == 0) &
       call ferror("wfn_read_molden_geometry","error reading number of atoms",faterr)
    if (allocated(at)) deallocate(at)
    allocate(at(n))

    ! read the geometry
    rewind(lu)
    do while(next_keyword())
       if (trim(lower(keyword)) == "atoms") exit
    end do
    if (len_trim(keyword) == 0) &
       call ferror("wfn_read_molden_geometry","error reading geometry",faterr)

    ! geometry header -> detect the units for the geometry
    lp = 1
    word1 = lgetword(line,lp)
    word2 = lgetword(line,lp)
    isang = (trim(lower(word2)) == "(angs)".or.trim(lower(word2)) == "(ang)") 

    ! read the atomic numbers and positions
    do i = 1, n
       ok = getline_raw(lu,line,.true.)
       read(line,*) fixword, idum, at(i)%z, at(i)%x
       if (isang) at(i)%x = at(i)%x / bohrtoa
       at(i)%name = nameguess(at(i)%z)
    end do

    call fclose(lu)

  contains

    function next_keyword()

      integer :: istart, iend
      logical :: next_keyword

      keyword = ""
      next_keyword = .false.

      do while(.true.)
         ok = getline_raw(lu,line,.false.)
         if (.not.ok) return
         if (index(lower(line),"[") > 0) exit
      end do
      next_keyword =.true.
      istart = index(lower(line),"[") + 1
      iend = index(lower(line),"]") - 1
      keyword = line(istart:iend)
      return

    end function next_keyword

  end subroutine wfn_read_molden_geometry

  !> Read the wavefunction from a wfn file
  subroutine wfn_read_wfn(file,f)
    use tools_io
    use types
    
    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f !< Output field

    integer :: luwfn, nat, iz, istat, i, j, num1, num2
    integer :: nalpha, ioc
    real*8 :: x(3), zreal, ene0, ene
    character*4 :: orbtyp
    character*2 :: dums, elem
    logical :: isfrac
    character*8 :: dum1

    f%useecp = .false.

    ! read number of atoms, primitives, orbitals
    luwfn = fopen_read(file)
    read (luwfn,*)
    read (luwfn,101) orbtyp, f%nmo, f%npri, nat
 
    ! atomic positions and numbers
    do i = 1, nat
       read(luwfn,106) elem, x, zreal
       iz = zatguess(elem)
       if (iz /= nint(zreal)) f%useecp = .true.
    end do
 
    ! center assignments, types of primitives
    if (allocated(f%icenter)) deallocate(f%icenter)
    allocate(f%icenter(f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfn','could not allocate memory for icenter',faterr)
    if (allocated(f%itype)) deallocate(f%itype)
    allocate(f%itype(f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfn','could not allocate memory for itype',faterr)
    read(luwfn,102) (f%icenter(i),i=1,f%npri)
    read(luwfn,102) (f%itype(i),i=1,f%npri)
    if (any(f%itype(1:f%npri) > 56)) then
       call ferror("wfn_read_wfn","primitive type not supported",faterr)
    endif

    ! primitive exponents
    if (allocated(f%e)) deallocate(f%e)
    allocate(f%e(f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfn','could not allocate memory for exponents',faterr)
    read(luwfn,103) (f%e(i),i=1,f%npri)
 
    ! deal with ecps
    dums=""
    do while (dums.ne."MO")
       read (luwfn,'(A2)') dums
    enddo
    backspace(luwfn)
 
    ! occupations and orbital coefficients
    if (allocated(f%occ)) deallocate(f%occ)
    allocate(f%occ(f%nmo),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfn','could not allocate memory for occupations',faterr)
    if (allocated(f%cmo)) deallocate(f%cmo)
    allocate(f%cmo(f%nmo,f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfn','could not allocate memory for orbital coefficients',faterr)
    isfrac = .false.
    num1 = 0
    num2 = 0
    ene0 = -1d30
    nalpha = -1
    do i = 1, f%nmo
       read(luwfn,104) f%occ(i), ene
       read(luwfn,105) (f%cmo(i,j),j=1,f%npri)
       ioc = nint(f%occ(i))
       if (abs(ioc-f%occ(i)) > 1d-10) then
          isfrac = .true.
       else if (ioc == 1) then
          num1 = num1 + 1
       else if (ioc == 2) then
          num2 = num2 + 1
       endif
       if (ene < ene0-1d-3) nalpha = i-1
       ene0 = ene
    end do
    read(luwfn,*) dum1
 
    ! figure out charge and multiplicity
    ! 0 - closed, 1 - open
    if (isfrac) then
       f%wfntyp = wfn_frac
       f%nalpha = 0
    else if (num1 == 0) then
       f%wfntyp = wfn_rhf
       f%nalpha = 0
    else if (num2 == 0) then
       f%wfntyp = wfn_uhf
       f%nalpha = nalpha
    else
       call ferror("wfn_read_wfn","restricted-open wfn files not supported",faterr)
    endif
    call fclose(luwfn)

    ! we are done
    f%init = .true.

101  format (4X,A4,10X,3(I5,15X))
102  format(20X,20I3)
103  format(10X,5E14.7)
104  format(35X,F12.7,15X,F12.6)
105  format(5(E16.8))
106  format(2X,A2,20X,3F12.8,10X,F5.1)

  end subroutine wfn_read_wfn

  !> Read the wavefunction from a wfx file
  subroutine wfn_read_wfx(file,f)
    use tools_io
    use types

    character*(*), intent(in) :: file !< Input file
    type(field), intent(out) :: f !< Output field

    integer :: luwfn, mult, ncore, istat, i
    character(len=:), allocatable :: line, line2
    logical :: ok

    f%useecp = .false.

    ! first pass
    luwfn = fopen_read(file)
    f%nmo = 0
    mult = 0
    ncore = 0
    f%npri = 0
    f%nedf = 0
    do while (getline_raw(luwfn,line))
       if (len_trim(line) < 1) exit
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Number of Occupied Molecular Orbitals>") then
             read (luwfn,*) f%nmo
          elseif (trim(line) == "<Electronic Spin Multiplicity>") then
             read (luwfn,*) mult
          elseif (trim(line) == "<Number of Core Electrons>") then
             read (luwfn,*) ncore
          elseif (trim(line) == "<Number of Primitives>") then
             read(luwfn,*) f%npri
          elseif (trim(line) == "<Number of EDF Primitives>") then
             read(luwfn,*) f%nedf
          endif
       endif
    enddo
    
    if (f%nmo == 0) call ferror("wfn_read_wfx","Number of Occupied Molecular Orbitals tag not found",faterr)
    if (mult == 0) call ferror("wfn_read_wfx","Electronic Spin Multiplicity tag not found",faterr)
    if (f%npri == 0) call ferror("wfn_read_wfx","Number of Primitives tag not found",faterr)
    if (ncore > 0) f%useecp = .true.

    ! allocate memory
    allocate(f%icenter(f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for icenter',faterr)
    allocate(f%itype(f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for itype',faterr)
    allocate(f%e(f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for exponents',faterr)
    allocate(f%occ(f%nmo),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for occupations',faterr)
    allocate(f%cmo(f%nmo,f%npri),stat=istat)
    if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for orbital coefficients',faterr)
    if (f%nedf > 0) then
       allocate(f%icenter_edf(f%nedf),stat=istat)
       if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for icenter',faterr)
       allocate(f%itype_edf(f%nedf),stat=istat)
       if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for itype',faterr)
       allocate(f%e_edf(f%nedf),stat=istat)
       if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for exponents',faterr)
       allocate(f%c_edf(f%npri),stat=istat)
       if (istat /= 0) call ferror('wfn_read_wfx','could not allocate memory for orbital coefficients',faterr)
    endif

    ! second pass
    rewind(luwfn)
    do while (.true.)
       read(luwfn,'(A)',end=20) line
       line2 = adjustl(line)
       line = line2
       if (line(1:1) == "<" .and. line(2:2) /= "/") then
          if (trim(line) == "<Primitive Centers>") then
             f%icenter = wfx_read_integers(luwfn,f%npri)
          elseif (trim(line) == "<Primitive Types>") then
             f%itype = wfx_read_integers(luwfn,f%npri)
             if (any(f%itype(1:f%npri) > 56)) &
                call ferror("wfn_read_wfx","primitive type not supported",faterr)
          elseif (trim(line) == "<Primitive Exponents>") then
             f%e = wfx_read_reals1(luwfn,f%npri)
          elseif (trim(line) == "<Molecular Orbital Occupation Numbers>") then
             f%occ = wfx_read_reals1(luwfn,f%nmo)
          elseif (trim(line) == "<Molecular Orbital Primitive Coefficients>") then
             read(luwfn,*)
             do i = 1, f%nmo
                read(luwfn,*)
                read(luwfn,*)
                f%cmo(i,:) = wfx_read_reals1(luwfn,f%npri)
             enddo
          elseif (trim(line) == "<EDF Primitive Centers>") then
             f%icenter_edf = wfx_read_integers(luwfn,f%nedf)
          elseif (trim(line) == "<EDF Primitive Types>") then
             f%itype_edf = wfx_read_integers(luwfn,f%nedf)
             if (any(f%itype_edf(1:f%nedf) > 56)) &
                call ferror("wfn_read_wfx","primitive type not supported",faterr)
          elseif (trim(line) == "<EDF Primitive Exponents>") then
             f%e_edf = wfx_read_reals1(luwfn,f%nedf)
          elseif (trim(line) == "<EDF Primitive Coefficients>") then
             f%c_edf = wfx_read_reals1(luwfn,f%nedf)
          endif
       endif
    enddo
20  continue

    ! wavefuntion type
    if (mult == 1) then
       f%wfntyp = wfn_rhf
       f%nalpha = 0
    else
       f%wfntyp = wfn_uhf
       f%nalpha = (f%nmo - mult + 1) / 2
    end if
    
    ! done
    f%init = .true.
    call fclose(luwfn)

  end subroutine wfn_read_wfx

  !> Register structural information.
  subroutine wfn_register_struct(ncel,atcel)
    use types

    integer, intent(in) :: ncel
    type(celatom), intent(in) :: atcel(:)

    integer :: i

    nat = ncel
    if (allocated(xat)) deallocate(xat)
    allocate(xat(3,nat))
    do i = 1, nat
       xat(:,i) = atcel(i)%r
    end do
    
  end subroutine wfn_register_struct

  !> Determine the density and derivatives at a given target point.
  !> xpos is in cartesian coordiantes and assume that the molecule
  !> has been displaced to the center of a big cube. Same transformation
  !> as in load xyz/wfn/wfx.
  subroutine wfn_rho2(f,xpos,nder,rho,grad,h,gkin,vir,stress,xmo)
    use tools_io
    use types
    use param

    type(field), intent(in) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    integer, intent(in) :: nder  !< Number of derivatives
    real*8, intent(out) :: rho !< Density
    real*8, intent(out) :: grad(3) !< Gradient
    real*8, intent(out) :: h(3,3) !< Hessian 
    real*8, intent(out) :: gkin !< G(r), kinetic energy density
    real*8, intent(out) :: vir !< Virial field
    real*8, intent(out) :: stress(3,3) !< Schrodinger stress tensor
    real*8, intent(out), optional :: xmo(f%nmo) !< Values of the MO

    integer, parameter :: imax(0:2) = (/1,4,10/)
    
    real*8 :: al, ex, xl(3,0:2), xl2
    integer :: ipri, iat, ityp, l(3), ix, i
    integer :: imo
    real*8 :: chi(f%npri,imax(nder)), maxc(f%npri), dd(3,nat), d2(nat)
    real*8 :: phi(f%nmo,imax(nder)), aocc
    real*8 :: hh(6)
    logical :: ldopri(f%npri,imax(nder))
    
    real*8, parameter :: cutoff_pri = 1d-15
    integer, parameter :: li(3,56) = reshape((/&
       0,0,0, & ! s
       1,0,0, 0,1,0, 0,0,1, & ! p
       2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1, & !d
       3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 0,2,1, &
       1,2,0, 1,0,2, 0,1,2, 1,1,1,& ! f
       4,0,0, 0,4,0, 0,0,4, 3,1,0, 3,0,1, 1,3,0, 0,3,1, 1,0,3,&
       0,1,3, 2,2,0, 2,0,2, 0,2,2, 2,1,1, 1,2,1, 1,1,2,& ! g
       0,0,5, 0,1,4, 0,2,3, 0,3,2, 0,4,1, 0,5,0, 1,0,4, 1,1,3,&
       1,2,2, 1,3,1, 1,4,0, 2,0,3, 2,1,2, 2,2,1, 2,3,0, 3,0,2,&
       3,1,1, 3,2,0, 4,0,1, 4,1,0, 5,0,0/),shape(li)) ! h

    ! identify the max coefficients
    maxc = 0d0
    do imo = 1, f%nmo
       do ipri = 1, f%npri
          maxc(ipri) = max(maxc(ipri),abs(f%cmo(imo,ipri)))
       enddo
    enddo

    ! calculate distances
    do iat = 1, nat
       dd(:,iat) = xpos - xat(:,iat)
       d2(iat) = dd(1,iat)*dd(1,iat)+dd(2,iat)*dd(2,iat)+dd(3,iat)*dd(3,iat)
    enddo

    ! primitive and derivatives 
    do ipri = 1, f%npri
       ityp = f%itype(ipri)
       iat = f%icenter(ipri)
       al = f%e(ipri)
       ex = exp(-al * d2(iat))

       l = li(1:3,ityp)
       do ix = 1, 3
          if (l(ix) == 0) then
             xl(ix,0) = 1d0
             xl(ix,1) = 0d0
             xl(ix,2) = 0d0
          else if (l(ix) == 1) then
             xl(ix,0) = dd(ix,iat)
             xl(ix,1) = 1d0
             xl(ix,2) = 0d0
          else if (l(ix) == 2) then
             xl(ix,0) = dd(ix,iat) * dd(ix,iat)
             xl(ix,1) = 2d0 * dd(ix,iat)
             xl(ix,2) = 2d0
          else if (l(ix) == 3) then
             xl(ix,0) = dd(ix,iat) * dd(ix,iat) * dd(ix,iat)
             xl(ix,1) = 3d0 * dd(ix,iat) * dd(ix,iat)
             xl(ix,2) = 6d0 * dd(ix,iat)
          else if (l(ix) == 4) then
             xl2 = dd(ix,iat) * dd(ix,iat)
             xl(ix,0) = xl2 * xl2
             xl(ix,1) = 4d0 * xl2 * dd(ix,iat)
             xl(ix,2) = 12d0 * xl2
          else if (l(ix) == 5) then
             xl2 = dd(ix,iat) * dd(ix,iat)
             xl(ix,0) = xl2 * xl2 * dd(ix,iat)
             xl(ix,1) = 5d0 * xl2 * xl2
             xl(ix,2) = 20d0 * xl2 * dd(ix,iat)
          else
             call ferror('wfn_rho2','power of L not supported',faterr)
          end if
       end do

       chi(ipri,1) = xl(1,0)*xl(2,0)*xl(3,0)*ex
       if (nder > 0) then
          chi(ipri,2) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * xl(2,0)*xl(3,0)*ex
          chi(ipri,3) = (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*xl(3,0)*ex
          chi(ipri,4) = (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(1,0)*xl(2,0)*ex
          if (nder > 1) then
             chi(ipri,5) = (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0) + 4*al*al*dd(1,iat)**(l(1)+2)) * xl(2,0)*xl(3,0)*ex
             chi(ipri,6) = (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0) + 4*al*al*dd(2,iat)**(l(2)+2)) * xl(3,0)*xl(1,0)*ex
             chi(ipri,7) = (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0) + 4*al*al*dd(3,iat)**(l(3)+2)) * xl(1,0)*xl(2,0)*ex
             chi(ipri,8) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(3,0)*ex
             chi(ipri,9) = (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(2,0)*ex
             chi(ipri,10)= (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*ex
          endif
       endif

       do ix = 1, imax(nder)
          ldopri(ipri,ix) = (abs(chi(ipri,ix))*maxc(ipri) > cutoff_pri)
       enddo
    enddo ! ipri = 1, npri

    ! build the MO avlues at the point
    phi = 0d0
    do ix = 1, imax(nder)
       do ipri = 1, f%npri
          if (.not.ldopri(ipri,ix)) cycle
          do imo = 1, f%nmo
             phi(imo,ix) = phi(imo,ix) + f%cmo(imo,ipri)*chi(ipri,ix)
          enddo
       enddo
    enddo

    ! contribution to the density, etc.
    rho = 0d0
    grad = 0d0
    hh = 0d0
    gkin = 0d0
    vir = 0d0
    stress = 0d0
    do imo = 1, f%nmo
       aocc = f%occ(imo) 
       rho = rho + aocc * phi(imo,1) * phi(imo,1)
       if (nder>0) then
          grad = grad + 2 * aocc * phi(imo,1) * phi(imo,2:4) 
          gkin = gkin + aocc * (phi(imo,2)*phi(imo,2) + phi(imo,3)*phi(imo,3) + phi(imo,4)*phi(imo,4))
       endif
       if (nder>1) then
          hh(1:3) = hh(1:3) + 2 * aocc * (phi(imo,1)*phi(imo,5:7)+phi(imo,2:4)**2)
          hh(4) = hh(4) + 2 * aocc * (phi(imo,1)*phi(imo,8)+phi(imo,2)*phi(imo,3))
          hh(5) = hh(5) + 2 * aocc * (phi(imo,1)*phi(imo,9)+phi(imo,2)*phi(imo,4))
          hh(6) = hh(6) + 2 * aocc * (phi(imo,1)*phi(imo,10)+phi(imo,3)*phi(imo,4))
          stress(1,1) = stress(1,1) + aocc * (phi(imo,1) * phi(imo,5)  - phi(imo,2)*phi(imo,2))
          stress(1,2) = stress(1,2) + aocc * (phi(imo,1) * phi(imo,8)  - phi(imo,2)*phi(imo,3))
          stress(1,3) = stress(1,3) + aocc * (phi(imo,1) * phi(imo,9)  - phi(imo,2)*phi(imo,4))
          stress(2,2) = stress(2,2) + aocc * (phi(imo,1) * phi(imo,6)  - phi(imo,3)*phi(imo,3))
          stress(2,3) = stress(2,3) + aocc * (phi(imo,1) * phi(imo,10) - phi(imo,3)*phi(imo,4))
          stress(3,3) = stress(3,3) + aocc * (phi(imo,1) * phi(imo,7)  - phi(imo,4)*phi(imo,4))
       endif
    enddo
    stress(2,1) = stress(1,2)
    stress(3,1) = stress(1,3)
    stress(3,2) = stress(2,3)
    stress = stress * 0.5d0
    gkin = 0.5d0 * gkin
    vir = stress(1,1)+stress(2,2)+stress(3,3)

    ! core contribution to the density
    if (f%nedf > 0) then
       do i = 1, f%nedf
          ityp = f%itype_edf(i)
          iat = f%icenter_edf(i)
          al = f%e_edf(i)
          ex = exp(-al * d2(iat))

          l = li(1:3,ityp)
          do ix = 1, 3
             if (l(ix) == 0) then
                xl(ix,0) = 1d0
                xl(ix,1) = 0d0
                xl(ix,2) = 0d0
             else if (l(ix) == 1) then
                xl(ix,0) = dd(ix,iat)
                xl(ix,1) = 1d0
                xl(ix,2) = 0d0
             else if (l(ix) == 2) then
                xl(ix,0) = dd(ix,iat) * dd(ix,iat)
                xl(ix,1) = 2d0 * dd(ix,iat)
                xl(ix,2) = 2d0
             else if (l(ix) == 3) then
                xl(ix,0) = dd(ix,iat) * dd(ix,iat) * dd(ix,iat)
                xl(ix,1) = 3d0 * dd(ix,iat) * dd(ix,iat)
                xl(ix,2) = 6d0 * dd(ix,iat)
             else if (l(ix) == 4) then
                xl2 = dd(ix,iat) * dd(ix,iat)
                xl(ix,0) = xl2 * xl2
                xl(ix,1) = 4d0 * xl2 * dd(ix,iat)
                xl(ix,2) = 12d0 * xl2
             else if (l(ix) == 5) then
                xl2 = dd(ix,iat) * dd(ix,iat)
                xl(ix,0) = xl2 * xl2 * dd(ix,iat)
                xl(ix,1) = 5d0 * xl2 * xl2
                xl(ix,2) = 20d0 * xl2 * dd(ix,iat)
             else
                call ferror('wfn_rho2','power of L not supported',faterr)
             end if
          end do
          
          rho = rho + f%c_edf(i) * xl(1,0)*xl(2,0)*xl(3,0)*ex
          if (nder > 0) then
             grad(1) = grad(1) + f%c_edf(i) * (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * xl(2,0)*xl(3,0)*ex
             grad(2) = grad(2) + f%c_edf(i) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*xl(3,0)*ex
             grad(3) = grad(3) + f%c_edf(i) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(1,0)*xl(2,0)*ex
             if (nder > 1) then
                hh(1) = hh(1) + f%c_edf(i) * (xl(1,2)-2*al*(2*l(1)+1)*xl(1,0) + 4*al*al*dd(1,iat)**(l(1)+2)) * xl(2,0)*xl(3,0)*ex
                hh(2) = hh(2) + f%c_edf(i) * (xl(2,2)-2*al*(2*l(2)+1)*xl(2,0) + 4*al*al*dd(2,iat)**(l(2)+2)) * xl(3,0)*xl(1,0)*ex
                hh(3) = hh(3) + f%c_edf(i) * (xl(3,2)-2*al*(2*l(3)+1)*xl(3,0) + 4*al*al*dd(3,iat)**(l(3)+2)) * xl(1,0)*xl(2,0)*ex
                hh(4) = hh(4) + f%c_edf(i) * (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(3,0)*ex
                hh(5) = hh(5) + f%c_edf(i) * (xl(1,1)-2*al*dd(1,iat)**(l(1)+1)) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * xl(2,0)*ex
                hh(6) = hh(6) + f%c_edf(i) * (xl(3,1)-2*al*dd(3,iat)**(l(3)+1)) * (xl(2,1)-2*al*dd(2,iat)**(l(2)+1)) * xl(1,0)*ex
             endif
          endif
       end do
    end if

    ! save the MO values
    if (present(xmo)) xmo = phi(1:f%nmo,1)

    ! re-order the hessian
    do i = 1, 3
       h(i,i) = hh(i)
    end do
    h(1,2) = hh(4)
    h(2,1) = h(1,2)
    h(1,3) = hh(5)
    h(3,1) = h(1,3)
    h(2,3) = hh(6)
    h(3,2) = h(2,3)

  end subroutine wfn_rho2

  !> Read a list of n integers from a logical unit
  function wfx_read_integers(lu,n) result(x)
    use tools_io
    integer, intent(in) :: lu, n
    integer :: x(n)

    integer :: kk, lp, idum
    character(len=:), allocatable :: line
    logical :: ok

    kk = 0
    lp = 1
    ok = getline_raw(lu,line,.true.)
    do while(.true.)
       if (.not.isinteger(idum,line,lp)) then
          lp = 1
          ok = getline_raw(lu,line)
          if (.not.ok .or. line(1:2) == "</") exit
       else
          kk = kk + 1
          if (kk > n) call ferror("wfx_read_integers","exceeded size of the array",faterr)
          x(kk) = idum
       endif
    enddo

  endfunction wfx_read_integers

  !> Read a list of n reals from a logical unit
  function wfx_read_reals1(lu,n) result(x)
    use tools_io
    integer, intent(in) :: lu, n
    real*8 :: x(n)

    integer :: kk, lp
    real*8 :: rdum
    character(len=:), allocatable :: line
    logical :: ok

    kk = 0
    lp = 1
    ok = getline_raw(lu,line,.true.)
    do while(.true.)
       if (.not.isreal(rdum,line,lp)) then
          lp = 1
          ok = getline_raw(lu,line)
          if (.not.ok .or. line(1:1) == "<") exit
       else
          kk = kk + 1
          if (kk > n) call ferror("wfx_read_reals1","exceeded size of the array",2)
          x(kk) = rdum
       endif
    enddo

  endfunction wfx_read_reals1

end module wfn_private
