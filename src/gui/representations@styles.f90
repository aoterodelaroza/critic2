! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Geometry-dependent display styles of the representations: atoms,
! molecules, bonds, labels, coordination polyhedra, symmetry elements.
submodule (representations) styles
  implicit none

contains

  !> Reset a symmetry-element style: recalculate the list of
  !> symmetry-element types of the system, keeping the previous
  !> visibility selection if the list has not changed size.
  module subroutine symelem_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    class(symelem_style), intent(inout) :: d
    type(representation), intent(in) :: r

    logical, allocatable :: shownold(:)

    ! reset the time and remember the previous selection
    d%timelastreset = glfwGetTime()
    if (allocated(d%shown)) call move_alloc(d%shown,shownold)
    d%isinit = .false.
    call d%se%end()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! recompute the element-type snapshot; the element positions are not kept
    ! here (they depend on the number of cells drawn)
    call sys(r%id)%c%list_symelems((/1,1,1/),.true.,d%se,typesonly=.true.)

    ! restore the previous visibility if the element count is unchanged,
    ! otherwise show all elements
    allocate(d%shown(d%se%ntype))
    if (allocated(shownold)) then
       if (size(shownold,1) == d%se%ntype) then
          d%shown = shownold
       else
          d%shown = .true.
       end if
    else
       d%shown = .true.
    end if
    d%isinit = .true.

  end subroutine symelem_style_reset

  !> Deallocate all arrays and end the symmetry-element style.
  module subroutine symelem_style_end(d)
    class(symelem_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    d%nop = 0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%iop)) deallocate(d%iop)
    call d%se%end()

  end subroutine symelem_style_end

  !> Reset all styles. Reset if itype = 0 (all), 1 (atom), 2 (bonds),
  !> 3 (labels), 4 (mol).
  module subroutine reset_all_styles(r,itype)
    class(representation), intent(inout) :: r
    integer, intent(in) :: itype

    if (itype == 0 .or. itype == 1) &
       call r%atoms%style%reset(r)
    if (itype == 0 .or. itype == 2) &
       call r%bonds%style%reset(r)
    if (itype == 0 .or. itype == 3) &
       call r%labels%style%reset(r)
    if (itype == 0 .or. itype == 4) &
       call r%mols%style%reset(r)
    if (itype == 0 .or. itype == 8) &
       call r%poly%style%reset(r)

  end subroutine reset_all_styles

  !> Reset atom style to the parameters and the contents of the system
  !> point at by representation r. Uses d%type to fill the arrays.
  module subroutine atom_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_ready, ok_system, atlisttype_species
    use gui_main, only: ColorElement
    use param, only: atmcov, atmvdw, jmlcol, jmlcol2
    class(atom_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, ispc, iz

    ! if not initialized, set type
    if (.not.d%isinit) d%type = atlisttype_species

    ! set the atom style to zero
    d%ntype = 0
    d%isinit = .false.
    if (allocated(d%rgb)) deallocate(d%rgb)
    if (allocated(d%rad)) deallocate(d%rad)

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! fill data
    d%ntype = sysc(r%id)%attype_number(d%type)
    allocate(d%rgb(3,d%ntype),d%rad(d%ntype))
    do i = 1, d%ntype
       ispc = sysc(r%id)%attype_species(d%type,i)
       iz = sys(r%id)%c%spc(ispc)%z

       ! color scheme or GUI colors
       if (r%atoms%color_type == 0) then
          d%rgb(:,i) = ColorElement(:,iz)
       elseif (r%atoms%color_type == 1) then
          d%rgb(:,i) = real(jmlcol(:,iz),c_float) / 255._c_float
       else
          d%rgb(:,i) = real(jmlcol2(:,iz),c_float) / 255._c_float
       end if

       ! scale covalent, vdw, or absolute value
       if (r%atoms%radii_type == 0) then
          d%rad(i) = r%atoms%radii_scale * atmcov(iz)
       elseif (r%atoms%radii_type == 1) then
          d%rad(i) = r%atoms%radii_scale * atmvdw(iz)
       else
          d%rad(i) = r%atoms%radii_value
       endif
    end do
    d%isinit = .true.

  end subroutine atom_style_reset

  !> Reset colors in an atom style to defaults.
  module subroutine atom_style_reset_colors(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_ready, ok_system
    use gui_main, only: ColorElement
    class(atom_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, ispc, iz

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    do i = 1, d%ntype
       ispc = sysc(r%id)%attype_species(d%type,i)
       iz = sys(r%id)%c%spc(ispc)%z
       d%rgb(:,i) = ColorElement(:,iz)
    end do

  end subroutine atom_style_reset_colors

  !> Deallocate all arrays and end the atom syle.
  module subroutine atom_style_end(d)
    class(atom_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    if (allocated(d%rgb)) deallocate(d%rgb)
    if (allocated(d%rad)) deallocate(d%rad)

  end subroutine atom_style_end

  !> Reset molecule style with default values. Use the information in
  !> representation r, or leave it empty if system is uninitalized.
  module subroutine mol_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    class(mol_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i

    ! set the atom style to zero
    d%ntype = 0
    d%isinit = .false.
    if (allocated(d%tint_rgb)) deallocate(d%tint_rgb)
    if (allocated(d%scale_rad)) deallocate(d%scale_rad)

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! fill
    d%ntype = sys(r%id)%c%nmol
    allocate(d%tint_rgb(3,d%ntype))
    allocate(d%scale_rad(d%ntype))
    do i = 1, sys(r%id)%c%nmol
       d%tint_rgb(:,i) = 1._c_float
       d%scale_rad(i) = 1d0
    end do
    d%isinit = .true.

  end subroutine mol_style_reset

  !> Deallocate all arrays and end the mol syle.
  module subroutine mol_style_end(d)
    class(mol_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    if (allocated(d%tint_rgb)) deallocate(d%tint_rgb)
    if (allocated(d%scale_rad)) deallocate(d%scale_rad)

  end subroutine mol_style_end

  !> Generate the neighbor stars from the data in the rij table using
  !> the geometry in system isys.
  module subroutine generate_neighstars(d,r)
    use systems, only: sys, sys_ready, ok_system
    class(bond_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    ! check all the info is available
    if (.not.d%isinit) return
    if (.not.ok_system(r%id,sys_ready)) return

    ! generate the new neighbor star using this representation's bonding
    ! parameters: per-species covalent radii and bond factor for non-metal
    ! bonds, bond delta for metal bonds (same criteria as the geometry window)
    call sys(r%id)%c%find_asterisms(d%nstar,atmrad=r%bonds%atmrad,bondfac=r%bonds%bfactor,&
       bonddelta=r%bonds%bdelta)

  end subroutine generate_neighstars

  !> Copy the neighbor stars from the given system.
  module subroutine copy_neighstars_from_system(d,isys)
    use systems, only: sys, sys_ready, ok_system
    class(bond_geom_style), intent(inout) :: d
    integer, intent(in) :: isys

    ! check all the info is available
    if (.not.d%isinit) return
    if (.not.ok_system(isys,sys_ready)) return

    ! copy the nstar
    if (allocated(sys(isys)%c%nstar)) &
       d%nstar = sys(isys)%c%nstar

  end subroutine copy_neighstars_from_system

  !> Reset bond style with default values, according to the given
  !> representation.
  module subroutine bond_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    class(bond_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, j, iz

    ! clear the bond style
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%nstar)) deallocate(d%nstar)
    d%use_sys_nstar = .true.

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return
    d%isinit = .true.

    ! fill temp options
    allocate(d%shown(sys(r%id)%c%nspc,sys(r%id)%c%nspc))
    d%shown = .true.

    ! fill data according to flavor
    if (r%flavor == repflavor_bonds_vdwcontacts) then
       ! van der waals contacts
       d%use_sys_nstar = .false.
       do i = 1, sys(r%id)%c%nspc
          if (sys(r%id)%c%spc(i)%z == 1) then
             d%shown(i,:) = .false.
             d%shown(:,i) = .false.
          end if
       end do
       call d%generate_neighstars(r)
    elseif (r%flavor == repflavor_bonds_hbonds) then
       ! hydrogen bonds
       d%use_sys_nstar = .false.
       d%shown = .false.
       do i = 1, sys(r%id)%c%nspc
          if (sys(r%id)%c%spc(i)%z /= 1) cycle
          do j = 1, sys(r%id)%c%nspc
             iz = sys(r%id)%c%spc(j)%z
             if (iz == 7 .or. iz == 8 .or. iz == 9 .or. iz == 16 .or. iz == 17) then
                d%shown(i,j) = .true.
                d%shown(j,i) = .true.
             end if
          end do
       end do
       call d%generate_neighstars(r)
    else
       ! other flavors are default (track system bonds)
       d%use_sys_nstar = .true.
       call d%copy_neighstars_from_system(r%id)
    end if

  end subroutine bond_style_reset

  !> Deallocate all arrays and end the bond syle.
  module subroutine bond_style_end(d)
    class(bond_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    d%use_sys_nstar = .true.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%nstar)) deallocate(d%nstar)

  end subroutine bond_style_end

  !> Reset label style with default values. Use the information in
  !> representation r.
  module subroutine label_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use tools_io, only: nameguess, string
    class(label_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i

    ! not initialized
    d%isinit = .false.
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%str)) deallocate(d%str)

    ! reset the time
    d%timelastreset = glfwGetTime()

    ! the critical point rows depend on the label type too
    call d%reset_cps(r)

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! fill according to the style
    d%isinit = .true.
    select case(r%labels%type)
    case (0,5,6)
       d%ntype = sys(r%id)%c%nspc
    case (2,3)
       d%ntype = sys(r%id)%c%ncel
    case (1,4,8)
       d%ntype = sys(r%id)%c%nneq
    case (7)
       d%ntype = sys(r%id)%c%nmol
    end select
    allocate(d%shown(d%ntype))
    allocate(d%str(d%ntype))

    ! fill shown, exclude hydrogens
    d%shown = .true.
    do i = 1, d%ntype
       select case(r%labels%type)
       case (0,5,6)
          if (sys(r%id)%c%spc(i)%z == 1) d%shown(i) = .false.
       case (2,3)
          if (sys(r%id)%c%spc(sys(r%id)%c%atcel(i)%is)%z == 1) d%shown(i) = .false.
       case (1,4,8)
          if (sys(r%id)%c%spc(sys(r%id)%c%at(i)%is)%z == 1) d%shown(i) = .false.
       end select
    end do

    ! fill text
    do i = 1, d%ntype
       if (r%labels%type == 0) then ! 0 = atomic symbol
          d%str(i) = trim(nameguess(sys(r%id)%c%spc(i)%z,.true.))
       elseif (r%labels%type == 1) then ! 1 = atom name
          d%str(i) = trim(sys(r%id)%c%at(i)%name)
       elseif (r%labels%type == 6) then ! 6 = Z
          d%str(i) = string(sys(r%id)%c%spc(i)%z)
       elseif (r%labels%type == 8) then ! 8 = wyckoff
          d%str(i) = string(sys(r%id)%c%at(i)%mult) // string(sys(r%id)%c%at(i)%wyc)
       else
          d%str(i) = string(i)
       end if
    end do

  end subroutine label_style_reset

  !> Reset the critical point rows of the label style from the
  !> non-nuclear CPs of field r%labels%ifield: one row per CP type
  !> (label type 0, text n/b/r/c), per symmetry-unique CP (1 = n/b/r/c,
  !> 4 = CP id, 8 = Wyckoff position), or per cell CP (2,3 = CP id). No
  !> rows for species, atomic number, or molecule labels. All rows
  !> start hidden.
  module subroutine label_style_reset_cps(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sys_ready, ok_system
    use tools_io, only: string
    class(label_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i
    character*1, allocatable :: wyc(:)

    d%ncp = 0
    if (allocated(d%cptyp)) deallocate(d%cptyp)
    if (allocated(d%cpshown)) deallocate(d%cpshown)
    if (allocated(d%cpstr)) deallocate(d%cpstr)
    d%timelastreset_cp = glfwGetTime()
    d%cpfield = r%labels%ifield
    if (.not.ok_system(r%id,sys_ready)) return
    if (.not.field_has_cps(r%id,d%cpfield)) return

    associate(c => sys(r%id)%c, f => sys(r%id)%f(d%cpfield))
      ! the rows, their CP types, and the text
      select case(r%labels%type)
      case (0)
         d%ncp = 4
      case (2,3)
         d%ncp = f%ncpcel - c%ncel
      case (1,4,8)
         d%ncp = f%ncp - c%nneq
      end select
      allocate(d%cptyp(d%ncp),d%cpshown(d%ncp),d%cpstr(d%ncp))
      d%cpshown = .false.
      do i = 1, d%ncp
         select case(r%labels%type)
         case (0)
            d%cptyp(i) = i - 1
            d%cpstr(i) = cps_letter(i-1)
         case (2,3)
            d%cptyp(i) = f%cp(f%cpcel(c%ncel+i)%idx)%typind
            d%cpstr(i) = string(c%ncel+i)
         case (1,4,8)
            d%cptyp(i) = f%cp(c%nneq+i)%typind
            if (r%labels%type == 4) then
               d%cpstr(i) = string(c%nneq+i)
            elseif (d%cptyp(i) >= 0 .and. d%cptyp(i) <= 3) then
               d%cpstr(i) = cps_letter(d%cptyp(i))
            else
               d%cpstr(i) = "?"
            end if
         end select
      end do

      ! Wyckoff positions of the symmetry-unique CPs
      if (r%labels%type == 8) then
         call cp_wyckoff(r%id,d%cpfield,wyc)
         if (allocated(wyc)) then
            do i = 1, d%ncp
               d%cpstr(i) = string(f%cp(c%nneq+i)%mult) // wyc(i)
            end do
         end if
      end if
    end associate

  end subroutine label_style_reset_cps

  !> Deallocate all arrays and end the label syle.
  module subroutine label_style_end(d)
    class(label_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%str)) deallocate(d%str)
    d%timelastreset_cp = 0d0
    d%cpfield = -1
    d%ncp = 0
    if (allocated(d%cptyp)) deallocate(d%cptyp)
    if (allocated(d%cpshown)) deallocate(d%cpshown)
    if (allocated(d%cpstr)) deallocate(d%cpstr)

  end subroutine label_style_end

  !> Allocate a coordination-polyhedra style with room for ntype center
  !> types and nspc species acting as corners.
  module subroutine coordpoly_style_alloc(d,ntype,nspc)
    use interfaces_glfw, only: glfwGetTime
    class(coordpoly_geom_style), intent(inout) :: d
    integer, intent(in) :: ntype
    integer, intent(in) :: nspc

    call d%end()
    d%timelastreset = glfwGetTime()
    d%ntype = ntype
    allocate(d%shown(ntype),d%corner(nspc,ntype),d%dmin(ntype),d%dmax(ntype))
    d%shown = .false.
    d%corner = .false.
    d%dmin = 0d0
    d%dmax = 0d0
    d%isinit = .true.

  end subroutine coordpoly_style_alloc

  !> Reset the coordination-polyhedra style to defaults from the
  !> system pointed at by representation r.
  module subroutine coordpoly_style_reset(d,r)
    use interfaces_glfw, only: glfwGetTime
    use systems, only: sys, sysc, sys_ready, ok_system, atlisttype_species
    use param, only: atmcov0
    use global, only: bondfactor_def
    class(coordpoly_geom_style), intent(inout) :: d
    type(representation), intent(in) :: r

    integer :: i, j, ispc, iz, jz, nspc
    real*8 :: dd
    logical, allocatable :: spccenter(:), spccorner(:)

    ! if not initialized, set type
    if (.not.d%isinit) d%type = atlisttype_species

    ! reset the style to zero and reset the time
    call d%end()
    d%timelastreset = glfwGetTime()

    ! check the system is sane
    if (.not.ok_system(r%id,sys_ready)) return

    ! which species act as polyhedron centers and which as corners
    call coordpoly_classify_species(r%id,spccenter,spccorner)

    ! fill the style, use the rcov sum times bondfactor as the distance cutoff
    nspc = sys(r%id)%c%nspc
    call d%alloc(sysc(r%id)%attype_number(d%type),nspc)
    do i = 1, d%ntype
       ispc = sysc(r%id)%attype_species(d%type,i)
       iz = sys(r%id)%c%spc(ispc)%z
       d%shown(i) = spccenter(ispc)
       do j = 1, nspc
          d%corner(j,i) = spccorner(j)
          if (d%corner(j,i)) then
             jz = sys(r%id)%c%spc(j)%z
             dd = (atmcov0(iz) + atmcov0(jz)) * bondfactor_def
             d%dmax(i) = max(d%dmax(i),dd)
          end if
       end do
    end do

  end subroutine coordpoly_style_reset

  !> Classify the species of system isys into the ones that can sit at
  !> the center of a coordination polyhedron (spccenter) and the ones
  !> that can act as its corners (spccorner); both are allocated to the
  !> number of species.
  module subroutine coordpoly_classify_species(isys,spccenter,spccorner)
    use systems, only: sys, sys_ready, ok_system
    use param, only: atmeneg, maxzat
    integer, intent(in) :: isys
    logical, allocatable, intent(inout) :: spccenter(:)
    logical, allocatable, intent(inout) :: spccorner(:)

    ! typical-anion species used as default corners: N, O, F, S, Cl, Br, I
    integer, parameter :: zanion(7) = (/7,8,9,16,17,35,53/)

    integer :: j, jz, nspc, navg, nvalid
    real*8 :: avgeneg
    logical :: useavg

    if (allocated(spccenter)) deallocate(spccenter)
    if (allocated(spccorner)) deallocate(spccorner)
    if (.not.ok_system(isys,sys_ready)) return

    ! count the number of species that can act as centers or corners
    nspc = sys(isys)%c%nspc
    allocate(spccenter(nspc),spccorner(nspc))
    spccenter = .false.
    spccorner = .false.
    nvalid = 0
    do j = 1, nspc
       jz = sys(isys)%c%spc(j)%z
       if (jz > 0 .and. jz <= maxzat) nvalid = nvalid + 1
    end do

    if (nvalid == 1) then
       ! single species: it acts as both center and corner
       do j = 1, nspc
          jz = sys(isys)%c%spc(j)%z
          if (jz > 0 .and. jz <= maxzat) then
             spccenter(j) = .true.
             spccorner(j) = .true.
          end if
       end do
    else
       ! If typical species that act as anions are present, use them
       ! as corners. Otherwise, use the electronegativity average to
       ! determine centers and corners.
       useavg = .true.
       do j = 1, nspc
          jz = sys(isys)%c%spc(j)%z
          if (jz > 0 .and. jz <= maxzat .and. any(zanion == jz)) then
             useavg = .false.
             exit
          end if
       end do

       if (.not.useavg) then
          ! corners are the anions, centers are everyone else
          do j = 1, nspc
             jz = sys(isys)%c%spc(j)%z
             if (jz > 0 .and. jz <= maxzat) then
                spccorner(j) = any(zanion == jz)
                spccenter(j) = .not.spccorner(j)
             end if
          end do
          ! all present species are anions: fall back to average
          if (.not.any(spccenter)) then
             useavg = .true.
             spccenter = .false.
             spccorner = .false.
          end if
       end if

       if (useavg) then
          ! average electronegativity over species with a defined value
          avgeneg = 0d0
          navg = 0
          do j = 1, nspc
             jz = sys(isys)%c%spc(j)%z
             if (jz > 0 .and. jz <= maxzat) then
                if (atmeneg(jz) > 0d0) then
                   avgeneg = avgeneg + atmeneg(jz)
                   navg = navg + 1
                end if
             end if
          end do
          if (navg > 0) avgeneg = avgeneg / navg
          do j = 1, nspc
             jz = sys(isys)%c%spc(j)%z
             if (jz > 0 .and. jz <= maxzat) then
                if (atmeneg(jz) > 0d0) then
                   spccenter(j) = atmeneg(jz) < avgeneg
                   spccorner(j) = atmeneg(jz) > avgeneg
                end if
             end if
          end do
       end if
    end if

  end subroutine coordpoly_classify_species

  !> Deallocate all arrays and end the coordination-polyhedra style.
  module subroutine coordpoly_style_end(d)
    class(coordpoly_geom_style), intent(inout) :: d

    d%isinit = .false.
    d%timelastreset = 0d0
    d%ntype = 0
    if (allocated(d%shown)) deallocate(d%shown)
    if (allocated(d%corner)) deallocate(d%corner)
    if (allocated(d%dmin)) deallocate(d%dmin)
    if (allocated(d%dmax)) deallocate(d%dmax)

  end subroutine coordpoly_style_end

end submodule styles
