! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Routines for reporting or exporting crystal/molecular structures
submodule (crystalmod) write
  implicit none

  !xx! private procedures

contains

  !> Write information about the crystal structure to the output. lcrys =
  !> information about the structure. lq = charges.
  module subroutine struct_report(c,lcrys,lq)
    use global, only: iunitname0, dunit0, iunit
    use tools_math, only: gcd
    use tools_io, only: uout, string, ioj_center, ioj_left, ioj_right
    use param, only: bohrtoa, maxzat, pi, atmass, pcamu, bohrtocm
    class(crystal), intent(inout) :: c
    logical, intent(in) :: lcrys
    logical, intent(in) :: lq

    integer :: i, j, k, iz, is
    integer :: nelec
    real*8 :: maxdv, xcm(3), x0(3), xlen(3), xang(3), xred(3,3)
    real*8 :: dens, mass, rnn2
    character(len=:), allocatable :: str1
    integer, allocatable :: nis(:)
    integer :: izp0

    character*1, parameter :: lvecname(3) = (/"a","b","c"/)

    if (lcrys) then
       ! Header
       if (.not.c%ismolecule) then
          write (uout,'("* Crystal structure")')
          write (uout,'("  From: ",A)') string(c%file)
          write (uout,'("  Lattice parameters (bohr): ",3(A,"  "))') &
             string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
          write (uout,'("  Lattice parameters (ang): ",3(A,"  "))') &
             string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
          write (uout,'("  Lattice angles (degrees): ",3(A,"  "))') &
             string(c%bb(1),'f',decimal=3), string(c%bb(2),'f',decimal=3), string(c%bb(3),'f',decimal=3)
       else
          write (uout,'("* Molecular structure")')
          write (uout,'("  From: ",A)') string(c%file)
          write (uout,'("  Encompassing cell dimensions (bohr): ",3(A,"  "))') &
             string(c%aa(1),'f',decimal=6), string(c%aa(2),'f',decimal=6), string(c%aa(3),'f',decimal=6)
          write (uout,'("  Encompassing cell dimensions (ang): ",3(A,"  "))') &
             string(c%aa(1)*bohrtoa,'f',decimal=6), string(c%aa(2)*bohrtoa,'f',decimal=6), string(c%aa(3)*bohrtoa,'f',decimal=6)
       endif

       ! Compute unit formula, and z
       allocate(nis(c%nspc))
       nis = 0
       do i = 1, c%nneq
          nis(c%at(i)%is) = nis(c%at(i)%is) + c%at(i)%mult
       end do
       maxdv = gcd(nis,c%nspc)
       write (uout,'("  Empirical formula: ",999(10(A,"(",A,") ")))') &
          (string(c%spc(i)%name), string(nint(nis(i)/maxdv)), i=1,c%nspc)
       deallocate(nis)
       if (.not.c%ismolecule) then
          write (uout,'("  Number of non-equivalent atoms in the unit cell: ",A)') string(c%nneq)
          write (uout,'("  Number of atoms in the unit cell: ",A)') string(c%ncel)
       else
          write (uout,'("  Number of atoms: ",A)') string(c%ncel)
       endif
       write (uout,'("  Number of atomic species: ",A)') string(c%nspc)
       nelec = 0
       mass = 0d0
       do i = 1, c%nneq
          iz = c%spc(c%at(i)%is)%z
          if (iz >= maxzat .or. iz <= 0) cycle
          nelec = nelec + iz * c%at(i)%mult
          mass = mass + atmass(iz) * c%at(i)%mult
       end do
       write (uout,'("  Number of electrons (with zero atomic charge): ",A)') string(nelec)
       if (.not.c%ismolecule) then
          write (uout,'("  Molar mass (amu, per unit cell): ",A)') string(mass,'f',decimal=3)
       else
          write (uout,'("  Molar mass (amu): ",A)') string(mass,'f',decimal=3)
       end if

       ! Cell volume and density, space group short report
       if (.not.c%ismolecule) then
          dens = (mass*pcamu) / (c%omega*bohrtocm**3)
          write (uout,'("  Cell volume (bohr^3): ",A)') string(c%omega,'f',decimal=5)
          write (uout,'("  Cell volume (ang^3): ",A)') string(c%omega * bohrtoa**3,'f',decimal=5)
          write (uout,'("  Density (g/cm^3): ",A)') string(dens,'f',decimal=5)
          ! space group, very short report
          if (c%havesym > 0 .and. c%spgavail) then
             write(uout,'("  Space group (H-M): ",A, " (",A,")")') &
                string(c%spg%international_symbol), string(c%spg%spacegroup_number)
          else
             write(uout,'("  Space group (H-M): ---")')
          end if
       end if

       write (uout,*)
    end if

    if (lq) then
       write (uout,'("+ List of atomic species: ")')
       write (uout,'("# spc = atomic species. Z = atomic number. name = atomic name (symbol).")')
       write (uout,'("# Q = charge.")')
       write (uout,'("# ",99(A," "))') string("spc",3,ioj_center), &
          string("Z",3,ioj_center), string("name",7,ioj_center),&
          string("Q",length=7,justify=ioj_center)
       do i = 1, c%nspc
          write (uout,'("  ",99(A," "))') string(i,3,ioj_center), &
             string(c%spc(i)%z,3,ioj_center), string(c%spc(i)%name,7,ioj_center),&
             string(c%spc(i)%qat,'f',length=7,decimal=4,justify=ioj_right)
       end do
       write (uout,*)
    end if

    if (lcrys) then
       ! List of atoms in crystallographic coordinates
       if (.not.c%ismolecule) then
          write (uout,'("+ List of non-equivalent atoms in the unit cell (cryst. coords.): ")')
          write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
          write (uout,'("# wyck = wyckoff position. name = atomic name (symbol). mult = multiplicity.")')
          write (uout,'("# Z = atomic number.")')

          write (uout,'("# ",99(A," "))') string("nat",3,ioj_center), &
             string("x",14,ioj_center), string("y",14,ioj_center),&
             string("z",14,ioj_center), string("spc",3,ioj_center), string("wyck",4,ioj_center), &
             string("name",7,ioj_center), string("mult",4,ioj_center), string("Z",3,ioj_center)
          do i=1, c%nneq
             is = c%at(i)%is
             write (uout,'("  ",99(A," "))') string(i,3,ioj_center),&
                (string(c%at(i)%x(j),'f',length=14,decimal=10,justify=3),j=1,3),&
                string(is,3,ioj_center), string(c%at(i)%mult,3,ioj_right) // c%at(i)%wyc, &
                string(c%spc(is)%name,7,ioj_center), &
                string(c%at(i)%mult,4,ioj_center), string(c%spc(is)%z,3,ioj_center)
          enddo
          write (uout,*)

          write (uout,'("+ List of atoms in the unit cell (cryst. coords.): ")')
          write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
          write (uout,'("# name = atomic name (symbol). Z = atomic number. nat = non-equivalent atom id.")')
          if (allocated(c%mol) .and. c%nmol > 0) then
             write (uout,'("# mol = molecular fragment.")')
             str1 = string("mol",3,ioj_center)
          else
             str1 = ""
          end if
          write (uout,'("# ",99(A," "))') string("at",3,ioj_center),&
             string("x",14,ioj_center), string("y",14,ioj_center),&
             string("z",14,ioj_center), string("spc",3,ioj_center), string("name",7,ioj_center),&
             string("Z",3,ioj_center), string("nat",3,ioj_center), str1
          do i=1,c%ncel
             is = c%atcel(i)%is
             if (allocated(c%mol) .and. c%nmol > 0) then
                str1 = string(c%idatcelmol(1,i),3,ioj_center)
             else
                str1 = ""
             end if
             write (uout,'("  ",99(A," "))') &
                string(i,3,ioj_center),&
                string(c%atcel(i)%x(1),'f',length=14,decimal=10,justify=3),&
                string(c%atcel(i)%x(2),'f',length=14,decimal=10,justify=3),&
                string(c%atcel(i)%x(3),'f',length=14,decimal=10,justify=3),&
                string(is,3,ioj_center),&
                string(c%spc(is)%name,7,ioj_center),&
                string(c%spc(is)%z,3,ioj_center),&
                string(c%atcel(i)%idx,3,ioj_center), str1
          enddo
          write (uout,*)

          write (uout,'("+ Lattice vectors (",A,")")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",A,": ",3(A," "))') lvecname(i),&
                (string(c%m_x2c(j,i)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,*)
       end if

       ! List of atoms in Cartesian coordinates
       write (uout,'("+ List of atoms in Cartesian coordinates (",A,"): ")') iunitname0(iunit)
       write (uout,'("# at = complete list atomic ID. xyz = Cartesian coordinates. spc = atomic species.")')
       write (uout,'("# name = atomic name (symbol). Z = atomic number. dnn = nearest-neighbor distance.")')
       if (allocated(c%mol) .and. c%nmol > 0) then
          write (uout,'("# nat = non-equivalent atom id. mol = molecular fragment.")')
          str1 = string("mol",3,ioj_center)
       else
          write (uout,'("# nat = non-equivalent atom id. ")')
          str1 = ""
       end if
       write (uout,'("# ",99(A," "))') string("at",3,ioj_center), &
          string("x",16,ioj_center), string("y",16,ioj_center),&
          string("z",16,ioj_center), string("spc",3,ioj_center), string("name",7,ioj_center),&
          string("Z",3,ioj_center), string("dnn",10,ioj_center), string("nat",3,ioj_center), str1
       do i=1,c%ncel
          is = c%atcel(i)%is
          if (allocated(c%mol) .and. c%nmol > 0) then
             str1 = string(c%idatcelmol(1,i),3,ioj_center)
          else
             str1 = ""
          end if
          rnn2 = c%get_rnn2(c%atcel(i)%idx)
          write (uout,'("  ",99(A," "))') &
             string(i,3,ioj_center),&
             (string((c%atcel(i)%r(j)+c%molx0(j))*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3),&
             string(is,3,ioj_center),string(c%spc(is)%name,7,ioj_center), string(c%spc(is)%z,3,ioj_center),&
             string(2d0*rnn2*dunit0(iunit),'f',length=10,decimal=4,justify=4), &
             string(c%atcel(i)%idx,3,ioj_center), str1
       enddo
       write (uout,*)

       ! Encompassing region for the molecule
       if (c%ismolecule) then
          write (uout,'("+ Limits of the molecular cell (in fractions of the unit cell).")')
          write (uout,'("# The part of the unit cell outside the molecular cell represents")')
          write (uout,'("# infinity (no CPs or gradient paths in it).")')
          write (uout,'("  x-axis: ",A," -> ",A)') trim(string(c%molborder(1),'f',10,4)), trim(string(1d0-c%molborder(1),'f',10,4))
          write (uout,'("  y-axis: ",A," -> ",A)') trim(string(c%molborder(2),'f',10,4)), trim(string(1d0-c%molborder(2),'f',10,4))
          write (uout,'("  z-axis: ",A," -> ",A)') trim(string(c%molborder(3),'f',10,4)), trim(string(1d0-c%molborder(3),'f',10,4))
          write (uout,*)
       end if

       ! Write symmetry operations
       if (.not.c%ismolecule) then
          write(uout,'("+ List of symmetry operations (",A,"):")') string(c%neqv)
          do k = 1, c%neqv
             write (uout,'("  Operation ",A,":")') string(k)
             write (uout,'(2("    ",4(A," ")/),"    ",4(A," "))') &
                ((string(c%rotm(i,j,k),'f',length=9,decimal=6,justify=3), j = 1, 4), i = 1, 3)
          enddo
          write (uout,*)

          call c%struct_report_symxyz(doaxes=.true.)

          write(uout,'("+ List of centering vectors (",A,"):")') string(c%ncv)
          do k = 1, c%ncv
             write (uout,'("  Vector ",A,": ",3(A," "))') string(k), &
                (string(c%cen(i,k),'f',length=9,decimal=6), i = 1, 3)
          enddo
          write (uout,*)

          call c%struct_report_symmetry()

          write (uout,'("+ Cartesian/crystallographic coordinate transformation matrices:")')
          write (uout,'("  A = car to crys (xcrys = A * xcar, ",A,"^-1)")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",3(A," "))') (string(c%m_c2x(i,j)/dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,'("  B = crys to car (xcar = B * xcrys, ",A,")")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",3(A," "))') (string(c%m_x2c(i,j)*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,'("  G = metric tensor (B''*B, ",A,"^2)")') iunitname0(iunit)
          do i = 1, 3
             write (uout,'("    ",3(A," "))') (string(c%gtensor(i,j)*dunit0(iunit)**2,'f',length=16,decimal=10,justify=5),j=1,3)
          end do
          write (uout,*)
       end if

       ! Discrete molecules, if available
       if (allocated(c%nstar) .and. allocated(c%mol) .and. c%nmol > 0) then
          write (uout,'("+ List of fragments in the system (",A,")")') string(c%nmol)
          write (uout,'("# Id = fragment ID. nat = number of atoms in fragment. C-o-m = center of mass (",A,").")')&
             iunitname0(iunit)
          write (uout,'("# Discrete = is this fragment finite?")')

          if (.not.c%ismolecule.and.c%ismol3d.and.allocated(c%idxmol)) then
             ! the table with the molecular equivalence integers
             write (uout,'("# idx = 0 (molecule is symmetry-unique), n>0 (equivalent to molecule n), -1 (symmetric to self)")')
             write (uout,'("# Id  nat           Center of mass            Discrete idx")')
             do i = 1, c%nmol
                if (c%ismolecule) then
                   xcm = (c%mol(i)%cmass()+c%molx0) * dunit0(iunit)
                else
                   xcm = c%c2x(c%mol(i)%cmass())
                end if
                write (uout,'(99("  ",A))') string(i,3,ioj_left), string(c%mol(i)%nat,4,ioj_left),&
                   (string(xcm(j),'f',10,6,3),j=1,3), string(c%mol(i)%discrete), string(c%idxmol(i),3,ioj_right)
             end do
          else
             ! the simple table with just the list of fragments
             write (uout,'("# Id  nat           Center of mass            Discrete")')
             do i = 1, c%nmol
                if (c%ismolecule) then
                   xcm = (c%mol(i)%cmass()+c%molx0) * dunit0(iunit)
                else
                   xcm = c%c2x(c%mol(i)%cmass())
                end if
                write (uout,'(99("  ",A))') string(i,3,ioj_left), string(c%mol(i)%nat,4,ioj_left),&
                   (string(xcm(j),'f',10,6,3),j=1,3), string(c%mol(i)%discrete)
             end do
          end if

          if (.not.c%ismolecule) then
             if (c%ismol3d .or. c%nlvac == 3) then
                write (uout,'(/"+ This is a molecular crystal.")')
                write (uout,'("  Number of molecules per cell (Z) = ",A)') string(c%nmol)
                izp0 = 0
                do i = 1, c%nmol
                   if (c%idxmol(i) < 0) then
                      izp0 = -1
                      exit
                   elseif (c%idxmol(i) == 0) then
                      izp0 = izp0 + 1
                   end if
                end do
                if (izp0 > 0) then
                   write (uout,'("  Number of molecules in the asymmetric unit (Z'') = ",A)') string(izp0)
                else
                   write (uout,'("  Number of molecules in the asymmetric unit (Z'') < 1")')
                end if

             else if (c%nlvac == 2) then
                write (uout,'(/"+ This is a 1D periodic (polymer) structure.")')
                write (uout,'("  Vacuum lattice vectors: (",2(A," "),A,"), (",2(A," "),A,")")') &
                   (string(c%lvac(j,1)),j=1,3), (string(c%lvac(j,2)),j=1,3)
                write (uout,'("  Connected lattice vectors: (",2(A," "),A,")")') &
                   (string(c%lcon(j,1)),j=1,3)
             else if (c%nlvac == 1) then
                write (uout,'(/"+ This is a 2D periodic (layered) structure.")')
                write (uout,'("  Vacuum lattice vectors: (",2(A," "),A,")")') &
                   (string(c%lvac(j,1)),j=1,3)
                write (uout,'("  Connected lattice vectors: (",2(A," "),A,"), (",2(A," "),A,")")') &
                   (string(c%lcon(j,1)),j=1,3), (string(c%lcon(j,2)),j=1,3)
             else
                write (uout,'(/"+ This is a 3D periodic structure.")')
             end if
          end if
          write (uout,*)
       end if

       ! Wigner-Seitz cell
       if (.not.c%ismolecule) then
          write (uout,'("+ Vertex of the WS cell in cryst. coords. (",A,")")') string(c%ws_nv)
          write (uout,'("# id = vertex ID. xyz = vertex cryst. coords. d = vertex distance to origin (",A,").")') iunitname0(iunit)
          write (uout,'(5("  ",A))') string("id",length=3,justify=ioj_right),&
             string("x",length=11,justify=ioj_center),&
             string("y",length=11,justify=ioj_center),&
             string("z",length=11,justify=ioj_center),&
             string("d ("//iunitname0(iunit)//")",length=14,justify=ioj_center)
          do i = 1, c%ws_nv
             x0 = c%x2c(c%ws_x(:,i))
             write (uout,'(5("  ",A))') string(i,length=3,justify=ioj_right), &
                (string(c%ws_x(j,i),'f',length=11,decimal=6,justify=4),j=1,3), &
                string(norm2(x0)*dunit0(iunit),'f',length=14,decimal=8,justify=4)
          enddo
          write (uout,*)

          write (uout,'("+ Faces of the WS cell (",A,")")') string(c%ws_nf)
          write (uout,'("# Face ID: vertexID1 vertexID2 ...")')
          do i = 1, c%ws_nf
             write (uout,'("  ",A,": ",999(A," "))') string(i,length=2,justify=ioj_right), &
                (string(c%ws_iside(j,i),length=2),j=1,c%ws_nside(i))
          end do
          write (uout,*)

          write (uout,'("+ Lattice vectors for the Wigner-Seitz neighbors")')
          write (uout,'("# FaceID: Voronoi lattice vector (cryst. coords.)")')
          do i = 1, c%ws_nf
             write (uout,'("  ",A,": ",99(A," "))') string(i,length=2,justify=ioj_right), &
                (string(c%ws_ineighx(j,i),length=2,justify=ioj_right),j=1,3)
          end do
          write (uout,*)

          write (uout,'("+ Lattice vectors for the Delaunay reduced cell (cryst. coords.)")')
          do i = 1, 3
             write (uout,'("  ",A,": ",99(A," "))') lvecname(i), &
                (string(nint(c%m_xr2x(j,i)),length=2,justify=ioj_right),j=1,3)
          end do

          do i = 1, 3
             x0 = c%m_xr2x(:,i)
             xred(:,i) = c%x2c(x0)
             xlen(i) = norm2(xred(:,i))
          end do
          xang(1) = acos(dot_product(xred(:,2),xred(:,3)) / xlen(2) / xlen(3)) * 180d0 / pi
          xang(2) = acos(dot_product(xred(:,1),xred(:,3)) / xlen(1) / xlen(3)) * 180d0 / pi
          xang(3) = acos(dot_product(xred(:,1),xred(:,2)) / xlen(1) / xlen(2)) * 180d0 / pi

          write (uout,'("  Delaunay reduced cell lengths: ",99(A," "))') &
             (string(xlen(j),'f',decimal=6,justify=ioj_right),j=1,3)
          write (uout,'("  Delaunay reduced cell angles: ",99(A," "))') &
             (string(xang(j),'f',decimal=3,justify=ioj_right),j=1,3)
          write (uout,*)

          write (uout,'("+ Is the cell orthogonal? ",L1)') c%isortho
          write (uout,'("+ Is the reduced cell orthogonal? ",L1/)') c%isortho_del
       end if
    end if

  end subroutine struct_report

  !> Write information about the crystal symmetry.
  module subroutine struct_report_symmetry(c)
    use tools_io, only: string, uout
    class(crystal), intent(in) :: c

    character(len=3) :: schpg
    integer :: holo, laue

    if (.not.c%ismolecule .and. c%havesym > 0) then
       write(uout,'("+ Crystal symmetry information")')
       if (c%spgavail) then
          write(uout,'("  Space group (Hermann-Mauguin): ",A, " (number ",A,")")') &
             string(c%spg%international_symbol), string(c%spg%spacegroup_number)
          write(uout,'("  Point group (Hermann-Mauguin): ",A)') string(c%spg%pointgroup_symbol)

          call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
          write(uout,'("  Point group (Schoenflies): ",A)') string(schpg)
          write(uout,'("  Holohedry: ",A)') string(holo_string(holo))
          write(uout,'("  Laue class: ",A)') string(laue_string(laue))
       else
          write(uout,'("  -- not available --")')
       end if
       write (uout,*)
    end if

  end subroutine struct_report_symmetry

  !> Write the list of symmetry operations to stdout, using
  !> crystallographic notation (if possible). If strfin is present,
  !> return the strings in that variable instead of writing them to
  !> uout.
  module subroutine struct_report_symxyz(c,strfin,doaxes)
    use tools_math, only: eig, det3
    use tools_io, only: uout, string, ioj_right
    use global, only: symprec
    use param, only: mlen, pi
    class(crystal), intent(in) :: c
    character(len=mlen), intent(out), optional :: strfin(c%neqv*c%ncv)
    logical, intent(in), optional :: doaxes

    real*8, parameter :: rfrac(25) = (/-12d0/12d0,-11d0/12d0,-10d0/12d0,&
       -9d0/12d0,-8d0/12d0,-7d0/12d0,-6d0/12d0,-5d0/12d0,-4d0/12d0,-3d0/12d0,&
       -2d0/12d0,-1d0/12d0,0d0/12d0,1d0/12d0,2d0/12d0,3d0/12d0,4d0/12d0,&
       5d0/12d0,6d0/12d0,7d0/12d0,8d0/12d0,9d0/12d0,10d0/12d0,11d0/12d0,12d0/12d0/)
    character*6, parameter :: sfrac(25) = (/"      ","-11/12","-5/6  ",&
       "-3/4  ","-2/3  ","-7/12 ","-1/2  ","-5/12 ","-1/3  ","-1/4  ","-1/6  ",&
       "-1/12 ","      ","1/12  ","1/6   ","1/4   ","1/3   ","5/12  ","1/2   ",&
       "7/12  ","2/3   ","3/4   ","5/6   ","11/12 ","      "/)
    character*1, parameter :: xyz(3) = (/"x","y","z"/)
    real*8, parameter :: eps = 1d-5

    logical :: ok, iszero, doax
    integer :: i1, i2, i, j, k, idx, rotnum
    character(len=mlen) :: strout(c%neqv*c%ncv)
    real*8 :: xtrans, rmat(3,3), eval(3), evali(3), rotaxis(3)
    real*8 :: trace, det, ang, ridx
    character(len=mlen), allocatable :: rotchar(:)

    doax = .false.
    if (present(doaxes)) doax = doaxes

    i = 0
    do i1 = 1, c%ncv
       do i2 = 1, c%neqv
          i = i + 1
          strout(i) = "<not found>"
       end do
    end do

    ! classify the rotations
    if (doax) then
       allocate(rotchar(c%neqv))
       do i2 = 1, c%neqv
          rmat = c%rotm(:,1:3,i2)
          trace = rmat(1,1)+rmat(2,2)+rmat(3,3)
          det = det3(rmat)

          if (abs(trace - 3d0) < 1d-5) then
             rotchar(i2) = " 1"
          elseif (abs(trace + 3d0) < 1d-5) then
             rotchar(i2) = "-1"
          else
             ! determine the angle of rotation
             ang = 0.5d0*(trace-det)
             if (abs(ang) > 1d0) ang = sign(1d0,ang)
             ang = acos(ang)
             if (abs(ang) < eps) then
                rotnum = 1
             else
                rotnum = nint((2d0*pi) / ang)
             end if
             if (det > 0d0) then
                rotchar(i2) = string(rotnum,2,ioj_right)
             else
                if (rotnum == 2) then
                   rotchar(i2) = " m"
                else
                   rotchar(i2) = string(-rotnum,2,ioj_right)
                end if
             end if

             ! determine the axis of rotation
             call eig(rmat,3,eval,evali)
             idx = 0
             do j = 1, 3
                if (abs(evali(j)) < eps .and. abs(eval(j)-det) < eps) then
                   idx = j
                   exit
                end if
             end do
             if (idx > 0) then
                rotaxis = rmat(:,idx)
                ! divide by the smallest non-zero element in the vector
                ridx = 1d40
                do j = 1, 3
                   if (abs(rotaxis(j)) > eps .and. abs(rotaxis(j)) < ridx) ridx = abs(rotaxis(j))
                end do
                rotaxis = rotaxis / ridx
             endif

             if (all(abs(rotaxis - nint(rotaxis)) < eps)) then
                rotchar(i2) = rotchar(i2)(1:2) // " [" // string(nint(rotaxis(1))) // "," //&
                   string(nint(rotaxis(2))) // "," // string(nint(rotaxis(3))) // "]"
             else
                rotchar(i2) = rotchar(i2)(1:2) // " [" // string(rotaxis(1),'f',decimal=1) // "," //&
                   string(rotaxis(2),'f',decimal=1) // "," // string(rotaxis(3),'f',decimal=1) // "]"
             end if

             rotaxis = c%x2c(rotaxis)
             rotaxis = rotaxis / norm2(rotaxis)
             rotchar(i2) = trim(rotchar(i2)) // "; [" // string(rotaxis(1),'f',decimal=3) // "," //&
                string(rotaxis(2),'f',decimal=3) // "," // string(rotaxis(3),'f',decimal=3) // "]"
          end if
       end do
    end if

    i = 0
    main: do i1 = 1, c%ncv
       loopi: do i2 = 1, c%neqv
          i = i + 1
          strout(i) = ""

          do j = 1, 3
             ! translation
             ok = .false.
             do k = 1, 25
                xtrans = c%rotm(j,4,i2)+c%cen(j,i1) - rfrac(k)
                xtrans = xtrans - nint(xtrans)
                if (abs(xtrans) < symprec) then
                   ok = .true.
                   strout(i) = trim(strout(i)) // sfrac(k)
                   iszero = (k == 13) .or. (k == 1) .or. (k == 25)
                   exit
                end if
             end do
             if (.not.ok) then
                strout(i) = "<not found>"
                cycle loopi
             end if

             ! rotation
             do k = 1, 3
                if (abs(c%rotm(j,k,i2) - 1d0) < symprec) then
                   if (iszero) then
                      strout(i) = trim(strout(i)) // xyz(k)
                   else
                      strout(i) = trim(strout(i)) // "+" // xyz(k)
                   end if
                   iszero = .false.
                elseif (abs(c%rotm(j,k,i2) + 1d0) < symprec) then
                   strout(i) = trim(strout(i)) // "-" // xyz(k)
                   iszero = .false.
                elseif (abs(c%rotm(j,k,i2)) > symprec) then
                   strout(i) = "<not found>"
                   cycle loopi
                end if
             end do

             ! the comma
             if (j < 3) &
                strout(i) = trim(strout(i)) // ","
          end do
          if (doax) then
             strout(i) = string(strout(i),30) // " ## " // trim(rotchar(i2))
          end if
       end do loopi
    end do main
    if (doax) deallocate(rotchar)

    if (present(strfin)) then
       strfin = strout
    else
       write(uout,'("+ List of symmetry operations in crystallographic notation:")')
       if (doax) &
          write(uout,'("# number: operation ... ## order of rotation [axis]; [axis in Cartesian]")')
       do k = 1, c%neqv*c%ncv
          write (uout,'("   ",A,": ",A)') string(k,length=3,justify=ioj_right), string(strout(k))
       enddo
       write (uout,*)
    end if

  end subroutine struct_report_symxyz

  !> Write the field info to a JSON object. The structure object to
  !> json with root p.
  module subroutine struct_write_json(c,json,p)
    use json_module, only: json_value, json_core
    class(crystal), intent(inout) :: c
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(inout) :: p

    type(json_value), pointer :: s, ap, arr

    integer :: i
    character(len=mlen), allocatable :: strfin(:)
    real*8 :: rnn2

    if (.not.c%isinit) return
    call json%create_object(s,'structure')
    call json%add(p,s)

    call json%add(s,'cell_lengths',c%aa)
    call json%add(s,'cell_angles',c%bb)
    call json%add(s,'cell_volume',c%omega)
    call json%add(s,'crys_to_cart_matrix',reshape(c%m_x2c,(/9/)))
    call json%add(s,'cart_to_crys_matrix',reshape(c%m_c2x,(/9/)))

    call json%add(s,'is_molecule',c%ismolecule)
    call json%add(s,'molecule_centering_vector',c%molx0)
    call json%add(s,'molecular_cell_border',c%molborder)
    call json%add(s,'periodicity',3-c%nlvac)

    call json%add(s,'number_of_species',c%nspc)
    call json%create_array(arr,'species')
    call json%add(s,arr)
    do i = 1, c%nspc
       call json%create_object(ap,'')
       call json%add(ap,'id',i)
       call json%add(ap,'name',trim(c%spc(i)%name))
       call json%add(ap,'atomic_number',c%spc(i)%z)
       call json%add(arr,ap)
       nullify(ap)
    end do
    nullify(arr)

    call json%add(s,'number_of_nonequivalent_atoms',c%nneq)
    call json%create_array(arr,'nonequivalent_atoms')
    call json%add(s,arr)
    do i = 1, c%nneq
       call json%create_object(ap,'')
       call json%add(ap,'id',i)
       call json%add(ap,'species',c%at(i)%is)
       call json%add(ap,'fractional_coordinates',c%at(i)%x(:))
       call json%add(ap,'cartesian_coordinates',c%at(i)%r(:))
       call json%add(ap,'multiplicity',c%at(i)%mult)
       call json%add(ap,'wyckoff_letter',c%at(i)%wyc)
       rnn2 = c%get_rnn2(i)
       call json%add(ap,'half_nn_distance',rnn2)
       call json%add(arr,ap)
       nullify(ap)
    end do
    nullify(arr)

    call json%add(s,'number_of_cell_atoms',c%ncel)
    call json%create_array(arr,'cell_atoms')
    call json%add(s,arr)
    do i = 1, c%ncel
       call json%create_object(ap,'')
       call json%add(ap,'id',i)
       call json%add(ap,'species',c%atcel(i)%is)
       call json%add(ap,'fractional_coordinates',c%atcel(i)%x(:))
       call json%add(ap,'cartesian_coordinates',c%atcel(i)%r(:))
       call json%add(ap,'nonequivalent_id',c%atcel(i)%idx)
       call json%add(ap,'symop_to_nneq',c%atcel(i)%ir)
       call json%add(ap,'centering_vector_to_nneq',c%atcel(i)%ic)
       call json%add(ap,'lattice_vector_to_nneq',c%atcel(i)%lvec)
       call json%add(arr,ap)
       nullify(ap)
    end do
    nullify(arr)

    if (.not.c%ismolecule) then
       call json%add(s,'have_symmetry',c%havesym > 0)
       if (c%spgavail) then
          call json%add(s,'space_group_hm',trim(c%spg%international_symbol))
          call json%add(s,'space_group_ita_number',c%spg%spacegroup_number)
       end if
       call json%add(s,'number_of_symops',c%neqv)

       allocate(strfin(c%neqv*c%ncv))
       call c%struct_report_symxyz(strfin)
       call json%create_array(arr,'symops')
       call json%add(s,arr)
       do i = 1, c%neqv
          call json%create_object(ap,'')
          call json%add(ap,'id',i)
          call json%add(ap,'operation',trim(strfin(i)))
          call json%add(ap,'rotation',reshape(c%rotm(1:3,1:3,i),(/9/)))
          call json%add(ap,'translation',c%rotm(:,4,i))
          call json%add(arr,ap)
          nullify(ap)
       end do
       nullify(arr)

       call json%add(s,'number_of_centering_vectors',c%ncv)
       call json%add(s,'centering_vectors',reshape(c%cen(1:3,1:c%ncv),(/3*c%ncv/)))
    endif
    call json%add(s,'number_of_molecular_fragments',c%nmol)
    nullify(s)

  end subroutine struct_write_json

  !> Write the structure to a file. Use the format derived from the
  !> extension of file and use default values for all options.
  module subroutine write_simple_driver(c,file,ti)
    use tools_io, only: lower, ferror, faterr, equal
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: wext, wext2, wroot
    integer :: idx

    wext = lower(file(index(file,'.',.true.)+1:))
    wroot = file(:index(file,'.',.true.)-1)

    if (equal(wext,'xyz').or.equal(wext,'gjf').or.equal(wext,'cml')) then
       call c%write_mol(file,wext,ti=ti)
    else if(equal(wext,'obj').or.equal(wext,'ply').or.equal(wext,'off')) then
       call c%write_3dmodel(file,wext,ti=ti)
    elseif (equal(wext,'gau')) then
       call c%write_gaussian(file,ti=ti)
    elseif (equal(wext,'in')) then
       idx = index(wroot,'.',.true.)
       if (idx > 0) then
          wext2 = lower(wroot(idx+1:))
          if (equal(wext2,'scf')) then
             call c%write_espresso(file,ti=ti)
          else
             idx = 0
          end if
       end if
       if (idx == 0) &
            call c%write_fhi(file,.true.,ti=ti)
    elseif (equal(wext,'pwi')) then
       call c%write_espresso(file,ti=ti)
    elseif (equal(wext,'poscar') .or. equal(wext,'contcar')) then
       call c%write_vasp(file,.false.,ti=ti)
    elseif (equal(wext,'abin')) then
       call c%write_abinit(file,ti=ti)
    elseif (equal(wext,'elk')) then
       call c%write_elk(file,ti=ti)
    elseif (equal(wext,'tess')) then
       call c%write_tessel(file,ti=ti)
    elseif (equal(wext,'incritic').or.equal(wext,'cri')) then
       call c%write_critic(file,ti=ti)
    elseif (equal(wext,'cif')) then
       call c%write_cif(file,.true.,ti=ti)
    elseif (equal(wext,'d12').or.equal(wext,'34')) then
       call c%write_d12(file,.true.,.false.,ti=ti)
    elseif (equal(wext,'res')) then
       call c%write_res(file,-1,ti=ti)
    elseif (equal(wext,'m')) then
       call c%write_escher(file,ti=ti)
    elseif (equal(wext,'db')) then
       call c%write_db(file,ti=ti)
    elseif (equal(wext,'gin')) then
       call c%write_gulp(file,ti=ti)
    elseif (equal(wext,'lammps')) then
       call c%write_lammps(file,ti=ti)
    elseif (equal(wext,'fdf')) then
       call c%write_siesta_fdf(file,ti=ti)
    elseif (equal(wext,'struct_in')) then
       call c%write_siesta_in(file,ti=ti)
    elseif (equal(wext,'hsd')) then
       call c%write_dftbp_hsd(file,ti=ti)
    elseif (equal(wext,'gen')) then
       call c%write_dftbp_gen(file,ti=ti)
    elseif (equal(wext,'pyscf')) then
       call c%write_pyscf(file,ti=ti)
    elseif (equal(wext,'fhi')) then
       call c%write_fhi(file,.true.,ti=ti)
    elseif (equal(wext,'frac')) then
       call c%write_tinkerfrac(file,ti=ti)
    elseif (equal(wext,'pdb')) then
       call c%write_pdb(file,ti=ti)
    else
       call ferror('struct_write','unrecognized file format',faterr)
    end if

  end subroutine write_simple_driver

  !> Write a xyz/gjf/cml file containing a finite piece of the crystal
  !> structure. fmt can be one of xyz, gjf, or cml. ix is the number
  !> of unit cells to plot.  If doborder is .true., add all atoms at
  !> the border. If onemotif is .true., write all molecules in the
  !> unit cell.  If molmotif is .true., complete molecules with atoms
  !> in adjacent cells. If docell, add sticks for the unit cell
  !> limits. If environ is true, write all molecules up to a distance
  !> renv (bohr) from the origin. If lnmer, partition the resulting
  !> list of molecules into n-mers, up to nth order. If rsph (bohr) is
  !> positive, then use all atoms in a sphere around xsph (cryst.). If
  !> rcub (bohr) is positive, use all atoms in a cube around xcub
  !> (cryst.). If luout is present, return the LU in that argument
  !> and do not close the file.
  module subroutine write_mol(c,file,fmt,ix0,doborder0,onemotif0,molmotif0,&
     environ0,renv0,lnmer0,nmer0,rsph0,xsph0,rcub0,xcub0,usenames0,luout,ti)
    use global, only: dunit0, iunit
    use tools_math, only: nchoosek, comb
    use tools_io, only: ferror, faterr, uout, string, ioj_left, string, ioj_right,&
       equal
    use types, only: realloc
    use fragmentmod, only: fragment, realloc_fragment
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in), optional :: ix0(3)
    logical, intent(in), optional :: doborder0, onemotif0, molmotif0, environ0
    real*8, intent(in), optional :: renv0
    logical, intent(in), optional :: lnmer0
    integer, intent(in), optional :: nmer0
    real*8, intent(in), optional :: rsph0, xsph0(3)
    real*8, intent(in), optional :: rcub0, xcub0(3)
    logical, intent(in), optional :: usenames0
    integer, intent(out), optional :: luout
    type(thread_info), intent(in), optional :: ti

    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: i, j, k, l, m, nmol, icel, lvec(3), ncm
    integer :: ncomb, nlimi, nlimj, icount
    integer, allocatable :: icomb(:), origmol(:)
    character(len=:), allocatable :: wroot, file0, aux
    logical :: doagain
    real*8, allocatable :: cmlist(:,:)
    real*8 :: xcm(3), dist
    integer :: ix(3)
    logical :: doborder, onemotif, molmotif, environ
    real*8 :: renv
    logical :: lnmer, usenames
    integer :: nmer
    real*8 :: rsph, xsph(3)
    real*8 :: rcub, xcub(3)

    ! set the default options
    ix = 1
    doborder = .false.
    onemotif = .false.
    molmotif = .false.
    environ = .false.
    renv = 0d0
    lnmer = .false.
    nmer = 1
    rsph = -1d0
    xsph = 0d0
    rcub = -1d0
    xcub = 0d0
    usenames = .false.
    if (present(ix0)) ix = ix0
    if (present(doborder0)) doborder = doborder0
    if (present(onemotif0)) onemotif = onemotif0
    if (present(molmotif0)) molmotif = molmotif0
    if (present(environ0)) environ = environ0
    if (present(renv0)) renv = renv0
    if (present(lnmer0)) lnmer = lnmer0
    if (present(nmer0)) nmer = nmer0
    if (present(rsph0)) rsph = rsph0
    if (present(xsph0)) xsph = xsph0
    if (present(rcub0)) rcub = rcub0
    if (present(xcub0)) xcub = xcub0
    if (present(usenames0)) usenames = usenames0

    ! determine the fragments
    if (onemotif) then
       call fr%merge_array(c%mol(1:c%nmol),.false.)
       allocate(fr0(c%nmol))
       fr0 = c%mol
       nmol = c%nmol
    elseif (environ) then
       ! calculate the centers of mass for all fragments in the molecular motif
       allocate(cmlist(3,c%nmol),fr0(c%nmol),origmol(c%nmol))
       nmol = c%nmol
       fr0 = c%mol
       do i = 1, c%nmol
          cmlist(:,i) = c%mol(i)%cmass()
          origmol(i) = i
       end do
       ncm = c%nmol

       doagain = .true.
       icel = 0
       do while(doagain)
          doagain = .false.
          icel = icel + 1
          do k = 1, 6
             if (k == 1 .or. k == 2) then
                nlimi = icel
                nlimj = icel
             elseif (k == 3 .or. k == 4) then
                nlimi = icel
                nlimj = icel-1
             else
                nlimi = icel-1
                nlimj = icel-1
             end if
             do i = -nlimi, nlimi
                do j = -nlimj, nlimj
                   lvec = icelcomb(k,i,j,icel)
                   do l = 1, c%nmol
                      xcm = c%c2x(cmlist(:,l)) + lvec
                      xcm = c%x2c(xcm)
                      dist = norm2(xcm)
                      if (dist <= renv) then
                         ncm = ncm + 1
                         if (ncm > size(cmlist,2)) then
                            call realloc(cmlist,3,2*ncm)
                            call realloc(origmol,2*ncm)
                            call realloc_fragment(fr0,2*ncm)
                         end if
                         cmlist(:,ncm) = xcm
                         fr0(ncm) = c%mol(l)
                         do m = 1, fr0(ncm)%nat
                            fr0(ncm)%at(m)%x = fr0(ncm)%at(m)%x + lvec
                            fr0(ncm)%at(m)%r = c%x2c(fr0(ncm)%at(m)%x)
                            fr0(ncm)%at(m)%lvec = fr0(ncm)%at(m)%lvec + lvec
                         end do
                         origmol(ncm) = l
                         doagain = .true.
                      end if
                   end do
                end do
             end do
          end do
       end do

       call realloc_fragment(fr0,ncm)
       call fr%merge_array(fr0,.false.)
       nmol = ncm
    else
       if (rcub > 0) then
          fr = c%listatoms_sphcub(rcub=rcub,xcub=xcub)
       elseif (rsph > 0) then
          fr = c%listatoms_sphcub(rsph=rsph,xsph=xsph)
       else
          fr = c%listatoms_cells(ix,doborder)
       endif
       if (molmotif) then
          call c%listmolecules(fr,nmol,fr0,isdiscrete)
          call fr%merge_array(fr0,.false.)
       end if
    end if

    if (lnmer .and..not.allocated(fr0)) &
       call ferror('write_mol','ONEMOTIF, MOLMOTIF, or ENVIRON are necessary with NMER',faterr)

    ! If environ, report the identities of all the molecules in the environment
    if (environ) then
       write (uout,'("+ List of fragments in the molecular environment")')
       write (uout,'("  Number of fragments: ",A)') string(nmol)
       write (uout,'("  Fragment number Id with nat atoms at center-of-mass comes from")')
       write (uout,'("  from fragment idmol in the Wigner-Seitz cell translated by")')
       write (uout,'("  lattice vector lvec.")')
       write (uout,'("# Id nat               Center of mass          idmol     lvec")')
       do i = 1, nmol
          if (c%ismolecule) then
             xcm = cmlist(:,i) * dunit0(iunit)
          else
             xcm = c%c2x(cmlist(:,i))
          end if
          write (uout,'(99("  ",A))') string(i,3,ioj_left), string(fr0(i)%nat,4,ioj_left),&
             (string(xcm(j),'f',10,6,3),j=1,3), string(origmol(i)), &
             (string(nint(fr0(i)%at(1)%x(j)-c%mol(origmol(i))%at(1)%x(j)),3,ioj_right),j=1,3)
       end do
       write (uout,*)
    end if

    ! if this is a molecule, translate to the proper origin
    if (c%ismolecule) then
       do i = 1, fr%nat
          fr%at(i)%r = fr%at(i)%r + c%molx0
       end do
       if (allocated(fr0)) then
          do i = 1, nmol
             do j = 1, fr0(i)%nat
                fr0(i)%at(j)%r = fr0(i)%at(j)%r + c%molx0
             end do
          end do
       end if
    end if

    if (.not.lnmer) then
       ! normal write
       call dowrite(file,fr,ti=ti)
    else
       wroot = file(:index(file,'.',.true.)-1)
       do i = 1, nmer
          if (i == 1) then
             if (nmer == 1) then
                nlimj = nmol
             else
                nlimj = c%nmol
             end if
             ! monomers
             do j = 1, nlimj
                file0 = trim(wroot) // "_" // string(j) // "." // fmt
                call dowrite(file0,fr0(j),ti=ti)
             end do
             write (uout,'("+ Written ",A," ",A,"-mers")') string(c%nmol), string(i)
          elseif (i == nmer) then
             ! n-mers
             allocate(icomb(i-1))
             icount = 0
             do l = 1, c%nmol
                ncomb = nchoosek(nmol,i-1)
                do j = 1, ncomb
                   call comb(nmol,i-1,j,icomb)
                   if (any(icomb == l)) cycle
                   file0 = trim(wroot) // "_" // string(l)
                   fr = fr0(l)
                   do k = 1, i-1
                      aux = trim(file0) // "_" // string(icomb(k))
                      file0 = aux
                      call fr%append(fr0(icomb(k)))
                   end do
                   aux = trim(file0) // "." // fmt
                   file0 = aux
                   call dowrite(file0,fr,ti=ti)
                   icount = icount + 1
                end do
             end do
             deallocate(icomb)
             write (uout,'("+ Written ",A," ",A,"-mers")') string(icount), string(i)
          else
             ! everything in between
             ncomb = nchoosek(nmol,i)
             allocate(icomb(i))
             icount = 0
             do j = 1, ncomb
                call comb(nmol,i,j,icomb)
                file0 = wroot
                call fr%init()
                do k = 1, i
                   aux = trim(file0) // "_" // string(icomb(k))
                   file0 = aux
                   call fr%append(fr0(icomb(k)))
                end do
                aux = trim(file0) // "." // fmt
                file0 = aux
                call dowrite(file0,fr,ti=ti)
                icount = icount + 1
             end do
             deallocate(icomb)
             write (uout,'("+ Written ",A," ",A,"-mers")') string(icount), string(i)
          end if
       end do
    end if

  contains
    subroutine dowrite(fileo,fro,ti)
      character*(*) :: fileo
      type(fragment) :: fro
      type(thread_info), intent(in), optional :: ti

      if (equal(fmt,"xyz")) then
         call fro%writexyz(fileo,usenames,ti=ti)
      elseif (equal(fmt,"gjf")) then
         call fro%writegjf(fileo,ti=ti)
      elseif (equal(fmt,"cml")) then
         if (c%ismolecule) then
            call fro%writecml(fileo,luout=luout,ti=ti)
         else
            call fro%writecml(fileo,c%m_x2c,luout=luout,ti=ti)
         end if
      else
         call ferror("write_mol","Unknown format",faterr)
      endif

    end subroutine dowrite
    function icelcomb(idx,i,j,icel)
      integer, intent(in) :: idx, i, j, icel
      integer :: icelcomb(3)

      icelcomb = 0
      if (idx == 1) then
         icelcomb = (/i,j,icel/)
      elseif (idx == 2) then
         icelcomb = (/i,j,-icel/)
      elseif (idx == 3) then
         icelcomb = (/i,icel,j/)
      elseif (idx == 4) then
         icelcomb = (/i,-icel,j/)
      elseif (idx == 5) then
         icelcomb = (/icel,i,j/)
      elseif (idx == 6) then
         icelcomb = (/-icel,i,j/)
      endif

    end function icelcomb
  end subroutine write_mol

  !> Write an obj/ply/off file containing the crystal structure. fmt
  !> can be one of obj, ply, or off. ix is the number of unit cells to
  !> plot. If doborder is .true., add all atoms at the border. If
  !> onemotif is .true., write all molecules in the unit cell. If
  !> molmotif is .true., complete molecules with atoms in adjacent
  !> cells. If docell, add sticks for the unit cell limits. If
  !> domolcell, add sticks tfor the molecular cell. If rsph (bohr) is
  !> positive, then use all atoms in a sphere around xsph (cryst.). If
  !> rcub (bohr) is positive, use all atoms in a cube around xcub
  !> (cryst.). If gr0 is present, then return the graphics handle and
  !> do not close the files.
  module subroutine write_3dmodel(c,file,fmt,ix0,doborder0,onemotif0,molmotif0,&
     docell0,domolcell0,rsph0,xsph0,rcub0,xcub0,gr0,ti)
    use graphics, only: grhandle
    use fragmentmod, only: fragment
    use tools_io, only: equal
    use param, only: maxzat, atmcov, jmlcol
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in), optional :: ix0(3)
    logical, intent(in), optional :: doborder0, onemotif0, molmotif0
    logical, intent(in), optional :: docell0, domolcell0
    real*8, intent(in), optional :: rsph0, xsph0(3)
    real*8, intent(in), optional :: rcub0, xcub0(3)
    type(grhandle), intent(out), optional :: gr0
    type(thread_info), intent(in), optional :: ti

    integer :: i, j
    real*8 :: d, xd(3), x0(3), x1(3), rr
    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: nmol
    type(grhandle) :: gr
    integer :: ix(3)
    logical :: doborder, onemotif, molmotif
    logical :: docell, domolcell
    real*8 :: rsph, xsph(3)
    real*8 :: rcub, xcub(3)

    real*8, parameter :: rfac = 1.4d0
    real*8, parameter :: x0cell(3,2,12) = reshape((/&
       0d0, 0d0, 0d0,   0d0, 0d0, 1d0,&
       0d0, 0d0, 0d0,   0d0, 1d0, 0d0,&
       0d0, 0d0, 0d0,   1d0, 0d0, 0d0,&
       1d0, 0d0, 0d0,   1d0, 1d0, 0d0,&
       1d0, 1d0, 0d0,   0d0, 1d0, 0d0,&
       0d0, 1d0, 0d0,   0d0, 1d0, 1d0,&
       0d0, 1d0, 1d0,   0d0, 0d0, 1d0,&
       0d0, 0d0, 1d0,   1d0, 0d0, 1d0,&
       1d0, 0d0, 1d0,   1d0, 0d0, 0d0,&
       1d0, 1d0, 1d0,   0d0, 1d0, 1d0,&
       1d0, 1d0, 1d0,   1d0, 0d0, 1d0,&
       1d0, 1d0, 1d0,   1d0, 1d0, 0d0/),shape(x0cell))

    ! set default values
    ix = 1
    doborder = .false.
    onemotif = .false.
    molmotif = .false.
    docell = .false.
    domolcell = .false.
    rsph = -1d0
    xsph = 0d0
    rcub = -1d0
    xcub = 0d0
    if (present(ix0)) ix = ix0
    if (present(doborder0)) doborder = doborder0
    if (present(onemotif0)) onemotif = onemotif0
    if (present(molmotif0)) molmotif = molmotif0
    if (present(docell0)) docell = docell0
    if (present(domolcell0)) domolcell = domolcell0
    if (present(rsph0)) rsph = rsph0
    if (present(xsph0)) xsph = xsph0
    if (present(rcub0)) rcub = rcub0
    if (present(xcub0)) xcub = xcub0

    ! open and get the atom list
    if (onemotif) then
       call fr%merge_array(c%mol(1:c%nmol),.false.)
    else
       if (rcub > 0) then
          fr = c%listatoms_sphcub(rcub=rcub,xcub=xcub)
       elseif (rsph > 0) then
          fr = c%listatoms_sphcub(rsph=rsph,xsph=xsph)
       else
          fr = c%listatoms_cells(ix,doborder)
       endif
       if (molmotif) then
          call c%listmolecules(fr,nmol,fr0,isdiscrete)
          call fr%merge_array(fr0,.false.)
       end if
    end if

    ! if this is a molecule, translate to the proper origin
    if (c%ismolecule) then
       do i = 1, fr%nat
          fr%at(i)%r = fr%at(i)%r + c%molx0
       end do
    end if

    call gr%open(fmt,file,ti=ti)

    ! add the balls
    do i = 1, fr%nat
       if (fr%spc(fr%at(i)%is)%z > maxzat) then
          rr = 0.21d0
       else
          rr = 0.6d0*atmcov(fr%spc(fr%at(i)%is)%z)
       endif
       call gr%ball(fr%at(i)%r,JMLcol(:,fr%spc(fr%at(i)%is)%z),rr)
    end do

    ! add the sticks
    do i = 1, fr%nat
       do j = i+1, fr%nat
          if (fr%spc(fr%at(i)%is)%z > maxzat .or. fr%spc(fr%at(j)%is)%z > maxzat) cycle
          xd = fr%at(i)%r - fr%at(j)%r
          d = norm2(xd)
          if (d < (atmcov(fr%spc(fr%at(i)%is)%z) + atmcov(fr%spc(fr%at(j)%is)%z)) * rfac) then
             xd = fr%at(i)%r + 0.5d0 * (fr%at(j)%r - fr%at(i)%r)
             call gr%stick(fr%at(i)%r,xd,JMLcol(:,fr%spc(fr%at(i)%is)%z),0.05d0)
             call gr%stick(fr%at(j)%r,xd,JMLcol(:,fr%spc(fr%at(j)%is)%z),0.05d0)
          end if
       end do
    end do

    ! add the cell
    if (docell) then
       do i = 1, 12
          x0 = c%x2c(x0cell(:,1,i)) + c%molx0
          x1 = c%x2c(x0cell(:,2,i)) + c%molx0
          call gr%stick(x0,x1,(/255,0,0/),0.03d0)
       end do
    end if

    ! add the molecular cell
    if (domolcell .and. c%ismolecule) then
       do i = 1, 12
          x0 = x0cell(:,1,i)
          x1 = x0cell(:,2,i)
          do j = 1, 3
             if (abs(x0(j)) < 1d-12) x0(j) = c%molborder(j)
             if (abs(x0(j)-1d0) < 1d-12) x0(j) = 1d0-c%molborder(j)
             if (abs(x1(j)) < 1d-12) x1(j) = c%molborder(j)
             if (abs(x1(j)-1d0) < 1d-12) x1(j) = 1d0-c%molborder(j)
          end do
          x0 = c%x2c(x0) + c%molx0
          x1 = c%x2c(x1) + c%molx0
          call gr%stick(x0,x1,(/0,0,255/),0.03d0)
       end do
    end if

    ! close or give the handles to the calling routine, cleanup
    if (present(gr0)) then
       gr0 = gr
    else
       call gr%close(ti=ti)
    end if

    if (allocated(fr0)) deallocate(fr0)

  end subroutine write_3dmodel

  !> Write a quantum espresso input template
  module subroutine write_espresso(c,file,rklength,ti)
    use tools_io, only: fopen_write, lower, fclose, string
    use param, only: atmass
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    real*8, intent(in), optional :: rklength
    type(thread_info), intent(in), optional :: ti

    integer :: i, lu, nk(3)
    real*8 :: rk

    rk = 40d0
    if (present(rklength)) rk = rklength

    lu = fopen_write(file,ti=ti)
    write (lu,'("&control")')
    write (lu,'(" title=''crystal'',")')
    write (lu,'(" prefix=''crystal'',")')
    write (lu,'(" pseudo_dir=''.'',")')
    write (lu,'(" calculation=''vc-relax'',")')
    write (lu,'("/")')
    write (lu,'("&system")')
    write (lu,'(" ibrav=0,")')
    write (lu,'(" nat=",A,",")') string(c%ncel)
    write (lu,'(" ntyp=",A,",")') string(c%nspc)
    write (lu,'(" ecutwfc=60.0,")')
    write (lu,'(" ecutrho=600.0,")')
    write (lu,'(" xdm=.true.,")')
    write (lu,'("/")')
    write (lu,'("&electrons"/" conv_thr = 1d-8,"/"/")')
    write (lu,'("&ions"/"/")')
    write (lu,'("&cell"/"/")')
    write (lu,'("ATOMIC_SPECIES")')
    do i = 1, c%nspc
       write (lu,'(A," ",F12.6," ",A,".UPF")') trim(c%spc(i)%name), atmass(c%spc(i)%z), trim(lower(c%spc(i)%name))
    end do
    write (lu,'(/"ATOMIC_POSITIONS crystal")')
    do i = 1, c%ncel
       write (lu,'(A,3(" ",F13.8," "))') trim(c%spc(c%atcel(i)%is)%name), c%atcel(i)%x
    end do

    call c%get_kpoints(rk,nk)
    write (lu,'(/"K_POINTS automatic"/3(A," ")" 1 1 1"/)') (string(nk(i)),i=1,3)

    write (lu,'("CELL_PARAMETERS bohr")')
    do i = 1, 3
       write (lu,'(3(F18.12," "))') c%m_x2c(:,i)
    end do
    call fclose(lu)

  end subroutine write_espresso

  !> Write a VASP POSCAR file. If verbose, write the atom sequence
  !> to the output. If append, append to an existing file.
  module subroutine write_vasp(c,file,verbose,append,ti)
    use tools_io, only: fopen_write, fopen_append, string, uout, fclose, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: verbose
    logical, intent(in), optional :: append
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: lbl1, lbl2, aux, auxname
    integer :: i, j, lu, ntyp
    logical :: append_, ok

    ! whether to append
    append_ = .false.
    if (present(append)) append_ = append

    ! Cell
    if (append_) then
       inquire(file=file,exist=ok)
       if (ok) then
          lu = fopen_append(file,ti=ti)
       else
          lu = fopen_write(file,ti=ti)
       end if
    else
       lu = fopen_write(file,ti=ti)
    end if
    write (lu,'("critic2 | ",A)') string(c%file)
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(3(F15.10," "))') c%m_x2c(:,i) * bohrtoa
    end do

    ! Number of atoms per type and Direct
    lbl1 = ""
    lbl2 = ""
    do i = 1, c%nspc
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       aux = lbl1 // " " // string(ntyp)
       lbl1 = aux
       auxname = lbl2 // " " // trim(nameguess(c%spc(i)%z,.true.))
       lbl2 = auxname
    end do
    write (lu,'(A)') lbl2
    write (lu,'(A)') lbl1
    write (lu,'("Direct")')

    ! Atomic positions
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'(3(F13.8," "))') c%atcel(j)%x
          end if
       end do
    end do
    call fclose(lu)

    if (verbose) &
       write (uout,'("+ Atom type sequence: ",999(A," "))') (string(c%spc(j)%name),j=1,c%nspc)

  end subroutine write_vasp

  !> Write an abinit input template
  module subroutine write_abinit(c,file,ti)
    use tools_io, only: fopen_write, string, fclose
    use param, only: pi
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: lbl1, aux
    integer :: ntyp
    real*8 :: aap(3), bbp(3), gpq(3,3)
    integer :: i, j, lu

    ! Find the lengths and angles of the cell
    gpq = matmul(transpose(c%m_x2c),c%m_x2c)
    do i = 1, 3
       aap(i) = sqrt(gpq(i,i))
    end do
    bbp(1)=acos(gpq(2,3)/sqrt(gpq(2,2)*gpq(3,3)))/pi*180d0
    bbp(2)=acos(gpq(1,3)/sqrt(gpq(1,1)*gpq(3,3)))/pi*180d0
    bbp(3)=acos(gpq(1,2)/sqrt(gpq(1,1)*gpq(2,2)))/pi*180d0

    ! Write input
    lu = fopen_write(file,ti=ti)
    write (lu,'("acell ",3(F14.10," "))') aap
    write (lu,'("angdeg ",3(F14.10," "))') bbp
    write (lu,'("ntypat ",I3)') c%nspc

    lbl1 = ""
    do i = 1, c%nspc
       aux = lbl1 // " " // string(c%spc(i)%z)
       lbl1 = aux
    end do
    write (lu,'("znucl ",A)') lbl1
    write (lu,'("natom ",I5)') c%ncel

    lbl1 = ""
    do i = 1, c%nspc
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       aux = lbl1 // " " // string(ntyp) // "*" // string(i)
       lbl1 = aux
    end do
    write (lu,'("typat ",A)') lbl1

    write (lu,'("xred ")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'(" ",3(F15.10," "))') c%atcel(j)%x
          end if
       end do
    end do
    write (lu,*)

    ! template for the rest
    write (lu,'("#Definition of the planewave basis set")')
    write (lu,'("ecut 15")')
    write (lu,'("")')
    write (lu,'("# k-grid")')
    write (lu,'("kptopt 1")')
    write (lu,'("nshiftk 4")')
    write (lu,'("shiftk  0.5 0.5 0.5  ")')
    write (lu,'("        0.5 0.0 0.0")')
    write (lu,'("        0.0 0.5 0.0")')
    write (lu,'("        0.0 0.0 0.5")')
    write (lu,'("ngkpt   1 1 1")')
    write (lu,'("")')
    write (lu,'("#Definition of the SCF procedure")')
    write (lu,'("nstep  20      ")')
    write (lu,'("toldfe 1.0d-10")')
    write (lu,'("diemac 12.0   ")')
    call fclose(lu)

  end subroutine write_abinit

  !> Write an elk input template
  module subroutine write_elk(c,file,ti)
    use tools_io, only: fopen_write, fclose
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: ntyp
    integer :: i, j, lu

    ! Write input
    lu = fopen_write(file,ti=ti)
    write (lu,'("tasks"/,"0"/)')
    write (lu,'("xctype"/,"20"/)')
    write (lu,'("avec")')
    do i = 1, 3
       write (lu,'("  ",3(F15.10," "))') c%m_x2c(:,i)
    end do
    write (lu,*)

    write (lu,'("sppath"/,"''./''"/)')

    write (lu,'("atoms")')
    write (lu,'("  ",I4)') c%nspc
    do i = 1, c%nspc
       write (lu,'("  ","''",A,".in''")') trim(c%spc(i)%name)
       ntyp = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) ntyp = ntyp + 1
       end do
       write (lu,'("  ",I3)') ntyp
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) then
             write (lu,'("  ",3(F14.10," "),"0.0 0.0 0.0")') c%atcel(j)%x
          end if
       end do
    end do
    write (lu,*)

    write (lu,'("ngridk"/,"  4 4 4"/)')
    write (lu,'("rgkmax"/,"  7.0"/)')
    write (lu,'("highq"/,"  .true."/)')
    call fclose(lu)

  end subroutine write_elk

  !> Write a Gaussian template input (periodic).
  module subroutine write_gaussian(c,file,ti)
    use tools_io, only: fopen_write, string, nameguess, ioj_left, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: wroot
    integer :: lu, i, j

    wroot = file(:index(file,'.',.true.)-1)

    lu = fopen_write(file,ti=ti)
    write (lu,'("%chk=",A,".chk")') wroot
    write (lu,'("%nprocs=8")')
    write (lu,'("%mem=2GB")')
    write (lu,'("#p pbepbe/sto-3g int(grid=ultrafine) pop=regular fmm=(print)")')
    write (lu,'("   iop1=timestamp iop(5/13=1,5/33=1,5/181=10,5/184=186)")')
    write (lu,'("   scf=(novaracc,noincfock,tight)")')
    write (lu,*)
    write (lu,'("title")')
    write (lu,*)
    write (lu,'("0 1")')
    do i = 1, c%ncel
       write (lu,'(99(A," "))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.),2,ioj_left),&
          (string(c%atcel(i)%r(j)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    do i = 1, 3
       write (lu,'(99(A," "))') string("Tv",2,ioj_left),&
          (string(c%m_x2c(j,i)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    write (lu,*)

    call fclose(lu)

  end subroutine write_gaussian

  !> Write a tessel input template
  module subroutine write_tessel(c,file,ti)
    use global, only: fileroot
    use tools_io, only: fopen_write, fclose
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i

    lu = fopen_write(file,ti=ti)

    write (lu,'("set camangle 75 -10 45")')
    write (lu,'("set background background {color rgb <1,1,1>}")')
    write (lu,'("set use_planes .false.")')
    write (lu,'("set ball_texture finish{specular 0.2 roughness 0.1 reflection 0.1}")')
    write (lu,'("set equalscale noscale")')
    write (lu,'("molecule")')
    write (lu,'("  crystal")')
    write (lu,'("    symmatrix seitz")')
    do i = 1, c%ncv
       write (lu,'("     ",A,3(F15.12," "))') "cen ",c%cen(:,i)
    end do
    write (lu,'("     ","#")')
    do i = 1, c%neqv
       write (lu,'("     ",3(F5.2," "),F15.12)') c%rotm(1,:,i)
       write (lu,'("     ",3(F5.2," "),F15.12)') c%rotm(2,:,i)
       write (lu,'("     ",3(F5.2," "),F15.12)') c%rotm(3,:,i)
       write (lu,'("     ","#")')
    end do
    write (lu,'("     ","endsymmatrix")')
    write (lu,'("     ",A,6(F12.8," "))') "cell", c%aa, c%bb
    write (lu,'("     ","crystalbox  -2.30 -2.30 -2.30 2.30 2.30 2.30")')
    write (lu,'("     ",A,6(F6.3," "))') "clippingbox ",-0.02,-0.02,-0.02,+1.02,+1.02,+1.02
    do i = 1, c%nneq
       write (lu,'("     ","neq ",3(F12.8," "),A10)') c%at(i)%x, trim(c%spc(c%at(i)%is)%name)
    end do
    write (lu,'("  endcrystal")')
    write (lu,'(A)') "  unitcell radius 0.01 rgb 1.0 0.5 0.5 many"
    write (lu,'("  molmotif allmaincell jmol")')
    write (lu,'("  off ",A,".off")') trim(adjustl(fileroot))
    write (lu,'("  vrml ",A,".wrl")') trim(adjustl(fileroot))
    write (lu,'("  povray ",A,".pov")') trim(adjustl(fileroot))
    write (lu,'("endmolecule")')
    write (lu,'("# run povray -D -UV +I",A,".pov +O",A,".png +W2000 +H2000 +A")') &
       trim(adjustl(fileroot)), trim(adjustl(fileroot))
    write (lu,'("end")')

    call fclose(lu)

  end subroutine write_tessel

  !> Write a critic2 input template
  module subroutine write_critic(c,file,ti)
    use tools_io, only: fopen_write, fclose, string
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)

    write (lu,'("crystal")')
    write (lu,'("  cell ",6(A," "))') (string(c%aa(i),'f',decimal=10),i=1,3),&
       (string(c%bb(i),'f',decimal=6),i=1,3)
    do i = 1, c%ncel
       write (lu,'("  neq ",3(A," "),A)') (string(c%atcel(i)%x(j),'f',decimal=10),j=1,3),&
          string(c%spc(c%atcel(i)%is)%name)
    end do
    write (lu,'("endcrystal")')
    write (lu,'("end")')
    call fclose(lu)

  end subroutine write_critic

  !> Write a simple cif file (filename = file) with the c crystal
  !> structure. If usesym0, write symmetry to the cif file; otherwise
  !> use P1.
  module subroutine write_cif(c,file,usesym0,ti)
    use global, only: fileroot, testing
    use tools_io, only: fopen_write, fclose, string, nameguess, deblank, nameguess,&
       ferror, faterr, ioj_left
    use param, only: bohrtoa, maxzat
    use tools_math, only: gcd
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: usesym0
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, iz, lu, idx, gcdz
    character(len=mlen), allocatable :: strfin(:)
    character*2 :: sym
    character*3 :: schpg
    character(len=:), allocatable :: str
    integer :: holo, laue, natmol
    logical :: usesym, doz
    integer :: datvalues(8)
    integer, allocatable :: atc(:,:), addlabel(:), spcuse(:)

    ! Hill order for chemical formula. First C, then H, then all the other
    ! elements in alphabetical order.
    integer, parameter :: hillord(maxzat) = &
       (/6,1,89,47,13,95,18,33,85,79,5,56,4,107,83,97,35,20,48,58,98,17,&
       96,112,27,24,55,29,105,110,66,68,99,63,9,26,114,100,87,31,64,32,2,72,80,&
       67,108,53,49,77,19,36,57,3,103,71,116,115,101,12,25,42,109,7,11,41,60,10,&
       113,28,102,93,8,118,76,15,91,82,46,61,84,59,78,94,88,37,75,104,111,45,86,&
       44,16,51,21,34,106,14,62,50,38,73,65,43,52,90,22,81,69,117,92,23,74,54,&
       39,70,30,40/)

    ! use symmetry?
    usesym = usesym0 .and. c%spgavail

    ! open output file
    lu = fopen_write(file,ti=ti)

    ! header (date and time mucks up testing)
    write (lu,'("data_",A)') string(deblank(fileroot))
    write (lu,'("_audit_creation_method ''critic2''")')
    if (.not.testing) then
       call date_and_time(values=datvalues)
       write (lu,'("_audit_creation_date ",A,"-",A,"-",A)') &
          (string(datvalues(i)),i=1,3)
    end if

    ! formula: count the number of element types
    allocate(atc(maxzat,c%nmol))
    doz = c%ismol3d
    atc = 0
    do i = 1, c%nmol
       do j = 1, c%mol(i)%nat
          if (c%mol(i)%spc(c%mol(i)%at(j)%is)%z > 0 .and. c%mol(i)%spc(c%mol(i)%at(j)%is)%z <= maxzat) then
             atc(c%mol(i)%spc(c%mol(i)%at(j)%is)%z,i) = atc(c%mol(i)%spc(c%mol(i)%at(j)%is)%z,i) + 1
          end if
       end do
       if (i > 1) then
          if (any(atc(:,i) - atc(:,1) /= 0)) then
             doz = .false.
          end if
       end if
    end do

    ! formula
    if (.not.doz) then
       ! not a molecular crystal or different types of molecules,
       ! collect all atc then calculate gdc of all non-zero
       ! numbers. atc(:,1) is now the formula unit and gcdz = Z.
       do i = 2, c%nmol
          atc(:,1) = atc(:,1) + atc(:,i)
       end do
       gcdz = -1
       do i = 1, maxzat
          if (atc(i,1) > 0) then
             if (gcdz < 0) then
                gcdz = atc(i,1)
             else
                gcdz = gcd(gcdz,atc(i,1))
             end if
          end if
       end do
       atc(:,1) = atc(:,1) / gcdz
    else
       ! a molecular crystal with always the same molecule: gcdz = Z
       natmol = sum(atc(:,1))
       if (abs(real(c%ncel,8)/real(natmol,8) - c%ncel/natmol) > 1d-10) &
          call ferror('write_cif','inconsistent number of atoms in fragment',faterr)
       gcdz = c%ncel / natmol
    end if

    ! formula: build the molecular formula/formula unit and write to cif
    str = ""
    do i = 1, maxzat
       idx = hillord(i)
       if (atc(idx,1) > 0) then
          sym = nameguess(idx,.true.)
          if (atc(idx,1) > 1) then
             str = str // trim(sym) // string(atc(idx,1)) // " "
          else
             str = str // trim(sym) //  " "
          end if
       end if
    end do
    str = trim(str)
    write (lu,'("_chemical_formula_sum ''",A,"''")') str
    write (lu,'("_cell_formula_units_Z ",A)') string(gcdz)
    deallocate(atc)

    ! cell dimensions
    write (lu,'("_cell_length_a ",F20.10)') c%aa(1)*bohrtoa
    write (lu,'("_cell_length_b ",F20.10)') c%aa(2)*bohrtoa
    write (lu,'("_cell_length_c ",F20.10)') c%aa(3)*bohrtoa
    write (lu,'("_cell_angle_alpha ",F14.4)') c%bb(1)
    write (lu,'("_cell_angle_beta ",F14.4)') c%bb(2)
    write (lu,'("_cell_angle_gamma ",F14.4)') c%bb(3)
    write (lu,'("_cell_volume ",F20.6)') c%omega * bohrtoa**3

    if (usesym) then
       allocate(strfin(c%neqv*c%ncv))
       call c%struct_report_symxyz(strfin)
       do i = 1, c%neqv*c%ncv
          if (index(strfin(i),"not found") > 0) then
             usesym = .false.
             exit
          end if
       end do
    end if

    ! write the symmetry, if applicable
    if (usesym) then
       call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
       write (lu,'("_space_group_crystal_system ",A)') string(holo_string(holo))
       write (lu,'("_space_group_IT_number ",A)') string(c%spg%spacegroup_number)
       write (lu,'("_space_group_name_H-M_alt ''",A,"''")') string(c%spg%international_symbol)
       write (lu,'("_space_group_name_Hall ''",A,"''")') string(c%spg%hall_symbol)

       write (lu,'("loop_")')
       write (lu,'("_symmetry_equiv_pos_site_id")')
       write (lu,'("_symmetry_equiv_pos_as_xyz")')
       do i = 1, c%neqv*c%ncv
          write (lu,'(" ",A," ''",A,"''")') string(i), string(strfin(i))
       end do
    else
       write (lu,'("_space_group_crystal_system triclinic")')
       write (lu,'("_space_group_IT_number 1")')
       write (lu,'("_space_group_name_H-M_alt ''P 1''")')
       write (lu,'("_space_group_name_Hall ''P 1''")')

       write (lu,'("loop_")')
       write (lu,'("_symmetry_equiv_pos_site_id")')
       write (lu,'("_symmetry_equiv_pos_as_xyz")')
       write (lu,'(" 1 ''x,y,z''")')
    end if
    if (allocated(strfin)) deallocate(strfin)

    ! calculate additional label for atom_site_label
    allocate(spcuse(c%nspc))
    spcuse = 0
    if (usesym) then
       allocate(addlabel(c%nneq))
       do i = 1, c%nneq
          spcuse(c%at(i)%is) = spcuse(c%at(i)%is) + 1
          addlabel(i) = spcuse(c%at(i)%is)
       end do
    else
       allocate(addlabel(c%ncel))
       do i = 1, c%ncel
          spcuse(c%atcel(i)%is) = spcuse(c%atcel(i)%is) + 1
          addlabel(i) = spcuse(c%atcel(i)%is)
       end do
    end if
    deallocate(spcuse)

    write (lu,'("loop_")')
    write (lu,'("_atom_site_label")')
    write (lu,'("_atom_site_type_symbol")')
    write (lu,'("_atom_site_fract_x")')
    write (lu,'("_atom_site_fract_y")')
    write (lu,'("_atom_site_fract_z")')
    if (usesym) then
       do i = 1, c%nneq
          iz = c%at(i)%is
          str = trim(c%spc(iz)%name) // string(addlabel(i))
          write (lu,'(5(A," "))') string(str,5,ioj_left),&
             string(nameguess(c%spc(iz)%z,.true.),5,ioj_left),&
             (string(c%at(i)%x(j),'f',decimal=14),j=1,3)
       end do
    else
       do i = 1, c%ncel
          iz = c%atcel(i)%is
          str = trim(c%spc(iz)%name) // string(addlabel(i))
          write (lu,'(5(A," "))') string(str,5,ioj_left),&
             string(nameguess(c%spc(iz)%z,.true.),5,ioj_left),&
             (string(c%atcel(i)%x(j),'f',decimal=14),j=1,3)
       end do
    end if
    deallocate(addlabel)
    call fclose(lu)

  end subroutine write_cif

  !> Write a simple d12 file
  module subroutine write_d12(c,file,dosym,doexternal,ti)
    use tools_io, only: fopen_write, fclose, string
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: dosym
    logical, intent(in) :: doexternal
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: file34
    character(len=3) :: schpg
    integer :: lu, holo, laue
    integer :: i, j, k, l, num, idang
    real*8 :: x(3)
    real*8 :: dum(3,3)
    logical :: need34
    type(crystal) :: caux

    need34 = .false.
    lu = fopen_write(file,ti=ti)
    write (lu,'("Title")')
    if (c%ismolecule) then
       write (lu,'("MOLECULE")')
       write (lu,'("1")')
       write (lu,'(A)') string(c%ncel)
       do i = 1, c%ncel
          write (lu,'(4(A," "))') string(c%spc(c%atcel(i)%is)%z), &
             (string((c%atcel(i)%r(j)+c%molx0(j))*bohrtoa,'f',15,8),j=1,3)
       end do
    elseif (doexternal) then
       write (lu,'("EXTERNAL")')
       need34 = .true.
    elseif (.not.dosym) then
       write (lu,'("CRYSTAL")')
       write (lu,'("0 0 0")')
       write (lu,'("1")')
       write (lu,'(6(A," "))') (string(c%aa(i)*bohrtoa,'f',15,8),i=1,3), (string(c%bb(j),'f',15,8),j=1,3)
       write (lu,'(A)') string(c%ncel)
       do i = 1, c%ncel
          write (lu,'(4(A," "))') string(c%spc(c%atcel(i)%is)%z), (string(c%atcel(i)%x(j),'f',15,8),j=1,3)
       end do
    else
       caux = c
       dum = caux%cell_standard(.false.,.false.,.true.)

       write (lu,'("CRYSTAL")')
       write (lu,'("0 0 0")')
       num = caux%spg%spacegroup_number
       write (lu,'(A)') string(num)
       if (num <= 2) then ! triclinic
          write (lu,'(6(A," "))') (string(caux%aa(i)*bohrtoa,'f',15,8),i=1,3), (string(caux%bb(j),'f',15,8),j=1,3)
       elseif (num <= 15) then ! monoclinic
          idang = maxloc(abs(caux%bb - 90d0),1)
          write (lu,'(6(A," "))') (string(caux%aa(i)*bohrtoa,'f',15,8),i=1,3), string(caux%bb(idang),'f',15,8)
       elseif (num <= 74) then ! orthorhombic
          write (lu,'(6(A," "))') (string(caux%aa(i)*bohrtoa,'f',15,8),i=1,3)
       elseif (num <= 142) then ! tetragonal
          write (lu,'(6(A," "))') string(caux%aa(1)*bohrtoa,'f',15,8), string(caux%aa(3)*bohrtoa,'f',15,8)
       elseif (num <= 167) then ! trigonal
          write (lu,'(6(A," "))') string(caux%aa(1)*bohrtoa,'f',15,8), string(caux%aa(3)*bohrtoa,'f',15,8)
       elseif (num <= 194) then ! hexagonal
          write (lu,'(6(A," "))') string(caux%aa(1)*bohrtoa,'f',15,8), string(caux%aa(3)*bohrtoa,'f',15,8)
       else ! cubic
          write (lu,'(6(A," "))') string(caux%aa(1)*bohrtoa,'f',15,8)
       end if
       write (lu,'(A)') string(caux%nneq)
       do i = 1, caux%nneq
          write (lu,'(4(A," "))') string(caux%spc(caux%at(i)%is)%z), (string(caux%at(i)%x(j),'f',15,8),j=1,3)
       end do
    end if
    write (lu,'("TESTGEOM")')
    write (lu,'("END")')
    write (lu,'("END")')
    write (lu,'("END")')
    call fclose(lu)

    if (need34) then
       file34 = file(:index(file,'.',.true.)-1) // ".fort.34"
       lu = fopen_write(file34,ti=ti)

       ! for a crystal, if symmetry is available
       ! header: dimensionality, centring type, and crystal holohedry
       call pointgroup_info(c%spg%pointgroup_symbol,schpg,holo,laue)
       write (lu,'("3 1 ",A)') string(holo-1)
       do i = 1, 3
          write (lu,'(3(A," "))') (string(c%m_x2c(j,i)*bohrtoa,'f',decimal=10),j=1,3)
       end do

       ! symmetry operations
       write (lu,'(A)') string(c%neqv*c%ncv)
       do i = 1, c%neqv
          do j = 1, c%ncv
             dum = transpose(matmul(matmul(c%m_x2c,c%rotm(1:3,1:3,i)),c%m_c2x))
             do k = 1, 3
                write (lu,'(3(" ",E19.12))') (dum(l,k),l=1,3)
             end do
             x = c%rotm(:,4,i)+c%cen(:,j)
             x = c%x2c(x)
             write (lu,'(3(" ",E19.12))') (x(l)*bohrtoa,l=1,3)
          end do
       end do

       ! atoms
       write (lu,'(A)') string(c%nneq)
       do i = 1, c%nneq
          write (lu,'(4(A," "))') string(c%spc(c%at(i)%is)%z), (string(c%at(i)%r(l)*bohrtoa,'f',decimal=10),l=1,3)
       end do

       call fclose(lu)
    end if

  end subroutine write_d12

  !> Write a shelx res file (filename = file) with the c crystal
  !> structure. If usesym0, write symmetry to the cif file; otherwise
  !> use P1. dosym = 0 (do not use symmetry), 1 (use symmetry), or
  !> -1 (use symmetry only if possible, do not emit warnings)
  module subroutine write_res(c,file,dosym,ti)
    use tools_io, only: fopen_write, fclose, string, ferror, warning, nameguess
    use tools_math, only: det3
    use param, only: bohrtoa, eye
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: dosym
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, lu, ilatt
    character(len=mlen), allocatable :: strfin(:)
    character(len=:), allocatable :: str
    logical :: usesym, ok3(3), ok2(2)
    real*8 :: dd

    real*8, parameter :: cen_i(3)  = (/0.5d0,0.5d0,0.5d0/)
    real*8, parameter :: cen_a(3)  = (/0.0d0,0.5d0,0.5d0/)
    real*8, parameter :: cen_b(3)  = (/0.5d0,0.0d0,0.5d0/)
    real*8, parameter :: cen_c(3)  = (/0.5d0,0.5d0,0.0d0/)
    real*8, parameter :: cen_r1(3) = (/2d0,1d0,1d0/) / 3d0
    real*8, parameter :: cen_r2(3) = (/1d0,2d0,2d0/) / 3d0
    real*8, parameter :: eps = 1d-5

    ! use symmetry?
    usesym = (dosym==1 .or. dosym==-1) .and. c%spgavail

10  continue

    ! open output file
    lu = fopen_write(file,ti=ti)

    ! header
    write (lu,'("TITL critic2 | ",A)') trim(c%file)
    write (lu,'("CELL 0.71073 ",6(A," "))') (string(c%aa(i)*bohrtoa,'f',12,8),i=1,3), &
       (string(c%bb(j),'f',10,6),j=1,3)
    if (usesym) then
       write (lu,'("ZERR ",A," 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001")') string(c%ncel/c%nneq)
    else
       write (lu,'("ZERR 1 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001")')
    end if
    if (usesym) then
       ! identify lattice type
       ilatt = 0
       if (c%ncv == 1) then
          ilatt = -1 ! P
       elseif (c%ncv == 2) then
          if (all(abs(c%cen(:,2) - cen_i - nint(c%cen(:,2) - cen_i)) < eps)) then
             ilatt = -2 ! I
          elseif (all(abs(c%cen(:,2) - cen_a - nint(c%cen(:,2) - cen_a)) < eps)) then
             ilatt = -5 ! A
          elseif (all(abs(c%cen(:,2) - cen_b - nint(c%cen(:,2) - cen_b)) < eps)) then
             ilatt = -6 ! B
          elseif (all(abs(c%cen(:,2) - cen_c - nint(c%cen(:,2) - cen_c)) < eps)) then
             ilatt = -7 ! C
          end if
       elseif (c%ncv == 3) then
          ok2 = .false.
          do i = 2, c%ncv
             if (.not.ok2(1)) ok2(1) = all(abs(c%cen(:,i) - cen_r1 - nint(c%cen(:,i) - cen_r1)) < eps)
             if (.not.ok2(2)) ok2(2) = all(abs(c%cen(:,i) - cen_r2 - nint(c%cen(:,i) - cen_r2)) < eps)
          end do
          if (all(ok2)) ilatt = -3
       elseif (c%ncv == 4) then
          ok3 = .false.
          do i = 2, c%ncv
             if (.not.ok3(1)) ok3(1) = all(abs(c%cen(:,i) - cen_a - nint(c%cen(:,i) - cen_a)) < eps)
             if (.not.ok3(2)) ok3(2) = all(abs(c%cen(:,i) - cen_b - nint(c%cen(:,i) - cen_b)) < eps)
             if (.not.ok3(3)) ok3(3) = all(abs(c%cen(:,i) - cen_c - nint(c%cen(:,i) - cen_c)) < eps)
          end do
          if (all(ok3)) ilatt = -4
       end if
       if (ilatt == 0) then
          if (dosym == 1) &
             call ferror('write_res','unknown set of centering vectors',warning)
          usesym = .false.
          call fclose(lu)
          goto 10
       end if

       ! identify centrosymmetry
       do i = 1, c%neqv
          if (all(abs(c%rotm(1:3,1:3,i) + eye) < eps) .and. all(abs(c%rotm(1:3,4,i)) < eps)) then
             ilatt = -ilatt
             exit
          end if
       end do
    else
       ilatt = -1
    end if
    write (lu,'("LATT ",A)') string(ilatt)

    if (usesym) then
       allocate(strfin(c%neqv*c%ncv))
       call c%struct_report_symxyz(strfin)
       do i = 2, c%neqv ! skip the identity
          dd = det3(c%rotm(1:3,1:3,i))
          if (dd > 0d0 .or. ilatt < 0) then
             if (index(strfin(i),"not found") > 0) then
                if (dosym == 1) &
                   call ferror('write_res','unknown set of centering vectors',warning)
                usesym = .false.
                call fclose(lu)
                goto 10
             end if
             write (lu,'("SYMM ",A)') trim(strfin(i))
          end if
       end do
    end if

    ! atomic species
    write (lu,'("SFAC ",999(A," "))') (trim(nameguess(c%spc(i)%z,.true.)),i=1,c%nspc)

    ! number of atoms of each type
    str = ""
    do i = 1, c%nspc
       str = str // " " // string(count(c%atcel(:)%is == i))
    end do
    write (lu,'("UNIT ",A)') str
    write (lu,'("FVAR 1.00")')

    ! list of atoms
    if (usesym) then
       do i = 1, c%nneq
          write (lu,'(999(A," "))') trim(c%spc(c%at(i)%is)%name) // string(i), string(c%at(i)%is), &
             (string(c%at(i)%x(j),'f',12,8),j=1,3), string(real(c%at(i)%mult,8)/(c%neqv*c%ncv),'f',12,8), &
             "0.05"
       end do
    else
       do i = 1, c%ncel
          write (lu,'(999(A," "))') trim(c%spc(c%atcel(i)%is)%name) // string(i), string(c%atcel(i)%is), &
             (string(c%atcel(i)%x(j),'f',12,8),j=1,3), "1.0", "0.05"
       end do
    end if

    ! close the file
    write (lu,'("END")')
    call fclose(lu)

  end subroutine write_res

  !> Write an escher octave script
  module subroutine write_escher(c,file,ti)
    use global, only: fileroot
    use tools_io, only: fopen_write, string, fclose
    use param, only: pi
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: lbl1, aux
    integer :: lu, i, n

    lu = fopen_write(file,ti=ti)

    ! count number of atoms per type
    write (lu,'("cr = struct();")')
    write (lu,'("cr.name = """,A,""";")') trim(adjustl(fileroot))
    write (lu,'("cr.a = [",1p,3(E22.14," "),"];")') c%aa
    write (lu,'("cr.b = [",1p,3(E22.14," "),"];")') c%bb * pi / 180d0
    write (lu,'("cr.nat = ",I6,";")') c%ncel
    write (lu,'("cr.ntyp = ",I6,";")') c%nspc
    write (lu,'("cr.r = [")')
    do i = 1, 3
       write (lu,'("  ",1p,3(E22.14," "))') c%m_x2c(:,i)
    end do
    write (lu,'("  ","];")')
    write (lu,'("cr.g = [")')
    do i = 1, 3
       write (lu,'("  ",1p,3(E22.14," "))') c%gtensor(:,i)
    end do
    write (lu,'("  ","];")')
    write (lu,'("cr.omega = ",1p,E22.14,";")') c%omega

    lbl1 = "cr.ztyp = ["
    do i = 1, c%nspc
       aux = lbl1 // " " // string(c%spc(i)%z)
       lbl1 = aux
    end do
    aux = lbl1 // "];"
    lbl1 = aux
    write (lu,'(A)') lbl1

    lbl1 = "cr.attyp = {"
    do i = 1, c%nspc
       if (i > 1) then
          aux = lbl1 // ","
          lbl1 = aux
       end if
       aux = lbl1 // '"' // string(c%spc(i)%name) // '"'
       lbl1 = aux
    end do
    aux = lbl1 // "};"
    lbl1 = aux
    write (lu,'(A)') lbl1

    lbl1 = "cr.typ = ["
    do i = 1, c%ncel
       aux = lbl1 // " " // string(c%atcel(i)%is)
       lbl1 = aux
    end do
    aux = lbl1 // "];"
    lbl1 = aux
    write (lu,'(A)') lbl1

    write (lu,'("cr.x = [")')
    n = 0
    do i = 1, c%ncel
       write (lu,'("  ",1p,3(E22.14," "))') c%atcel(i)%x
    end do
    write (lu,'("  ];")')

    call fclose(lu)

  end subroutine write_escher

  !> Write a db file for the dcp automatic input generator
  module subroutine write_db(c,file,ti)
    use tools_io, only: fopen_write, string, fclose, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)
    write (lu,'("type crystal_energy")')
    write (lu,'("kpts 4")')
    write (lu,'("crys")')
    write (lu,'(6(A," "))') (string(c%aa(i)*bohrtoa,'f',18,10),i=1,3), (string(c%bb(i),'f',18,10),i=1,3)
    do i = 1, c%ncel
       write (lu,'(A," ",1p,3(A," "))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.)), &
          (string(c%atcel(i)%x(j),'f',18,10),j=1,3)
    end do
    write (lu,'("end")')
    call fclose(lu)

  end subroutine write_db

  !> Write a gulp input script
  module subroutine write_gulp(c,file,ti)
    use tools_io, only: fopen_write, nameguess, fclose, string
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)
    write (lu,'("eem")')
    write (lu,'("cell ",6(A," "))') (string(c%aa(j) * bohrtoa,'f',13,9),j=1,3), &
       (string(c%bb(j),'f',10,5),j=1,3)
    write (lu,'("fractional")')
    do i = 1, c%ncel
       write (lu,'(A5," ",3(A," "))') trim(c%spc(c%atcel(i)%is)%name),&
          (string(c%atcel(i)%x(j),'f',15,9),j=1,3)
    end do

    call fclose(lu)

  end subroutine write_gulp

  !> Write a lammps data file
  module subroutine write_lammps(c,file,ti)
    use tools_io, only: fopen_write, ferror, faterr, fclose
    use tools_math, only: m_x2c_from_cellpar
    use param, only: bohrtoa, atmass
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, l
    integer :: lu
    real*8 :: rnew(3,3)

    lu = fopen_write(file,ti=ti)

    ! header
    write (lu,'("LAMMPS data file created by critic2. (experimental)",/)')
    write (lu,'(I9," atoms")') c%ncel
    write (lu,'(I9," atom types")') c%nspc
    write (lu,*)

    ! metrics of the cell --> this needs more testing
    rnew = m_x2c_from_cellpar(c%aa,c%bb)
    if (abs(c%m_x2c(1,2)) > 1d-12 .or. abs(c%m_x2c(1,3)) > 1d-12 .or.&
       abs(c%m_x2c(2,3)) > 1d-12) then
       call ferror ('write_lammps','non-orthogonal cells not implemented',faterr)
    end if
    write (lu,'(2(F18.10," ")," xlo xhi")') 0d0, c%m_x2c(1,1)*bohrtoa
    write (lu,'(2(F18.10," ")," ylo yhi")') 0d0, c%m_x2c(2,2)*bohrtoa
    write (lu,'(2(F18.10," ")," zlo zhi")') 0d0, c%m_x2c(3,3)*bohrtoa
    write (lu,'(3(F18.10," ")," xy xz yz")') 0d0, 0d0, 0d0
    write (lu,*)

    write (lu,'("Masses"/)')
    do i = 1, c%nspc
       write (lu,'(I3," ",F10.4)') i, atmass(c%spc(i)%z)
    end do
    write (lu,*)

    write (lu,'("Atoms"/)')
    l = 0
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is /= i) cycle
          l = l + 1
          write (lu,'(I7," ",I3," ",F4.1,3(F15.8," "))') l, i, 0d0, c%atcel(j)%r*bohrtoa
       end do
    end do

    call fclose(lu)

  end subroutine write_lammps

  !> Write a siesta fdf data file
  module subroutine write_siesta_fdf(c,file,ti)
    use tools_io, only: fopen_write, nameguess, lower, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: i, j
    integer :: lu

    lu = fopen_write(file,ti=ti)

    ! header
    write (lu,'("# fdf file created by critic2.",/)')
    write (lu,'("SystemName crystal")')
    write (lu,'("SystemLabel crystal")')
    write (lu,*)

    write (lu,'("NumberOfSpecies ",I3)') c%nspc
    write (lu,'("NumberOfAtoms ", I6)') c%ncel
    write (lu,'("%block Chemical_Species_Label")')
    do i = 1, c%nspc
       write (lu,'(I3,I3," ",A2)') i, c%spc(i)%z, lower(nameguess(c%spc(i)%z,.true.))
    end do
    write (lu,'("%endblock Chemical_Species_Label")')
    write (lu,*)

    write (lu,'("LatticeConstant 1.0 ang")')
    write (lu,'("%block LatticeParameters")')
    write (lu,'(3(F16.10," "),3(F16.8," "))') c%aa*bohrtoa, c%bb
    write (lu,'("%endblock LatticeParameters")')
    write (lu,'("AtomicCoordinatesFormat Fractional")')
    write (lu,'("%block AtomicCoordinatesAndAtomicSpecies")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) &
             write (lu,'(3(F18.12," "),I3)') c%atcel(j)%x, i
       end do
    end do
    write (lu,'("%endblock AtomicCoordinatesAndAtomicSpecies")')
    write (lu,*)

    write (lu,'("XC.functional GGA")')
    write (lu,'("XC.authors PBE")')
    write (lu,'("SpinPolarized .false.")')
    write (lu,'("MaxSCFIterations 100")')
    write (lu,'("MeshCutoff 100. Ry")')
    write (lu,'("DM.NumberPulay 3")')
    write (lu,*)

    write (lu,'("PAO.BasisSize DZP")')
    write (lu,*)

    write (lu,'("kgrid_cutoff 10.0 ang")')
    write (lu,*)

    write (lu,'("ElectronicTemperature 5 K")')
    write (lu,*)

    write (lu,'("# options")')
    write (lu,'("LongOutput")')
    write (lu,'("SaveRho")')
    write (lu,'("SaveBaderCharge")')
    write (lu,'("DM.UseSaveDM")')
    write (lu,'("WriteDenchar")')
    write (lu,'("WriteCoorXmol")')

    call fclose(lu)

  end subroutine write_siesta_fdf

  !> Write a siesta STRUCT_IN data file
  module subroutine write_siesta_in(c,file,ti)
    use tools_io, only: fopen_write, uout, nameguess, string, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, j, k

    lu = fopen_write(file,ti=ti)

    ! lattice vectors
    r = transpose(c%m_x2c) * bohrtoa
    do i = 1, 3
       write (lu,'(3(F20.12," "))') r(i,:)
    end do

    ! atoms
    write (lu,*) c%ncel
    j = 0
    do i = 1, c%nspc
       do k = 1, c%ncel
          if (c%atcel(k)%is == i) &
             write (lu,'(I3," ",I3," ",3(F20.12," "))') i, c%spc(i)%z, c%atcel(k)%x
       end do
    end do

    call fclose(lu)

    ! Write the chemical species block to standard output
    write (uout,'("%block Chemical_Species_Label")')
    do i = 1, c%nspc
       write (uout,'(3("  ",A))') string(i), string(c%spc(i)%z), &
          string(nameguess(c%spc(i)%z,.true.))
    end do
    write (uout,'("%endblock Chemical_Species_Label")')
    write (uout,*)

  end subroutine write_siesta_in

  !> Write a DFTB+ human-friendly structured data format (hsd) file
  module subroutine write_dftbp_hsd(c,file,ti)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: maxzat0
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    real*8, parameter :: hderiv(maxzat0) = (/&
     -0.1857d0,      0.d0,      0.d0,    0.d0,      0.d0, -0.1492d0,& ! 1:6   (H-C)
     -0.1535d0, -0.1575d0, -0.1623d0,    0.d0, -0.0454d0,   -0.02d0,& ! 7:12  (N-Mg)
          0.d0,      0.d0,   -0.14d0, -0.11d0, -0.0697d0,      0.d0,& ! 13:18 (Al-Ar)
     -0.0339d0, -0.0340d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 19:24 (K-Cr)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,   -0.03d0,& ! 25:30 (Mn-Zn)
          0.d0,      0.d0,      0.d0,    0.d0, -0.0573d0,      0.d0,& ! 31:36 (Ga-Kr)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 37:42 (Rb-Mo)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 43:48 (Tc-Cd)
          0.d0,      0.d0,      0.d0,    0.d0, -0.0433d0,      0.d0,& ! 49:54 (In-Xe)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 55:60 (Cs-Nd)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 61:66 (Pm-Dy)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 67:72 (Ho-Hf)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 73:78 (Ta-Pt)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 79:84 (Au-Po)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 85:90 (At-Th)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 91:96 (Pa-Cm)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 97:102 (Bk-No)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 103:108 (Lr-Hs)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0,      0.d0,& ! 109:114 (Mt-Uuq)
          0.d0,      0.d0,      0.d0,    0.d0,&                       ! 115:118 (Uup-Uuo)
          0.d0,      0.d0,      0.d0,    0.d0,      0.d0/)            ! 119:123

    character*1, parameter :: maxang(maxzat0) = (/&
       "s", "x", "x", "x", "x", "p",& ! 1:6   (H-C)
       "p", "p", "p", "x", "p", "p",& ! 7:12  (N-Mg)
       "x", "x", "d", "d", "d", "x",& ! 13:18 (Al-Ar)
       "p", "p", "x", "x", "x", "x",& ! 19:24 (K-Cr)
       "x", "x", "x", "x", "x", "d",& ! 25:30 (Mn-Zn)
       "x", "x", "x", "x", "d", "x",& ! 31:36 (Ga-Kr)
       "x", "x", "x", "x", "x", "x",& ! 37:42 (Rb-Mo)
       "x", "x", "x", "x", "x", "x",& ! 43:48 (Tc-Cd)
       "x", "x", "x", "x", "d", "x",& ! 49:54 (In-Xe)
       "x", "x", "x", "x", "x", "x",& ! 55:60 (Cs-Nd)
       "x", "x", "x", "x", "x", "x",& ! 61:66 (Pm-Dy)
       "x", "x", "x", "x", "x", "x",& ! 67:72 (Ho-Hf)
       "x", "x", "x", "x", "x", "x",& ! 73:78 (Ta-Pt)
       "x", "x", "x", "x", "x", "x",& ! 79:84 (Au-Po)
       "x", "x", "x", "x", "x", "x",& ! 85:90 (At-Th)
       "x", "x", "x", "x", "x", "x",& ! 91:96 (Pa-Cm)
       "x", "x", "x", "x", "x", "x",& ! 97:102 (Bk-No)
       "x", "x", "x", "x", "x", "x",& ! 103:108 (Lr-Hs)
       "x", "x", "x", "x", "x", "x",& ! 109:114 (Mt-Uuq)
       "x", "x", "x", "x",&           ! 115:118 (Uup-Uuo)
       "x", "x", "x", "x", "x"/)      ! 119:123

    integer :: lu, i

    lu = fopen_write(file,ti=ti)
    write (lu,'("Geometry = GenFormat {")')
    call c%write_dftbp_gen(file,lu,ti=ti)
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("Driver = ConjugateGradient {")')
    write(lu,'("       MovedAtoms = 1:-1")')
    write(lu,'("       MaxForceComponent = 1e-5")')
    write(lu,'("       MaxSteps = 3000")')
    write(lu,'("       LatticeOpt = Yes")')
    write(lu,'("       OutputPrefix = ""geo_end""")')
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("Hamiltonian = DFTB{")')
    write(lu,'("  ThirdOrderFull = Yes")')
    write(lu,'("  SCC = Yes")')
    write(lu,'("  SCCTolerance = 1e-7")')
    write(lu,'("  MaxSCCIterations = 125")')
    write(lu,'("  MaxAngularMomentum = {")')
    do i = 1, c%nspc
       write (lu,'("    ",A," = ",A)') string(nameguess(c%spc(i)%z,.true.)), &
          string(maxang(c%spc(i)%z))
    end do
    write(lu,'("  }")')
    write(lu,'("  SlaterKosterFiles = Type2FileNames {")')
    write(lu,'("    Prefix = ""xxx""")')
    write(lu,'("    Separator = ""-""")')
    write(lu,'("    Suffix = "".skf""")')
    write(lu,'("    LowerCaseTypeName = No")')
    write(lu,'("  }")')
    if (.not.c%ismolecule) then
       write(lu,'("  KPointsAndWeights = SupercellFolding {")')
       write(lu,'("    4 0 0 ")')
       write(lu,'("    0 4 0")')
       write(lu,'("    0 0 4")')
       write(lu,'("    0.5 0.5 0.5")')
       write(lu,'("  }")')
    end if
    write(lu,'("  DampXH = Yes")')
    write(lu,'("  DampXHExponent = 4.2")')
    write(lu,'("  HubbardDerivs {")')
    do i = 1, c%nspc
       write (lu,'("    ",A," = ",A)') string(nameguess(c%spc(i)%z,.true.)), &
          string(hderiv(c%spc(i)%z),'f',decimal=4)
    end do
    write(lu,'("  }")')
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("Options {")')
    write(lu,'("  WriteDetailedXML = Yes")')
    write(lu,'("}")')
    write(lu,'("")')
    write(lu,'("ParserOptions {")')
    write(lu,'("  ParserVersion = 4")')
    write(lu,'("}")')
    write(lu,'("")')
    call fclose(lu)

  end subroutine write_dftbp_hsd

  !> Write a DFTB+ human-friendly gen structure file
  module subroutine write_dftbp_gen(c,file,lu0,ti)
    use tools_io, only: fopen_write, nameguess, string, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in), optional :: lu0
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, k
    real*8 :: r(3,3)
    character(len=:), allocatable :: strtyp, aux

    ! open file
    if (present(lu0)) then
       lu = lu0
    else
       lu = fopen_write(file,ti=ti)
    end if

    ! atom types
    strtyp = ""
    do i = 1, c%nspc
       aux = strtyp // " " // string(nameguess(c%spc(i)%z,.true.))
       strtyp = aux
    end do

    if (c%ismolecule) then
       ! molecule
       write (lu,'(A," C")') string(c%ncel)
       write (lu,'(A)') strtyp

       ! Cartesian coordinates
       do i = 1, c%nspc
          do k = 1, c%ncel
             if (c%atcel(k)%is == i) &
                write (lu,'(99(A," "))') string(k), string(i), &
                (string(c%atcel(k)%r(j)*bohrtoa,'f',20,12),j=1,3)
          end do
       end do
    else
       ! crystal
       write (lu,'(A," F")') string(c%ncel)
       write (lu,'(A)') strtyp

       ! fractional coordinates
       do i = 1, c%nspc
          do k = 1, c%ncel
             if (c%atcel(k)%is == i) &
                write (lu,'(99(A," "))') string(k), string(i), &
                (string(c%atcel(k)%x(j),'f',20,12),j=1,3)
          end do
       end do

       ! lattice vectors
       r = c%m_x2c * bohrtoa
       write (lu,'(3(A," "))') (string(0d0,'f',20,12),j=1,3)
       do i = 1, 3
          write (lu,'(3(A," "))') (string(r(j,i),'f',20,12),j=1,3)
       end do
    endif

    ! close file
    if (.not.present(lu0)) call fclose(lu)

  end subroutine write_dftbp_gen

  !> Write the crystal structure in pyscf format (python script)
  module subroutine write_pyscf(c,file,ti)
    use tools_io, only: fopen_write, fclose, string, nameguess
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j

    lu = fopen_write(file,ti=ti)
    if (c%ismolecule) then
       write (lu,'("from pyscf import gto")')
       write (lu,'("")')
       write (lu,'("mol = gto.Mole()")')
       write (lu,'("mol.atom = ''''''")')
       do i = 1, c%ncel
          write (lu,'(A,3(" ",A))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.)),&
             (string(c%atcel(i)%r(j),'f',18,10),j=1,3)
       end do
       write (lu,'("''''''")')
       write (lu,'("mol.basis = ''aug-cc-pvtz''")')
       write (lu,'("mol.verbose = 4")')
       write (lu,'("mol.spin = 0")')
       write (lu,'("mol.charge = 0")')
       write (lu,'("mol.build()")')
    else
       write (lu,'("from pyscf.pbc import gto")')
       write (lu,'("")')
       write (lu,'("cell = gto.Cell()")')
       write (lu,'("cell.atom = ''''''")')
       do i = 1, c%ncel
          write (lu,'(A,3(" ",A))') string(nameguess(c%spc(c%atcel(i)%is)%z,.true.)),&
             (string(c%atcel(i)%r(j),'f',18,10),j=1,3)
       end do
       write (lu,'("''''''")')
       write (lu,'("cell.a = ''''''")')
       do i = 1, 3
          write (lu,'("  ",1p,3(E22.14," "))') c%m_x2c(:,i)
       end do
       write (lu,'("''''''")')

       write (lu,'("cell.unit = ''Bohr''")')
       write (lu,'("cell.basis = ''gth-szv''")')
       write (lu,'("cell.pseudo = ''gth-pade''")')
       write (lu,'("cell.verbose = 4")')
       write (lu,'("cell.build()")')
    end if
    call fclose(lu)

  end subroutine write_pyscf

  !> Write the crystal or molecualr structure in FHIaims geometry.in
  !> format.  If frac = .true., use atom_frac instead of atom if c is
  !> a crystal.
  module subroutine write_fhi(c,file,frac,rklength,ti)
    use tools_io, only: fopen_write, fclose, string, nameguess
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: frac
    real*8, intent(in), optional :: rklength
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, nk(3)
    character*2 :: name

    nk = 0
    if (present(rklength)) &
       call c%get_kpoints(rklength,nk)

    lu = fopen_write(file,ti=ti)
    write (lu,'("## FHIaims input file generated by critic2.")')
    if (.not.c%ismolecule) then
       do i = 1, 3
          write (lu,'("lattice_vector ",3(A," "))') (string(c%m_x2c(j,i) * bohrtoa,'f',18,10),j=1,3)
       end do
    end if
    do i = 1, c%ncel
       name = nameguess(c%spc(c%atcel(i)%is)%z,.true.)
       if (frac .and. .not.c%ismolecule) then
          write (lu,'("atom_frac ",4(" ",A))') (string(c%atcel(i)%x(j),'f',18,10),j=1,3),&
             name
       else
          write (lu,'("atom ",4(" ",A))') (string(c%atcel(i)%r(j)*bohrtoa,'f',18,10),j=1,3),&
             name
       end if
    end do
    call fclose(lu)

    if (present(rklength) .and. all(nk > 0)) then
       lu = fopen_write(trim(file) // "_control",ti=ti)
       write (lu,'("xc b86bpbe")')
       write (lu,'("spin none")')
       write (lu,'("charge 0")')
       write (lu,'("output_level MD_light")')
       write (lu,'("xdm lightdenser")')
       write (lu,'("relativistic atomic_zora scalar")')
       write (lu,'("")')
       write (lu,'("sc_accuracy_rho 1e-7")')
       write (lu,'("")')
       write (lu,'("relax_geometry bfgs 0.005")')
       write (lu,'("relax_unit_cell full")')
       write (lu,'("max_relaxation_steps 200")')
       write (lu,'("")')
       write (lu,'("k_grid",3(" ",A))') (string(nk(i)),i=1,3)
       call fclose(lu)
    end if

  end subroutine write_fhi

  !> Write a TINKER frac file with tiny FF parametrization.
  module subroutine write_tinkerfrac(c,file,ti)
    use tools_io, only: fopen_write, string, fclose, nameguess, ioj_right, ioj_left, ferror, faterr,&
       warning
    use tools, only: tiny_atom_type
    use types, only: realloc
    use param, only: bohrtoa
    use fragmentmod, only: fragment
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, nn, lu
    integer :: iz, nneig, ityp
    type(fragment) :: fr
    integer, allocatable :: imap(:), ineig(:)
    logical :: dowarn

    ! open file
    lu = fopen_write(file,ti=ti)

    ! line 1: number of atoms and title (the spaces at the end are there for the USPEX reader)
    write (lu,'(A,"          ")') string(c%ncel)

    ! line 2: cell parameters (the spaces at the end are there for the USPEX reader)
    write (lu,'(6(A," "),"          ")') (string(c%aa(i)*bohrtoa,'f',decimal=8),i=1,3), &
       (string(c%bb(i),'f',decimal=8),i=1,3)

    ! merge all molecules into the onemotif
    call fr%merge_array(c%mol(1:c%nmol),.false.)

    ! build the atom map
    allocate(imap(c%ncel))
    imap = 0
    nn = 0
    do i = 1, fr%nat
       imap(fr%at(i)%cidx) = i
    end do
    if (any(imap == 0)) &
       call ferror ('write_tinkerfrac',&
       'some atoms are not mapped into molecular fragments; perhaps not a molecular crystal?',faterr)

    ! do the atom typing and write the atom list
    dowarn = .false.
    allocate(ineig(10))
    do i = 1, fr%nat
       ! type of atom
       iz = fr%spc(fr%at(i)%is)%z
       nneig = c%nstar(fr%at(i)%cidx)%ncon
       ityp = tiny_atom_type(iz,nneig)
       dowarn = dowarn .or. (ityp == 0)

       ! neighbor list
       if (nneig > size(ineig,1)) call realloc(ineig,nneig)
       do j = 1, nneig
          ineig(j) = imap(c%nstar(fr%at(i)%cidx)%idcon(j))
       end do

       write (lu,'(999(A," "))') string(i,5,ioj_right), &
          string(nameguess(iz,.true.),5,ioj_left), &
          (string(fr%at(i)%x(j),'f',12,8,ioj_right),j=1,3), string(ityp),&
          (string(ineig(j)),j=1,nneig)
    end do
    deallocate(imap,ineig)

    ! close and clean up
    if (dowarn) &
       call ferror('write_tinkerfrac','some atom types could not be assigned',warning)
    call fclose(lu)

  end subroutine write_tinkerfrac

  !> Write a PDB file.
  module subroutine write_pdb(c,file,ti)
    use tools_io, only: fopen_write, string, fclose, upper, ioj_right, nameguess
    use tools_math, only: gcd
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti

    integer :: lu
    integer, allocatable :: nis(:)
    integer :: i
    real*8 :: maxdv, x(3)
    character(len=:), allocatable :: str
    character*2 :: atsym
    character*11 :: spg

    ! open file
    lu = fopen_write(file,ti=ti)
    write (lu,'("HEADER    UNKNOWN FUNCTION                        01-JAN-01   9Z9Z              ")')
    write (lu,'("TITLE     PDB FILE GENERATED BY CRITIC2, IGNORE ALL ENTRIES EXCEPT              ")')
    write (lu,'("TITLE    2 HETATM, CRYST1, SCALEn                                               ")')
    write (lu,'("COMPND    MOL_ID: 1                                                             ")')
    write (lu,'("SOURCE    MOL_ID: 1                                                             ")')
    write (lu,'("KEYWDS    UNKNOWN                                                               ")')
    write (lu,'("EXPDTA    X-RAY DIFFRACTION                                                     ")')
    write (lu,'("AUTHOR    CRITIC2                                                               ")')
    write (lu,'("REVDAT   1   01-JAN-01 9Z9Z    0                                                ")')
    write (lu,'("REMARK   2                                                                      ")')
    write (lu,'("REMARK   2 RESOLUTION.    2.35 ANGSTROMS.                                       ")')
    write (lu,'("REMARK   3 REFINEMENT.                                                          ")')
    write (lu,'("HET    UNL  A   11  ",A,"     UNKNOWN                                           ")') &
       string(c%ncel,5)
    write (lu,'("HETNAM     UNL ",A,"          ")') string(upper(c%file),55)

    ! empirical formula
    allocate(nis(c%nspc))
    nis = 0
    do i = 1, c%nneq
       nis(c%at(i)%is) = nis(c%at(i)%is) + c%at(i)%mult
    end do
    maxdv = gcd(nis,c%nspc)
    str = ""
    do i = 1, c%nspc
       str = str // string(c%spc(i)%name) // "(" // string(nint(nis(i)/maxdv)) // ") "
    end do
    str = trim(adjustl(str))
    write (lu,'("FORMUL   1  UNL    ",A51,"          ")') str

    ! lattice parameters
    if (.not.c%ismolecule) then
       if (c%havesym > 0) then
          spg = c%spg%international_symbol
       else
          spg = "P1"
       end if
       write (lu,'("CRYST1",3(F9.3),3(F7.2)," ",2(A),"          ")') c%aa * bohrtoa, c%bb,&
          string(spg,11), string(c%nmol,4)
       write (lu,'("SCALE1    ",3(F10.6),"     ",F10.5,"                         ")') c%m_c2x(1,:) / bohrtoa, 0d0
       write (lu,'("SCALE2    ",3(F10.6),"     ",F10.5,"                         ")') c%m_c2x(2,:) / bohrtoa, 0d0
       write (lu,'("SCALE3    ",3(F10.6),"     ",F10.5,"                         ")') c%m_c2x(3,:) / bohrtoa, 0d0
    end if

    ! coordinates
    do i = 1, c%ncel
       if (c%ismolecule) then
          x = (c%atcel(i)%r + c%molx0) * bohrtoa
       else
          x = c%atcel(i)%r * bohrtoa
       end if
       atsym = adjustr(nameguess(c%spc(c%atcel(i)%is)%z,.true.))
       write (lu,'("HETATM",A," ",A," UNL A    1   ",3(F8.3),2(F6.2),"          ",A2,"  ")') &
          string(i,5,ioj_right), string(c%spc(c%atcel(i)%is)%name,4,ioj_right),&
          x, 1d0, 1d0, atsym
    end do
    call fclose(lu)

  end subroutine write_pdb

  !> Write a CASTEP cell input template
  module subroutine write_castep_cell(c,file,rklength,ti)
    use tools_io, only: fopen_write, fclose, string, nameguess, upper
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    real*8, intent(in), optional :: rklength
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, lu, nk(3), iz
    real*8 :: rk
    character*2 :: sym

    rk = 40d0
    if (present(rklength)) rk = rklength

    lu = fopen_write(file,ti=ti)
    write (lu,'("## CASTEP input file generated by critic2.")')
    write (lu,*)
    write (lu,'("%block lattice_abc")')
    write (lu,'("ang")')
    write (lu,'("  ",3(" ",A))') (string(c%aa(i) * bohrtoa,'f',decimal=10),i=1,3)
    write (lu,'("  ",3(" ",A))') (string(c%bb(i),'f',decimal=5),i=1,3)
    write (lu,'("%endblock lattice_abc")')
    write (lu,*)

    write (lu,'("%block positions_frac")')
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       sym = nameguess(iz,.true.)
       sym(1:1) = upper(sym(1:1))
       write (lu,'("  ",4(" ",A))') sym, (string(c%atcel(i)%x(j),'f',decimal=10),j=1,3)
    end do
    write (lu,'("%endblock positions_frac")')
    write (lu,*)

    call c%get_kpoints(rk,nk)
    write (lu,'("kpoint_mp_grid ",3(A," ")/)') (string(nk(i)),i=1,3)

    call fclose(lu)

  end subroutine write_castep_cell

  !> Write a grid to a cube file. The input is the crystal (c),
  !> the grid in 3D array form (g), the filename (file), and whether
  !> to write the whole cube or only the header (onlyheader). If xd0
  !> is given, use it as the metric of the cube; otherwise, use the
  !> unit cell. If x00 is given, use it as the origin of the cube
  !> (in bohr). Otherwise, use the crystal's molx0.
  module subroutine writegrid_cube(c,g,file,onlyheader,binary,xd0,x00,ishift0,ti)
    use global, only: precisecube
    use tools_io, only: fopen_write, fclose
    use param, only: eye
    class(crystal), intent(in) :: c
    real*8, intent(in), allocatable :: g(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: onlyheader
    logical, intent(in) :: binary
    real*8, intent(in), optional :: xd0(3,3)
    real*8, intent(in), optional :: x00(3)
    integer, intent(in), optional :: ishift0(3)
    type(thread_info), intent(in), optional :: ti

    integer :: n(3), i, ix, iy, lu, ishift(3)
    integer :: iix, iiy, iiz
    real*8 :: xd(3,3), x0(3), rshift(3)

    ! grid size
    if (onlyheader.and..not.allocated(g)) then
       n = 1
    else
       do i = 1, 3
          n(i) = size(g,i)
       end do
    end if

    ! process optional arguments
    ishift = 0
    if (present(ishift0)) then
       ishift = ishift0
    end if
    rshift = real(ishift,8) / n
    rshift = c%x2c(rshift)
    if (present(xd0)) then
       xd = xd0
    else
       xd = eye
       do i = 1, 3
          xd(:,i) = c%x2c(xd(:,i))
          xd(:,i) = xd(:,i) / real(n(i),8)
       end do
    endif
    if (present(x00)) then
       x0 = x00
    else
       x0 = c%molx0
    endif

    ! write
    if (binary) then
       lu = fopen_write(file,form="unformatted",ti=ti)
       write(lu) c%ncel, x0
       write(lu) n, xd
       do i = 1, c%ncel
          write(lu) c%spc(c%atcel(i)%is)%z, 0d0, c%atcel(i)%r(:) + c%molx0 - rshift
       end do
       write (lu) g
    else
       lu = fopen_write(file,ti=ti)
       write(lu,'("critic2-cube")')
       write(lu,'("critic2-cube")')
       if (precisecube) then
          write(lu,'(I5," ",3(E22.14E3," "))') c%ncel, x0
          do i = 1, 3
             write(lu,'(I5," ",3(E22.14E3," "))') n(i), xd(:,i)
          end do
          do i = 1, c%ncel
             write(lu,'(I4,F5.1,3(E22.14E3," "))') c%spc(c%atcel(i)%is)%z, 0d0, &
                c%atcel(i)%r(:) + c%molx0 - rshift
          end do
       else
          write(lu,'(I5,3(F12.6))') c%ncel, x0
          do i = 1, 3
             write(lu,'(I5,3(F12.6))') n(i), xd(:,i)
          end do
          do i = 1, c%ncel
             write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') c%spc(c%atcel(i)%is)%z, 0d0, &
                c%atcel(i)%r(:) + c%molx0 - rshift
          end do
       end if
       if (.not.onlyheader) then
          do iix = 0, n(1)-1
             ix = modulo(iix + ishift(1),n(1)) + 1
             do iiy = 0, n(2)-1
                iy = modulo(iiy + ishift(2),n(2)) + 1
                if (precisecube) then
                   write (lu,'(6(" ",E22.14E3))') (g(ix,iy,modulo(iiz+ishift(3),n(3))+1),iiz=0,n(3)-1)
                else
                   write (lu,'(1p,6(" ",E12.5E3))') (g(ix,iy,modulo(iiz+ishift(3),n(3))+1),iiz=0,n(3)-1)
                end if
             enddo
          enddo
       end if
    end if
    call fclose(lu)

  end subroutine writegrid_cube

  !> Write a grid to a VASP CHGCAR file. The input is the crystal (c),
  !> the grid in 3D array form (g), the filename (file), and whether
  !> to write the whole cube or only the header (onlyheader).
  module subroutine writegrid_vasp(c,g,file,onlyheader,ishift0,ti)
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    real*8, intent(in), allocatable :: g(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: onlyheader
    integer, intent(in), optional :: ishift0(3)
    type(thread_info), intent(in), optional :: ti

    integer :: n(3), i, j, ix, iy, iz, lu, ishift(3), nat
    character(len=:), allocatable :: line0, aux
    real*8 :: xshift(3)

    ! grid size
    if (onlyheader.and..not.allocated(g)) then
       n = 1
    else
       do i = 1, 3
          n(i) = size(g,i)
       end do
    end if

    ! process optional arguments
    ishift = 0
    if (present(ishift0)) then
       ishift = ishift0
    end if
    xshift = real(ishift,8) / n

    lu = fopen_write(file,ti=ti)
    write (lu,'("CHGCAR generated by critic2")')
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(1p,3(E22.14," "))') c%m_x2c(:,i) * bohrtoa
    end do

    ! species line
    line0 = ""
    do i = 1, c%nspc
       aux = line0 // " " // string(nameguess(c%spc(i)%z,.true.))
       line0 = aux
    end do
    write (lu,'(A)') line0

    ! number of atoms line
    line0 = ""
    do i = 1, c%nspc
       nat = 0
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) nat = nat + 1
       end do
       aux = line0 // " " // string(nat)
       line0 = aux
    end do
    write (lu,'(A)') line0

    ! positions
    write (lu,'("Direct")')
    do i = 1, c%nspc
       do j = 1, c%ncel
          if (c%atcel(j)%is == i) &
             write (lu,'(1p,3(E22.14," "))') c%atcel(j)%x - xshift
       end do
    end do
    write (lu,*)
    write (lu,'(3(I5," "))') n
    if (.not.onlyheader) then
       write (lu,'(5(" ",e22.14))') &
          (((g(modulo(ix+ishift(1),n(1))+1,modulo(iy+ishift(2),n(2))+1,modulo(iz+ishift(3),n(3))+1)*c%omega,&
          ix=0,n(1)-1),iy=0,n(2)-1),iz=0,n(3)-1)
    end if
    call fclose(lu)

  end subroutine writegrid_vasp

  !> Write a grid to a xsf file. The input is the crystal (c), the
  !> grid in 3D array form (g), the filename (file), and whether to
  !> write the whole xsf or only the structure (onlyheader).
  module subroutine writegrid_xsf(c,g,file,onlyheader,ishift0,ti)
    use tools_io, only: fopen_write, fclose, string, ferror, faterr
    use param, only: bohrtoa
    class(crystal), intent(in) :: c
    real*8, intent(in), allocatable :: g(:,:,:)
    character*(*), intent(in) :: file
    logical, intent(in) :: onlyheader
    integer, intent(in), optional :: ishift0(3)
    type(thread_info), intent(in), optional :: ti

    real*8 :: x(3), rshift(3)
    integer :: n(3), i, j, ix, iy, iz, lu, ishift(3)

    if (c%ismolecule) then
       call ferror('writegrid_xsf','molecular xsf files not supported yet',faterr)
    end if

    ! grid size
    if (onlyheader .and..not.allocated(g)) then
       n = 1
    else
       do i = 1, 3
          n(i) = size(g,i)
       end do
    end if

    ! process optional arguments
    ishift = 0
    if (present(ishift0)) then
       ishift = ishift0
    end if
    rshift = real(ishift,8) / n
    rshift = c%x2c(rshift)

    lu = fopen_write(file,ti=ti)
    write (lu,'("## xsf file generated by critic2")')
    write (lu,'("CRYSTAL")')
    write (lu,'("PRIMVEC")')
    do i = 1, 3
       write (lu,'(1p,3(E22.14," "))') c%m_x2c(:,i) * bohrtoa
    end do
    write (lu,'("CONVVEC")')
    do i = 1, 3
       write (lu,'(1p,3(E22.14," "))') c%m_x2c(:,i) * bohrtoa
    end do

    write (lu,'("PRIMCOORD")')
    write (lu,'(A," 1")') string(c%ncel)
    do j = 1, c%ncel
       write (lu,'(A," ",1p,3(E22.14," "))') string(c%spc(c%atcel(j)%is)%z), &
          (c%atcel(j)%r - rshift) * bohrtoa
    end do
    if (.not.onlyheader) then
       write (lu,'("BEGIN_BLOCK_DATAGRID3D")')
       write (lu,'("critic2:grid")')
       write (lu,'("DATAGRID_3D_GRID")')
       write (lu,'(3(A," "))') (string(n(j)+1),j=1,3)
       x = 0d0
       write (lu,'(1p,3(" ",e22.14))') x
       do i = 1, 3
          x = 0d0
          x(i) = 1d0
          x = c%x2c(x) * bohrtoa
          write (lu,'(1p,3(" ",e22.14))') x
       end do
       write (lu,'(1p,5(" ",e22.14))') (((g(modulo(ix+ishift(1),n(1))+1,modulo(iy+ishift(2),n(2))+1,&
          modulo(iz+ishift(3),n(3))+1),ix=0,n(1)),iy=0,n(2)),iz=0,n(3))
       write (lu,'("END_DATAGRID_3D")')
       write (lu,'("END_BLOCK_DATAGRID3D")')
    end if

    call fclose(lu)

  end subroutine writegrid_xsf

end submodule write
