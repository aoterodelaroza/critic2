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

! The struct_read_qein and qe_latgen routines in this module were
! adapted from
! Quantum ESPRESSO, version 4.3.2.
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! Routines that write the crystal structure to the output in several formats
module struct_writers
  implicit none

  private
  public :: struct_write_mol
  public :: struct_write_3dmodel
  public :: struct_write_espresso
  public :: struct_write_vasp
  public :: struct_write_abinit
  public :: struct_write_elk
  public :: struct_write_gaussian
  public :: struct_write_tessel
  public :: struct_write_critic
  public :: struct_write_cif
  public :: struct_write_d12
  public :: struct_write_escher
  public :: struct_write_gulp
  public :: struct_write_lammps
  public :: struct_write_siesta_fdf
  public :: struct_write_siesta_in
  public :: struct_write_dftbp_hsd
  public :: struct_write_dftbp_gen

contains

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
  subroutine struct_write_mol(c,file,fmt,ix,doborder,onemotif,molmotif,&
     environ,renv,lnmer,nmer,rsph,xsph,rcub,xcub,luout)
    use fragmentmod, only: fragment_merge_array, fragment_cmass, fragment_init
    use struct_basic, only: crystal
    use global, only: dunit0, iunit
    use graphics, only: writecml, writexyz, writegjf
    use tools_math, only: norm, nchoosek, comb
    use tools_io, only: ferror, faterr, uout, string, ioj_left, string, ioj_right,&
       equal
    use types, only: fragment, realloc
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in) :: ix(3)
    logical, intent(in) :: doborder, onemotif, molmotif, environ
    real*8, intent(in) :: renv
    logical, intent(in) :: lnmer
    integer, intent(in) :: nmer
    real*8, intent(in) :: rsph, xsph(3)
    real*8, intent(in) :: rcub, xcub(3)
    integer, intent(out), optional :: luout

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

    ! calculate the motifs if necessary
    if (onemotif .or. environ) &
       call c%checkflags(.false.,ast0=.true.)

    ! determine the fragments
    if (onemotif) then
       fr = fragment_merge_array(c%mol(1:c%nmol))
       allocate(fr0(c%nmol))
       fr0 = c%mol
       nmol = c%nmol
    elseif (environ) then
       ! calculate the centers of mass for all fragments in the molecular motif
       allocate(cmlist(3,c%nmol),fr0(c%nmol),origmol(c%nmol))
       nmol = c%nmol
       fr0 = c%mol
       do i = 1, c%nmol
          cmlist(:,i) = fragment_cmass(c%mol(i))
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
                      dist = norm(xcm)
                      if (dist <= renv) then
                         ncm = ncm + 1
                         if (ncm > size(cmlist,2)) then
                            call realloc(cmlist,3,2*ncm)
                            call realloc(origmol,2*ncm)
                            call realloc(fr0,2*ncm)
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
       fr = fragment_merge_array(fr0)
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
          fr = fragment_merge_array(fr0)
       end if
    end if

    if (lnmer .and..not.allocated(fr0)) &
       call ferror('struct_write_mol','ONEMOTIF, MOLMOTIF, or ENVIRON are necessary with NMER',faterr)

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
          write (uout,'(99(2X,A))') string(i,3,ioj_left), string(fr0(i)%nat,4,ioj_left),&
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
       call dowrite(file,fr)
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
                call dowrite(file0,fr0(j))
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
                      fr = fragment_merge_array((/fr,fr0(icomb(k))/))
                   end do
                   aux = trim(file0) // "." // fmt
                   file0 = aux
                   call dowrite(file0,fr)
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
                call fragment_init(fr)
                do k = 1, i
                   aux = trim(file0) // "_" // string(icomb(k))
                   file0 = aux
                   fr = fragment_merge_array((/fr,fr0(icomb(k))/))
                end do
                aux = trim(file0) // "." // fmt
                file0 = aux
                call dowrite(file0,fr)
                icount = icount + 1
             end do
             deallocate(icomb)
             write (uout,'("+ Written ",A," ",A,"-mers")') string(icount), string(i)
          end if
       end do
    end if

  contains
    subroutine dowrite(fileo,fro)
      character*(*) :: fileo
      type(fragment) :: fro
      
      if (equal(fmt,"xyz")) then
         call writexyz(fileo,fro)
      elseif (equal(fmt,"gjf")) then
         call writegjf(fileo,fro)
      elseif (equal(fmt,"cml")) then
         if (c%ismolecule) then
            call writecml(fileo,fro,luout=luout)
         else
            call writecml(fileo,fro,c%crys2car,luout=luout)
         end if
      else
         call ferror("struct_write_mol","Unknown format",faterr)
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
  end subroutine struct_write_mol

  !> Write an obj/ply/off file containing the crystal structure. fmt
  !> can be one of obj, ply, or off. ix is the number of unit cells to
  !> plot. If doborder is .true., add all atoms at the border. If
  !> onemotif is .true., write all molecules in the unit cell. If
  !> molmotif is .true., complete molecules with atoms in adjacent
  !> cells. If docell, add sticks for the unit cell
  !> limits. If domolcell, add sticks tfor the molecular cell. If rsph
  !> (bohr) is positive, then use all atoms in a sphere around xsph
  !> (cryst.). If rcub (bohr) is positive, use all atoms in a cube
  !> around xcub (cryst.). If lu0 and lumtl0 are present, then return
  !> the logical units for the obj and the mtl files and do not close
  !> the files.
  subroutine struct_write_3dmodel(c,file,fmt,ix,doborder,onemotif,molmotif,&
     docell,domolcell,rsph,xsph,rcub,xcub,lu0,lumtl0)
    use graphics, only: graphics_open, graphics_ball, graphics_stick, graphics_close
    use fragmentmod, only: fragment_merge_array
    use struct_basic, only: crystal
    use tools_math, only: norm
    use types, only: fragment
    use tools_io, only: equal
    use param, only: maxzat, atmcov, jmlcol
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in) :: ix(3)
    logical, intent(in) :: doborder, onemotif, molmotif
    logical, intent(in) :: docell, domolcell
    real*8, intent(in) :: rsph, xsph(3)
    real*8, intent(in) :: rcub, xcub(3)
    integer, intent(out), optional :: lu0, lumtl0

    integer :: lu, lumtl
    integer :: i, j
    real*8 :: d, xd(3), x0(3), x1(3), rr
    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: nmol

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

    ! open and get the atom list
    if (onemotif) then
       fr = fragment_merge_array(c%mol(1:c%nmol))
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
          fr = fragment_merge_array(fr0)
       end if
    end if

    ! if this is a molecule, translate to the proper origin
    if (c%ismolecule) then
       do i = 1, fr%nat
          fr%at(i)%r = fr%at(i)%r + c%molx0
       end do
    end if

    lumtl = 0
    call graphics_open(fmt,file,lu,lumtl)

    ! add the balls
    do i = 1, fr%nat
       if (fr%at(i)%z > maxzat) then
          rr = 0.21d0
       else
          rr = 0.6d0*atmcov(fr%at(i)%z)
       endif
       call graphics_ball(fmt,lu,fr%at(i)%r,JMLcol(:,fr%at(i)%z),rr)
    end do

    ! add the sticks
    do i = 1, fr%nat
       do j = i+1, fr%nat
          if (fr%at(i)%z > maxzat .or. fr%at(j)%z > maxzat) cycle
          xd = fr%at(i)%r - fr%at(j)%r
          d = norm(xd)
          if (d < (atmcov(fr%at(i)%z) + atmcov(fr%at(j)%z)) * rfac) then
             xd = fr%at(i)%r + 0.5d0 * (fr%at(j)%r - fr%at(i)%r)
             call graphics_stick(fmt,lu,fr%at(i)%r,xd,JMLcol(:,fr%at(i)%z),0.05d0)
             call graphics_stick(fmt,lu,fr%at(j)%r,xd,JMLcol(:,fr%at(j)%z),0.05d0)
          end if
       end do
    end do

    ! add the cell
    if (docell) then
       do i = 1, 12
          x0 = c%x2c(x0cell(:,1,i)) + c%molx0
          x1 = c%x2c(x0cell(:,2,i)) + c%molx0
          call graphics_stick(fmt,lu,x0,x1,(/255,0,0/),0.03d0)
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
          call graphics_stick(fmt,lu,x0,x1,(/0,0,255/),0.03d0)
       end do
    end if

    ! close or give the handles to the calling routine, cleanup
    if (present(lumtl0) .and. present(lu0)) then
       lu0 = lu
       lumtl0 = lumtl
    else
       call graphics_close(fmt,lu,lumtl)
    end if

    if (allocated(fr0)) deallocate(fr0)

  end subroutine struct_write_3dmodel

  !> Write a quantum espresso input template
  subroutine struct_write_espresso(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, nameguess, lower, fclose
    use param, only: maxzat0, atmass

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    character(len=:), allocatable :: lbl1
    integer :: i, lu
    logical :: ztyp(maxzat0)

    ztyp = .false.
    do i = 1, c%nneq
       if (c%at(i)%z > 0 .and. c%at(i)%z <= maxzat0) then
          ztyp(c%at(i)%z) = .true.
       end if
    end do

    lu = fopen_write(file)
    write (lu,'("&control")')
    write (lu,'(" title=''crystal'',")')
    write (lu,'(" prefix=''crystal'',")')
    write (lu,'(" pseudo_dir=''.'',")')
    write (lu,'(" calculation=''vc-relax'',")')
    write (lu,'("/")')
    write (lu,'("&system")')
    write (lu,'(" ibrav=0,")')
    write (lu,'(" celldm(1)=1.0,")')
    write (lu,'(" nat=",I6,",")') c%ncel
    write (lu,'(" ntyp=",I3,",")') count(ztyp)
    write (lu,'(" ecutwfc=60.0,")')
    write (lu,'(" ecutrho=600.0,")')
    write (lu,'(" xdm=.true.,")')
    write (lu,'(" xdm_a1=0.4073,")')
    write (lu,'(" xdm_a2=2.4150,")')
    write (lu,'("/")')
    write (lu,'("&electrons"/" conv_thr = 1d-8,"/"/")')
    write (lu,'("&ions"/"/")')
    write (lu,'("&cell"/"/")')
    write (lu,'("ATOMIC_SPECIES")')
    do i = 1, size(ztyp)
       if (ztyp(i)) then
          lbl1 = lower(nameguess(i,.true.))
          write (lu,'(A2,X,F12.6,X,A,".UPF")') lbl1(1:2), atmass(i), trim(lbl1)
       end if
    end do
    write (lu,'(/"ATOMIC_POSITIONS crystal")')
    do i = 1, c%ncel
       write (lu,'(A2,3(X,F13.8,X))') lower(nameguess(c%at(c%atcel(i)%idx)%z,.true.)), c%atcel(i)%x
    end do
    write (lu,'(/"K_POINTS automatic"/"2 2 2 1 1 1"/)')
    write (lu,'("CELL_PARAMETERS cubic")')
    do i = 1, 3
       write (lu,'(3(F18.12,X))') c%crys2car(:,i)
    end do
    call fclose(lu)

  end subroutine struct_write_espresso

  !> Write a VASP POSCAR template
  subroutine struct_write_vasp(file,c,verbose)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, nameguess, string, uout, fclose
    use param, only: bohrtoa, maxzat0

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    logical, intent(in) :: verbose

    character(len=:), allocatable :: lbl1
    integer :: i, j, lu
    integer :: ntyp(maxzat0)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, c%ncel
       ntyp(c%at(c%atcel(i)%idx)%z) = ntyp(c%at(c%atcel(i)%idx)%z) + 1
    end do

    ! Cell
    lu = fopen_write(file)
    write (lu,'("crystal")')
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(3(F15.10,X))') c%crys2car(:,i) * bohrtoa
    end do

    ! Number of atoms per type and Direct
    lbl1 = ""
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          lbl1 = lbl1 // " " // string(ntyp(i))
       end if
    end do
    write (lu,'(A)') lbl1
    write (lu,'("Direct")')

    ! Atomic positions
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          do j = 1, c%ncel
             if (c%at(c%atcel(j)%idx)%z == i) then
                write (lu,'(3(F13.8,X))') c%atcel(j)%x
             end if
          end do
       end if
    end do
    call fclose(lu)

    lbl1 = ""
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          lbl1 = lbl1 // " " // string(nameguess(i,.true.))
       end if
    end do

    if (verbose) &
       write (uout,'("+ Atom type sequence: ",A)') lbl1

  end subroutine struct_write_vasp

  !> Write an abinit input template
  subroutine struct_write_abinit(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, string, fclose
    use param, only: pi, maxzat0

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    character(len=:), allocatable :: lbl1
    integer :: ntyp(maxzat0), iz
    real*8 :: aap(3), bbp(3), gpq(3,3)
    integer :: i, j, lu

    ! count number of atoms per type
    ntyp = 0
    do i = 1, c%ncel
       ntyp(c%at(c%atcel(i)%idx)%z) = ntyp(c%at(c%atcel(i)%idx)%z) + 1
    end do

    ! Find the lengths and angles of the cell
    gpq = matmul(transpose(c%crys2car),c%crys2car)
    do i = 1, 3
       aap(i) = sqrt(gpq(i,i))
    end do
    bbp(1)=acos(gpq(2,3)/sqrt(gpq(2,2)*gpq(3,3)))/pi*180d0
    bbp(2)=acos(gpq(1,3)/sqrt(gpq(1,1)*gpq(3,3)))/pi*180d0
    bbp(3)=acos(gpq(1,2)/sqrt(gpq(1,1)*gpq(2,2)))/pi*180d0

    ! Write input
    lu = fopen_write(file)
    write (lu,'("acell ",3(F14.10,X))') aap
    write (lu,'("angdeg ",3(F14.10,X))') bbp
    write (lu,'("ntypat ",I3)') count(ntyp > 0)

    lbl1 = ""
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          lbl1 = lbl1 // " " // string(i)
       end if
    end do
    write (lu,'("znucl ",A)') lbl1
    write (lu,'("natom ",I5)') c%ncel

    lbl1 = ""
    iz = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          iz = iz + 1
          lbl1 = lbl1 // " " // string(ntyp(i)) // "*" // string(iz)
       end if
    end do
    write (lu,'("typat ",A)') lbl1

    write (lu,'("xred ")')
    do i = 1, size(ntyp)
       do j = 1, c%ncel
          if (c%at(c%atcel(j)%idx)%z == i) then
             write (lu,'(X,3(F15.10,X))') c%atcel(j)%x
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

  end subroutine struct_write_abinit

  !> Write an elk input template
  subroutine struct_write_elk(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, nameguess, fclose
    use param, only: maxzat0
    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: ntyp(maxzat0)
    integer :: i, j, lu

    ! count number of atoms per type
    ntyp = 0
    do i = 1, c%ncel
       ntyp(c%at(c%atcel(i)%idx)%z) = ntyp(c%at(c%atcel(i)%idx)%z) + 1
    end do

    ! Write input
    lu = fopen_write(file)
    write (lu,'("tasks"/,"0"/)')
    write (lu,'("xctype"/,"20"/)')
    write (lu,'("avec")')
    do i = 1, 3
       write (lu,'(2X,3(F15.10,X))') c%crys2car(:,i)
    end do
    write (lu,*)

    write (lu,'("sppath"/,"''./''"/)')

    write (lu,'("atoms")')
    write (lu,'(2X,I4)') count(ntyp > 0)
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          write (lu,'(2X,"''",A,".in''")') trim(nameguess(i,.true.))
          write (lu,'(2X,I3)') ntyp(i)
          do j = 1, c%ncel
             if (c%at(c%atcel(j)%idx)%z == i) then
                write (lu,'(2X,3(F14.10,X),"0.0 0.0 0.0")') c%atcel(j)%x
             end if
          end do
       end if
    end do
    write (lu,*)

    write (lu,'("ngridk"/,"  4 4 4"/)')
    write (lu,'("rgkmax"/,"  7.0"/)')
    call fclose(lu)

  end subroutine struct_write_elk

  !> Write a Gaussian template input (periodic).
  subroutine struct_write_gaussian(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, string, nameguess, ioj_left, fclose
    use param, only: bohrtoa
    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    character(len=:), allocatable :: wroot
    integer :: lu, i, j

    wroot = file(:index(file,'.',.true.)-1)

    lu = fopen_write(file)
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
       write (lu,'(99(A,X))') string(nameguess(c%at(c%atcel(i)%idx)%z,.true.),2,ioj_left),&
          (string(c%atcel(i)%r(j)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    do i = 1, 3
       write (lu,'(99(A,X))') string("Tv",2,ioj_left),&
          (string(c%crys2car(j,i)*bohrtoa,'f',14,8,ioj_left),j=1,3)
    end do
    write (lu,*)

    call fclose(lu)

  end subroutine struct_write_gaussian

  !> Write a tessel input template
  subroutine struct_write_tessel(file,c)
    use struct_basic, only: crystal
    use global, only: fileroot
    use tools_io, only: fopen_write, fclose

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: lu, i

    lu = fopen_write(file)

    write (lu,'("set camangle 75 -10 45")')
    write (lu,'("set background background {color rgb <1,1,1>}")')
    write (lu,'("set use_planes .false.")')
    write (lu,'("set ball_texture finish{specular 0.2 roughness 0.1 reflection 0.1}")')
    write (lu,'("set equalscale noscale")')
    write (lu,'("molecule")')
    write (lu,'("  crystal")')
    write (lu,'("    symmatrix seitz")')
    do i = 1, c%ncv
       write (lu,'(5X,A,3(F15.12,X))') "cen ",c%cen(:,i)
    end do
    write (lu,'(5X,"#")')
    do i = 1, c%neqv
       write (lu,'(5X,3(F5.2,X),F15.12)') c%rotm(1,:,i)
       write (lu,'(5X,3(F5.2,X),F15.12)') c%rotm(2,:,i)
       write (lu,'(5X,3(F5.2,X),F15.12)') c%rotm(3,:,i)
       write (lu,'(5X,"#")')
    end do
    write (lu,'(5X,"endsymmatrix")')
    write (lu,'(5X,A,6(F12.8,X))') "cell", c%aa, c%bb
    write (lu,'(5X,"crystalbox  -2.30 -2.30 -2.30 2.30 2.30 2.30")')
    write (lu,'(5X,A,6(F6.3,X))') "clippingbox ",-0.02,-0.02,-0.02,+1.02,+1.02,+1.02
    do i = 1, c%nneq
       write (lu,'(5X,"neq ",3(F12.8," "),A10)') c%at(i)%x,trim(c%at(i)%name)
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

  end subroutine struct_write_tessel

  !> Write a critic2 input template
  subroutine struct_write_critic(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, fclose

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: lu, i

    lu = fopen_write(file)

    write (lu,'("crystal")')
    write (lu,'("  cell ",3(F15.11,X),3(F9.5,X))') c%aa, c%bb
    do i = 1, c%ncel
       write (lu,'("  neq ",3(F12.8," "),A10)') c%atcel(i)%x,&
          trim(c%at(c%atcel(i)%idx)%name)
    end do
    write (lu,'("endcrystal")')
    write (lu,'("end")')
    call fclose(lu)

  end subroutine struct_write_critic

  !> Write a simple cif file
  subroutine struct_write_cif(file,c)
    use struct_basic, only: crystal
    use global, only: fileroot
    use tools_io, only: fopen_write, fclose, string, nameguess
    use param, only: bohrtoa

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: i, iz, lu

    lu = fopen_write(file)

    write (lu,'("data_",A)') string(fileroot)
    write (lu,'("_cell_volume ",F20.6)') c%omega * bohrtoa**3
    write (lu,'("_symmetry_space_group_name_H-M ''P 1''")');
    write (lu,'("_symmetry_Int_Tables_number 1")');
    write (lu,'("loop_")');
    write (lu,'("_symmetry_equiv_pos_site_id")');
    write (lu,'("_symmetry_equiv_pos_as_xyz")');
    write (lu,'("1 x,y,z")');
    write (lu,'("_cell_length_a ",F20.10)') c%aa(1)*bohrtoa
    write (lu,'("_cell_length_b ",F20.10)') c%aa(2)*bohrtoa
    write (lu,'("_cell_length_c ",F20.10)') c%aa(3)*bohrtoa
    write (lu,'("_cell_angle_alpha ",F14.4)') c%bb(1)
    write (lu,'("_cell_angle_beta ",F14.4)') c%bb(2)
    write (lu,'("_cell_angle_gamma ",F14.4)') c%bb(3)
    write (lu,'("_cell_formula_units_Z 1")')
    write (lu,'("loop_")');
    write (lu,'("_atom_site_label")');
    write (lu,'("_atom_site_type_symbol")');
    write (lu,'("_atom_site_fract_x")');
    write (lu,'("_atom_site_fract_y")');
    write (lu,'("_atom_site_fract_z")');
    do i = 1, c%ncel
       iz = c%atcel(i)%idx
       write (lu,'(A5,X,A3,X,3(F20.14,X))') c%at(iz)%name, &
          nameguess(c%at(iz)%z,.true.), c%atcel(i)%x
    end do
    call fclose(lu)

  end subroutine struct_write_cif

  !> Write a simple cif file
  subroutine struct_write_d12(file,c,dosym)
    use struct_basic, only: crystal, pointgroup_info, holo_unk, holo_tric,&
       holo_mono, holo_ortho, holo_tetra, holo_trig, holo_hex, holo_cub
    use tools_io, only: fopen_write, fclose, string, ferror, faterr
    use global, only: symprec
    use param, only: bohrtoa, pi

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    logical, intent(in) :: dosym

    integer :: holo, laue, i, j, nmin
    character(len=3) :: schpg
    type(crystal) :: nc
    integer :: lu, spgnum, irhomb, count90, count120
    real*8 :: xmin(6)
    logical :: ok

    ! we need symmetry for this
    if (dosym) then
       nc = c
       call nc%cell_standard(.false.,.false.,.false.)
       call pointgroup_info(nc%spg%pointgroup_symbol,schpg,holo,laue)
       xmin = 0d0
       irhomb = 0
       if (holo == holo_unk) then
          call ferror("struct_write_d12","unknown holohedry",faterr)
       elseif (holo == holo_tric) then
          nmin = 6
          xmin(1:3) = nc%aa * bohrtoa
          xmin(4:6) = nc%bb
       elseif (holo == holo_mono) then
          nmin = 4
          xmin(1:3) = nc%aa * bohrtoa
          ok = .false.
          do i = 1, 3
             if (abs(nc%bb(i) - 90d0) > symprec) then
                xmin(4) = nc%bb(i)
                ok = .true.
                exit
             end if
          end do
          if (.not.ok) &
             xmin(4) = 90d0
       elseif (holo == holo_ortho) then
          nmin = 3
          xmin(1:3) = nc%aa * bohrtoa
       elseif (holo == holo_tetra) then
          nmin = 2
          xmin(1) = nc%aa(1) * bohrtoa
          xmin(2) = nc%aa(3) * bohrtoa
       elseif (holo == holo_trig) then
          count90 = count(abs(nc%bb - 90d0) < 1d-1)
          count120 = count(abs(sin(nc%bb * pi / 180d0) - sqrt(3d0)/2d0) < 1d-2)
          nmin = 2
          xmin(1) = nc%aa(1) * bohrtoa
          if (count90 == 2 .and. count120 == 1) then
             ! hexagonal axes
             xmin(2) = nc%aa(3) * bohrtoa
             irhomb = 0
          else
             ! rhombohedral axes
             xmin(2) = nc%bb(1)
             irhomb = 1
          end if
       elseif (holo == holo_hex) then
          nmin = 2
          xmin(1) = nc%aa(1) * bohrtoa
          xmin(2) = nc%aa(3) * bohrtoa
       elseif (holo == holo_cub) then
          nmin = 1
          xmin(1) = nc%aa(1) * bohrtoa
       end if
       spgnum = nc%spg%spacegroup_number
    else
       spgnum = 1
       irhomb = 0
       nmin = 6
       xmin(1:3) = c%aa * bohrtoa
       xmin(4:6) = c%bb
    end if

    lu = fopen_write(file)
    write (lu,'("Title")')
    write (lu,'("CRYSTAL")')
    write (lu,'("0 ",A," 0")') string(irhomb)
    write (lu,'(A)') string(spgnum)
    write (lu,'(6(A,X))') (string(xmin(i),'f',15,8),i=1,nmin)
    if (dosym) then
       write (lu,'(A)') string(nc%nneq)
       do i = 1, nc%nneq
          write (lu,'(4(A,X))') string(nc%at(i)%z), (string(nc%at(i)%x(j),'f',15,8),j=1,3)
       end do
    else
       write (lu,'(A)') string(c%ncel)
       do i = 1, c%ncel
          write (lu,'(4(A,X))') string(c%at(c%atcel(i)%idx)%z), (string(c%atcel(i)%x(j),'f',15,8),j=1,3)
       end do
    end if
    write (lu,'("SETPRINT")')
    write (lu,'("1")')
    write (lu,'("3 1")')
    write (lu,'("END")')
    write (lu,'("xx basis xx")')
    write (lu,'("99 0")')
    write (lu,'("END")')
    write (lu,'("SHRINK")')
    write (lu,'("4 4")')
    write (lu,'("TOLDEE")')
    write (lu,'("7")')
    write (lu,'("END")')

    call fclose(lu)

  end subroutine struct_write_d12

  !> Write an escher octave script
  subroutine struct_write_escher(file,c)
    use struct_basic, only: crystal
    use global, only: fileroot
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: pi, maxzat0

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    character(len=:), allocatable :: lbl1
    integer :: lu, i, j, n
    integer :: ntyp(maxzat0)

    lu = fopen_write(file)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, c%ncel
       ntyp(c%at(c%atcel(i)%idx)%z) = ntyp(c%at(c%atcel(i)%idx)%z) + 1
    end do

    write (lu,'("cr = struct();")')
    write (lu,'("cr.name = """,A,""";")') trim(adjustl(fileroot))
    write (lu,'("cr.a = [",1p,3(E22.14,X),"];")') c%aa
    write (lu,'("cr.b = [",1p,3(E22.14,X),"];")') c%bb * pi / 180d0
    write (lu,'("cr.nat = ",I6,";")') c%ncel
    write (lu,'("cr.ntyp = ",I6,";")') count(ntyp > 0)
    write (lu,'("cr.r = [")')
    do i = 1, 3
       write (lu,'(2X,1p,3(E22.14,X))') c%crys2car(:,i)
    end do
    write (lu,'(2X,"];")')
    write (lu,'("cr.g = [")')
    do i = 1, 3
       write (lu,'(2X,1p,3(E22.14,X))') c%gtensor(:,i)
    end do
    write (lu,'(2X,"];")')
    write (lu,'("cr.omega = ",1p,E22.14,";")') c%omega

    lbl1 = "cr.ztyp = ["
    n = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          lbl1 = lbl1 // " " // string(i)
       end if
    end do
    lbl1 = lbl1 // "];"
    write (lu,'(A)') lbl1

    lbl1 = "cr.attyp = {"
    n = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          n = n + 1
          if (n > 1) lbl1 = lbl1 // ","
          lbl1 = lbl1 // '"' // string(nameguess(i,.true.)) // '"'
       end if
    end do
    lbl1 = lbl1 // "};"
    write (lu,'(A)') lbl1

    lbl1 = "cr.typ = ["
    n = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          n = n + 1
          do j = 1, ntyp(i)
             lbl1 = lbl1 // " " // string(n)
          end do
       end if
    end do
    lbl1 = lbl1 // "];"
    write (lu,'(A)') lbl1

    write (lu,'("cr.x = [")')
    n = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          do j = 1, c%ncel
             if (c%at(c%atcel(j)%idx)%z == i) then
                write (lu,'(2X,1p,3(E22.14,X))') c%atcel(j)%x
             endif
          end do
       end if
    end do
    write (lu,'("  ];")')

    call fclose(lu)

  end subroutine struct_write_escher

  !> Write a gulp input script
  subroutine struct_write_gulp(file,c,dodreiding)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, faterr, ferror, nameguess, fclose
    use tools_math, only: norm
    use param, only: bohrtoa, atmcov, pi

    character*(*), intent(in) :: file
    type(crystal), intent(inout) :: c
    logical :: dodreiding

    integer, parameter :: maxneigh = 20
    real*8, parameter :: rfac = 1.4d0
    real*8, parameter :: hbmin = 1.6 / bohrtoa
    real*8, parameter :: hbmax = 3.2 / bohrtoa

    character*(5) :: lbl
    integer :: lu, i, j, iz, jz, n, k, kz
    integer :: idx
    integer :: nneigh(maxneigh), ineigh(maxneigh,c%nneq)
    real*8 :: dneigh(maxneigh,c%nneq), dhb(maxneigh,c%nneq), avgang(c%nneq)
    integer :: nhb(maxneigh), ihb(maxneigh,c%nneq)
    real*8 :: d, x1(3), x2(3), ang
    logical :: ok, isat

    ! check that we have an environment
    call c%checkflags(.true.,env0=.true.)

    lu = fopen_write(file)
    if (.not. dodreiding) then
       write (lu,'("eem")')
       write (lu,'("cell ",3(F13.9,X),3(F10.5,X))') c%aa * bohrtoa, c%bb
       write (lu,'("fractional")')
       do i = 1, c%ncel
          write (lu,'(A5,X,3(F15.9,X))') trim(c%at(c%atcel(i)%idx)%name),&
             c%atcel(i)%x
       end do
    else
       ! calculate bonded neighbors
       nneigh = 0
       nhb = 0
       do i = 1, c%nneq
          iz = c%at(i)%z
          n = 0
          ! determine covalent bonds
          do j = 1, c%nenv
             jz = c%at(c%atenv(j)%idx)%z
             d = sqrt(dot_product(c%atenv(j)%r-c%at(i)%r,c%atenv(j)%r-c%at(i)%r))
             if (d < 1d-10) cycle
             if (d < (atmcov(iz) + atmcov(jz)) * rfac) then
                n = n + 1
                if (n > maxneigh) call ferror("struct_write_gulp","too many neighbors",faterr)
                ineigh(n,i) = j
                dneigh(n,i) = d
                continue
             endif
          end do
          nneigh(i) = n
          ! determine average angles with bonded neighbors
          avgang(i) = 0d0
          n = 0
          do j = 1, nneigh(i)
             do k = j+1, nneigh(i)
                n = n + 1
                x1 = c%atenv(ineigh(j,i))%r - c%at(i)%r
                x2 = c%atenv(ineigh(k,i))%r - c%at(i)%r
                ang = abs(acos(dot_product(x1,x2) / norm(x1) / norm(x2)) * 180d0 / pi)
                avgang(i) = avgang(i) + ang
             end do
          end do
          if (n > 0) avgang(i) = avgang(i) / n

          ! determine hydrogen bonds, only for hydrogen
          if (iz == 1) then
             n = 0
             do j = 1, c%nenv
                jz = c%at(c%atenv(j)%idx)%z
                ! only with N, O, and S
                if (jz==7 .or. jz==8 .or. jz==9 .or. jz==16 .or. jz==17 .or. jz==35 .or. jz==53) then
                   d = sqrt(dot_product(c%atenv(j)%r-c%at(i)%r,c%atenv(j)%r-c%at(i)%r))
                   ! only in the correct distance range
                   if (d > hbmin .and. d < hbmax) then
                      ! only if the angles to all other neighbor atoms is more than 145
                      ok = .true.
                      do k = 1, nneigh(i)
                         x1 = c%atenv(ineigh(k,i))%r - c%at(i)%r
                         x2 = c%atenv(j)%r - c%at(i)%r
                         ang = abs(acos(dot_product(x1,x2) / norm(x1) / norm(x2)) * 180d0 / pi)
                         kz = c%at(c%atenv(ineigh(k,i))%idx)%z
                         isat = (kz==7 .or. kz==8 .or. kz==9 .or. kz==16 .or. kz==17 .or. kz==35 .or. kz==53)
                         if (ang < 145d0 .or..not.isat) then
                            ok = .false.
                            exit
                         end if
                      end do
                      ! oh, sure, fine... you're a hydrogen bond
                      if (ok) then
                         n = n + 1
                         ihb(n,i) = j
                         dhb(n,i) = d
                      endif
                   end if
                end if
             end do
             nhb(i) = n
          endif
       end do

       write (lu,'("eem")')
       write (lu,'("cell ",3(F13.9,X),3(F10.5,X))') c%aa * bohrtoa, c%bb
       write (lu,'("fractional")')
       do i = 1, c%ncel
          idx = c%atcel(i)%idx
          iz = c%at(idx)%z
          ang = avgang(idx)
          ! the first two letters is the atomic symbol
          lbl = adjustl(nameguess(iz))
          ! hydrogen types: H_ (normal), H___A (hydrogen-bonded to N, O, or S), H___b (bridging)
          if (iz == 1 .and. nneigh(idx) > 1) lbl = "H___b"
          if (iz == 1 .and. nhb(idx) > 0) lbl = "H___A"
          ! boron: sp3 (109.47 angles) and sp2 (120 angles)
          if (iz == 5 .and. abs(ang-109.47d0)<abs(ang-120d0)) lbl = "B_3"
          if (iz == 5 .and. abs(ang-109.47d0)>abs(ang-120d0)) lbl = "B_2"
          ! carbon: sp3 (109.47 angles), sp2 (120 angles), sp (180 angles)
          if (iz == 6 .and. abs(ang-109.47d0)<abs(ang-120d0) .and. abs(ang-109.47d0)<abs(ang-180d0)) lbl = "C_3"
          if (iz == 6 .and. abs(ang-120d0)<abs(ang-109.47d0) .and. abs(ang-120d0)<abs(ang-180d0)) lbl = "C_2"
          if (iz == 6 .and. abs(ang-180d0)<abs(ang-109.47d0) .and. abs(ang-180d0)<abs(ang-120d0)) lbl = "C_1"
          ! nitrogen: sp3 (109.47 angles), sp2 (120 angles), sp (180 angles)
          if (iz == 7 .and. abs(ang-109.47d0)<abs(ang-120d0) .and. abs(ang-109.47d0)<abs(ang-180d0)) lbl = "N_3"
          if (iz == 7 .and. abs(ang-120d0)<abs(ang-109.47d0) .and. abs(ang-120d0)<abs(ang-180d0)) lbl = "N_2"
          if (iz == 7 .and. (abs(ang-180d0)<abs(ang-109.47d0) .and. abs(ang-180d0)<abs(ang-120d0) .or. nneigh(idx) == 1)) lbl = "N_1"
          ! oxygen: sp3 (109.47 angles), sp2 (120 angles), sp (180 angles)
          if (iz == 8 .and. abs(ang-109.47d0)<abs(ang-120d0) .and. abs(ang-109.47d0)<abs(ang-180d0)) lbl = "O_3"
          if (iz == 8 .and. (abs(ang-120d0)<abs(ang-109.47d0) .and. abs(ang-120d0)<abs(ang-180d0) .or. nneigh(idx) == 1)) lbl = "O_2"
          if (iz == 8 .and. abs(ang-180d0)<abs(ang-109.47d0) .and. abs(ang-180d0)<abs(ang-120d0)) lbl = "O_1"
          ! Al, Si, P, S, Ga, Ge, As, Se, In, Sn, Sb, Te -> only sp3 is known
          if (iz == 13 .or. iz == 14 .or. iz == 15 .or. iz == 16 .or. iz == 31 .or. iz == 32 .or. iz == 33 .or. iz == 34 .or.&
             iz == 49 .or. iz == 50 .or. iz == 51 .or. iz == 52) then
             lbl(3:3) = "3"
          end if
          write (lu,'(A5,X,3(F15.9,X))') adjustl(trim(lbl)), c%atcel(i)%x
       end do
    end if
    call fclose(lu)

  end subroutine struct_write_gulp

  !> Write a lammps data file
  subroutine struct_write_lammps(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, ferror, faterr, fclose
    use tools_math, only: crys2car_from_cellpar
    use param, only: bohrtoa, maxzat0, atmass

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: i, j, k, l
    integer :: ntyp(maxzat0), lu
    real*8 :: rnew(3,3)

    lu = fopen_write(file)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, c%ncel
       ntyp(c%at(c%atcel(i)%idx)%z) = ntyp(c%at(c%atcel(i)%idx)%z) + 1
    end do

    ! header
    write (lu,'("LAMMPS data file created by critic2. (experimental)",/)')
    write (lu,'(I9," atoms")') c%ncel
    write (lu,'(I9," atom types")') count(ntyp > 0)
    write (lu,*)

    ! metrics of the cell --> this needs more testing
    rnew = crys2car_from_cellpar(c%aa,c%bb)
    if (abs(c%crys2car(1,2)) > 1d-12 .or. abs(c%crys2car(1,3)) > 1d-12 .or.&
       abs(c%crys2car(2,3)) > 1d-12) then
       call ferror ('struct_write_lammps','non-orthogonal cells not implemented',faterr)
    end if
    write (lu,'(2(F18.10,X)," xlo xhi")') 0d0, c%crys2car(1,1)*bohrtoa
    write (lu,'(2(F18.10,X)," ylo yhi")') 0d0, c%crys2car(2,2)*bohrtoa
    write (lu,'(2(F18.10,X)," zlo zhi")') 0d0, c%crys2car(3,3)*bohrtoa
    write (lu,'(3(F18.10,X)," xy xz yz")') 0d0, 0d0, 0d0
    write (lu,*)

    write (lu,'("Masses"/)')
    j = 0
    do i = 1, maxzat0
       if (ntyp(i) > 0) then
          j = j + 1
          write (lu,'(I3,X,F10.4)') j, atmass(i)
       end if
    end do
    write (lu,*)

    write (lu,'("Atoms"/)')
    k = 0
    l = 0
    do i = 1, maxzat0
       if (ntyp(i) == 0) cycle
       k = k + 1
       do j = 1, c%ncel
          if (c%at(c%atcel(j)%idx)%z /= i) cycle
          l = l + 1
          write (lu,'(I7,X,I3,X,F4.1,3(F15.8,X))') l, k, 0d0, c%atcel(j)%r*bohrtoa
       end do
    end do

    call fclose(lu)

  end subroutine struct_write_lammps

  !> Write a siesta fdf data file
  subroutine struct_write_siesta_fdf(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, nameguess, lower, fclose
    use param, only: bohrtoa, maxzat0
    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: i, j, k
    integer :: ntyp(maxzat0), lu, nspecies

    lu = fopen_write(file)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, c%ncel
       ntyp(c%at(c%atcel(i)%idx)%z) = ntyp(c%at(c%atcel(i)%idx)%z) + 1
    end do
    nspecies = count(ntyp > 0)

    ! header
    write (lu,'("# fdf file created by critic2.",/)')
    write (lu,'("SystemName crystal")') 
    write (lu,'("SystemLabel crystal")') 
    write (lu,*)

    write (lu,'("NumberOfSpecies ",I3)') nspecies
    write (lu,'("NumberOfAtoms ", I6)') c%ncel
    write (lu,'("%block Chemical_Species_Label")') 
    j = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          j = j + 1
          write (lu,'(I3,I3,X,A2)') j, i, lower(nameguess(i,.true.))
       end if
    end do
    write (lu,'("%endblock Chemical_Species_Label")') 
    write (lu,*)

    write (lu,'("LatticeConstant 1.0 ang")')
    write (lu,'("%block LatticeParameters")')
    write (lu,'(3(F16.10,X),3(F16.8,X))') c%aa*bohrtoa, c%bb
    write (lu,'("%endblock LatticeParameters")')
    write (lu,'("AtomicCoordinatesFormat Fractional")')
    write (lu,'("%block AtomicCoordinatesAndAtomicSpecies")')
    k = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          k = k + 1
          do j = 1, c%ncel
             if (c%at(c%atcel(j)%idx)%z == i) then
                write (lu,'(3(F18.12,X),I3)') c%atcel(j)%x, k
             endif
          end do
       end if
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

  end subroutine struct_write_siesta_fdf

  !> Write a siesta STRUCT_IN data file
  subroutine struct_write_siesta_in(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, uout, nameguess, string, fclose
    use param, only: bohrtoa, maxzat0

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, j, k, ntyp(maxzat0), nspecies

    lu = fopen_write(file)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, c%ncel
       ntyp(c%at(c%atcel(i)%idx)%z) = ntyp(c%at(c%atcel(i)%idx)%z) + 1
    end do
    nspecies = count(ntyp > 0)

    ! lattice vectors
    r = transpose(c%crys2car) * bohrtoa
    do i = 1, 3
       write (lu,'(3(F20.12,X))') r(i,:)
    end do

    ! atoms
    write (lu,*) c%ncel
    j = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          j = j + 1
          do k = 1, c%ncel
             if (c%at(c%atcel(k)%idx)%z == i) then
                write (lu,'(I3,X,I3,X,3(F20.12,X))') j, i, c%atcel(k)%x
             end if
          end do
       end if
    end do

    call fclose(lu)

    ! Write the chemical species block to standard output
    write (uout,'("%block Chemical_Species_Label")') 
    j = 0
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          j = j + 1
          write (uout,'(3(2X,A))') string(j), string(i), &
             string(nameguess(i,.true.))
       end if
    end do
    write (uout,'("%endblock Chemical_Species_Label")') 
    write (uout,*)

  end subroutine struct_write_siesta_in

  !> Write a DFTB+ human-friendly structured data format (hsd) file
  subroutine struct_write_dftbp_hsd(file,c)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, string, nameguess, fclose
    use param, only: maxzat0

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    logical :: ltyp(maxzat0)

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

    lu = fopen_write(file)
    write (lu,'("Geometry = GenFormat {")')
    call struct_write_dftbp_gen(file,c,lu,ltyp)
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
    do i = 1, size(ltyp)
       if (ltyp(i)) then
          write (lu,'(4X,A," = ",A)') string(nameguess(i,.true.)), &
             string(maxang(i))
       end if
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
    do i = 1, size(ltyp)
       if (ltyp(i)) then
          write (lu,'(4X,A," = ",A)') string(nameguess(i,.true.)), &
             string(hderiv(i),'f',decimal=4)
       end if
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

  end subroutine struct_write_dftbp_hsd

  !> Write a DFTB+ human-friendly gen structure file
  subroutine struct_write_dftbp_gen(file,c,lu0,ltyp0)
    use struct_basic, only: crystal
    use tools_io, only: fopen_write, nameguess, string, fclose
    use param, only: bohrtoa, maxzat0

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    integer, intent(in), optional :: lu0
    logical, intent(out), optional :: ltyp0(maxzat0)

    integer :: lu, nspecies, n, nt, i, j, k
    logical :: ltyp(maxzat0)
    real*8 :: r(3,3)
    character(len=:), allocatable :: strtyp

    ! count atoms types
    ltyp = .false.
    do i = 1, c%ncel
       ltyp(c%at(c%atcel(i)%idx)%z) = .true.
    end do
    nspecies = count(ltyp)
    if (present(ltyp0)) ltyp0 = ltyp

    ! open file
    if (present(lu0)) then
       lu = lu0
    else
       lu = fopen_write(file)
    end if

    ! atom types
    strtyp = ""
    do i = 1, size(ltyp)
       if (ltyp(i)) then
          strtyp = strtyp // " " // string(nameguess(i,.true.))
       end if
    end do

    if (c%ismolecule) then
       ! molecule
       write (lu,'(A," C")') string(c%ncel)
       write (lu,'(A)') strtyp

       ! Cartesian coordinates
       n = 0
       nt = 0
       do i = 1, size(ltyp)
          if (ltyp(i)) then
             nt = nt + 1
             do k = 1, c%ncel
                if (c%at(c%atcel(k)%idx)%z == i) then
                   n = n + 1
                   write (lu,'(99(A,X))') string(n), string(nt), &
                      (string(c%atcel(k)%r(j)*bohrtoa,'f',20,12),j=1,3)
                end if
             end do
          end if
       end do
    else
       ! crystal
       write (lu,'(A," F")') string(c%ncel)
       write (lu,'(A)') strtyp

       ! fractional coordinates
       n = 0
       nt = 0
       do i = 1, size(ltyp)
          if (ltyp(i)) then
             nt = nt + 1
             do k = 1, c%ncel
                if (c%at(c%atcel(k)%idx)%z == i) then
                   n = n + 1
                   write (lu,'(99(A,X))') string(n), string(nt), &
                      (string(c%atcel(k)%x(j),'f',20,12),j=1,3)
                end if
             end do
          end if
       end do

       ! lattice vectors
       r = c%crys2car * bohrtoa
       write (lu,'(3(A,X))') (string(0d0,'f',20,12),j=1,3)
       do i = 1, 3
          write (lu,'(3(A,X))') (string(r(j,i),'f',20,12),j=1,3)
       end do
    endif

    ! close file
    if (.not.present(lu0)) call fclose(lu)

  end subroutine struct_write_dftbp_gen

end module struct_writers
