! Copyright (c) 2015 Alberto Otero de la Roza
! <alberto@fluor.quimica.uniovi.es>,
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
  public :: struct_write_escher
  public :: struct_write_gulp
  public :: struct_write_lammps
  public :: struct_write_siesta_fdf
  public :: struct_write_siesta_in
  private :: cell_for_espresso
  private :: crystal_system

  ! crystal system flags
  integer, parameter :: csys_cub = 1
  integer, parameter :: csys_hex = 2
  integer, parameter :: csys_trig = 3
  integer, parameter :: csys_tetr = 4
  integer, parameter :: csys_orth = 5
  integer, parameter :: csys_mono = 6
  integer, parameter :: csys_tric = 7

  ! transformations to primitive cell for different centering types
  real*8, parameter :: toc_real(3,3,9) = reshape((/ 1.d0, 0.d0, 0.d0,&
     & 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0,& ! P (1)
     1.d0, 0.d0, 0.d0, 0.d0, .5d0, .5d0, 0.d0,-.5d0, .5d0,& ! A (2)
     .5d0, 0.d0,-.5d0, 0.d0, 1.d0, 0.d0, .5d0, 0.d0, .5d0,& ! B (3)
     .5d0, .5d0, 0.d0, -.5d0, .5d0, 0.d0, 0.d0, 0.d0, 1.d0,& ! C (4)
    -.5d0, .5d0, .5d0, .5d0,-.5d0, .5d0, .5d0, .5d0,-.5d0,& ! I (5)
     0.d0, .5d0, .5d0, .5d0, 0.d0, .5d0, .5d0, .5d0, 0.d0,& ! F (6)
     1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0,& ! R (7)
     2d0/3d0,  1d0/3d0,  1d0/3d0, -1d0/3d0,  1d0/3d0,  1d0/3d0, -1d0&
     &/3d0, -2d0/3d0,  1d0/3d0,& ! R+ (8)
    -2d0/3d0, -1d0/3d0,  1d0/3d0, 1d0/3d0, -1d0/3d0,  1d0/3d0, 1d0&
    &/3d0,  2d0/3d0,  1d0/3d0 & ! R- (9)
     /),shape(toc_real)) !< transformation to primitive real cell
  real*8, parameter :: toc_rec(3,3,9) = reshape((/ 1.d0, 0.d0, 0.d0,&
     & 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0,& ! P (1)
     1.d0, 0.d0, 0.d0, 0.d0, 1.d0,-1.d0, 0.d0, 1.d0, 1.d0,& ! A (2)
     1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0,& ! B (3)
     1.d0,-1.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0,& ! C (4)
     0.d0, 1.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0, 0.d0,& ! I (5)
    -1.d0, 1.d0, 1.d0, 1.d0,-1.d0, 1.d0, 1.d0, 1.d0,-1.d0,& ! F (6)
     1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0,& ! R (7)
     1.d0,-1.d0, 0.d0, 0.d0, 1.d0,-1.d0, 1.d0, 1.d0, 1.d0,& ! R+ (8)
    -1.d0, 1.d0, 0.d0, 0.d0,-1.d0, 1.d0, 1.d0, 1.d0, 1.d0 & ! R- (9)
     /),shape(toc_rec)) !< transformation to primitive reciprocal cell

contains

  !> Write an xyz/gjf file containing a finite piece of the crystal
  !> structure. fmt can be one of xyz or gjf. ix is the number of
  !> unit cells to plot.  If doborder is .true., add all atoms at the
  !> border. If molmotif is .true., complete molecules with atoms in
  !> adjacent cells. If docell, add sticks for the unit cell limits. If
  !> rsph (bohr) is positive, then use all atoms in a sphere around xsph
  !> (cryst.). If rcub (bohr) is positive, use all atoms in a cube
  !> around xcub (cryst.). 
  subroutine struct_write_mol(c,file,fmt,ix,doborder,molmotif,doburst,dopairs,rsph,xsph,rcub,xcub)
    use fragmentmod
    use struct_basic
    use graphics
    use tools_io
    use types
    use param

    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in) :: ix(3)
    logical, intent(in) :: doborder, molmotif, doburst, dopairs
    real*8, intent(in) :: rsph, xsph(3)
    real*8, intent(in) :: rcub, xcub(3)

    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: i, j, nmol
    character(len=:), allocatable :: wroot, file0

    if (rcub > 0) then
       fr = c%listatoms_sphcub(rcub=rcub,xcub=xcub)
    elseif (rsph > 0) then
       fr = c%listatoms_sphcub(rsph=rsph,xsph=xsph)
    elseif (doburst.or.dopairs) then
       fr = c%listatoms_cells(ix,doborder)
       call c%listmolecules(fr,nmol,fr0,isdiscrete)
    elseif (molmotif) then
       fr = c%listatoms_cells(ix,doborder)
       call c%listmolecules(fr,nmol,fr0,isdiscrete)
       fr = fragment_merge_array(fr0)
    else
       fr = c%listatoms_cells(ix,doborder)
    endif

    if (.not.doburst.and..not.dopairs) then
       if (equal(fmt,"xyz")) then
          call writexyz(file,fr)
       elseif (equal(fmt,"gjf")) then
          call writegjf(file,fr)
       else
          call ferror("struct_write_mol","Unknown format",faterr)
       endif
    else
       if (doburst) then
          wroot = file(:index(file,'.',.true.)-1)
          do i = 1, nmol
             file0 = wroot // "_" // string(i) // "." // fmt
             if (equal(fmt,"xyz")) then
                call writexyz(file0,fr0(i))
             elseif (equal(fmt,"gjf")) then
                call writegjf(file0,fr0(i))
             else
                call ferror("struct_write_mol","Unknown format",faterr)
             endif
          end do
       end if
       if (dopairs) then
          wroot = file(:index(file,'.',.true.)-1)
          do i = 1, nmol
             do j = i+1, nmol
                file0 = wroot // "_" // string(i) // "_" // string(j) // "." // fmt
                fr = fragment_merge_array((/fr0(i),fr0(j)/))
                if (equal(fmt,"xyz")) then
                   call writexyz(file0,fr)
                elseif (equal(fmt,"gjf")) then
                   call writegjf(file0,fr)
                else
                   call ferror("struct_write_mol","Unknown format",faterr)
                endif
             end do
          end do
       end if
    end if
  end subroutine struct_write_mol

  !> Write an obj file containing the crystal structure. fmt can be
  !> one of obj, ply, or off. ix is the number of unit cells to plot.
  !> If doborder is .true., add all atoms at the border. If molmotif is
  !> .true., complete molecules with atoms in adjacent cells. If
  !> docell, add sticks for the unit cell limits. If rsph (bohr) is positive,
  !> then use all atoms in a sphere around xsph (cryst.). If rcub (bohr) is
  !> positive, use all atoms in a cube around xcub (cryst.). If lu0 and
  !> lumtl0 are present, then return the logical units for the obj and
  !> the mtl files and do not close the files.
  subroutine struct_write_3dmodel(c,file,fmt,ix,doborder,molmotif,doburst,&
     docell,domolcell,rsph,xsph,rcub,xcub,lu0,lumtl0)
    use graphics
    use struct_basic
    use tools_math
    use types
    use tools_io
    use param

    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*3, intent(in) :: fmt
    integer, intent(in) :: ix(3)
    logical, intent(in) :: doborder, molmotif, doburst, docell, domolcell
    real*8, intent(in) :: rsph, xsph(3)
    real*8, intent(in) :: rcub, xcub(3)
    integer, intent(out), optional :: lu0, lumtl0

    integer :: lu, lumtl
    integer :: i, j
    real*8 :: d, xd(3), x0(3), x1(3)
    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    integer :: k, nmol
    character(len=:), allocatable :: wroot, file0

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
    nmol = 1
    if (rcub > 0) then
       fr = c%listatoms_sphcub(rcub=rcub,xcub=xcub)
    elseif (rsph > 0) then
       fr = c%listatoms_sphcub(rsph=rsph,xsph=xsph)
    elseif (doburst .or. molmotif) then
       fr = c%listatoms_cells(ix,doborder)
       call c%listmolecules(fr,nmol,fr0,isdiscrete)
    else
       fr = c%listatoms_cells(ix,doborder)
    endif

    if (.not.doburst) then
       file0 = file
       if (.not.molmotif) then
          allocate(fr0(1))
          fr0(1) = fr
       end if
    else
       wroot = file(:index(file,'.',.true.)-1)
    end if

    if (doburst) then
       file0 = wroot // "_" // string(k) // "." // fmt
    end if

    lumtl = 0
    if (equal(fmt,"obj")) then
       call obj_open(file0,lu,lumtl)
    elseif (equal(fmt,"ply")) then
       call ply_open(file0,lu)
    elseif (equal(fmt,"off")) then
       call off_open(file0,lu)
    endif

    do k = 1, nmol
       ! add the balls
       do i = 1, fr0(k)%nat
          if (equal(fmt,"obj")) then
             call obj_ball(lu,fr0(k)%at(i)%r,JMLcol(:,fr0(k)%at(i)%z),0.6d0*atmcov(fr0(k)%at(i)%z))
          elseif (equal(fmt,"ply")) then
             call ply_ball(lu,fr0(k)%at(i)%r,JMLcol(:,fr0(k)%at(i)%z),0.6d0*atmcov(fr0(k)%at(i)%z))
          elseif (equal(fmt,"off")) then
             call off_ball(lu,fr0(k)%at(i)%r,JMLcol(:,fr0(k)%at(i)%z),0.6d0*atmcov(fr0(k)%at(i)%z))
          end if
       end do

       ! add the sticks
       do i = 1, fr0(k)%nat
          do j = i+1, fr0(k)%nat
             xd = fr0(k)%at(i)%r - fr0(k)%at(j)%r
             d = norm(xd)
             if (d < (atmcov(fr0(k)%at(i)%z) + atmcov(fr0(k)%at(j)%z)) * rfac) then
                xd = fr0(k)%at(i)%r + 0.5d0 * (fr0(k)%at(j)%r - fr0(k)%at(i)%r)
                if (equal(fmt,"obj")) then
                   call obj_stick(lu,fr0(k)%at(i)%r,xd,JMLcol(:,fr0(k)%at(i)%z),0.05d0)
                   call obj_stick(lu,fr0(k)%at(j)%r,xd,JMLcol(:,fr0(k)%at(j)%z),0.05d0)
                elseif (equal(fmt,"ply")) then
                   call ply_stick(lu,fr0(k)%at(i)%r,xd,JMLcol(:,fr0(k)%at(i)%z),0.05d0)
                   call ply_stick(lu,fr0(k)%at(j)%r,xd,JMLcol(:,fr0(k)%at(j)%z),0.05d0)
                elseif (equal(fmt,"off")) then
                   call off_stick(lu,fr0(k)%at(i)%r,xd,JMLcol(:,fr0(k)%at(i)%z),0.05d0)
                   call off_stick(lu,fr0(k)%at(j)%r,xd,JMLcol(:,fr0(k)%at(j)%z),0.05d0)
                end if
             end if
          end do
       end do

       ! add the cell
       if (docell) then
          do i = 1, 12
             x0 = c%x2c(x0cell(:,1,i))
             x1 = c%x2c(x0cell(:,2,i))
             if (equal(fmt,"obj")) then
                call obj_stick(lu,x0,x1,(/255,0,0/),0.03d0)
             elseif (equal(fmt,"ply")) then
                call ply_stick(lu,x0,x1,(/255,0,0/),0.03d0)
             elseif (equal(fmt,"off")) then
                call off_stick(lu,x0,x1,(/255,0,0/),0.03d0)
             end if
          end do
       end if
    end do

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
          x0 = c%x2c(x0)
          x1 = c%x2c(x1)
          if (equal(fmt,"obj")) then
             call obj_stick(lu,x0,x1,(/0,0,255/),0.03d0)
          elseif (equal(fmt,"ply")) then
             call ply_stick(lu,x0,x1,(/0,0,255/),0.03d0)
          elseif (equal(fmt,"off")) then
             call off_stick(lu,x0,x1,(/0,0,255/),0.03d0)
          end if
       end do
    end if

    ! close or give the handles to the calling routine, cleanup
    if (present(lumtl0) .and. present(lu0)) then
       lu0 = lu
       lumtl0 = lumtl
    else
       if (equal(fmt,"obj")) then
          call obj_close(lu,lumtl)
       elseif (equal(fmt,"ply")) then
          call ply_close(lu)
       elseif (equal(fmt,"off")) then
          call off_close(lu)
       end if
    end if

    if (allocated(fr0)) deallocate(fr0)

  end subroutine struct_write_3dmodel

  !> Write a quantum espresso input template
  subroutine struct_write_espresso(file,c,doprim)
    use struct_basic
    use tools_io
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    logical, intent(in) :: doprim

    character(len=:), allocatable :: lbl1
    integer :: i, lu
    integer :: ncelq, zq(c%ncel)
    real*8 :: x2cq(3,3), xq(3,c%ncel)
    logical :: ishex, ztyp(120)

    call cell_for_espresso(c,ncelq,x2cq,zq,xq,ishex,doprim)

    ztyp = .false.
    do i = 1, c%nneq
       if (.not.ztyp(c%at(i)%z)) ztyp(c%at(i)%z) = .true.
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
    write (lu,'(" nat=",I6,",")') ncelq
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
    do i = 1, ncelq
       write (lu,'(A2,3(X,F13.8,X))') lower(nameguess(zq(i),.true.)), xq(:,i)
    end do
    write (lu,'(/"K_POINTS automatic"/"2 2 2 1 1 1"/)')
    if (ishex) then
       write (lu,'("CELL_PARAMETERS hexagonal")')
    else
       write (lu,'("CELL_PARAMETERS cubic")')
    end if
    do i = 1, 3
       write (lu,'(3(F18.12,X))') x2cq(i,:)
    end do
    call fclose(lu)

  end subroutine struct_write_espresso

  !> Write a VASP POSCAR template
  subroutine struct_write_vasp(file,c,doprim,verbose)
    use struct_basic
    use tools_io
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    logical, intent(in) :: doprim, verbose

    character(len=:), allocatable :: lbl1
    integer :: i, j, lu
    integer :: ncelq, zq(c%ncel), ntyp(120)
    real*8 :: x2cq(3,3), xq(3,c%ncel)
    logical :: ishex

    ! Use the same cell as espresso.
    call cell_for_espresso(c,ncelq,x2cq,zq,xq,ishex,doprim)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, ncelq
       ntyp(zq(i)) = ntyp(zq(i)) + 1
    end do

    ! Cell
    lu = fopen_write(file)
    write (lu,'("crystal")')
    write (lu,'("1.0")')
    do i = 1, 3
       write (lu,'(3(F15.10,X))') x2cq(i,:) * bohrtoa
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
          do j = 1, ncelq
             if (zq(j) == i) then
                write (lu,'(3(F13.8,X))') xq(:,j)
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
  subroutine struct_write_abinit(file,c,doprim)
    use struct_basic
    use tools_io
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    logical, intent(in) :: doprim

    character(len=:), allocatable :: lbl1
    integer :: ncelq, zq(c%ncel), ntyp(120), iz
    real*8 :: x2cq(3,3), xq(3,c%ncel)
    logical :: ishex
    real*8 :: aap(3), bbp(3), gpq(3,3)
    integer :: i, j, lu

    ! Use the same cell as espresso.
    call cell_for_espresso(c,ncelq,x2cq,zq,xq,ishex,doprim)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, ncelq
       ntyp(zq(i)) = ntyp(zq(i)) + 1
    end do

    ! Find the lengths and angles of the primitive cell
    gpq = matmul(x2cq,transpose(x2cq))
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
    write (lu,'("natom ",I5)') ncelq

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
       do j = 1, ncelq
          if (zq(j) == i) then
             write (lu,'(X,3(F15.10,X))') xq(:,j)
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
  subroutine struct_write_elk(file,c,doprim)
    use struct_basic
    use tools_io
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    logical, intent(in) :: doprim

    integer :: ncelq, zq(c%ncel)
    real*8 :: x2cq(3,3), xq(3,c%ncel)
    logical :: ishex
    integer :: ntyp(100)
    integer :: i, j, lu

    ! Use the same cell as espresso.
    call cell_for_espresso(c,ncelq,x2cq,zq,xq,ishex,doprim)

    ! count number of atoms per type
    ntyp = 0
    do i = 1, ncelq
       ntyp(zq(i)) = ntyp(zq(i)) + 1
    end do

    ! Write input
    lu = fopen_write(file)
    write (lu,'("tasks"/,"0"/)')
    write (lu,'("xctype"/,"20"/)')
    write (lu,'("avec")')
    do i = 1, 3
       write (lu,'(2X,3(F15.10,X))') x2cq(i,:)
    end do
    write (lu,*)

    write (lu,'("sppath"/,"''./''"/)')

    write (lu,'("atoms")')
    write (lu,'(2X,I4)') count(ntyp > 0)
    do i = 1, size(ntyp)
       if (ntyp(i) > 0) then
          write (lu,'(2X,"''",A,".in''")') trim(nameguess(i,.true.))
          write (lu,'(2X,I3)') ntyp(i)
          do j = 1, ncelq
             if (zq(j) == i) then
                write (lu,'(2X,3(F14.10,X),"0.0 0.0 0.0")') xq(:,j)
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
  subroutine struct_write_gaussian(file,c,doprim)
    use struct_basic
    use tools_io
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
    logical, intent(in) :: doprim

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
    use struct_basic
    use global
    use tools_io
    use param

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
    use struct_basic
    use tools_io
    use param

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
    use struct_basic
    use tools_io
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: i, iz, lu

    lu = fopen_write(file)

    write (lu,'("data_default")')
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

  !> Write an escher octave script
  subroutine struct_write_escher(file,c)
    use struct_basic
    use global
    use tools_io
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    character(len=:), allocatable :: lbl1
    integer :: lu, i, j, n
    integer :: ntyp(120)

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
    use struct_basic
    use global
    use tools_io
    use tools_math
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c
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
    use struct_basic
    use global
    use tools_io
    use tools_math
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: i, j, k, l
    integer :: ntyp(100), lu
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
    do i = 1, 100
       if (ntyp(i) > 0) then
          j = j + 1
          write (lu,'(I3,X,F10.4)') j, atmass(i)
       end if
    end do
    write (lu,*)

    write (lu,'("Atoms"/)')
    k = 0
    l = 0
    do i = 1, 100
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
    use struct_basic
    use global
    use tools_io
    use tools_math
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: i, j, k
    integer :: ntyp(100), lu, nspecies

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
    use struct_basic
    use global
    use tools_io
    use tools_math
    use param

    character*(*), intent(in) :: file
    type(crystal), intent(in) :: c

    integer :: lu
    real*8 :: r(3,3)
    integer :: i, j, k, ntyp(100), nspecies

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

  !> Build a cell that is compliant with QE input rules, and perhaps
  !> reduce to a primitive cell.
  subroutine cell_for_espresso(c,ncelq,x2cq,zq,xq,ishex,lprim)
    use struct_basic
    use tools_math
    use tools_io
    use param

    type(crystal), intent(in) :: c
    integer, intent(out) :: ncelq
    real*8, intent(out) :: x2cq(3,3)
    integer, intent(out) :: zq(c%ncel)
    real*8, intent(out) :: xq(3,c%ncel)
    logical, intent(out) :: ishex
    logical, intent(in) :: lprim

    integer :: i, j
    integer :: csys
    real*8 :: r(3,3), rinv(3,3), gprim(3,3), rootp, d0(3)
    real*8 :: aap(3), ccp(3), ss3, tx, ty, tz, x0(3), dist2
    real*8 :: ss(3), cc(3)
    integer :: doprim, n0
    logical :: dotrig

    ! it defaults to the critic2 cell
    ncelq = c%ncel
    x2cq = transpose(c%crys2car)
    do i = 1, c%ncel
       zq(i) = c%at(c%atcel(i)%idx)%z
       xq(:,i) = c%atcel(i)%x
    end do
    ishex = .false.
    dotrig = .false.

    ! the user says no primitive
    if (.not.lprim) return

    ! consider each crystal system and centering
    csys = crystal_system(c)

    doprim = 0
    select case (csys)
    case (csys_cub)
       if (c%lcent == 1) then
          ! primitive
          x2cq = transpose(c%crys2car)
       else if (c%lcent == 6) then
          ! F-centered
          x2cq(1,:) = (/-c%aa(1)/2d0,         0d0, c%aa(3)/2d0 /)
          x2cq(2,:) = (/         0d0, c%aa(2)/2d0, c%aa(3)/2d0 /)
          x2cq(3,:) = (/-c%aa(1)/2d0, c%aa(2)/2d0,         0d0 /)
          doprim = 1
       else if (c%lcent == 5) then
          ! I-centered
          x2cq(1,:) = (/ c%aa(1)/2d0,  c%aa(2)/2d0,  c%aa(3)/2d0 /)
          x2cq(2,:) = (/-c%aa(1)/2d0,  c%aa(2)/2d0,  c%aa(3)/2d0 /)
          x2cq(3,:) = (/-c%aa(1)/2d0, -c%aa(2)/2d0,  c%aa(3)/2d0 /)
          doprim = 1
       else
          doprim = 2
       end if
    case (csys_hex,csys_trig)
       ishex = .true.
       if (c%lcent == 1) then
          x2cq = transpose(c%crys2car)
       elseif (c%ncv == 3) then
          doprim = 2
          dotrig = .true.
       end if
    case (csys_tetr)
       if (c%lcent == 1) then
          ! primitive
          x2cq(1,:) = (/ c%aa(1),     0d0,     0d0 /)
          x2cq(2,:) = (/     0d0, c%aa(2),     0d0 /)
          x2cq(3,:) = (/     0d0,     0d0, c%aa(3) /)
       elseif (c%lcent == 5) then
          ! I-centered
          x2cq(1,:) = (/ c%aa(1)/2d0, -c%aa(2)/2d0, c%aa(3)/2d0 /)
          x2cq(2,:) = (/ c%aa(1)/2d0,  c%aa(2)/2d0, c%aa(3)/2d0 /)
          x2cq(3,:) = (/-c%aa(1)/2d0, -c%aa(2)/2d0, c%aa(3)/2d0 /)
          doprim = 1
       else
          doprim = 2
       end if
    case (csys_orth)
       if (c%lcent == 1) then
          ! primitive
          x2cq(1,:) = (/ c%aa(1),     0d0,     0d0 /)
          x2cq(2,:) = (/     0d0, c%aa(2),     0d0 /)
          x2cq(3,:) = (/     0d0,     0d0, c%aa(3) /)
       elseif (c%lcent == 4) then
          ! C-centered
          x2cq(1,:) = (/ c%aa(1)/2d0, c%aa(2)/2d0,     0d0 /)
          x2cq(2,:) = (/-c%aa(1)/2d0, c%aa(2)/2d0,     0d0 /)
          x2cq(3,:) = (/         0d0,         0d0, c%aa(3) /)
          doprim = 1
       elseif (c%lcent == 6) then
          ! F-centered
          x2cq(1,:) = (/ c%aa(1)/2d0,         0d0, c%aa(3)/2d0 /)
          x2cq(2,:) = (/ c%aa(1)/2d0, c%aa(2)/2d0,         0d0 /)
          x2cq(3,:) = (/         0d0, c%aa(2)/2d0, c%aa(3)/2d0 /)
          doprim = 1
       elseif (c%lcent == 5) then
          ! I-centered
          x2cq(1,:) = (/ c%aa(1)/2d0,  c%aa(2)/2d0,  c%aa(3)/2d0 /)
          x2cq(2,:) = (/-c%aa(1)/2d0,  c%aa(2)/2d0,  c%aa(3)/2d0 /)
          x2cq(3,:) = (/-c%aa(1)/2d0, -c%aa(2)/2d0,  c%aa(3)/2d0 /)
          doprim = 1
       else
          doprim = 2
       end if
    case (csys_mono)
       if (c%lcent == 1) then
          if (abs(cc(3)) > 1d-5) then
             ! primitive unique axis c
             x2cq(1,:) = (/ c%aa(1),   0d0, 0d0 /)
             x2cq(2,:) = (/ c%aa(2)*cc(3), c%aa(2)*ss(3), 0d0/)
             x2cq(3,:) = (/ 0d0, 0d0, c%aa(3) /)
          else if (abs(cc(2)) > 1d-5) then
             ! primitive unique axis b
             x2cq(1,:) = (/ c%aa(1),   0d0, 0d0 /)
             x2cq(2,:) = (/   0d0, c%aa(2), 0d0/)
             x2cq(3,:) = (/ c%aa(3)*cc(2), 0d0, c%aa(3)*ss(2) /)
          else
             ! primitive unique axis a
             doprim = 2
          end if
       else if (c%lcent == 3) then
          if (abs(cc(3)) > 1d-5) then
             ! B centering, unique axis c
             x2cq(1,:) = (/ c%aa(1)/2d0, 0d0, -c%aa(3)/2d0 /)
             x2cq(2,:) = (/ c%aa(2)*cc(3), c%aa(2)*ss(3), 0d0/)
             x2cq(3,:) = (/ c%aa(1)/2d0, 0d0, c%aa(3)/2d0 /)
             doprim = 1
          else
             doprim = 2
          end if
       else
          ! some other centering and setting
          doprim = 2
       end if
    case (csys_tric)
       ! qe likes cholesky for triclinic
       doprim = 2
    case default
       call ferror("cell_for_espresso","unknown crystal system",faterr)
    end select

    if (doprim == 1) then
       ! use qe-suggested transformation to the primitive
       r = x2cq
       do i = 1,3
          r(i,:) = c%c2x(r(i,:))
       end do
       r = transpose(r)
       rinv = matinv(r)
    else if (doprim == 2) then
       ! then the setting was not considered in the manual, use my transformation
       r = toc_real(:,:,c%lcent)
       rinv = toc_rec(:,:,c%lcent)
       ! special case -> trigonal
       if (c%lcent == 7 .and. c%ncv == 3) then
          if (abs(c%cen(3,2)-2d0/3d0)<1d-10.and.abs(c%cen(3,3)-2d0/3d0)<1d-10) then
             ! r-
             r = toc_real(:,:,9)
             rinv = toc_rec(:,:,9)
          else
             ! r+
             r = toc_real(:,:,8)
             rinv = toc_rec(:,:,8)
          endif
       endif
    end if

    if (doprim > 0) then
       ! transform the cell
       gprim = matmul(transpose(r),matmul(c%gtensor,r))
       do i = 1, 3
          aap(i) = sqrt(gprim(i,i))
       end do
       ccp(1) = gprim(2,3)/aap(2)/aap(3)
       ccp(2) = gprim(1,3)/aap(1)/aap(3)
       ccp(3) = gprim(1,2)/aap(1)/aap(2)
       ss3 = sqrt(1d0-ccp(3)**2)

       ! transform the atoms
       n0 = 0
       main: do i = 1, ncelq
          x0 = matmul(rinv,xq(:,i))
          x0 = x0 - floor(x0)
          do j = 1, n0
             d0 = x0 - xq(:,j)
             call c%shortest(d0,dist2)
             if (dist2 < 1d-10) cycle main
          end do
          n0 = n0 + 1
          xq(:,n0) = x0
          zq(n0) = zq(i)
       end do main
       ncelq = n0

       if (dotrig) then
          tx = sqrt((1-ccp(3))/2)
          ty = sqrt((1-ccp(3))/6)
          tz = sqrt((1+2*ccp(3))/3)
          x2cq(1,:) = (/tx,-ty,tz/) * aap(1)
          x2cq(2,:) = (/0d0,2d0*ty,tz/) * aap(1)
          x2cq(3,:) = (/-tx,-ty,tz/) * aap(1)
       else if (doprim > 1) then
          ! poor man's version of the crystal system determination for qe
          ! we have a primitive now
          if (all(abs(ccp) < 1d-5)) then
             ! cubic, tetragonal, orthorhombic
             x2cq = 0d0
             do i = 1, 3
                x2cq(i,i) = aap(i)
             end do
          elseif (abs(ccp(1))<1d-5 .and. abs(ccp(2))<1d-5 .and. abs(ccp(3)+0.5d0)<1d-5) then
             ! hexagonal/trigonal
             x2cq(1,:) = (/1d0, 0d0, 0d0/) * aap(1)
             x2cq(2,:) = (/-0.5d0, sqrt(3d0)/2d0, 0d0/) * aap(1)
             x2cq(3,:) = (/0d0, 0d0, aap(3)/aap(1)/) * aap(1)
          elseif (abs(ccp(3))>1d-5 .and. abs(ccp(1))<1d-5 .and. abs(ccp(2))<1d-5) then
             ! monoclinic c
             x2cq(1,:) = (/ aap(1),   0d0, 0d0 /)
             x2cq(2,:) = (/ aap(2)*ccp(3), aap(2)*sqrt(1-ccp(3)**2), 0d0/)
             x2cq(3,:) = (/ 0d0, 0d0, c%aa(3) /)
          elseif (abs(ccp(2))>1d-5 .and. abs(ccp(1))<1d-5 .and. abs(ccp(3))<1d-5) then
             ! monoclinic b
             x2cq(1,:) = (/ aap(1),   0d0, 0d0 /)
             x2cq(2,:) = (/   0d0, aap(2), 0d0/)
             x2cq(3,:) = (/ aap(3)*ccp(2), 0d0, aap(3)*sqrt(1-ccp(2)**2) /)
          else
             ! triclinic
             rootp=sqrt(1d0-ccp(1)*ccp(1)-ccp(2)*ccp(2)-ccp(3)*ccp(3)+2d0*ccp(1)*ccp(2)*ccp(3))
             x2cq = 0d0
             x2cq(1,1) = aap(1)
             x2cq(1,2) = aap(2)*ccp(3)
             x2cq(1,3) = aap(3)*ccp(2)
             x2cq(2,2) = aap(2)*ss3
             x2cq(2,3) = aap(3)*(ccp(1)-ccp(3)*ccp(2))/ss3
             x2cq(3,3) = aap(3)*rootp/ss3
             x2cq = transpose(x2cq)
          end if
       end if
    end if

  end subroutine cell_for_espresso

  !> Find the crystal system from the Laue group
  function crystal_system(c) result(csys)
    use struct_basic
    use tools_math
    type(crystal), intent(in) :: c
    integer :: csys

    csys = 0
    if (c%lauec == 1) then
       csys = csys_tric
    elseif (c%lauec == 2) then
       csys = csys_mono
    elseif (c%lauec == 3) then
       csys = csys_orth
    elseif (c%lauec == 4 .or. c%lauec == 5) then
       csys = csys_tetr
    elseif (c%lauec == 6 .or. c%lauec == 7) then
       csys = csys_trig
    elseif (c%lauec == 8 .or. c%lauec == 9) then
       csys = csys_hex
    elseif (c%lauec == 10 .or. c%lauec == 11) then
       csys = csys_cub
    endif

  end function crystal_system

end module struct_writers
