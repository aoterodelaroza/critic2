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

! Routines for reading and handling molecular and crystal vibrations.
submodule (crystalmod) vibrationsmod
  implicit none

  ! abramowitz, stegun, numerical differentiation formulas, table 25.2 (1st and 2nd derivatives)
  real*8, parameter :: der_coef1(3) = (/-1d0,  0d0, 1d0/)
  real*8, parameter :: der_coef2(3) = (/ 1d0, -2d0, 1d0/)

  ! frequency conversion factor: sqrt(Hartree/bohr^2/amu) to cm-1
  ! sqrt(Hartree / bohr^2 / amu) / 1e12 -> THz
  ! sqrt(4.3597447222071e-18) / 5.29177210903e-11 / sqrt(1.66053906660e-27) /  2 / pi / c / 100
  real*8, parameter :: freqfactor = 5140.487143715827d0

  ! other conversion factors
  real*8, parameter :: cminv_to_kJmol = 1.196265656955833d-02 ! c * h * NA / 10
  real*8, parameter :: cminv_to_K = 1.438776959983815d0 ! c * h * 100 / kb
  real*8, parameter :: cminv_to_hartree = 4.5563352529120d-8 * 100d0 ! CODATA2018
  real*8, parameter :: amu_to_me = 1.66053906660e-27 / 9.1093837015e-31 ! CODATA2018
  real*8, parameter :: cminv_to_angfreq_au = 4.556335252903557d-06 ! 100 * c * a0 * sqrt(me / Ha) * 2 * pi, using CODATA2018 values

  ! Force-constant units of the codes that can generate a phonopy
  ! FORCE_CONSTANTS file, and the factor that converts them to the
  ! internal units (Hartree/bohr^2).
  integer, parameter :: fc2_nunit = 6
  character*20, parameter :: fc2_unitname(fc2_nunit) = (/&
     "eV/ang^2            ", "eV/(ang*bohr)       ", "Ry/bohr^2           ",&
     "mRy/bohr^2          ", "Hartree/bohr^2      ", "Hartree/(ang*bohr)  "/)
  integer, parameter :: fc2_ngen = 22
  character*10, parameter :: fc2_genname(fc2_ngen) = (/&
     "vasp      ", "aims      ", "fhiaims   ", "lammps    ", "pwmat     ",&
     "crystal   ", "castep    ", "alamode   ",&
     "abinit    ", "siesta    ", "abacus    ",&
     "qe        ", "pwscf     ", "espresso  ",&
     "wien2k    ", "wien      ",&
     "elk       ", "dftb+     ", "dftbp     ", "turbomole ", "fleur     ",&
     "cp2k      "/)
  integer, parameter :: fc2_genunit(fc2_ngen) = (/&
     1, 1, 1, 1, 1,&
     1, 1, 1,&
     2, 2, 2,&
     3, 3, 3,&
     4, 4,&
     5, 5, 5, 5, 5,&
     6/)

  ! tolerances used when reading force constants
  real*8, parameter :: fc2_epsint = 1d-3 ! integrality of the supercell matrix
  real*8, parameter :: fc2_epsmetric = 1d-4 ! relative consistency of the supercell metric
  real*8, parameter :: fc2_epspos = 1d-2 ! atom identification distance (bohr)
  real*8, parameter :: fc2_epssym = 1d-2 ! relative tolerance for the transpose-symmetry check

  !xx! private procedures
  ! subroutine vibrations_detect_format(file,ivformat)
  ! subroutine read_matdyn_modes(v,c,file,ivformat,errmsg,ti)
  ! subroutine read_qe_dyn(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_ascii(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_yaml(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_hdf5(v,c,file,errmsg)
  ! subroutine read_phonopy_fc2(v,c,file,sline,errmsg,ti)
  ! function fc2_generator_index(sline)
  ! subroutine fc2_lattice_points(smat,nlat,madj,lvec,lkey,errmsg)
  ! function fc2_ilat(nlat,lkey,madj,lv)
  ! function fc2_pure_translation(c,t,nlat,lvec,lkey,madj,perm,dev)
  ! subroutine read_crystal_out(v,c,file,errmsg,ti)
  ! subroutine read_gaussian_log(v,c,file,errmsg,ti)
  ! subroutine read_gaussian_fchk(v,c,file,errmsg,ti)
  ! subroutine read_castep_phonon(v,c,file,errmsg,ti)

contains

  !> Terminate a vibrations object
  module subroutine vibrations_end(v,keepfc2,keepvibs)
    use param, only: ivformat_unknown
    class(vibrations), intent(inout) :: v
    logical, intent(in), optional :: keepfc2
    logical, intent(in), optional :: keepvibs

    logical :: keepfc2_, keepvibs_

    keepfc2_ = .false.
    if (present(keepfc2)) keepfc2_ = keepfc2
    keepvibs_ = .false.
    if (present(keepvibs)) keepvibs_ = keepvibs

    if (.not.keepvibs_) then
       v%hasvibs = .false.
       v%file = ""
       v%ivformat = ivformat_unknown
       v%nqpt = 0
       if (allocated(v%qpt)) deallocate(v%qpt)
       if (allocated(v%freq)) deallocate(v%freq)
       if (allocated(v%vec)) deallocate(v%vec)
       v%nfreq = 0
    end if
    if (.not.keepfc2_) then
       v%hasfc2 = .false.
       v%fc2_file = ""
       v%fc2_smat = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
       v%fc2_madj = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
       v%fc2_nlat = 0
       v%fc2_nsat = 0
       v%fc2_ncel = 0
       v%fc2_iscompact = .false.
       v%fc2_acoustic = -1
       v%fc2_vs_center = -1d0
       v%fc2_vs_delta = -1d0
       if (allocated(v%fc2)) deallocate(v%fc2)
       if (allocated(v%fc2_lvec)) deallocate(v%fc2_lvec)
       if (allocated(v%fc2_lkey)) deallocate(v%fc2_lkey)
    end if

  end subroutine vibrations_end

  !> Read a file (file) containing the vibrational information for
  !> this system and add the q-points to c%vib. The vibration data
  !> format is detected from the file extension if ivformat is
  !> isformat_unknown or format ivformat is used otherwise. Returns
  !> non-zero errmsg if error.
  module subroutine vibrations_read_file(v,c,file,sline,ivformat,errmsg,ti)
    use param, only: ivformat_unknown, ivformat_matdynmodes, ivformat_matdyneig,&
       ivformat_qedyn, ivformat_phonopy_ascii, ivformat_phonopy_yaml, ivformat_phonopy_hdf5,&
       ivformat_crystal_out, ivformat_gaussian_log, ivformat_gaussian_fchk,&
       ivformat_phonopy_fc2, ivformat_castep_phonon
    use types, only: realloc
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file, sline
    integer, intent(in) :: ivformat
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: ivf

    ! detect the format. A generator keyword at the start of the option
    ! string means this is a phonopy FORCE_CONSTANTS file, whatever the
    ! file happens to be called.
    if (ivformat == ivformat_unknown) then
       if (fc2_generator_index(sline) > 0) then
          ivf = ivformat_phonopy_fc2
       else
          call vibrations_detect_format(file,ivf)
          if (ivf == ivformat_unknown) then
             errmsg = "Unknown vibration file format: " // trim(file)
             return
          end if
       end if
    else
       ivf = ivformat
    end if

    ! Clear only the half of the object this format provides, so that
    ! force constants and frequencies/modes from different files can
    ! coexist (see the hasfc2/hasvibs comment in the type definition).
    if (ivf == ivformat_phonopy_fc2) then
       call v%end(keepvibs=.true.)
    else
       call v%end(keepfc2=.true.)
    end if

    ! read the vibrations
    if (ivf == ivformat_matdynmodes .or. ivf == ivformat_matdyneig) then
       call read_matdyn_modes(v,c,file,ivf,errmsg,ti)
    elseif (ivf == ivformat_qedyn) then
       call read_qe_dyn(v,c,file,errmsg,ti)
    elseif (ivf == ivformat_phonopy_ascii) then
       call read_phonopy_ascii(v,c,file,errmsg,ti)
    elseif (ivf == ivformat_phonopy_yaml) then
       call read_phonopy_yaml(v,c,file,errmsg,ti)
    elseif (ivf == ivformat_phonopy_hdf5) then
       call read_phonopy_hdf5(v,c,file,errmsg)
    elseif (ivf == ivformat_phonopy_fc2) then
       call read_phonopy_fc2(v,c,file,sline,errmsg,ti)
    elseif (ivf == ivformat_crystal_out) then
       call read_crystal_out(v,c,file,errmsg,ti)
    elseif (ivf == ivformat_gaussian_log) then
       call read_gaussian_log(v,c,file,errmsg,ti)
    elseif (ivf == ivformat_gaussian_fchk) then
       call read_gaussian_fchk(v,c,file,errmsg,ti)
    elseif (ivf == ivformat_castep_phonon) then
       call read_castep_phonon(v,c,file,errmsg,ti)
    else
       errmsg = "Unknown vibration file format: " // trim(file)
       return
    end if
    if (len(errmsg) > 0) return
    if (v%hasvibs) then
       if (v%nfreq == 0) then
          errmsg = "No frequencies found in file: " // trim(file)
          return
       end if
       if (v%nqpt == 0) then
          errmsg = "No q-points found in file: " // trim(file)
          return
       end if
    end if

  end subroutine vibrations_read_file

  !> Print a summary of the vibrational info.
  module subroutine vibrations_print_summary(v)
    use tools_io, only: uout, string, ioj_right
    use param, only: ivformat_matdynmodes, ivformat_matdyneig, ivformat_qedyn,&
       ivformat_phonopy_ascii, ivformat_phonopy_yaml, ivformat_phonopy_hdf5,&
       ivformat_crystal_out, ivformat_gaussian_log, ivformat_castep_phonon,&
       ivformat_gaussian_fchk
    class(vibrations), intent(inout) :: v

    integer :: i, j

    write (uout,'("+ PRINT summary")')

    ! dynamical matrices
    if (v%hasvibs) then
       write (uout,*)
       write (uout,'("# Dynamical matrices available")')
       write (uout,'("  File: ",A)') trim(v%file)
       if (v%ivformat == ivformat_matdynmodes) then
          write (uout,'("  Format: QE matdyn.modes")')
       elseif (v%ivformat == ivformat_matdyneig) then
          write (uout,'("  Format: QE matdyn.eig")')
       elseif (v%ivformat == ivformat_qedyn) then
          write (uout,'("  Format: QE ph.x dyn files")')
       elseif (v%ivformat == ivformat_phonopy_ascii) then
          write (uout,'("  Format: phonopy ascii file")')
       elseif (v%ivformat == ivformat_phonopy_yaml) then
          write (uout,'("  Format: phonopy yaml file")')
       elseif (v%ivformat == ivformat_phonopy_hdf5) then
          write (uout,'("  Format: phonopy hdf5 file")')
       elseif (v%ivformat == ivformat_crystal_out) then
          write (uout,'("  Format: crystal output file")')
       elseif (v%ivformat == ivformat_gaussian_log) then
          write (uout,'("  Format: Gaussian output file")')
       elseif (v%ivformat == ivformat_gaussian_fchk) then
          write (uout,'("  Format: Gaussian formatted checkpoint file")')
       elseif (v%ivformat == ivformat_castep_phonon) then
          write (uout,'("  Format: CASTEP phonon file")')
       end if

       write (uout,*)
       if (v%nqpt == 0) then
          write (uout,'("# No q-points known")')
       else
          write (uout,'("# List of q-points: ",A)') string(v%nqpt)
          write (uout,'("# id           x              y              z")')
          do i = 1, v%nqpt
             write (uout,'(4(A," "))') string(i,5), (string(v%qpt(j,i),'f',14,8,ioj_right),j=1,3)
          end do
       end if
       write (uout,'("# Number of frequencies at each q-point: ",A)') string(v%nfreq)
    else
       write (uout,*)
       write (uout,'("# Dynamical matrices NOT available")')
    end if

    if (v%hasfc2) then
       write (uout,*)
       write (uout,'("# Second-order force constants available")')
       write (uout,'("  File: ",A)') trim(v%fc2_file)
       if (v%fc2_iscompact) then
          write (uout,'("  Format: phonopy FORCE_CONSTANTS (compact)")')
       else
          write (uout,'("  Format: phonopy FORCE_CONSTANTS (full)")')
       end if
       write (uout,'("  Supercell matrix (rows = supercell lattice vectors):")')
       do i = 1, 3
          write (uout,'("    ",3(A," "))') (string(v%fc2_smat(i,j),4,ioj_right),j=1,3)
       end do
       write (uout,'("  Number of lattice points in the supercell: ",A)') string(v%fc2_nlat)
       write (uout,'("  Number of atoms in the cell: ",A)') string(v%fc2_ncel)
       write (uout,'("  Number of atoms in the supercell: ",A)') string(v%fc2_nsat)
       write (uout,'("  Size of the force constant matrix (Mb): ",A)') &
          string(72d0*v%fc2_ncel*v%fc2_nsat/1024d0**2,'f',10,2)
       write (uout,'("  Units: Hartree/bohr^2")')
       write (uout,'("  Note: the Cartesian components of the force constants are assumed to be")')
       write (uout,'("        in the same frame as the current structure.")')
    else
       write (uout,*)
       write (uout,'("# Second-order force constants NOT available")')
    end if

  end subroutine vibrations_print_summary

  !> Print information about the stored FC2 to the standard output.
  !> The structural info comes from crystal structures c.
  module subroutine vibrations_print_fc2(v,c,disteps,fc2eps,environ)
    use tools_io, only: ferror, faterr
    use tools_io, only: uout, ioj_center, ioj_right, string
    use global, only: iunit, iunitname0, dunit0
    use tools, only: mergesort
    use param, only: icrd_crys
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in), optional :: disteps, fc2eps
    logical, intent(in), optional :: environ

    !! NOT REIMPLEMENTED YET. This routine assumed the old fc2(3,3,ncel,ncel)
    !! layout, in which both indices ran over the atoms of the cell. The body
    !! is kept commented out below until it is rewritten for the cell-compact
    !! fc2(3,3,ncel,nsat) array (second index over the supercell).
    call ferror('vibrations_print_fc2','not yet reimplemented for the cell-compact force constants',faterr)

!!    integer :: i, j, k, idum, jdum, id
!!    integer :: npair, cidx, nidx
!!    real*8 :: maxdisteps, dist, disteps_, fc2eps_, up2d, xx(3), norm
!!    logical :: isintra, environ_
!!    real*8, allocatable :: dista(:)
!!    integer, allocatable :: idx(:), i2(:,:), pair(:,:)
!!    logical, allocatable :: isdone(:), isdonej(:)
!!    integer :: nat
!!    integer, allocatable :: eid(:)
!!    real*8, allocatable :: distr(:)
!!    integer, allocatable :: lvec(:,:)
!!
!!    ! header and checks
!!    write (uout,'("+ PRINT 2nd-order force constant information (FC2)")')
!!    if (.not.v%hasfc2.or..not.allocated(v%fc2)) then
!!       write (uout,'("!! No FC2 info has been loaded, exiting.")')
!!       return
!!    end if
!!
!!    ! process options
!!    disteps_ = huge(1d0)
!!    fc2eps_ = 0d0
!!    environ_ = .false.
!!    if (present(disteps)) then
!!       disteps_ = disteps
!!       if (disteps_ < huge(1d0)) &
!!          write (uout,'("# Cutoff distance for pairs = ",A," ",A)') &
!!             string(disteps_,'f',decimal=5), iunitname0(iunit)
!!    end if
!!    if (present(fc2eps)) then
!!       fc2eps_ = fc2eps
!!       if (fc2eps_ > 0d0) &
!!          write (uout,'("# Cutoff for FC2 values = ",A)') string(fc2eps_,'e',decimal=5)
!!    end if
!!    if (present(environ)) environ_ = environ
!!
!!    ! allocate information about atom pairs, and sort by distance
!!    npair = c%ncel*(c%ncel + 1) / 2
!!    allocate(dista(npair),idx(npair),i2(2,npair))
!!    allocate(pair(2,npair))
!!    k = 0
!!    do i = 1, c%ncel
!!       do j = i, c%ncel
!!          k = k + 1
!!          dista(k) = c%eql_distance(c%atcel(i)%x, c%atcel(j)%x)
!!          idx(k) = k
!!          i2(1,k) = i
!!          i2(2,k) = j
!!       end do
!!    end do
!!    maxdisteps = maxval(dista)
!!    call mergesort(dista,idx,1,npair)
!!
!!    if (.not.environ) then
!!       ! print the whole fc2 matrix
!!       write (uout,'("# Id1,Id2 = complete cell IDs. At1,At2 = atomic symbols. isintra = is intramolecular?")')
!!       write (uout,'("# norm = Frobenius norm FC2. xx,... = FC2 components")')
!!       write (uout,'("# All FC2 in Hartree/bohr^2.")')
!!       write (uout,'("# Id1 At1  Id2 At2  dist(bohr) isintra FC2: norm      xx              xy              xz   &
!!          &           yx              yy              yz              zx              zy              zz   ")')
!!       k = 0
!!       do idum = 1, c%ncel
!!          do jdum = idum, c%ncel
!!             k = k + 1
!!             i = i2(1,idx(k))
!!             j = i2(2,idx(k))
!!             dist = dista(idx(k))
!!             isintra = (c%idatcelmol(1, i) == c%idatcelmol(1, j))
!!             norm = sqrt(sum(v%fc2(:,:,i,j)**2))
!!             if (dist <= disteps_ .and. norm >= fc2eps_) then
!!                write (uout,'(99(A," "))') string(i, 5, ioj_center), string(c%spc(c%atcel(i)%is)%name,2),&
!!                   string(j, 5, ioj_center), string(c%spc(c%atcel(j)%is)%name,2),&
!!                   string(dist*dunit0(iunit),'f',12,6,ioj_right), string(isintra), string(norm,'e',15,8,ioj_right),&
!!                   string(v%fc2(1,1,i,j),'e',15,8,ioj_right), string(v%fc2(1,2,i,j),'e',15,8,ioj_right),&
!!                   string(v%fc2(1,3,i,j),'e',15,8,ioj_right), string(v%fc2(2,1,i,j),'e',15,8,ioj_right),&
!!                   string(v%fc2(2,2,i,j),'e',15,8,ioj_right), string(v%fc2(2,3,i,j),'e',15,8,ioj_right),&
!!                   string(v%fc2(3,1,i,j),'e',15,8,ioj_right), string(v%fc2(3,2,i,j),'e',15,8,ioj_right),&
!!                   string(v%fc2(3,3,i,j),'e',15,8,ioj_right)
!!             end if
!!          end do
!!       end do
!!    else
!!       ! print the FC2 in environments of nneq atoms
!!       allocate(isdone(c%nneq),isdonej(c%ncel))
!!       isdone = .false.
!!       up2d = min(disteps_,maxdisteps+1d-4)
!!
!!       ! run over non-equivalent atoms, i is the complete list index
!!       write (uout,*)
!!       do i = 1, c%ncel
!!          id = c%atcel(i)%idx
!!          if (isdone(id)) cycle
!!          isdone(id) = .true.
!!
!!          ! find the atomic environment
!!          call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.true.,nat,eid,distr,lvec,up2d=up2d)
!!
!!          write (uout,'("+ Environment of atom ",A," (spc=",A,", nid=",A,") at ",3(A," "))') &
!!             string(i), string(c%spc(c%at(id)%is)%name), string(id), &
!!             (string(c%atcel(i)%x(j),'f',length=10,decimal=6),j=1,3)
!!          write (uout,'("# Up to distance (",A,"): ",A)') iunitname0(iunit), string(up2d*dunit0(iunit),'f',length=10,decimal=6)
!!          write (uout,'("# Number of atoms in the environment: ",A)') string(nat)
!!          write (uout,'("# nid = non-equivalent list atomic ID. id = complete list ID plus lattice vector (lvec).")')
!!          write (uout,'("# name = atomic name. dist = distance. isintra = is intramolecular?")')
!!          write (uout,'("# norm = Frobenius norm FC2. xx,... = FC2 components")')
!!          write (uout,'("# All FC2 in Hartree/bohr^2.")')
!!          write (uout,'("#nid   id      lvec     name  dist(",A,")  isintra FC2: norm        xx              &
!!             &xy              xz              yx              yy              yz              &
!!             &zx              zy              zz")') iunitname0(iunit)
!!          isdonej = .false.
!!          do j = 1, nat
!!             cidx = eid(j)
!!             if (isdonej(cidx)) cycle
!!             isdonej(cidx) = .true.
!!
!!             nidx = c%atcel(cidx)%idx
!!             if (c%ismolecule) then
!!                xx = (c%atcel(cidx)%r + c%molx0) * dunit0(iunit)
!!             else
!!                xx = c%atcel(cidx)%x + lvec(:,j)
!!             end if
!!
!!             isintra = (c%idatcelmol(1,i) == c%idatcelmol(1,cidx))
!!             norm = sqrt(sum(v%fc2(:,:,i,cidx)**2))
!!             if (distr(j) <= disteps_ .and. norm >= fc2eps_) then
!!                write (uout,'("  ",2(A," "),"(",A," ",A," ",A,")",99(" ",A))') string(nidx,4,ioj_center), string(cidx,4,ioj_center),&
!!                   (string(lvec(k,j),2,ioj_right),k=1,3), string(c%spc(c%atcel(cidx)%is)%name,7,ioj_center),&
!!                   string(distr(j)*dunit0(iunit),'f',12,6,4), string(isintra), string(norm,'e',15,8,ioj_right),&
!!                   string(v%fc2(1,1,i,cidx),'e',15,8,ioj_right), string(v%fc2(1,2,i,cidx),'e',15,8,ioj_right),&
!!                   string(v%fc2(1,3,i,cidx),'e',15,8,ioj_right), string(v%fc2(2,1,i,cidx),'e',15,8,ioj_right),&
!!                   string(v%fc2(2,2,i,cidx),'e',15,8,ioj_right), string(v%fc2(2,3,i,cidx),'e',15,8,ioj_right),&
!!                   string(v%fc2(3,1,i,cidx),'e',15,8,ioj_right), string(v%fc2(3,2,i,cidx),'e',15,8,ioj_right),&
!!                   string(v%fc2(3,3,i,cidx),'e',15,8,ioj_right)
!!             end if
!!          end do
!!          write (uout,*)
!!       end do
!!    end if
!!
  end subroutine vibrations_print_fc2

  !> Print frequency information about a particular q-point with given
  !> id.
  module subroutine vibrations_print_freq(v,id)
    use tools_io, only: uout, ferror, faterr, string, ioj_right
    class(vibrations), intent(inout) :: v
    integer, intent(in) :: id

    integer :: i, j

    ! header
    write (uout,'("+ PRINT frequency info for a given q-point")')

    ! check the q-point id
    if (id <= 0 .or. id > v%nqpt) &
       call ferror('vibrations_print_freq','q-point index out of range',faterr)
    write (uout,'("# q-point (",A,", cryst. coords.) = ",3(A," "))') string(id), &
       (string(v%qpt(j,id),'f',decimal=8,justify=ioj_right),j=1,3)

    ! write the info
    write (uout,'("# List of frequencies: ",A)') string(v%nfreq)
    write (uout,'("#Id      Frequency(cm^-1)")')
    do i = 1, v%nfreq
       write (uout,'(4(A," "))') string(i,5), string(v%freq(i,id),'f',decimal=8,justify=ioj_right)
    end do

  end subroutine vibrations_print_freq

  !> Print information about a particular eigenvector with band index
  !> ifreq and q-point index idq. If cartesian, convert to Cartesian
  !> coordinates by dividing by sqrt(mass); otherwise, mass-weighted
  !> coordinates (orthonormal eigenvecs).
  module subroutine vibrations_print_eigenvector(v,c,ifreq,idq,cartesian)
    use tools_io, only: string, ferror, faterr, ioj_right, ioj_center, uout
    use param, only: atmass
    class(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    integer, intent(in) :: ifreq, idq
    logical, intent(in) :: cartesian

    integer :: i, j, k, nidx
    character*3, parameter :: dir(3) = (/" x "," y "," z "/)
    complex*16 :: evec

    ! header
    write (uout,'("+ PRINT eigenvector for a given q-point and band index")')

    ! check the q-point idq
    if (idq <= 0 .or. idq > v%nqpt) &
       call ferror('vibrations_print_eigenvector','q-point index out of range',faterr)
    if (ifreq <= 0 .or. ifreq > v%nfreq) &
       call ferror('vibrations_print_eigenvector','frequency index out of range',faterr)
    write (uout,'("# q-point (",A,", cryst. coords.) = ",3(A," "))') string(idq),&
       (string(v%qpt(j,idq),'f',decimal=8,justify=ioj_right),j=1,3)
    write (uout,'("# frequency (",A,") = ",A," cm^-1")') string(ifreq),&
       string(v%freq(ifreq,idq),'f',decimal=8,justify=ioj_right)
    write (uout,'("# nid = non-equivalent list atomic ID. id = complete list ID plus lattice vector (lvec).")')
    write (uout,'("# name = atomic name. dir = direction. evec = eigenvector, real and imaginary parts.")')
    write (uout,'("#nid   id     name  ---      cryst. coords.     ---  dir   evec(real)         evec(imag) ")')
    do i = 1, c%ncel
       do k = 1, 3
          nidx = c%atcel(i)%idx
          evec = v%vec(k,i,ifreq,idq)
          if (cartesian) evec = evec / sqrt(atmass(c%spc(c%atcel(i)%is)%z))
          write (uout,'("  ",99(A," "))') string(nidx,4,ioj_center), string(i,4,ioj_center),&
             string(c%spc(c%atcel(i)%is)%name,7,ioj_center),&
             (string(c%atcel(i)%x(j),'f',length=10,decimal=6),j=1,3),&
             dir(k), string(real(evec,8),'e',18,10,ioj_right),&
             string(aimag(evec),'e',18,10,ioj_right)
       end do
    end do

  end subroutine vibrations_print_eigenvector

  !> Modify the current FC2 by imposing acoustic sum rules and symmetrize the FC2:
  !>   sum_i phi_(ia,jb) = 0
  !>   sum_j phi_(ia,jb) = 0
  !>   phi_(ia,jb) = phi(jb,ia)
  module subroutine vibrations_apply_acoustic(v,c,verbose)
    use tools_io, only: uout, string
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    logical, intent(in), optional :: verbose

    integer :: i, j, ia, is
    real*8 :: summ(3,3), drift
    logical :: verbose_

    ! return if no FC2 is available
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    ! optional parameters
    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose

    ! Enforce sum_j Phi(i,j) = 0 exactly, by moving the whole residual
    ! into the self term (phonopy's set_translational_symmetry_compact_fc,
    ! c/phonopy.c). The other direction, sum_i Phi(i,j) = 0, then follows
    ! from the translational symmetry of the full array.
    drift = 0d0
    do ia = 1, c%ncel
       is = (ia-1)*v%fc2_nlat + 1
       summ = 0d0
       do j = 1, v%fc2_nsat
          if (j == is) cycle
          summ = summ + v%fc2(:,:,ia,j)
       end do
       drift = max(drift,maxval(abs(summ + v%fc2(:,:,ia,is))))
       do i = 1, 3
          do j = 1, 3
             v%fc2(i,j,ia,is) = -0.5d0 * (summ(i,j) + summ(j,i))
          end do
       end do
    end do

    if (verbose_) then
       write (uout,'("+ FC2 ACOUSTIC_SUM_RULES: apply the acoustic sum rule to the force constants")')
       write (uout,'("  Sum rule violation before (Hartree/bohr^2) = ",A)') string(drift,'e',12,4)
    end if

  end subroutine vibrations_apply_acoustic

  !> Write the FC2 to a phonopy FORCE_CONSTANTS file, in the compact
  !> format and in the units of the given generator (default: qe,
  !> Ry/bohr^2). If no file is given or if file is empty, use
  !> FORCE_CONSTANTS. The array written is exactly phonopy's compact
  !> force constants for a primitive cell equal to the current cell, so
  !> it can be compared directly against phonopy's own output.
  module subroutine vibrations_write_fc2(v,c,file,full,verbose)
    use tools_io, only: fopen_write, fclose, uout, string, ioj_right
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(in), optional :: file
    logical, intent(in), optional :: full
    logical, intent(in), optional :: verbose

    integer :: lu, ia, ja, il, jl, kl, is, js, ks, k, igen, nrow
    real*8 :: fac
    logical :: verbose_, full_
    character(len=:), allocatable :: file_

    ! return if no FC2 is available
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    ! optional parameters
    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose
    full_ = .false.
    if (present(full)) full_ = full
    file_ = "FORCE_CONSTANTS"
    if (present(file)) then
       if (len(file) > 0) file_ = file
    end if

    ! written in QE units (Ry/bohr^2), like the previous version of this routine
    igen = 3
    fac = 2d0

    if (full_) then
       nrow = v%fc2_nsat
    else
       nrow = c%ncel
    end if

    lu = fopen_write(file_)
    write (lu,'(A,X,A)') string(nrow), string(v%fc2_nsat)
    do ia = 1, c%ncel
       do il = 1, v%fc2_nlat
          is = (ia-1)*v%fc2_nlat + il

          ! in the compact format only the images at the origin are written
          if (.not.full_ .and. il > 1) cycle

          do js = 1, v%fc2_nsat
             ! Phi(ia at L_il, ja at L_jl) = Phi(ia at 0, ja at L_jl - L_il),
             ! by translational invariance; in the compact case L_il = 0 and
             ! this is just ks = js.
             ja = (js-1)/v%fc2_nlat + 1
             jl = js - (ja-1)*v%fc2_nlat
             kl = fc2_ilat(v%fc2_nlat,v%fc2_lkey,v%fc2_madj,v%fc2_lvec(:,jl)-v%fc2_lvec(:,il))
             ks = (ja-1)*v%fc2_nlat + kl

             write (lu,'(A,X,A)') string(is), string(js)
             write (lu,'(3(A,X))') (string(fac*v%fc2(1,k,ia,ks),'f',22,15,ioj_right),k=1,3)
             write (lu,'(3(A,X))') (string(fac*v%fc2(2,k,ia,ks),'f',22,15,ioj_right),k=1,3)
             write (lu,'(3(A,X))') (string(fac*v%fc2(3,k,ia,ks),'f',22,15,ioj_right),k=1,3)
          end do
       end do
    end do
    call fclose(lu)

    if (verbose_) then
       if (full_) then
          write (uout,'("+ Written FC2 file (",A,", full format): ",A)') &
             trim(fc2_unitname(igen)), file_
       else
          write (uout,'("+ Written FC2 file (",A,", compact format): ",A)') &
             trim(fc2_unitname(igen)), file_
       end if
    end if

  end subroutine vibrations_write_fc2

  !> Numerical checks on the force constants, meant to catch a wrong
  !> atom ordering or a wrong unit factor right after a read. Reports
  !> to standard output if verbose. Returns non-zero errmsg only if the
  !> force constants are certainly wrong.
  module subroutine vibrations_check_fc2(v,c,verbose,errmsg)
    use tools_math, only: eigsym
    use tools_io, only: uout, string, ioj_right
    use param, only: atmass
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    logical, intent(in) :: verbose
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: ia, ja, il, jl, is, js, ier
    integer :: lv(3)
    real*8 :: fmax, dsym, dasr, summ(3,3), blk(3,3), eval(3), emin, omin, omax, om

    errmsg = ""
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    fmax = maxval(abs(v%fc2))
    if (fmax <= 0d0) then
       errmsg = "The force constants read are all zero"
       return
    end if

    ! Transpose symmetry: Phi(ia at 0, ja at L) = Phi(ja at 0, ia at -L)^t.
    ! This is the sharpest test available: it is scale-free (so a wrong
    ! unit factor does not affect it) but any permutation of the
    ! supercell atoms breaks it by an amount of order one.
    dsym = 0d0
    do ia = 1, c%ncel
       do js = 1, v%fc2_nsat
          ja = (js-1)/v%fc2_nlat + 1
          jl = js - (ja-1)*v%fc2_nlat
          lv = -v%fc2_lvec(:,jl)
          il = fc2_ilat(v%fc2_nlat,v%fc2_lkey,v%fc2_madj,lv)
          if (il == 0) cycle
          is = (ia-1)*v%fc2_nlat + il
          dsym = max(dsym,maxval(abs(v%fc2(:,:,ia,js) - transpose(v%fc2(:,:,ja,is)))))
       end do
    end do
    dsym = dsym / fmax

    ! Acoustic sum rule drift. Blind to a permutation of the columns
    ! (it sums over all of them), but it catches a truncated read, a
    ! wrong supercell size or a badly completed compact file.
    dasr = 0d0
    do ia = 1, c%ncel
       summ = 0d0
       do js = 1, v%fc2_nsat
          summ = summ + v%fc2(:,:,ia,js)
       end do
       dasr = max(dasr,maxval(abs(summ)))
    end do
    dasr = dasr / fmax

    ! Self blocks: symmetric and positive definite for a stable
    ! structure. Also gives an order-of-magnitude frequency, which is
    ! an unmistakable check on the unit conversion factor.
    emin = huge(1d0)
    omin = huge(1d0)
    omax = 0d0
    do ia = 1, c%ncel
       is = (ia-1)*v%fc2_nlat + 1
       blk = 0.5d0 * (v%fc2(:,:,ia,is) + transpose(v%fc2(:,:,ia,is)))
       om = blk(1,1) + blk(2,2) + blk(3,3)
       call eigsym(blk,3,eval,ier)
       if (ier == 0) emin = min(emin,eval(1))
       om = sqrt(max(om,0d0) / 3d0 / atmass(c%spc(c%atcel(ia)%is)%z)) * freqfactor
       omin = min(omin,om)
       omax = max(omax,om)
    end do

    if (verbose) then
       write (uout,'("+ Checks on the force constants")')
       write (uout,'("  Largest |FC2| element (Hartree/bohr^2) = ",A)') string(fmax,'e',12,4)
       write (uout,'("  Transpose symmetry violation (relative) = ",A)') string(dsym,'e',12,4)
       write (uout,'("  Acoustic sum rule violation (relative) = ",A)') string(dasr,'e',12,4)
       write (uout,'("  Lowest eigenvalue of the self blocks (Hartree/bohr^2) = ",A)') string(emin,'e',12,4)
       write (uout,'("  Estimated frequency range (cm^-1) = ",A," to ",A)') &
          string(omin,'f',10,2), string(omax,'f',10,2)
       if (dsym > fc2_epssym) &
          write (uout,'("  WARNING: the transpose symmetry is badly violated; the atom ordering,")')
       if (dsym > fc2_epssym) &
          write (uout,'("           the supercell or the structure are probably wrong")')
       if (emin < 0d0) &
          write (uout,'("  WARNING: a self block is not positive definite")')
    end if

    if (dsym > 1d0) &
       errmsg = "The force constants violate the transpose symmetry Phi(i,j) = Phi(j,i)^t by " //&
          string(dsym,'e',12,4) // " (relative). The atom ordering, the supercell, or the structure &
          &do not correspond to the file"

  end subroutine vibrations_check_fc2

  !> Calculate frequencies and eigenvectors from the FC2 for a single
  !> q (fractional coordiantes in reciprocal space); adds the
  !> resulting frequencies and eigenvectors to the v type.  If the
  !> optional arguments freqo or veco are present, return the
  !> frequencies and/or eigenvectors in these variables instead of
  !> overwriting the vibrational info in v. Frequencies in output are
  !> in ascending order.
  module subroutine vibrations_calculate_q(v,c,q,freqo,veco)
    use tools_io, only: ferror, faterr
    use param, only: icrd_crys, atmass, tpi, img
    use types, only: realloc
    use tools_math, only: eigherm
    use tools_io, only: ferror, faterr
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: q(3)
    real*8, intent(inout), allocatable, optional :: freqo(:)
    complex*16, intent(inout), allocatable, optional :: veco(:,:)

    !! NOT REIMPLEMENTED YET. This routine assumed the old fc2(3,3,ncel,ncel)
    !! layout, in which both indices ran over the atoms of the cell. The body
    !! is kept commented out below until it is rewritten for the cell-compact
    !! fc2(3,3,ncel,nsat) array (second index over the supercell).
    call ferror('vibrations_calculate_q','not yet reimplemented for the cell-compact force constants',faterr)

!!    complex*16 :: dm_local(3,3), phase
!!    real*8 :: sqrt_ij
!!    integer :: nat
!!    real*8 :: x1(3), x2(3), maxdisteps
!!    integer :: i, jj, j, ier
!!    integer, allocatable :: nida(:), lvec(:,:), mult(:)
!!    real*8, allocatable:: eval(:), dist(:), freq(:), rused(:)
!!    complex*16, allocatable :: dm(:,:)
!!    logical :: varoutput
!!
!!    real*8, parameter :: epsgen = 1d-4
!!    real*8, parameter :: epsneigh = 1d-5
!!
!!    ! return if no FC2 is available
!!    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return
!!
!!    ! process input arguments
!!    varoutput = present(freqo) .or. present(veco)
!!
!!    ! maximum distance of all pairs of atoms in the unit cell
!!    maxdisteps = 0d0
!!    do i = 1, c%ncel
!!       do j = i+1, c%ncel
!!          maxdisteps = max(maxdisteps,c%eql_distance(c%atcel(i)%x, c%atcel(j)%x))
!!       end do
!!    end do
!!
!!    ! build the dynamical matrix for this q-point (Parlinski's recipe)
!!    allocate(dm(c%ncel*3,c%ncel*3),rused(c%ncel),mult(c%ncel))
!!    dm = 0d0
!!    do i = 1, c%ncel
!!       rused = huge(1d0)
!!       mult = 0
!!       call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.true.,nat,nida,dist=dist,lvec=lvec,up2d=maxdisteps+epsgen)
!!       do jj = 1, nat
!!          j = nida(jj)
!!          if (dist(jj) <= rused(j) + epsneigh) then
!!             x1 = c%atcel(i)%x
!!             x2 = c%atcel(j)%x + lvec(:,jj)
!!
!!             sqrt_ij = sqrt(atmass(c%spc(c%atcel(i)%is)%z) * atmass(c%spc(c%atcel(j)%is)%z))
!!             phase = exp(-tpi * img * dot_product(q,x1 - x2))
!!
!!             dm_local = v%fc2(:,:,i,j) * phase / sqrt_ij
!!             dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) = dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) + dm_local
!!
!!             rused(j) = min(rused(j),dist(jj))
!!             mult(j) = mult(j) + 1
!!          end if
!!       end do
!!
!!       ! apply the multiplicities
!!       do j = 1, c%ncel
!!          dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) = dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) / mult(j)
!!       end do
!!    end do
!!
!!    ! impose hermiticity
!!    dm = (dm + transpose(conjg(dm)))/2
!!
!!    ! diagonalize
!!    allocate(eval(c%ncel*3),freq(c%ncel*3))
!!    call eigherm(dm,c%ncel*3,eval,ier)
!!    if (ier /= 0) &
!!       call ferror('vibrations_calculate_q','Error in diagonalization',faterr)
!!
!!    ! calculate the frequencies
!!    freq = sign(sqrt(abs(eval)) * freqfactor,eval)
!!
!!    ! save output
!!    if (varoutput) then
!!       if (present(freqo)) freqo = freq
!!       if (present(veco)) veco = dm
!!    else
!!       ! save the info
!!       v%nqpt = v%nqpt + 1
!!       if (.not.allocated(v%qpt)) allocate(v%qpt(3,v%nqpt))
!!       if (.not.allocated(v%freq)) then
!!          v%nfreq = 3 * c%ncel
!!          allocate(v%freq(3*c%ncel,v%nqpt))
!!       end if
!!       if (.not.allocated(v%vec)) allocate(v%vec(3,c%ncel,3*c%ncel,v%nqpt))
!!       if (v%nqpt > size(v%qpt,2)) call realloc(v%qpt,3,2*v%nqpt)
!!       if (v%nqpt > size(v%freq,2)) call realloc(v%freq,size(v%freq,1),2*v%nqpt)
!!       if (v%nqpt > size(v%vec,4)) call realloc(v%vec,size(v%vec,1),size(v%vec,2),size(v%vec,3),2*v%nqpt)
!!       v%qpt(:,v%nqpt) = q
!!       v%freq(:,v%nqpt) = freq
!!       do i = 1, 3*c%ncel
!!          v%vec(:,:,i,v%nqpt) = reshape(dm(:,i),(/3,c%ncel/))
!!       end do
!!       v%hasvibs = .true.
!!    end if
!!
  end subroutine vibrations_calculate_q

  !> Calculate sound velocities in the reciprocal space direction
  !> given by q (Cartesian coordinates). q cannot be zero, and it is
  !> normalized on input. c = crystal structure. vs = output sound
  !> velocities in increasing order. Automatically applies acoustic
  !> sum rules and changes the FC2. Sound velocities are in ascending
  !> order and in units of m/s.
  module subroutine vibrations_calculate_vs(v,c,q,vs)
    use tools_io, only: ferror, faterr
    use tools_io, only: ferror, faterr
    use param, only: cm1tohz, bohrtom
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: q(3)
    real*8, intent(out) :: vs(3)

    !! NOT REIMPLEMENTED YET. This routine assumed the old fc2(3,3,ncel,ncel)
    !! layout, in which both indices ran over the atoms of the cell. The body
    !! is kept commented out below until it is rewritten for the cell-compact
    !! fc2(3,3,ncel,nsat) array (second index over the supercell).
    call ferror('vibrations_calculate_vs','not yet reimplemented for the cell-compact force constants',faterr)

!!    real*8 :: normq, qn(3), qthis(3) !, h, fd0(3), fdn, ratio, fdiff(3)
!!    integer :: i, j
!!    real*8, allocatable :: freq(:)
!!    real*8 :: fcur(3,3)
!!
!!    real*8, parameter :: epsnq = 1d-10 ! cutoff for zero q-point
!!
!!    ! return if no FC2 is available
!!    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return
!!
!!    ! normalize q
!!    normq = norm2(q)
!!    if (normq < epsnq) &
!!       call ferror('vibrations_calculate_vs','zero-length q-point',faterr)
!!    qn = q / normq
!!
!!    ! prepare for the vs calculation if not already done
!!    if (any(v%fc2_acoustic < 0)) &
!!       call v%calculate_vs_prepare(c,q,vs,.true.)
!!
!!    ! calculate the frequencies for the finite difference
!!    do i = -1, 1
!!       qthis = (v%fc2_vs_center + i * v%fc2_vs_delta) * qn
!!       qthis = c%rc2rx(qthis)
!!       call v%calculate_q(c,qthis,freqo=freq)
!!       fcur(i+2,:) = freq(v%fc2_acoustic)
!!    end do
!!
!!    ! first calculation of vs
!!    vs = 0d0
!!    do j = 1, 3
!!       vs(j) = sum(fcur(:,j) * der_coef1) / (2d0 * v%fc2_vs_delta)
!!    end do
!!    vs = vs * cm1tohz * bohrtom
!!
  end subroutine vibrations_calculate_vs

  !> Prepare for the calculation of sound velocities by locating the
  !> acoustic branches and finding the reciprocal distance range
  !> where the slope is most easily calculated. For the preparation,
  !> use reciprocal space direction given by q (Cartesian
  !> coordinates). q cannot be zero, and it is normalized on input. c
  !> = crystal structure. vs = output sound velocities in increasing
  !> order. Automatically applies acoustic sum rules and changes the
  !> FC2. Sound velocities are in ascending order and in units of m/s.
  !> If verbose, write progress to standard output.
  module subroutine vibrations_calculate_vs_prepare(v,c,q,vs,verbose)
    use tools_io, only: ferror, faterr
    use tools_io, only: ferror, faterr, uout, string
    use param, only: cm1tohz, bohrtom
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: q(3)
    real*8, intent(out) :: vs(3)
    logical, intent(in) :: verbose

    !! NOT REIMPLEMENTED YET. This routine assumed the old fc2(3,3,ncel,ncel)
    !! layout, in which both indices ran over the atoms of the cell. The body
    !! is kept commented out below until it is rewritten for the cell-compact
    !! fc2(3,3,ncel,nsat) array (second index over the supercell).
    call ferror('vibrations_calculate_vs_prepare','not yet reimplemented for the cell-compact force constants',faterr)

!!    real*8 :: normq, qn(3), qthis(3)!, h, fd0(3), fdn, ratio, fdiff(3)
!!    integer :: i, j, idx(3), iminlast
!!    real*8, allocatable :: freq(:)!, ff(:,:)
!!    real*8 :: fcur(3,3), f2, f2min
!!
!!    real*8, parameter :: ac_thr = 0.01d0 ! threshold for acoustic ident, cm-1
!!    real*8, parameter :: epsnq = 1d-10 ! cutoff for zero q-point
!!    real*8, parameter :: hstep = 1d-4 ! step for the bracketing
!!    integer, parameter :: isafesteps = 10 ! number of steps after minimum before exiting
!!
!!    ! return if no FC2 is available
!!    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return
!!
!!    ! header
!!    if (verbose) &
!!       write (uout,'("+ Preparing for sound velocity calculation")')
!!
!!    ! normalize q
!!    normq = norm2(q)
!!    if (normq < epsnq) &
!!       call ferror('vibrations_calculate_vs','zero-length q-point',faterr)
!!    qn = q / normq
!!    if (verbose) then
!!       qthis = c%rc2rx(qn)
!!       qthis = qthis / norm2(qthis)
!!       write (uout,'("  Direction (rec-cryst): ",3(A," "))') (string(qthis(i),'f',10,5),i=1,3)
!!    end if
!!
!!    ! apply acoustic sum rules
!!    if (verbose) &
!!       write (uout,'("  Applying acoustic sum rules...")')
!!    call v%apply_acoustic(c)
!!
!!    ! calculate the frequencies at gamma, locate the acoustic branches
!!    qthis = 0d0
!!    call v%calculate_q(c,qthis,freqo=freq)
!!    idx = 0
!!    j = 0
!!    do i = 1, size(freq,1)
!!       if (abs(freq(i)) < ac_thr) then
!!          j = j + 1
!!          if (j > 3) &
!!             call ferror('vibrations_calculate_vs','too many zero frequencies',faterr)
!!          idx(j) = i
!!       end if
!!    end do
!!    if (any(idx == 0)) &
!!       call ferror('vibrations_calculate_vs','too few zero frequencies',faterr)
!!    if (verbose) then
!!       write (uout,'("  Acoustic branches: ",2(A," (",A,"), "),A," (",A,")")') &
!!          string(idx(1)), string(freq(idx(1)),'f',decimal=8),&
!!          string(idx(2)), string(freq(idx(2)),'f',decimal=8),&
!!          string(idx(3)), string(freq(idx(3)),'f',decimal=8)
!!    end if
!!    v%fc2_acoustic = idx
!!
!!    ! initialize the vs bracketing
!!    if (verbose) &
!!       write (uout,'("# h(bohr-1)  ---- frequencies (cm-1)    ----  ----- sound velocities (m/s) -----   curvature (au)")')
!!    do i = 1, 3
!!       qthis = i * hstep * qn
!!       qthis = c%rc2rx(qthis)
!!       call v%calculate_q(c,qthis,freqo=freq)
!!       fcur(i,:) = freq(idx)
!!       if (verbose .and. i < 3) &
!!          write (uout,'("  "4(A," "))') string(i * hstep,'e',10,3), (string(freq(idx(j)),'e',10,3),j=1,3)
!!    end do
!!
!!    ! first calculation of vs
!!    vs = 0d0
!!    f2 = 0d0
!!    do j = 1, 3
!!       vs(j) = sum(fcur(:,j) * der_coef1) / (2d0 * hstep)
!!       f2 = f2 + abs(sum(fcur(:,j) * der_coef2) / (hstep * hstep))
!!    end do
!!    vs = vs * cm1tohz * bohrtom
!!    if (verbose) then
!!       write (uout,'("  "8(A," ")," (best)")') string(3 * hstep,'e',10,3),&
!!          (string(freq(idx(j)),'e',10,3),j=1,3),&
!!          (string(vs(j),'f',12,3),j=1,3), string(f2,'f',15,6)
!!    end if
!!    f2min = f2
!!    iminlast = 3
!!
!!    ! continue taking steps until we find a minimum of f2
!!    i = 3
!!    do while (.true.)
!!       i = i + 1
!!       qthis = i * hstep * qn
!!       qthis = c%rc2rx(qthis)
!!       call v%calculate_q(c,qthis,freqo=freq)
!!       fcur(1:2,:) = fcur(2:3,:)
!!       fcur(3,:) = freq(idx)
!!
!!       vs = 0d0
!!       f2 = 0d0
!!       do j = 1, 3
!!          vs(j) = sum(fcur(:,j) * der_coef1) / (2d0 * hstep)
!!          f2 = f2 + abs(sum(fcur(:,j) * der_coef2) / (hstep * hstep))
!!       end do
!!       vs = vs * cm1tohz * bohrtom
!!       if (verbose) then
!!          if (f2 < f2min) then
!!             write (uout,'("  "8(A," ")," (best)")') string(i * hstep,'e',10,3),&
!!                (string(freq(idx(j)),'e',10,3),j=1,3),&
!!                (string(vs(j),'f',12,3),j=1,3), string(f2,'f',15,6)
!!          else
!!             write (uout,'("  "8(A," "))') string(i * hstep,'e',10,3), (string(freq(idx(j)),'e',10,3),j=1,3),&
!!                (string(vs(j),'f',12,3),j=1,3), string(f2,'f',15,6)
!!          end if
!!       end if
!!       if (f2 < f2min) then
!!          f2min = f2
!!          iminlast = i
!!       end if
!!
!!       if (iminlast > 0 .and. (i-iminlast) > isafesteps) exit
!!    end do
!!
!!    ! save the info and final message
!!    v%fc2_vs_delta = hstep
!!    v%fc2_vs_center = (iminlast-1) * hstep
!!    if (verbose) then
!!       write (uout,'("# Finite difference distances (bohr-1): ",3(A," ")/)') &
!!          string((iminlast-2) * hstep,'e',10,3), string((iminlast-1) * hstep,'e',10,3),&
!!          string(iminlast * hstep,'e',10,3)
!!    end if
!!
  end subroutine vibrations_calculate_vs_prepare

  !> Calculate thermodynamic properties at temperature T using the
  !> vibrational frequencies in v. The routine assumes the frequencies
  !> have all equal weight (i.e. it is a mesh)
  module subroutine vibrations_calculate_thermo(v,t,zpe,fvib,svib,cv)
    use tools_io, only: ferror, faterr
    class(vibrations), intent(inout) :: v
    real*8, intent(in) :: t
    real*8, intent(out) :: zpe, fvib, svib, cv

    integer :: i, j, ncount
    real*8 :: nu, x, y, nut, nue, rt, ff, l1mx, nutdiv
    real*8 :: ym1
    real*8, parameter :: cutoff_frequency = 1d0 ! cutoff frequency, cm-1

    real*8, parameter :: R = 8.314462618d0 ! molar gas constant (kB/NA, J/K/mol)
    real*8, parameter :: small1 = 50000d0 * cminv_to_K / huge(1d0) ! protection against zerodiv in nu/(kB*T)
    real*8, parameter :: small2 = 0.5d0 * log(huge(1d0)) ! protection against overflow in exp(nu/kB*T)**2
    real*8, parameter :: small3 = 1d0 / sqrt(huge(1d0)) ! protection against 1/(e^(nu/kt)-1)^2 overflow

    ! calculate the thermodynamic properties
    zpe = 0d0
    fvib = 0d0
    svib = 0d0
    cv = 0d0
    ncount = 0
    do i = 1, v%nqpt
       do j = 1, v%nfreq
          ! prepare and cycle if this is a low frequency
          nu = v%freq(j,i)
          nut = nu * cminv_to_K      ! frequency in K
          nue = nu * cminv_to_kJmol ! frequency in kJ/mol
          rt = R / 1000d0 * t        ! RT in kJ/mol
          if (nu < cutoff_frequency) cycle
          ncount = ncount + 1

          ! some quantities we will need, with protection
          if (t < small1) then
             x = 0d0
          else
             nutdiv = nut / t
             x = exp(-nutdiv)
          end if
          l1mx = log(1d0 - x)

          ! fvib and zpe calculation
          zpe = zpe + 0.5d0 * nue
          fvib = fvib + 0.5d0 * nue + rt * log(1d0 - x)

          ! svib
          if (t > small1) then
             svib = svib + R * (-l1mx + nutdiv * x / (1d0 - x))
          end if

          ! cv
          if (nutdiv <= small2) then
             y = exp(nutdiv)
             ym1 = y - 1
             if (ym1 >= small3) then
                cv = cv + R * nutdiv * nutdiv * y / ((y-1) * (y-1))
             end if
          end if
       end do
    end do
    if (ncount == 0) &
       call ferror('vibrations_calculate_thermo','no frequencies available for THERMO',faterr)

    ! renormalize after taking out the zero frequencies
    ff = real(v%nfreq,8) / real(ncount,8)
    zpe = zpe * ff
    fvib = fvib * ff
    svib = svib * ff
    cv = cv * ff

  end subroutine vibrations_calculate_thermo

  !> Nullify the FC2 elements where the atoms are farther apart
  !> than dist.
  module subroutine vibrations_trim_fc2(v,c,dist,verbose)
    use tools_io, only: ferror, faterr
    use tools_io, only: uout, string
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: dist
    logical, intent(in), optional :: verbose

    !! NOT REIMPLEMENTED YET. This routine assumed the old fc2(3,3,ncel,ncel)
    !! layout, in which both indices ran over the atoms of the cell. The body
    !! is kept commented out below until it is rewritten for the cell-compact
    !! fc2(3,3,ncel,nsat) array (second index over the supercell).
    call ferror('vibrations_trim_fc2','not yet reimplemented for the cell-compact force constants',faterr)

!!    logical :: verbose_
!!    integer :: i, j, n, nmax
!!    real*8 :: dthis, perc
!!
!!    ! return if no FC2 is available
!!    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return
!!
!!    ! optional parameters
!!    verbose_ = .false.
!!    if (present(verbose)) verbose_ = verbose
!!
!!    ! nullify
!!    n = 0
!!    do i = 1, c%ncel
!!       do j = i+1, c%ncel
!!          dthis = c%eql_distance(c%atcel(i)%x, c%atcel(j)%x)
!!          if (dthis >= dist) then
!!             v%fc2(:,:,i,j) = 0d0
!!             v%fc2(:,:,j,i) = 0d0
!!             n = n + 1
!!          end if
!!       end do
!!    end do
!!
!!    ! final message
!!    if (verbose_) then
!!       nmax = 3 * 3 * c%ncel * v%fc2_nsat
!!       perc = real(2 * n * 3 * 3,8) / real(nmax,8) * 100d0
!!       write (uout,'("+ FC2 TRIM: nullify all FC components beyond a distance")')
!!       write (uout,'("  Number of terms made zero = ",A," out of ",A," (",A,"%)")') &
!!          string(n), string(nmax), string(perc,'f',decimal=2)
!!    end if
!!
  end subroutine vibrations_trim_fc2

  !> Nullify the FC2 elements whose absolute value is below eps0.
  module subroutine vibrations_zero_fc2(v,c,eps0,verbose)
    use tools_io, only: uout, string
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: eps0
    logical, intent(in), optional :: verbose

    logical :: verbose_
    integer :: i, j, k, l, n, nmax
    real*8 :: perc

    ! return if no FC2 is available
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    ! optional parameters
    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose

    ! nullify
    n = 0
    do j = 1, v%fc2_nsat
       do i = 1, c%ncel
          do l = 1, 3
             do k = 1, 3
                if (abs(v%fc2(k,l,i,j)) < eps0) then
                   v%fc2(k,l,i,j) = 0d0
                   n = n + 1
                end if
             end do
          end do
       end do
    end do

    ! final message
    if (verbose_) then
       nmax = 3 * 3 * c%ncel * v%fc2_nsat
       perc = real(n,8) / real(nmax,8) * 100d0
       write (uout,'("+ FC2 ZERO: nullify all FC components with abs lower than a value")')
       write (uout,'("  Number of terms made zero = ",A," out of ",A," (",A,"%)")') &
          string(n), string(nmax), string(perc,'f',decimal=2)
    end if

  end subroutine vibrations_zero_fc2

  !> Generate a rattled structure based on the vibrational frequencies
  !> (at Gamma) in v and the molecular/crystal structure in c using
  !> temperature temp. The phonon rattled structure samples normal
  !> modes according to a Boltzmann distribution based on the given
  !> temperature. Returns the seed for the new structure in
  !> seed.
  !> Sources: internal-documents/phonon_rattle,
  !> https://arxiv.org/pdf/2210.03531 , Messiah (chap. 12). Equivalent
  !> to the implementation in the alamode code
  !> (https://github.com/ttadano/alamode/).
  module subroutine vibrations_phonon_rattle(v,c,temp,seed)
    use crystalseedmod, only: crystalseed
    use tools, only: mergesort
    use tools_math, only: gauss_random
    use param, only: atmass, pi, hbar, clight, pcamu, bohrtom
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: temp
    type(crystalseed), intent(out) :: seed

    real*8, allocatable :: xat(:,:)
    real*8 :: ff, fterm, xn
    real*8 :: qsq
    integer :: i, j, iqpt, ini
    integer, allocatable :: idx(:)

    real*8, parameter :: fcap = 250d0 ! lower bound for the frequencies (cm-1)

    if (c%ismolecule) then
       ! use the first and only q-point
       iqpt = 1
       if (.not.v%hasvibs.or..not.allocated(v%freq).or..not.allocated(v%qpt).or.&
          .not.allocated(v%vec)) return
    else
       ! try to find the gamma q-point
       iqpt = -1
       if (v%hasvibs.and.allocated(v%freq).and.allocated(v%qpt).and.allocated(v%vec)) then
          do i = 1, v%nqpt
             if (all(abs(v%qpt(:,i)) < 1d-5)) then
                iqpt = i
                exit
             end if
          end do
       end if

       ! if no gamma q-point is available, check if we have the FC2
       if (iqpt < 0 .and. v%hasfc2 .and. allocated(v%fc2)) then
          ! get the frequencies (cm-1) and eigenvectors at Gamma
          call v%calculate_q(c,(/0d0,0d0,0d0/))
          iqpt = v%nqpt
       end if

       ! could not figure out how to obtain the gamma frequencies and eigenvectors
       if (iqpt < 0) return
    end if

    ! copy the seed from the given system to the output seed
    call c%makeseed(seed,.false.)

    ! prepare the atomic positions array in Cartesian coordinates
    allocate(xat(3,c%ncel))
    do j = 1, c%ncel
       if (c%ismolecule) then
          xat(:,j) = c%atcel(j)%r + c%molx0
       else
          xat(:,j) = c%atcel(j)%r
       end if
    end do

    ! sort the frequencies by absolute value
    allocate(idx(v%nfreq))
    do i = 1, v%nfreq
       idx(i) = i
    end do
    if (.not.c%ismolecule) then
       call mergesort(abs(v%freq(:,iqpt)),idx,1,v%nfreq)
       ini = 4
    else
       ini = 1
    end if

    ! run over all modes; skip the first three
    do i = ini, v%nfreq
       ! cap frequencies at some low value
       ff = max(fcap,abs(v%freq(idx(i),iqpt)))

       ! calculate the <q^2> according to the Boltzmann distribution
       qsq = 0.5d0 * &
          hbar / (2d0 * pi * clight * ff * 100d0) * (1e3 / pcamu) / bohrtom**2 / &
          tanh(0.5d0 * ff * cminv_to_K / temp)

       ! sample normal distribution
       xn = sqrt(qsq) * gauss_random()

       do j = 1, seed%nat
          fterm = xn / sqrt(atmass(seed%spc(seed%is(j))%z))
          xat(:,j) = xat(:,j) + fterm * real(v%vec(:,j,idx(i),iqpt),8)
       end do
    end do

    ! copy the new atomic positions into the seed
    if (c%ismolecule) seed%useabr = 0
    do j = 1, seed%nat
       if (c%ismolecule) then
          seed%x(:,j) = xat(:,j)
       else
          seed%x(:,j) = c%c2x(xat(:,j))
       end if
    end do

  end subroutine vibrations_phonon_rattle

  !xx! private procedures

  !> Detect the format for a file containing molecular or
  !> crystal vibrations.
  subroutine vibrations_detect_format(file,ivformat)
    use param, only: ivformat_unknown, ivformat_matdynmodes, ivformat_matdyneig,&
       ivformat_qedyn, ivformat_phonopy_ascii, ivformat_phonopy_yaml,&
       ivformat_phonopy_hdf5, ivformat_crystal_out, ivformat_gaussian_log,&
       ivformat_gaussian_fchk, ivformat_phonopy_fc2, ivformat_castep_phonon, dirsep
    use tools_io, only: equal, lower
    character*(*), intent(in) :: file
    integer, intent(out) :: ivformat

    character(len=:), allocatable :: basename, wextdot, wextdot2, wext_
    integer :: idx

    basename = file(index(file,dirsep,.true.)+1:)
    wext_ = basename(index(basename,'_',.true.)+1:)
    idx = index(basename,'.',.true.)
    wextdot = basename(idx+1:)
    if (idx > 0) then
       wextdot2 = basename(index(basename(1:idx-1),'.',.true.)+1:)
    else
       wextdot2 = ""
    end if

    if (equal(lower(wextdot),'modes')) then
       ivformat = ivformat_matdynmodes
    elseif (equal(lower(wextdot),'eig')) then
       ivformat = ivformat_matdyneig
    elseif (equal(lower(wextdot),'ascii')) then
       ivformat = ivformat_phonopy_ascii
    elseif (equal(lower(wextdot),'yaml')) then
       ivformat = ivformat_phonopy_yaml
    elseif (equal(lower(wextdot),'hdf5')) then
       ivformat = ivformat_phonopy_hdf5
    elseif (equal(lower(wextdot),'out')) then
       ivformat = ivformat_crystal_out
    elseif (equal(lower(wextdot),'log')) then
       ivformat = ivformat_gaussian_log
    elseif (equal(lower(wextdot),'fchk')) then
       ivformat = ivformat_gaussian_fchk
    elseif (equal(lower(wextdot),'phonon')) then
       ivformat = ivformat_castep_phonon
    elseif (isdynfile(wextdot)) then
       ivformat = ivformat_qedyn
    elseif (index(basename,'FORCE_CONSTANTS') > 0) then
       ivformat = ivformat_phonopy_fc2
    else
       ivformat = ivformat_unknown
    endif

  contains
    function isdynfile(str)
      logical :: isdynfile
      character*(*), intent(in) :: str

      isdynfile = .false.
      if (len(str) < 3) return
      if (str(1:3) /= "dyn") return
      if (len(str) > 3) then
         if (str(4:4) == "0") return
      end if
      isdynfile = .true.

    end function isdynfile

  end subroutine vibrations_detect_format

  !> Read vibration data with Quantum ESPRESSO matdyn.modes and
  !> matdyn.eig format from a file, and return it in vib. If error,
  !> return non-zero errmsg.
  subroutine read_matdyn_modes(v,c,file,ivformat,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalseedmod, only: read_alat_from_qeout
    use param, only: atmass, isformat_r_qeout, ivformat_matdynmodes
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: ivformat
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, nline, nline1, nat, idx
    logical :: ok
    character(len=:), allocatable :: line
    integer :: iqpt, ifreq, iat, iz
    real*8 :: xdum(6), alat
    ! integer :: jfreq ! checking normalization
    ! complex*16 :: summ

    ! initialize
    errmsg = ""
    lu = -1

    ! read the alat from the crystal source file
    if (c%isformat /= isformat_r_qeout) then
       errmsg = "Error reading alat: the crystal structure must be a QE output"
       goto 999
    end if
    call read_alat_from_qeout(c%file,alat,errmsg,ti)
    if (len_trim(errmsg) > 0) goto 999

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat

    ! first pass: determine nqpt and nfreq
    nline = 0
    v%nqpt = 0
    v%nfreq = 0
    do while (getline_raw(lu,line,.false.))
       nline = nline + 1
       if (len(line) >= 4) then
          if (line(1:4) == " q =") v%nqpt = v%nqpt + 1
       end if
       if (v%nqpt == 1) then
          if (len(line) >= 9) then
             if (line(1:9) == "     freq") then
                v%nfreq = v%nfreq + 1
                if (v%nfreq == 1) then
                   nline1 = nline
                elseif (v%nfreq == 2) then
                   nat = nline - nline1 - 1
                end if
             end if
          end if
       end if
    end do
    if (v%nfreq == 0) then
       errmsg = "Found no frequencies in file: " // trim(file)
       goto 999
    end if
    if (v%nqpt == 0) then
       errmsg = "Found no q-points in file: "  // trim(file)
       goto 999
    end if
    if (nat /= c%ncel) then
       errmsg = "Number of atomic displacements inconsistent with number of atoms in crystal structure"
       goto 999
    end if

    ! allocate arrays
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,v%nqpt))
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(v%nfreq,v%nqpt))
    if (allocated(v%vec)) deallocate(v%vec)
    allocate(v%vec(3,c%ncel,v%nfreq,v%nqpt))

    ! second pass: read the information
    errmsg = "Error reading modes file: " // trim(file)
    rewind(lu)
    v%nqpt = 0
    do while (getline_raw(lu,line,.false.))
       ok = .false.
       if (len(line) >= 4) ok = (line(1:4) == " q =")
       if (ok) then
          v%nqpt = v%nqpt + 1
          read (line(5:),*,end=999,err=999) v%qpt(:,v%nqpt)
          v%qpt(:,v%nqpt) = c%rc2rx(v%qpt(:,v%nqpt)) / alat

          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          do ifreq = 1, v%nfreq
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) goto 999
             idx = index(line,'=',back=.true.)
             read(line(idx+1:),*) v%freq(ifreq,v%nqpt)
             do iat = 1, c%ncel
                ok = getline_raw(lu,line,.false.)
                if (.not.ok) goto 999
                idx = index(line,'(')
                read(line(idx+1:),*,end=999,err=999) xdum
                v%vec(1,iat,ifreq,v%nqpt) = cmplx(xdum(1),xdum(2),8)
                v%vec(2,iat,ifreq,v%nqpt) = cmplx(xdum(3),xdum(4),8)
                v%vec(3,iat,ifreq,v%nqpt) = cmplx(xdum(5),xdum(6),8)
             end do
          end do
       end if
    end do

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    if (ivformat == ivformat_matdynmodes) then
       do iat = 1, c%ncel
          iz = c%spc(c%atcel(iat)%is)%z
          v%vec(:,iat,:,:) = v%vec(:,iat,:,:) * sqrt(atmass(iz))
       end do
    end if

    ! normalize
    do iqpt = 1, v%nqpt
       do ifreq = 1, v%nfreq
          v%vec(:,:,ifreq,iqpt) = v%vec(:,:,ifreq,iqpt) / &
             sqrt(sum(v%vec(:,:,ifreq,iqpt)*conjg(v%vec(:,:,ifreq,iqpt))))
       end do
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, v%nqpt
    !    do ifreq = 1, v%nfreq
    !       do jfreq = 1, v%nfreq
    !          summ = sum(v%vec(:,:,ifreq,iqpt)*conjg(v%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do
    ! stop 1

    ! wrap up
    errmsg = ""
    v%hasvibs = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_matdyn_modes

  !> Read vibration data with Quantum ESPRESSO *.dyn* file format, and
  !> return it in vib. If error, return non-zero errmsg.
  subroutine read_qe_dyn(v,c,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalseedmod, only: read_alat_from_qeout
    use param, only: atmass, isformat_r_qeout, ivformat_qedyn
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok
    character(len=:), allocatable :: line
    integer :: iz, i, j, lu
    real*8 :: xdum(6), alat
    ! integer :: jfreq, iqpt, ifreq ! checking normalization
    ! complex*16 :: summ

    ! initialize
    errmsg = "Error reading dyn file: " // trim(file)
    lu = -1

    ! read the alat from the crystal source file
    if (c%isformat /= isformat_r_qeout) then
       errmsg = "Error reading alat: the crystal structure must be a QE output"
       goto 999
    end if
    call read_alat_from_qeout(c%file,alat,errmsg,ti)
    if (len_trim(errmsg) > 0) goto 999

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_qedyn

    ! advance to the header
    v%nqpt = 1
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,v%nqpt))
    do while (getline_raw(lu,line,.false.))
       if (len(line) >= 39) then
          if (line(1:39) == "     Diagonalizing the dynamical matrix") exit
       end if
    end do

    ! read qpt coordinates
    ok = getline_raw(lu,line,.false.)
    ok = ok .and. getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    if (len(line) < 11) goto 999
    read (line(11:),*,end=999,err=999) v%qpt(:,1)
    v%qpt(:,1) = c%rc2rx(v%qpt(:,1)) / alat

    ! allocate frequencies
    v%nfreq = 3 * c%ncel
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(v%nfreq,v%nqpt))
    if (allocated(v%vec)) deallocate(v%vec)
    allocate(v%vec(3,c%ncel,v%nfreq,v%nqpt))

    ! read the frequencies
    ok = getline_raw(lu,line,.false.)
    ok = ok .and. getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    do i = 1, v%nfreq
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999
       if (len(line) < 42) goto 999
       read(line(43:),*,end=999,err=999) v%freq(i,1)

       do j = 1, c%ncel
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          if (len(line) < 2) goto 999
          read(line(3:),*,end=999,err=999) xdum
          v%vec(1,j,i,1) = cmplx(xdum(1),xdum(2),8)
          v%vec(2,j,i,1) = cmplx(xdum(3),xdum(4),8)
          v%vec(3,j,i,1) = cmplx(xdum(5),xdum(6),8)
       end do
    end do
    ok = getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    if (index(line,'****') == 0) then
       errmsg = "Number of atomic displacements inconsistent with number of atoms in crystal structure"
       goto 999
    end if

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       v%vec(:,i,:,:) = v%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, v%nfreq
       v%vec(:,:,i,1) = v%vec(:,:,i,1) / &
          sqrt(sum(v%vec(:,:,i,1)*conjg(v%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, v%nqpt
    !    do ifreq = 1, v%nfreq
    !       do jfreq = 1, v%nfreq
    !          summ = sum(v%vec(:,:,ifreq,iqpt)*conjg(v%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    v%hasvibs = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_qe_dyn

  !> Read vibration data from a phonopy ascii file (ANIME keyword),
  !> and return it in vib. If error, return non-zero errmsg.
  subroutine read_phonopy_ascii(v,c,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw
    use param, only: atmass, ivformat_phonopy_ascii, cm1tothz
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok
    character(len=:), allocatable :: line
    integer :: lu, idx, i, ifreq, iz
    real*8 :: xdum(6)
    ! integer :: jfreq, iqpt ! checking normalization
    ! complex*16 :: summ

    ! initialize
    errmsg = "Error reading ascii file: " // trim(file)
    lu = -1

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_phonopy_ascii

    ! allocate and initialize
    v%nqpt = 1
    v%nfreq = 3 * c%ncel
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,v%nqpt))
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(v%nfreq,v%nqpt))
    if (allocated(v%vec)) deallocate(v%vec)
    allocate(v%vec(3,c%ncel,v%nfreq,v%nqpt))

    ! advance to the header
    do while (getline_raw(lu,line,.false.))
       if (len(line) >= 10) then
          if (line(1:10) == "#metaData:") exit
       end if
    end do

    ! get the q-points
    idx = index(line,"[")
    line = line(idx+1:)
    call strip(line)
    read (line,*,err=999,end=999) v%qpt(:,1), v%freq(1,1)

    ! get the rest of the info
    do ifreq = 1, v%nfreq
       ! read the eigenvector
       do i = 1, c%ncel
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          call strip(line)
          read (line,*,err=999,end=999) xdum
          v%vec(1,i,ifreq,1) = cmplx(xdum(1),xdum(4),8)
          v%vec(2,i,ifreq,1) = cmplx(xdum(2),xdum(5),8)
          v%vec(3,i,ifreq,1) = cmplx(xdum(3),xdum(6),8)
       end do

       ! advance to the next metadata
       if (ifreq < v%nfreq) then
          ok = getline_raw(lu,line,.false.)
          ok = ok .and. getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          idx = index(line,"[")
          line = line(idx+1:)
          call strip(line)
          read (line,*,err=999,end=999) xdum(1:3), v%freq(ifreq+1,1)
       end if
    end do

    ! THz to cm-1
    v%freq = v%freq / cm1tothz

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       v%vec(:,i,:,:) = v%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, v%nfreq
       v%vec(:,:,i,1) = v%vec(:,:,i,1) / &
          sqrt(sum(v%vec(:,:,i,1)*conjg(v%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, v%nqpt
    !    do ifreq = 1, v%nfreq
    !       do jfreq = 1, v%nfreq
    !          summ = sum(v%vec(:,:,ifreq,iqpt)*conjg(v%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    v%hasvibs = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  contains
    subroutine strip(str)
      character*(*), intent(inout) :: str

      integer :: i

      do i = 1, len(str)
         if (str(i:i) == ";" .or. str(i:i) == "#" .or. str(i:i) == "\") &
            str(i:i) = " "
      end do

    end subroutine strip

  end subroutine read_phonopy_ascii

  !> Read vibration data from a phonopy yaml file (MESH/WRITE_MESH or
  !> QPOINTS, plus EIGENVECTORS keywords), and return it in vib. If
  !> error, return non-zero errmsg. This implements a custom-made
  !> parser for the yaml file, to avoid having to require libyaml or
  !> incorporating a proper YAML reader into the code. If more yaml
  !> files come along, this will change.
  subroutine read_phonopy_yaml(v,c,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw, string
    use param, only: ivformat_phonopy_yaml, cm1tothz
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok, readvec
    character(len=:), allocatable :: line
    integer :: lu, idx, i, j, ifreq
    real*8 :: xdum(2), norm
    ! integer :: jfreq, iqpt ! checking normalization
    ! complex*16 :: summ

    real*8, parameter :: eps_error = 1d-100

    ! initialize
    errmsg = "Error reading yaml file: " // trim(file)
    lu = -1

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_phonopy_yaml

    ! allocate and initialize
    v%nqpt = 0
    v%nfreq = 3 * c%ncel
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,10))
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(v%nfreq,10))
    if (allocated(v%vec)) deallocate(v%vec)
    allocate(v%vec(3,c%ncel,v%nfreq,10))

    ! advance to the header
    do while (getline_raw(lu,line,.false.))
       if (index(line,"phonon:") > 0) exit
    end do

    !!! Reader implemented to fit this YAML structure:
    ! phonon:
    ! - q-position: [    0.2500000,    0.0000000,    0.0000000 ]
    !   distance_from_gamma:  0.021997038
    !   weight: 2
    !   band:
    !   - # 1
    !     frequency:     0.3378353291
    !     eigenvector:
    !     - # atom 1
    !       - [  0.00451257038677,  0.00000000000000 ]
    !       - [ -0.11674136018427,  0.04305724458615 ]
    readvec = .false.
    do while (getline_raw(lu,line,.false.))
       if (index(line,"q-position") > 0) then
          ! a new q-point, increase nqpt and reallocate
          v%nqpt = v%nqpt + 1
          ifreq = 0
          if (v%nqpt > size(v%qpt,2)) then
             call realloc(v%qpt,3,2*v%nqpt)
             call realloc(v%freq,v%nfreq,2*v%nqpt)
             call realloc(v%vec,3,c%ncel,v%nfreq,2*v%nqpt)
          end if

          ! a new q-point
          idx = index(line,":")
          line = line(idx+1:)
          call strip(line)
          read (line,*,err=999,end=999) v%qpt(:,v%nqpt)

       elseif (index(line,"frequency") > 0) then
          ! a new frequency (band)
          ifreq = ifreq + 1
          if (ifreq > v%nfreq) then
             errmsg = "Inconsistent number of frequencies"
             goto 999
          end if
          idx = index(line,":")
          line = line(idx+1:)
          read (line,*,err=999,end=999) v%freq(ifreq,v%nqpt)

       elseif (index(line,"eigenvector") > 0) then
          readvec = .true.
          ! an eigenvector
          do i = 1, c%ncel
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) then
                errmsg = "error reading eigenvector " // string(i)
                goto 999
             end if
             do j = 1, 3
                ok = getline_raw(lu,line,.false.)
                if (.not.ok) then
                   errmsg = "error reading eigenvector " // string(i) // ", coord " // string(j)
                   goto 999
                end if
                idx = index(line,"[")
                line = line(idx+1:)
                call strip(line)
                read (line,*,err=999,end=999) xdum
                v%vec(j,i,ifreq,v%nqpt) = cmplx(xdum(1),xdum(2),8)
             end do
          end do
       end if
    end do
    call realloc(v%qpt,3,v%nqpt)
    call realloc(v%freq,v%nfreq,v%nqpt)
    call realloc(v%vec,3,c%ncel,v%nfreq,v%nqpt)
    if (.not.readvec) then
       errmsg = "the yaml file contains no eigenvector information"
       goto 999
    end if

    ! THz to cm-1
    v%freq = v%freq / cm1tothz

    ! normalize
    do i = 1, v%nfreq
       norm = sqrt(real(sum(v%vec(:,:,i,1)*conjg(v%vec(:,:,i,1))),8))
       if (norm < eps_error) then
          errmsg = "zero-length eigenvector found (the yaml file has eigenvector?)"
          goto 999
       end if
       v%vec(:,:,i,1) = v%vec(:,:,i,1) / norm
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, v%nqpt
    !    do ifreq = 1, v%nfreq
    !       do jfreq = 1, v%nfreq
    !          summ = sum(v%vec(:,:,ifreq,iqpt)*conjg(v%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    v%hasvibs = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  contains
    subroutine strip(str)
      character*(*), intent(inout) :: str

      integer :: i

      do i = 1, len(str)
         if (str(i:i) == "[" .or. str(i:i) == "," .or. str(i:i) == "]") &
            str(i:i) = " "
      end do

    end subroutine strip

  end subroutine read_phonopy_yaml

  !> Read vibration data from a phonopy hdf5 file (MESH/WRITE_MESH or
  !> QPOINTS, plus EIGENVECTORS keywords), and return it in v. If
  !> error, return non-zero errmsg.
  subroutine read_phonopy_hdf5(v,c,file,errmsg)
#ifdef HAVE_HDF5
    use iso_c_binding, only: c_loc, c_ptr
    use hdf5
#endif
    use param, only: ivformat_phonopy_hdf5, cm1tothz
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg

#ifdef HAVE_HDF5
    integer(HID_T) :: fid, dsetid, spaceid, ctypeid
    integer(HSIZE_T) :: dim(2), maxdim(2)
    integer :: i, ier, ifreq, iqpt
    complex*16, allocatable, target :: vaux(:,:,:,:)
    type(c_ptr) :: cptr
    ! integer :: jfreq ! checking normalization
    ! complex*16 :: summ

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_phonopy_hdf5

    ! open file
    call h5fopen_f(file,H5F_ACC_RDONLY_F, fid, ier)
    if (ier < 0) then
       errmsg = "error opening hdf5 file: " // trim(file)
       return
    end if

    ! read qpoints
    call h5dopen_f(fid,"qpoint",dsetid,ier)
    if (ier < 0) goto 999
    call h5dget_space_f(dsetid,spaceid,ier)
    if (ier < 0) goto 999
    CALL h5sget_simple_extent_dims_f(spaceid,dim,maxdim,ier)
    if (ier < 0) goto 999
    v%nqpt = int(dim(2))
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(dim(1),dim(2)))
    call h5dread_f(dsetid,H5T_NATIVE_DOUBLE,v%qpt,dim,ier)
    if (ier < 0) goto 999

    ! read frequencies
    call h5dopen_f(fid,"frequency",dsetid,ier)
    if (ier < 0) goto 999
    call h5dget_space_f(dsetid,spaceid,ier)
    if (ier < 0) goto 999
    CALL h5sget_simple_extent_dims_f(spaceid,dim,maxdim,ier)
    if (ier < 0) goto 999
    v%nfreq = int(dim(1))
    if (3*c%ncel /= v%nfreq) then
       errmsg = "Inconsistent number of frequencies: " // trim(file)
       goto 999
    end if
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(dim(1),dim(2)))
    call h5dread_f(dsetid,H5T_NATIVE_DOUBLE,v%freq,dim,ier)
    if (ier < 0) goto 999

    ! THz to cm-1
    v%freq = v%freq / cm1tothz

    ! read eigenvectors
    call h5dopen_f(fid,"eigenvector",dsetid,ier)
    if (ier < 0) goto 999

    ! compound type for eigenvectors
    CALL H5Tcreate_f(H5T_COMPOUND_F,16_8,ctypeid,ier)
    if (ier < 0) goto 999
    CALL H5Tinsert_f(ctypeid,"r",0_8,H5T_NATIVE_DOUBLE,ier)
    if (ier < 0) goto 999
    CALL H5Tinsert_f(ctypeid,"i",8_8,H5T_NATIVE_DOUBLE,ier)
    if (ier < 0) goto 999

    ! rearrange data = it is written in a weird, weird order
    allocate(vaux(v%nfreq,3,c%ncel,v%nqpt))
    cptr = c_loc(vaux)
    call h5dread_f(dsetid,ctypeid,cptr,ier)
    if (ier < 0) goto 999
    allocate(v%vec(3,c%ncel,v%nfreq,v%nqpt))
    do iqpt = 1, v%nqpt
       do ifreq = 1, v%nfreq
          do i = 1, c%ncel
             v%vec(:,i,ifreq,iqpt) = vaux(ifreq,:,i,iqpt)
          end do
       end do
    end do
    deallocate(vaux)

    ! wrap up
    v%hasvibs = .true.
    errmsg = ""
999 continue
    call h5close_f(ier)

#else
    errmsg = "Critic2 must be compiled against the HDF5 library to read phonopy hdf5 files"
#endif

  end subroutine read_phonopy_hdf5

  !> Read the force constants from a phonopy FORCE_CONSTANTS file and
  !> populate v. Use the structure information in the crystal
  !> structure c. If error, return non-zero errmsg.
  !> Read the 2nd-order force constants from a phonopy FORCE_CONSTANTS
  !> file (text format) and return them in v, for the crystal structure
  !> c, which must be the unit cell phonopy was run with. sline carries
  !> the options from the VIBRATIONS LOAD line: a mandatory generator
  !> keyword (the code that wrote the file, which fixes the units), an
  !> optional supercell specification (three integers for a diagonal
  !> supercell, nine integers for a general one in phonopy DIM order, or
  !> the name of a file containing the supercell; default is no
  !> supercell), and the optional keyword ASR to enforce the acoustic
  !> sum rule after reading. If error, return non-zero errmsg.
  !>
  !> The file carries no atomic positions, so the supercell atom order
  !> is reconstructed following phonopy (structure/cells.py): atom-major
  !> (all lattice images of cell atom 1, then those of cell atom 2, and
  !> so on) with the lattice points enumerated over the surrounding
  !> frame of the supercell matrix, first axis fastest, keeping the
  !> first occurrence of each class modulo the supercell lattice. A
  !> compact file (fewer rows than supercell atoms) gives the force
  !> constants only for the atoms of phonopy's primitive cell; the rows
  !> for the remaining atoms of the cell are generated with the pure
  !> translations of the crystal, which is all that is needed (a
  !> rotation would not commute with the supercell lattice).
  subroutine read_phonopy_fc2(v,c,file,sline,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword, getword, equal,&
       isinteger, string
    use tools_math, only: matinv
    use param, only: hartoev, bohrtoa
    type(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file, sline
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: word, line, scfile
    integer :: lu, lp, lp0, i, j, ia, il, ip, jp, is, js, iap, ilp, irow, idum
    integer :: igen, ndim, nlat, nsat, nrow, n1, n2, ierr
    integer :: smat(3,3), madj(3,3), idim(9), lv(3)
    real*8 :: fc2factor, rmat(3,3), gsc(3,3), guc(3,3), gchk(3,3), dev
    real*8 :: t(3), xsc(3), dx(3), rminv(3,3)
    logical :: iscompact, haveseed, doasr, found
    integer, allocatable :: lvec(:,:), lkey(:,:), p2s(:), s2row(:), perm(:)
    logical, allocatable :: coldone(:)
    real*8, allocatable :: fcrow(:,:,:,:)
    type(crystal) :: sc
    type(crystalseed) :: seed

    ! initialize
    errmsg = "Error reading FORCE_CONSTANTS file: " // trim(file)
    lu = -1
    scfile = ""
    if (c%ismolecule) then
       errmsg = "FORCE_CONSTANTS (phonopy) files can only be read for crystals"
       goto 999
    end if

    ! generator keyword: mandatory, gives the units of the file
    lp = 1
    word = lgetword(sline,lp)
    igen = fc2_generator_index(sline)
    if (igen == 0) then
       errmsg = "A generator keyword (vasp, qe, aims, alamode, ...) is required to read FORCE_CONSTANTS"
       goto 999
    end if
    select case (fc2_genunit(igen))
    case (1)
       fc2factor = bohrtoa**2 / hartoev ! eV/ang^2
    case (2)
       fc2factor = bohrtoa / hartoev ! eV/(ang*bohr)
    case (3)
       fc2factor = 0.5d0 ! Ry/bohr^2
    case (4)
       fc2factor = 0.5d-3 ! mRy/bohr^2
    case (5)
       fc2factor = 1d0 ! Hartree/bohr^2
    case (6)
       fc2factor = bohrtoa ! Hartree/(ang*bohr)
    end select

    ! supercell specification and options
    ndim = 0
    haveseed = .false.
    doasr = .false.
    do while (.true.)
       lp0 = lp
       if (isinteger(idum,sline,lp)) then
          ndim = ndim + 1
          if (ndim > 9) then
             errmsg = "Too many integers in the supercell specification (VIBRATIONS LOAD)"
             goto 999
          end if
          idim(ndim) = idum
          cycle
       end if
       word = lgetword(sline,lp)
       if (len_trim(word) == 0) exit
       if (equal(word,"asr").or.equal(word,"acoustic").or.equal(word,"acoustic_sum_rules")) then
          doasr = .true.
       elseif (haveseed) then
          errmsg = "Unknown keyword in VIBRATIONS LOAD (FORCE_CONSTANTS): " // trim(word)
          goto 999
       else
          ! the supercell, read from a structure file
          lp = lp0
          scfile = getword(sline,lp)
          call seed%read_any_file(scfile,-1,errmsg,ti=ti)
          if (len_trim(errmsg) > 0) goto 999
          call sc%struct_new(seed,errmsg,noenv=.true.,ti=ti)
          if (len_trim(errmsg) > 0) goto 999
          if (sc%ismolecule) then
             errmsg = "The supercell file in VIBRATIONS LOAD is a molecule"
             goto 999
          end if
          haveseed = .true.
       end if
    end do

    ! build the supercell matrix (rows = supercell lattice vectors, in
    ! units of the cell vectors; this is phonopy's convention)
    smat = 0
    do i = 1, 3
       smat(i,i) = 1
    end do
    if (haveseed) then
       if (ndim > 0) then
          errmsg = "Give either the supercell matrix or a supercell file in VIBRATIONS LOAD, not both"
          goto 999
       end if

       ! columns of m_c2x*sc%m_x2c are the supercell vectors in cell coordinates
       rmat = transpose(matmul(c%m_c2x,sc%m_x2c))
       smat = nint(rmat)
       dev = maxval(abs(rmat - real(smat,8)))
       if (dev > fc2_epsint) then
          errmsg = "The cell in " // trim(scfile) // " is not an integer supercell of the current &
             &structure (deviation " // string(dev,'e',12,4) // ")"
          goto 999
       end if

       ! the product above assumes both cells use the same Cartesian
       ! frame; the metric identity G_sc = smat*G_uc*smat^t does not
       guc = matmul(transpose(c%m_x2c),c%m_x2c)
       gsc = matmul(transpose(sc%m_x2c),sc%m_x2c)
       gchk = matmul(real(smat,8),matmul(guc,transpose(real(smat,8))))
       dev = maxval(abs(gchk - gsc)) / max(maxval(abs(gsc)),1d-10)
       if (dev > fc2_epsmetric) then
          errmsg = "The supercell matrix deduced from " // trim(scfile) // " is inconsistent with the &
             &cell metric (relative deviation " // string(dev,'e',12,4) // "); the two structures may &
             &be given in different Cartesian frames"
          goto 999
       end if
    elseif (ndim == 3) then
       smat = 0
       do i = 1, 3
          smat(i,i) = idim(i)
       end do
    elseif (ndim == 9) then
       ! Same order as phonopy's DIM keyword. Note that phonopy reshapes
       ! those nine integers row-wise into its supercell matrix M, but
       ! then builds the supercell with lattice vectors given by the rows
       ! of transpose(M) (structure/cells.py, _trim_cell with the old
       ! style, which is the default). smat is the latter, so the nine
       ! integers fill it column-wise.
       smat = reshape(idim,(/3,3/))
    elseif (ndim /= 0) then
       errmsg = "The supercell specification needs 3 or 9 integers (found " // string(ndim) // ")"
       goto 999
    end if

    ! lattice points of the supercell, in phonopy order
    call fc2_lattice_points(smat,nlat,madj,lvec,lkey,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    nsat = c%ncel * nlat

    ! open the file and read the header
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if
    errmsg = "Error reading the header of FORCE_CONSTANTS file: " // trim(file)
    if (.not.getline_raw(lu,line,.false.)) goto 999
    lp = 1
    if (.not.isinteger(n1,line,lp)) goto 999
    if (.not.isinteger(n2,line,lp)) n2 = n1

    if (n1 == n2 .and. n1 == nsat) then
       iscompact = .false.
       nrow = nsat
    elseif (n1 < n2 .and. n2 == nsat) then
       iscompact = .true.
       nrow = n1
    else
       errmsg = "Inconsistent number of atoms in " // trim(file) // ": the file has (" // string(n1) //&
          "," // string(n2) // ") but the structure and the supercell give " // string(c%ncel) // "*" //&
          string(nlat) // " = " // string(nsat) // " supercell atoms"
       goto 999
    end if

    ! read the force constants
    allocate(fcrow(3,3,nrow,nsat),p2s(nrow),s2row(nsat),coldone(nsat),stat=ierr)
    if (ierr /= 0) then
       errmsg = "Could not allocate " // string(nint(72d0*nrow*nsat/1024d0**2)) //&
          " Mb to read the force constants"
       goto 999
    end if
    s2row = 0
    errmsg = "Error reading the force constants from file: " // trim(file)
    do irow = 1, nrow
       coldone = .false.
       do jp = 1, nsat
          read (lu,*,err=999,end=999) i, j
          if (i < 1 .or. i > nsat .or. j < 1 .or. j > nsat) then
             errmsg = "Atom index out of range in FORCE_CONSTANTS block " // string(irow) // "," // string(jp)
             goto 999
          end if
          if (jp == 1) then
             p2s(irow) = i
             if (s2row(i) /= 0) then
                errmsg = "Repeated row atom index " // string(i) // " in FORCE_CONSTANTS"
                goto 999
             end if
             s2row(i) = irow
          elseif (i /= p2s(irow)) then
             errmsg = "Inconsistent row atom index in FORCE_CONSTANTS block " // string(irow) // "," // string(jp)
             goto 999
          end if
          if (coldone(j)) then
             errmsg = "Repeated column atom index " // string(j) // " in FORCE_CONSTANTS row " // string(irow)
             goto 999
          end if
          coldone(j) = .true.
          read (lu,*,err=999,end=999) rmat(1,:)
          read (lu,*,err=999,end=999) rmat(2,:)
          read (lu,*,err=999,end=999) rmat(3,:)
          fcrow(:,:,irow,j) = rmat
       end do
       if (.not.all(coldone)) then
          errmsg = "Incomplete set of columns in FORCE_CONSTANTS row " // string(irow)
          goto 999
       end if
       if (.not.iscompact .and. p2s(irow) /= irow) then
          errmsg = "Full FORCE_CONSTANTS file with a non-identity row map, at row " // string(irow)
          goto 999
       end if
    end do
    call fclose(lu)
    lu = -1
    fcrow = fcrow * fc2factor

    ! build the cell-compact array: rows for all the atoms of the cell
    if (allocated(v%fc2)) deallocate(v%fc2)
    allocate(v%fc2(3,3,c%ncel,nsat),perm(nsat),stat=ierr)
    if (ierr /= 0) then
       errmsg = "Could not allocate " // string(nint(72d0*c%ncel*nsat/1024d0**2)) //&
          " Mb for the force constants"
       goto 999
    end if
    v%fc2 = 0d0
    do ia = 1, c%ncel
       ! the source is the image of cell atom ia at the origin
       is = (ia-1)*nlat + 1
       irow = s2row(is)
       if (irow > 0) then
          v%fc2(:,:,ia,:) = fcrow(:,:,irow,:)
          cycle
       end if

       ! not in the file: find a pure translation carrying it to an atom that is
       found = .false.
       do ip = 1, nrow
          js = p2s(ip)
          iap = (js-1)/nlat + 1
          ilp = js - (iap-1)*nlat
          if (c%spc(c%atcel(iap)%is)%z /= c%spc(c%atcel(ia)%is)%z) cycle
          t = c%atcel(iap)%x + real(lvec(:,ilp),8) - c%atcel(ia)%x
          if (.not.fc2_pure_translation(c,t,nlat,lvec,lkey,madj,perm,dev)) cycle
          if (perm(is) /= js) cycle
          found = .true.
          exit
       end do
       if (.not.found) then
          errmsg = "No pure translation carries atom " // string(ia) // " onto any of the atoms in the &
             &compact FORCE_CONSTANTS file. The current structure may not be the unit cell used by &
             &phonopy, or phonopy used incompatible primitive axes; rerun it with &
             &FULL_FORCE_CONSTANTS = .TRUE."
          goto 999
       end if
       do js = 1, nsat
          v%fc2(:,:,ia,js) = fcrow(:,:,ip,perm(js))
       end do
    end do

    ! if a supercell file was given, check the atom ordering against it.
    ! This validates the lattice-point table, the atom-major convention
    ! and the atom order of the current structure, all at once.
    if (haveseed) then
       if (seed%nat /= nsat) then
          errmsg = "The supercell file " // trim(scfile) // " has " // string(seed%nat) //&
             " atoms, but " // string(nsat) // " were expected"
          goto 999
       end if
       rminv = real(smat,8)
       call matinv(rminv,3,ierr)
       if (ierr /= 0) then
          errmsg = "Singular supercell matrix"
          goto 999
       end if
       do ia = 1, c%ncel
          do il = 1, nlat
             is = (ia-1)*nlat + il
             xsc = matmul(c%atcel(ia)%x + real(lvec(:,il),8),rminv)
             if (seed%spc(seed%is(is))%z /= c%spc(c%atcel(ia)%is)%z) then
                errmsg = "Species mismatch at atom " // string(is) // " of the supercell file " //&
                   trim(scfile) // ": its atom ordering is not phonopy's, or the current structure is &
                   &not the phonopy unit cell"
                goto 999
             end if
             dx = xsc - seed%x(:,is)
             dx = dx - nint(dx)
             dev = norm2(matmul(sc%m_x2c,dx))
             if (dev > fc2_epspos) then
                errmsg = "Atom " // string(is) // " of the supercell file " // trim(scfile) //&
                   " is " // string(dev,'e',12,4) // " bohr away from its expected phonopy position; &
                   &check the atom order of the current structure"
                goto 999
             end if
          end do
       end do
    end if

    ! wrap up
    v%fc2_file = file
    v%hasfc2 = .true.
    v%fc2_smat = smat
    v%fc2_madj = madj
    v%fc2_nlat = nlat
    v%fc2_nsat = nsat
    v%fc2_ncel = c%ncel
    v%fc2_iscompact = iscompact
    if (allocated(v%fc2_lvec)) deallocate(v%fc2_lvec)
    if (allocated(v%fc2_lkey)) deallocate(v%fc2_lkey)
    call move_alloc(lvec,v%fc2_lvec)
    call move_alloc(lkey,v%fc2_lkey)
    v%fc2_acoustic = -1
    v%fc2_vs_center = -1d0
    v%fc2_vs_delta = -1d0

    ! apply the acoustic sum rule, if requested
    if (doasr) &
       call v%apply_acoustic(c,.false.)

    errmsg = ""
    return

999 continue
    if (lu > 0) call fclose(lu)
    call v%end(keepvibs=.true.)

  end subroutine read_phonopy_fc2

  !> Lattice points of the supercell defined by the integer supercell
  !> matrix smat (rows = supercell lattice vectors in units of the cell
  !> vectors), in phonopy's order (structure/cells.py): the surrounding
  !> frame of the matrix is swept with the first axis fastest and the
  !> first occurrence of each class modulo the supercell lattice is
  !> kept. Returns the number of lattice points (nlat = det(smat)), the
  !> integer adjugate of smat (madj, with smat*madj = nlat*I), the
  !> lattice points in cell fractional coordinates (lvec) and their
  !> residue keys modulo the supercell lattice (lkey). Two lattice
  !> points are equivalent modulo the supercell lattice if and only if
  !> they have the same key, and the key is exact integer arithmetic.
  subroutine fc2_lattice_points(smat,nlat,madj,lvec,lkey,errmsg)
    use tools_io, only: string
    integer, intent(in) :: smat(3,3)
    integer, intent(out) :: nlat
    integer, intent(out) :: madj(3,3)
    integer, allocatable, intent(inout) :: lvec(:,:), lkey(:,:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, k, ia, ib, ic, nl
    integer :: cmin(3), cmax(3), frame(3), vv(3), ll(3), kk(3)
    logical :: found

    errmsg = ""

    ! determinant and adjugate, in exact integer arithmetic
    nlat = smat(1,1)*(smat(2,2)*smat(3,3)-smat(2,3)*smat(3,2))&
       - smat(1,2)*(smat(2,1)*smat(3,3)-smat(2,3)*smat(3,1))&
       + smat(1,3)*(smat(2,1)*smat(3,2)-smat(2,2)*smat(3,1))
    if (nlat <= 0) then
       errmsg = "The supercell matrix must have positive determinant (found " // string(nlat) // ")"
       return
    end if
    madj(1,1) = smat(2,2)*smat(3,3) - smat(2,3)*smat(3,2)
    madj(1,2) = -smat(1,2)*smat(3,3) + smat(1,3)*smat(3,2)
    madj(1,3) = smat(1,2)*smat(2,3) - smat(1,3)*smat(2,2)
    madj(2,1) = -smat(2,1)*smat(3,3) + smat(2,3)*smat(3,1)
    madj(2,2) = smat(1,1)*smat(3,3) - smat(1,3)*smat(3,1)
    madj(2,3) = -smat(1,1)*smat(2,3) + smat(1,3)*smat(2,1)
    madj(3,1) = smat(2,1)*smat(3,2) - smat(2,2)*smat(3,1)
    madj(3,2) = -smat(1,1)*smat(3,2) + smat(1,2)*smat(3,1)
    madj(3,3) = smat(1,1)*smat(2,2) - smat(1,2)*smat(2,1)

    ! the surrounding frame: bounding box of the eight corners of the supercell
    cmin = 0
    cmax = 0
    do k = 0, 7
       vv = 0
       if (btest(k,0)) vv = vv + smat(1,:)
       if (btest(k,1)) vv = vv + smat(2,:)
       if (btest(k,2)) vv = vv + smat(3,:)
       cmin = min(cmin,vv)
       cmax = max(cmax,vv)
    end do
    frame = cmax - cmin

    ! sweep the frame with the first axis fastest, keeping first occurrences
    if (allocated(lvec)) deallocate(lvec)
    if (allocated(lkey)) deallocate(lkey)
    allocate(lvec(3,nlat),lkey(3,nlat))
    nl = 0
    loopc: do ic = 0, frame(3)-1
       do ib = 0, frame(2)-1
          do ia = 0, frame(1)-1
             ll = (/ia,ib,ic/)
             kk = modulo(matmul(ll,madj),nlat)
             found = .false.
             do i = 1, nl
                if (all(kk == lkey(:,i))) then
                   found = .true.
                   exit
                end if
             end do
             if (found) cycle
             nl = nl + 1
             lvec(:,nl) = ll
             lkey(:,nl) = kk
             if (nl == nlat) exit loopc
          end do
       end do
    end do loopc
    if (nl /= nlat) then
       errmsg = "Could only generate " // string(nl) // " of the " // string(nlat) // " supercell &
          &lattice points; the supercell matrix may be invalid"
       return
    end if

  end subroutine fc2_lattice_points

  !> Index of the generator (the code that wrote the file, which fixes
  !> the units of the force constants) named by the first word of sline.
  !> Zero if the first word is not a known generator.
  function fc2_generator_index(sline) result(igen)
    use tools_io, only: lgetword, equal
    character*(*), intent(in) :: sline
    integer :: igen

    integer :: i, lp
    character(len=:), allocatable :: word

    lp = 1
    word = lgetword(sline,lp)
    igen = 0
    do i = 1, fc2_ngen
       if (equal(word,trim(fc2_genname(i)))) then
          igen = i
          return
       end if
    end do

  end function fc2_generator_index

  !> Index in the lattice-point table of the point equivalent to lv
  !> (integer vector in cell fractional coordinates) modulo the
  !> supercell lattice. Zero if not found, which cannot happen for a
  !> complete table.
  function fc2_ilat(nlat,lkey,madj,lv) result(il)
    integer, intent(in) :: nlat
    integer, intent(in) :: lkey(:,:)
    integer, intent(in) :: madj(3,3)
    integer, intent(in) :: lv(3)
    integer :: il

    integer :: i, kk(3)

    kk = modulo(matmul(lv,madj),nlat)
    do i = 1, nlat
       if (all(kk == lkey(:,i))) then
          il = i
          return
       end if
    end do
    il = 0

  end function fc2_ilat

  !> Is the vector t (cell fractional coordinates) a pure translation of
  !> the crystal c? Checking the atoms of the cell is enough: if t maps
  !> the cell onto itself then the induced map of the supercell atoms,
  !> (ia,L) -> (ja,L+n), is automatically a permutation, because adding
  !> a fixed integer vector commutes with the reduction modulo the
  !> supercell lattice. If t is a translation, return that permutation
  !> of the supercell atoms (perm) and the largest position mismatch
  !> found (dev, bohr).
  function fc2_pure_translation(c,t,nlat,lvec,lkey,madj,perm,dev) result(ok)
    use param, only: icrd_crys
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: t(3)
    integer, intent(in) :: nlat
    integer, intent(in) :: lvec(:,:), lkey(:,:)
    integer, intent(in) :: madj(3,3)
    integer, intent(out) :: perm(:)
    real*8, intent(out) :: dev
    logical :: ok

    integer :: ia, ja, il, jl, is
    integer :: lv(3), at2(c%ncel), lv2(3,c%ncel)
    real*8 :: dd

    ok = .false.
    dev = 0d0
    do ia = 1, c%ncel
       ja = c%identify_atom(c%atcel(ia)%x + t,icrd_crys,lvec=lv,dist=dd,distmax=fc2_epspos)
       if (ja == 0) return
       if (dd > fc2_epspos) return
       if (c%spc(c%atcel(ja)%is)%z /= c%spc(c%atcel(ia)%is)%z) return
       at2(ia) = ja
       lv2(:,ia) = lv
       dev = max(dev,dd)
    end do

    ! the permutation of the supercell atoms induced by t
    do ia = 1, c%ncel
       ja = at2(ia)
       do il = 1, nlat
          is = (ia-1)*nlat + il
          jl = fc2_ilat(nlat,lkey,madj,lv2(:,ia) + lvec(:,il))
          if (jl == 0) return
          perm(is) = (ja-1)*nlat + jl
       end do
    end do
    ok = .true.

  end function fc2_pure_translation

  !> Read vibration data from a crystal output file, and return it in
  !> vib. If error, return non-zero errmsg.
  subroutine read_crystal_out(v,c,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw
    use param, only: ivformat_crystal_out, atmass
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok
    character(len=:), allocatable :: line
    integer :: lu, i, j, k, n0, n1, iz
    real*8 :: xdum(6)
    ! integer :: ifreq, jfreq, iqpt ! checking normalization
    ! complex*16 :: summ

    ! initialize
    errmsg = "Error reading output file: " // trim(file)
    lu = -1

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_crystal_out

    ! allocate and initialize -> only gamma point
    v%nqpt = 1
    v%nfreq = 3 * c%ncel
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,1))
    v%qpt = 0d0
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(v%nfreq,1))
    if (allocated(v%vec)) deallocate(v%vec)
    allocate(v%vec(3,c%ncel,v%nfreq,1))

    ! advance to the header
    do while (getline_raw(lu,line,.false.))
       if (index(line,"NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES (IN BOHR)") > 0) exit
    end do
    ok = getline_raw(lu,line,.false.)
    if (.not.ok) goto 999

    ! read the info
    do i = 1, v%nfreq, 6
       n0 = i
       n1 = min(i+5,v%nfreq)

       ! read the frequency line
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999
       line = line(14:)
       read (line,*,end=999,err=999) v%freq(n0:n1,1)
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999

       ! read the eigenvector
       do j = 1, c%ncel
          do k = 1, 3
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) goto 999
             line = line(14:)
             read (line,*,end=999,err=999) xdum(1:n1-n0+1)
             v%vec(k,j,n0:n1,1) = xdum(1:n1-n0+1)
          end do
       end do
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999
    end do

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       v%vec(:,i,:,:) = v%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, v%nfreq
       v%vec(:,:,i,1) = v%vec(:,:,i,1) / &
          sqrt(sum(v%vec(:,:,i,1)*conjg(v%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, v%nqpt
    !    do ifreq = 1, v%nfreq
    !       do jfreq = 1, v%nfreq
    !          summ = sum(v%vec(:,:,ifreq,iqpt)*conjg(v%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do
    ! stop 1

    ! wrap up
    v%hasvibs = .true.
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_crystal_out

  !> Read vibration data from a Gaussian output file, and return it in
  !> v. If error, return non-zero errmsg.
  subroutine read_gaussian_log(v,c,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw, isreal
    use tools_math, only: rmsd_walker
    use param, only: ivformat_gaussian_log, atmass, bohrtoa, img
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok, okl(3), have_std
    character(len=:), allocatable :: line
    integer :: lu, i, j, nread, iz, nn, lp, idum(2), idum3(3)
    real*8 :: xdum, xread(9), rot(3,3), xstd(3,c%ncel), xinp(3,c%ncel)
    ! integer :: jfreq ! checking normalization
    ! complex*16 :: summ

    ! initialize
    errmsg = "Error reading vibrations from file: " // trim(file)
    lu = -1

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_gaussian_log

    ! allocate and initialize -> only gamma point
    v%nqpt = 1
    v%nfreq = 3 * c%ncel
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,1))
    v%qpt = 0d0
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(v%nfreq,1))
    if (allocated(v%vec)) deallocate(v%vec)
    allocate(v%vec(3,c%ncel,v%nfreq,1))

    ! The normal-mode displacements below are printed in Gaussian's standard
    ! orientation, but the stored geometry is the input orientation. Capture the
    ! last standard-orientation block so we can rotate the modes into the input
    ! frame afterwards (see below).
    have_std = .false.
    nread = 0
    do while (.true.)
       ! advance to the next frequencies header
       ok = .false.
       do while (getline_raw(lu,line,.false.))
          if (index(line," Harmonic frequencies (cm**-1)") > 0) nread = 0
          if (index(line,"Standard orientation:") > 0) then
             ! skip the 4 header/separator lines, then read the atomic coordinates
             do i = 1, 4
                if (.not.getline_raw(lu,line,.false.)) goto 999
             end do
             do i = 1, c%ncel
                if (.not.getline_raw(lu,line,.false.)) goto 999
                read (line,*,end=999,err=999) idum3, xstd(:,i)
             end do
             xstd = xstd / bohrtoa
             have_std = .true.
          end if
          if (index(line," Frequencies --") > 0) then
             ok = .true.
             exit
          end if
       end do
       if (.not.ok) exit

       ! read the frequencies
       line = line(16:)
       lp = 1
       do i = 1, 3
          okl(i) = isreal(xdum,line,lp)
          if (okl(i)) v%freq(nread+i,1) = xdum
       end do
       nn = count(okl)

       ! skip to the displacements
       do while(.true.)
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          if (index(line," Atom ") > 0) exit
       end do

       ! read the frequencies
       do i = 1, c%ncel
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          read (line,*,end=999,err=999) idum, xread(1:3*nn)
          do j = 1, nn
             v%vec(:,i,nread+j,1) = xread(3*(j-1)+1:3*j)
          end do
       end do
       nread = nread + nn
    end do
    v%nfreq = nread
    call realloc(v%freq,v%nfreq,v%nqpt)
    call realloc(v%vec,3,c%ncel,v%nfreq,v%nqpt)

    ! The mode vectors are in the standard orientation but the
    ! geometry is in the input orientation. Rotate the modes into the
    ! input frame using the rotation that superimposes the
    ! standard-orientation geometry onto the stored one.  Skipped when
    ! there is no standard orientation block (e.g. nosymm), in which
    ! case the modes are already in the input frame.
    if (have_std) then
       do i = 1, c%ncel
          xinp(:,i) = c%atcel(i)%r
       end do
       xdum = rmsd_walker(xstd,xinp,rot)
       do i = 1, c%ncel
          do j = 1, v%nfreq
             v%vec(:,i,j,1) = matmul(rot,real(v%vec(:,i,j,1),8)) + &
                img * matmul(rot,aimag(v%vec(:,i,j,1)))
          end do
       end do
    end if

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       v%vec(:,i,:,:) = v%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, v%nfreq
       v%vec(:,:,i,1) = v%vec(:,:,i,1) / &
          sqrt(sum(v%vec(:,:,i,1)*conjg(v%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do ifreq = 1, v%nfreq
    !    do jfreq = 1, v%nfreq
    !       summ = sum(v%vec(:,:,ifreq,1)*conjg(v%vec(:,:,jfreq,1)))
    !       write (*,*) ifreq, jfreq, summ
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    v%hasvibs = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_gaussian_log

  !> Read vibration data from a Gaussian formatted checkpoint file,
  !> and return it in v. If error, return non-zero errmsg.
  subroutine read_gaussian_fchk(v,c,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw, isreal
    use param, only: ivformat_gaussian_fchk, atmass
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok
    character(len=:), allocatable :: line
    integer :: lu, i, j, k, nread, iz, nn, ll, n0, n1, nnm
    real*8 :: xdum(5)
    character*1 :: ch
    ! integer :: ifreq, jfreq ! checking normalization
    ! complex*16 :: summ

    ! initialize
    errmsg = "Error reading vibrations from file: " // trim(file)
    lu = -1

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_gaussian_fchk

    ! allocate and initialize -> only gamma point
    v%nqpt = 1
    v%nfreq = 0
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,1))
    v%qpt = 0d0
    if (allocated(v%freq)) deallocate(v%freq)
    if (allocated(v%vec)) deallocate(v%vec)

    do while (getline_raw(lu,line,.false.))
       ll = len(line)

       ! read the number of frequencies
       if (ll >= 22) then
          if (line(1:22) == "Number of Normal Modes") then
             line = line(23:)
             read (line,*,end=999,err=999) ch, v%nfreq
             allocate(v%freq(v%nfreq,1))
             allocate(v%vec(3,c%ncel,v%nfreq,1))
          end if
       end if

       ! read frequencies
       if (ll >= 6) then
          if (line(1:6) == "Vib-E2") then
             do i = 1, v%nfreq, 5
                n0 = i
                n1 = min(i+4,v%nfreq)
                ok = getline_raw(lu,line,.false.)
                if (.not.ok) goto 999
                read (line,*,end=999,err=999) v%freq(n0:n1,1)
             end do
          end if
       end if

       ! read modes
       if (ll >= 9) then
          if (line(1:9) == "Vib-Modes") then

             nn = -1
             do i = 1, v%nfreq
                do j = 1, c%ncel
                   do k = 1, 3
                      nn = nn + 1
                      nnm = mod(nn,5) + 1
                      if (mod(nn,5) == 0) then
                         nread = min(v%nfreq * c%ncel * 3 - nn,5)
                         ok = getline_raw(lu,line,.false.)
                         if (.not.ok) goto 999
                         read (line,*,end=999,err=999) xdum(1:nread)
                      end if
                      v%vec(k,j,i,1) = xdum(nnm)
                   end do
                end do
             end do

          end if
       end if
    end do

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       v%vec(:,i,:,:) = v%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, v%nfreq
       v%vec(:,:,i,1) = v%vec(:,:,i,1) / &
          sqrt(sum(v%vec(:,:,i,1)*conjg(v%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do ifreq = 1, v%nfreq
    !    do jfreq = 1, v%nfreq
    !       summ = sum(v%vec(:,:,ifreq,1)*conjg(v%vec(:,:,jfreq,1)))
    !       write (*,*) ifreq, jfreq, summ
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    v%hasvibs = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_gaussian_fchk

  !> Read vibration data from a phonopy ascii file (ANIME keyword),
  !> and return it in vib. If error, return non-zero errmsg.
  subroutine read_castep_phonon(v,c,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw
    use types, only: realloc
    use param, only: ivformat_castep_phonon
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, idum, ifreq, iat
    character(len=:), allocatable :: line
    character*20 :: aux
    real*8 :: xx(6)
    ! integer :: jfreq, iqpt ! checking normalization
    ! complex*16 :: summ

    ! initialize
    errmsg = "Error reading ascii file: " // trim(file)
    lu = -1

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    v%file = file
    v%ivformat = ivformat_castep_phonon

    ! allocate and initialize
    v%nqpt = 0
    v%nfreq = 3 * c%ncel
    if (allocated(v%qpt)) deallocate(v%qpt)
    allocate(v%qpt(3,10))
    if (allocated(v%freq)) deallocate(v%freq)
    allocate(v%freq(v%nfreq,10))
    if (allocated(v%vec)) deallocate(v%vec)
    allocate(v%vec(3,c%ncel,v%nfreq,10))

    ! advance to the header
    do while (getline_raw(lu,line,.false.))
       if (len(line) >= 10) then
          if (line(1:10) == "     q-pt=") then
             v%nqpt = v%nqpt + 1
             if (v%nqpt > size(v%qpt,2)) then
                call realloc(v%qpt,3,2*v%nqpt)
                call realloc(v%freq,v%nfreq,2*v%nqpt)
                call realloc(v%vec,3,c%ncel,v%nfreq,2*v%nqpt)
             end if
             read (line,*,err=999,end=999) aux, idum, v%qpt(:,v%nqpt)
             do i = 1, v%nfreq
                read (lu,*,err=999,end=999) idum, v%freq(i,v%nqpt)
             end do
             read (lu,*)
             read (lu,*)
             do ifreq = 1, v%nfreq
                do iat = 1, c%ncel
                   read (lu,*,err=999,end=999) idum, idum, xx
                   v%vec(1,iat,ifreq,v%nqpt) = cmplx(xx(1),xx(2),8)
                   v%vec(2,iat,ifreq,v%nqpt) = cmplx(xx(3),xx(4),8)
                   v%vec(3,iat,ifreq,v%nqpt) = cmplx(xx(5),xx(6),8)
                end do
             end do
          end if
       end if
    end do
    if (v%nqpt == 0) then
       errmsg = "No q-points found"
       goto 999
    end if
    call realloc(v%qpt,3,v%nqpt)
    call realloc(v%freq,v%nfreq,v%nqpt)
    call realloc(v%vec,3,c%ncel,v%nfreq,v%nqpt)

    ! normalize
    do i = 1, v%nfreq
       v%vec(:,:,i,1) = v%vec(:,:,i,1) / sqrt(sum(v%vec(:,:,i,1)*conjg(v%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, v%nqpt
    !    do ifreq = 1, v%nfreq
    !       do jfreq = 1, v%nfreq
    !          summ = sum(v%vec(:,:,ifreq,iqpt)*conjg(v%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do
    ! stop 1

    ! wrap up
    errmsg = ""
    v%hasvibs = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_castep_phonon

end submodule vibrationsmod
