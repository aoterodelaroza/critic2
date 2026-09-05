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
!
! The routines that deal with the second-order force constants (reading,
! calculating, creating displacements, building the dynamical matrix, etc.)
! have been adapted from phonopy (version 2.38.2), by Atsushi Togo.
! https://github.com/phonopy/phonopy
! Phonopy is Copyright (c) 2014-2024, Phonopy. All rights reserved.
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
  character(len=*), parameter :: fc2_unitlist = "eV/ang^2, eV/(ang*bohr), Ry/bohr^2, &
     &mRy/bohr^2, Hartree/bohr^2, Hartree/(ang*bohr)"
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
  real*8, parameter :: fc2_epsdisp = 1d-4 ! displacement mismatch in a force file (bohr)
  real*8, parameter :: fc2_epscond = 1d-10 ! smallest eigenvalue ratio of the least-squares matrix
  real*8, parameter :: fc2_epssvec = 1d-4 ! same-length tolerance for the shortest supercell images (bohr)
  real*8, parameter :: fc2_epsgeo = 1d-8 ! geometry change that invalidates the force constants
  real*8, parameter :: fc2_epsdataset = 1d-6 ! structure mismatch with the displacement dataset
  real*8, parameter :: thermo_epszero = 1d-2 ! floor of the THERMO cutoff: |nu| below this is a numerical zero (cm^-1)
  real*8, parameter :: thermo_epsimag = 1d0 ! a mode below -this is imaginary, not a numerical zero (cm^-1)

  ! Displacement dataset file: everything create_forces needs to
  ! interpret the forces of the displaced supercells, written by
  ! create_displacements so that the two keywords can be run in
  ! different critic2 sessions without the user having to remember the
  ! supercell and the displacement length. The reference structure is
  ! recorded as well, so a re-relaxed or reordered cell is caught.
  type fc2_dataset
     integer :: version = 0 ! version of the dataset file format
     integer :: smat(3,3) = 0 ! supercell matrix (rows = supercell lattice vectors)
     real*8 :: dist = 0d0 ! displacement length (bohr)
     character(len=mlen) :: template = "" ! file name template of the displaced structures
     integer :: ndisp = 0 ! number of displacements
     integer :: nat = 0 ! number of atoms in the reference cell
     real*8 :: m_x2c(3,3) = 0d0 ! reference crystal-to-Cartesian matrix (bohr)
     integer, allocatable :: z(:) ! (nat) atomic numbers of the reference cell
     real*8, allocatable :: x(:,:) ! (3,nat) fractional coordinates of the reference cell
     integer, allocatable :: datom(:) ! (ndisp) displaced supercell atom
     integer, allocatable :: ddir(:,:) ! (3,ndisp) displacement direction (supercell fractional)
     character(len=mlen), allocatable :: fname(:) ! (ndisp) file written for each displacement
  end type fc2_dataset

  !xx! private procedures
  ! subroutine vibrations_detect_format(file,ivformat)
  ! subroutine read_matdyn_modes(v,c,file,ivformat,errmsg,ti)
  ! subroutine read_qe_dyn(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_ascii(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_yaml(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_hdf5(v,c,file,errmsg)
  ! subroutine read_phonopy_fc2(v,c,file,sline,errmsg,ti)
  ! function fc2_unit_index(sline,lp)
  ! function fc2_unit_factor(iunit)
  ! subroutine disp_default_template(c,template,template_,iwf)
  ! subroutine fc2_disp_setup(c,smat,dist,verbose,sc,seed,nlat,nsat,madj,lvec,lkey,nindep,indep,ndisp,datom,ddir,errmsg,ti)
  ! subroutine fc2_compatible_ops(c,smat,nlat,madj,nkeep,ikeep,rp)
  ! function fc2_dispvec(sc,ddir,dist)
  ! function fc2_expand_star(template,i,npad)
  ! subroutine fc2_atom_perm(c,io,icv,perm,ok)
  ! subroutine fc2_check_disp(file,sc,seed,iat,dexp,jmax,dmax,errmsg,ti)
  ! subroutine fc2_write_dataset(c,file,smat,dist,template,ndisp,datom,ddir,fname,errmsg,ti)
  ! subroutine fc2_read_dataset(file,ds,errmsg,ti)
  ! subroutine fc2_check_dataset(c,ds,smat,dist,ndisp,datom,ddir,errmsg)
  ! function fc2_smatstr(m)
  ! subroutine fc2_output_template(template,otemplate,errmsg)
  ! subroutine fc2_smat_from_cell(c,scfile,smat,errmsg,ti)
  ! subroutine thermo_sum(freq,nf,nq,t,cutoff,zpe,fvib,svib,cv,nint,ntot,nimag)
  ! subroutine dos_run(c,freq,nf,nq,file,nk,qshift,sigma,npts,verbose,errmsg)
  ! subroutine dos_gaussian(freq,nf,nq,sigma,fmin,step,npts,dos)
  ! subroutine dos_tetrahedra(c,freq,nf,nq,nk,fmin,step,npts,dos)
  ! function tetra_nstates(om,e)
  ! subroutine fc2_read_forces(file,nat,f,errmsg,ti)
  ! subroutine fc2_lattice_points(smat,nlat,madj,lvec,lkey,errmsg)
  ! subroutine fc2_site_directions(nsym,irot,nd,dsel)
  ! function fc2_scpos(x,lv,madj,nlat)
  ! function fc2_ilat(nlat,lkey,madj,lv)
  ! function fc2_pure_translation(c,t,nlat,lvec,lkey,madj,perm,dev)
  ! subroutine fc2_stamp_geometry(v,c)
  ! subroutine fc2_build_svec(v,c,errmsg)
  ! subroutine fc2_check_current(v,c,errmsg)
  ! subroutine read_crystal_out(v,c,file,errmsg,ti)
  ! subroutine read_gaussian_log(v,c,file,errmsg,ti)
  ! subroutine read_gaussian_fchk(v,c,file,errmsg,ti)
  ! subroutine read_castep_phonon(v,c,file,errmsg,ti)

contains

  !> Terminate a vibrations object.
  pure module subroutine vibrations_end(v,keepfc2,keepvibs)
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
       if (allocated(v%fc2_x)) deallocate(v%fc2_x)
       v%fc2_m = 0d0
       if (allocated(v%fc2_svec)) deallocate(v%fc2_svec)
       if (allocated(v%fc2_sptr)) deallocate(v%fc2_sptr)
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

    ! Detect the format from the file name. Force constants are read
    ! only through VIBRATIONS LOAD_FC2, which passes the format
    ! explicitly, because the file carries neither units nor positions
    ! and both have to be given on the keyword line.
    if (ivformat == ivformat_unknown) then
       call vibrations_detect_format(file,ivf)
       if (ivf == ivformat_unknown) then
          errmsg = "Unknown vibration file format: " // trim(file)
          return
       elseif (ivf == ivformat_phonopy_fc2) then
          errmsg = "Force constants are read with VIBRATIONS LOAD_FC2 (which also takes their units &
             &and the supercell), not with VIBRATIONS LOAD"
          return
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
       if (len_trim(v%file) > 0) then
          write (uout,'("  File: ",A)') trim(v%file)
       else
          write (uout,'("  Source: calculated from the force constants (CALCQ)")')
       end if
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
          write (uout,'("  Format: compact")')
       else
          write (uout,'("  Format: full")')
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

  !> Supercell matrix from the integers given in the input: ndim = 1
  !> (an n x n x n supercell), 3 (a diagonal supercell) or 9 (a
  !> general transformation). By default the nine integers are in the
  !> NEWCELL order: they fill, column by column, the matrix whose
  !> columns are the new lattice vectors in cell coordinates. If
  !> phonopy is true they are in the order of phonopy's DIM tag
  !> instead. Returns smat, rows = supercell lattice vectors in units
  !> of the cell vectors. Non-zero errmsg if error.
  module subroutine supercell_matrix_from_ints(ndim,idim,phonopy,smat,flipped,errmsg)
    use tools_math, only: idet3
    use tools_io, only: string
    integer, intent(in) :: ndim
    integer, intent(in) :: idim(:)
    logical, intent(in) :: phonopy
    integer, intent(out) :: smat(3,3)
    logical, intent(out) :: flipped
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i

    errmsg = ""
    flipped = .false.
    smat = 0
    if (ndim == 1 .or. ndim == 3) then
       do i = 1, 3
          smat(i,i) = idim(min(i,ndim))
       end do
    elseif (ndim == 9) then
       if (phonopy) then
          ! phonopy reshapes the nine integers row-wise into its matrix M
          ! and builds the supercell with the rows of transpose(M)
          smat = reshape(idim(1:9),(/3,3/))
       else
          ! NEWCELL: the columns of reshape(idim) are the new lattice vectors
          smat = transpose(reshape(idim(1:9),(/3,3/)))
       end if
    else
       errmsg = "The supercell needs 1, 3 or 9 integers (found " // string(ndim) // ")"
       return
    end if
    if (idet3(smat) == 0) then
       errmsg = "The supercell matrix is singular"
    elseif (idet3(smat) < 0) then
       smat = -smat
       flipped = .true.
    end if

  end subroutine supercell_matrix_from_ints

  !> Build everything the finite-difference force-constant calculation
  !> needs from the structure c and the supercell matrix smat: the
  !> supercell (sc), its lattice-point table (nlat, madj, lvec, lkey), the
  !> independent atoms (nindep, indep) and the list of displacements
  !> (ndisp, datom = supercell atom displaced, ddir = direction in
  !> supercell fractional coordinates). dist is the displacement length
  !> (bohr). If error, return non-zero errmsg.
  !>
  !> This routine was adapted from phonopy, by A. Togo.
  subroutine fc2_disp_setup(c,smat,dist,verbose,sc,seed,nlat,nsat,madj,lvec,lkey,&
     nindep,indep,ndisp,datom,ddir,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use global, only: symprec
    use tools_io, only: uout, string, ioj_left, ioj_right, ferror, warning
    use tools_math, only: idet3
    use param, only: bohrtoa
    type(crystal), intent(inout) :: c
    integer, intent(in) :: smat(3,3)
    real*8, intent(in) :: dist
    logical, intent(in) :: verbose
    type(crystal), intent(inout) :: sc
    type(crystalseed), intent(inout) :: seed
    integer, intent(out) :: nlat, nsat, madj(3,3)
    integer, allocatable, intent(inout) :: lvec(:,:), lkey(:,:)
    integer, intent(out) :: nindep
    integer, allocatable, intent(inout) :: indep(:)
    integer, intent(out) :: ndisp
    integer, allocatable, intent(inout) :: datom(:), ddir(:,:)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character*3 :: pgsymb
    integer :: nopfull, ia, il, is, i, k, nsym, nd, nkeep
    integer :: irot(3,3,48), dsel(3,6), ikeep(48), rp(3,3,48)
    real*8 :: rotm(3,3,48), tt(3)

    errmsg = ""
    if (c%ismolecule) then
       errmsg = "the finite-difference force constants can only be used with crystals"
       return
    end if
    if (dist <= 0d0) then
       errmsg = "The displacement length must be positive"
       return
    end if

    ! the supercell, in phonopy order (see read_phonopy_fc2)
    call fc2_lattice_points(smat,nlat,madj,lvec,lkey,errmsg)
    if (len_trim(errmsg) > 0) return
    nsat = c%ncel * nlat

    call seed%end()
    seed%isused = .true.
    seed%file = c%file
    seed%name = "supercell"
    seed%isformat = c%isformat
    seed%nspc = c%nspc
    allocate(seed%spc(c%nspc))
    seed%spc = c%spc
    seed%useabr = 2
    seed%m_x2c = matmul(c%m_x2c,transpose(real(smat,8)))
    seed%nat = nsat
    allocate(seed%x(3,nsat),seed%is(nsat),seed%atname(nsat))
    do ia = 1, c%ncel
       do il = 1, nlat
          is = (ia-1)*nlat + il
          seed%x(:,is) = fc2_scpos(c%atcel(ia)%x,lvec(:,il),madj,nlat)
          seed%is(is) = c%atcel(ia)%is
          seed%atname(is) = c%spc(c%atcel(ia)%is)%name
       end do
    end do

    ! The symmetry of the supercell is the current symmetry of the
    ! crystal (as loaded, or set with SYM), not re-detected: the
    ! operations compatible with the supercell lattice
    ! (fc2_compatible_ops), their translations mapped to the supercell
    ! basis, and as centering vectors the lattice points of the cell
    ! inside the supercell combined with the centering vectors of the
    ! cell. The kept operations are closed under products and inverses,
    ! so they form a group, and struct_new builds the orbits of the
    ! periodic supercell from them.
    nopfull = c%neqv * c%ncv * nlat
    call fc2_compatible_ops(c,smat,nlat,madj,nkeep,ikeep,rp)
    seed%havesym = 1
    seed%findsym = 0
    seed%checkrepeats = .false.
    seed%neqlist = .false.
    if (allocated(seed%rotm)) deallocate(seed%rotm)
    if (allocated(seed%cen)) deallocate(seed%cen)
    allocate(seed%rotm(3,4,nkeep),seed%cen(3,c%ncv*nlat))
    seed%neqv = nkeep
    do k = 1, nkeep
       seed%rotm(1:3,1:3,k) = rp(:,:,k)
       tt = fc2_scpos(c%rotm(1:3,4,ikeep(k)),(/0,0,0/),madj,nlat)
       seed%rotm(1:3,4,k) = tt - floor(tt)
    end do
    seed%ncv = 0
    do k = 1, c%ncv
       do il = 1, nlat
          seed%ncv = seed%ncv + 1
          tt = fc2_scpos(c%cen(:,k),lvec(:,il),madj,nlat)
          seed%cen(:,seed%ncv) = tt - floor(tt)
       end do
    end do
    call sc%struct_new(seed,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) return
    if (c%neqv == 1 .and. c%ncv == 1) &
       call ferror('fc2_disp_setup','the crystal has no symmetry, so every atom is displaced along six &
          &directions; use SYM to find it before this step',warning)

    ! independent atoms: the first atom of each orbit (reduceatoms
    ! creates the orbits in cell order, so this list is increasing)
    nindep = sc%nneq
    if (allocated(indep)) deallocate(indep)
    allocate(indep(nindep))
    indep = 0
    do i = 1, sc%ncel
       if (indep(sc%atcel(i)%idx) == 0) indep(sc%atcel(i)%idx) = i
    end do

    if (verbose) then
       write (uout,'("+ Supercell matrix (rows = supercell lattice vectors):")')
       do i = 1, 3
          write (uout,'("    ",3(A," "))') (string(smat(i,k),4,ioj_right),k=1,3)
       end do
       write (uout,'("  Number of lattice points in the supercell: ",A)') string(nlat)
       write (uout,'("  Number of atoms in the supercell: ",A)') string(nsat)
       write (uout,'("  Symmetry operations of the crystal in the supercell: ",A)') string(nopfull)
       write (uout,'("  Of those, compatible with the supercell lattice (used here): ",A)') &
          string(sc%neqv*sc%ncv)
       rotm(:,:,1:sc%neqv) = sc%rotm(1:3,1:3,1:sc%neqv)
       write (uout,'("  Point group of the symmetry compatible with the supercell: ",A)') &
          trim(pointgroup_symbol(sc%neqv,rotm(:,:,1:sc%neqv)))
       write (uout,'("  Displacement length (bohr): ",A," (",A," ang)")') &
          string(dist,'f',10,6), string(dist*bohrtoa,'f',10,6)
       write (uout,'("  Number of independent atoms: ",A)') string(nindep)
       write (uout,'("# List of independent atoms")')
       write (uout,'("# id     atom  site symmetry  nsite  ndisp")')
    end if

    ! displacements for each independent atom
    ndisp = 0
    if (allocated(datom)) deallocate(datom)
    if (allocated(ddir)) deallocate(ddir)
    allocate(datom(6*nindep),ddir(3,6*nindep))
    do i = 1, nindep
       is = indep(i)
       pgsymb = sc%sitesymm(sc%atcel(is)%x,max(fc2_epspos,symprec),nsym,rotm)
       irot(:,:,1:nsym) = nint(rotm(:,:,1:nsym))

       call fc2_site_directions(nsym,irot(:,:,1:nsym),nd,dsel)
       datom(ndisp+1:ndisp+nd) = is
       ddir(:,ndisp+1:ndisp+nd) = dsel(:,1:nd)
       ndisp = ndisp + nd
       if (verbose) &
          write (uout,'(2X,A,X,A,X,A,X,A,X,A)') string(i,4), string(is,8), string(pgsymb,14,ioj_left),&
             string(nsym,6), string(nd,6)
    end do

  end subroutine fc2_disp_setup

  !> phonopy's choice of displacement directions for a site whose
  !> site-symmetry rotations are the nsym integer matrices irot (in
  !> the supercell basis): the nd directions dsel(:,1:nd) (integer,
  !> fractional coordinates of the supercell), including the opposite
  !> of each direction when no site operation reverses it. Every
  !> displacement count in critic2 (CREATE_DISPLACEMENTS, NEWCELL NICE)
  !> goes through this routine. Note that the result depends on the
  !> supercell basis, not only on the lattice, because the candidate
  !> directions are fixed integer triples in that basis.
  pure subroutine fc2_site_directions(nsym,irot,nd,dsel)
    integer, intent(in) :: nsym
    integer, intent(in) :: irot(3,3,nsym)
    integer, intent(out) :: nd
    integer, intent(out) :: dsel(3,6)

    ! phonopy's directions_diag, tried in this order
    integer, parameter :: ndirs = 13
    integer, parameter :: dirs(3,ndirs) = reshape((/&
       1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,0,1, 0,1,1, 1,-1,0, 1,0,-1, 0,1,-1,&
       1,1,1, 1,1,-1, 1,-1,1, -1,1,1/),(/3,ndirs/))

    integer :: id, jd, k, k1, k2, nsel, isel(3)
    integer :: rd(3,nsym,ndirs)
    logical :: isminus

    ! images of all the candidate directions under the site symmetry
    do id = 1, ndirs
       do k = 1, nsym
          rd(:,k,id) = matmul(irot(:,:,k),dirs(:,id))
       end do
    end do

    ! one direction whose images span space (needs two operations
    ! besides the identity)
    nsel = 0
    if (nsym >= 3) then
       one: do id = 1, ndirs
          do k1 = 1, nsym-1
             do k2 = k1+1, nsym
                if (itriple(dirs(:,id),rd(:,k1,id),rd(:,k2,id)) /= 0) then
                   nsel = 1
                   isel(1) = id
                   exit one
                end if
             end do
          end do
       end do one
    end if

    ! two directions (needs one operation besides the identity)
    if (nsel == 0 .and. nsym >= 2) then
       two: do id = 1, ndirs
          do k1 = 1, nsym
             do jd = 1, ndirs
                if (itriple(dirs(:,id),rd(:,k1,id),dirs(:,jd)) /= 0) then
                   nsel = 2
                   isel(1) = id
                   isel(2) = jd
                   exit two
                end if
             end do
          end do
       end do two
    end if

    ! three directions: the axes
    if (nsel == 0) then
       nsel = 3
       isel = (/1,2,3/)
    end if

    ! add the opposite displacement if no site operation reverses it
    nd = 0
    do k1 = 1, nsel
       id = isel(k1)
       nd = nd + 1
       dsel(:,nd) = dirs(:,id)
       isminus = .true.
       do k = 1, nsym
          if (all(rd(:,k,id) == -dirs(:,id))) then
             isminus = .false.
             exit
          end if
       end do
       if (isminus) then
          nd = nd + 1
          dsel(:,nd) = -dirs(:,id)
       end if
    end do

  contains
    ! determinant of the matrix with columns a, b, c (no temporaries)
    pure function itriple(a,b,c)
      integer, intent(in) :: a(3), b(3), c(3)
      integer :: itriple

      itriple = a(1) * (b(2)*c(3) - b(3)*c(2)) + a(2) * (b(3)*c(1) - b(1)*c(3)) + &
         a(3) * (b(1)*c(2) - b(2)*c(1))

    end function itriple
  end subroutine fc2_site_directions

  !> The symmetry operations of the crystal compatible with the
  !> supercell lattice given by smat (rows = supercell vectors in cell
  !> units): those whose rotation is an integer matrix in the supercell
  !> basis, R' = S^-T R S^T with S^-1 = adj(S)/n. Returns the size n of
  !> the supercell (nlat) and adj(S) (madj), and the nkeep compatible
  !> operations: their indices in c%rotm (ikeep) and their rotations in
  !> the supercell basis (rp). The kept set is closed under products
  !> and inverses (integer unimodular matrices), so it is a group.
  subroutine fc2_compatible_ops(c,smat,nlat,madj,nkeep,ikeep,rp)
    use tools_math, only: idet3, iadj3
    type(crystal), intent(in) :: c
    integer, intent(in) :: smat(3,3)
    integer, intent(out) :: nlat, madj(3,3), nkeep
    integer, intent(out) :: ikeep(48), rp(3,3,48)

    integer :: k
    real*8 :: rr(3,3)

    nlat = idet3(smat)
    madj = iadj3(smat)
    nkeep = 0
    do k = 1, c%neqv
       rr = matmul(transpose(real(madj,8)),matmul(c%rotm(1:3,1:3,k),transpose(real(smat,8)))) / real(nlat,8)
       if (any(abs(rr - nint(rr)) > 1d-6)) cycle
       nkeep = nkeep + 1
       ikeep(nkeep) = k
       rp(:,:,nkeep) = nint(rr)
    end do

  end subroutine fc2_compatible_ops

  !> Write the supercell and the displaced supercells needed to
  !> compute the second-order force constants by finite
  !> differences. smat is the supercell matrix. template is the file
  !> name pattern, with the character * replaced by the zero-padded
  !> index of the displacement (0 = the undisplaced supercell); if
  !> empty, the name and format are chosen from the source file of the
  !> structure. rklength is passed to the writers that generate a
  !> k-point grid (QE, FHIaims, CASTEP); since the supercell is larger
  !> than the cell, the same rklength gives a coarser grid than in the
  !> cell. If error, return non-zero errmsg.
  !>
  !> This routine was adapted from phonopy, by A. Togo.
  module subroutine create_displacements(c,smat0,dist,template,dataset,scfile,verbose,errmsg,ti,rklength)
    use crystalseedmod, only: crystalseed
    use tools_io, only: uout, string, ioj_right
    use param, only: isformat_w_unknown, isformat_w_vasp, isformat_w_abinit,&
       isformat_w_elk, isformat_w_siesta_struct, isformat_w_dftbp_hsd
    class(crystal), intent(inout) :: c
    integer, intent(in) :: smat0(3,3)
    real*8, intent(in) :: dist
    character*(*), intent(in) :: template
    character*(*), intent(in) :: dataset
    character*(*), intent(in) :: scfile
    logical, intent(in) :: verbose
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti
    real*8, intent(in), optional :: rklength

    character(len=:), allocatable :: template_, fname, otemplate, errmsg2
    character(len=mlen), allocatable :: fnames(:)
    integer :: iwf, nlat, nsat, madj(3,3), is, i, k, nindep, ndisp, npad, nk(3), smat(3,3)
    real*8 :: y(3), dc(3)
    logical :: seen(c%nspc)
    integer, allocatable :: lvec(:,:), lkey(:,:), indep(:), datom(:), ddir(:,:)
    type(crystalseed) :: seed
    type(crystal) :: sc, scd

    errmsg = ""

    ! the supercell: from a structure file if one was given
    smat = smat0
    if (len_trim(scfile) > 0) then
       call fc2_smat_from_cell(c,scfile,smat,errmsg,ti)
       if (len_trim(errmsg) > 0) return
       if (verbose) &
          write (uout,'("+ Supercell read from: ",A)') trim(scfile)
    end if

    ! file names and format
    call disp_default_template(c,template,template_,iwf)
    if (index(template_,'*') == 0) then
       errmsg = "The template for the displaced structures must contain a * character"
       return
    end if
    if (iwf == isformat_w_unknown) then
       errmsg = "Unknown file extension in the template for the displaced structures: " // template_
       return
    end if

    ! Some formats write the atoms grouped by species; the file order would
    ! then differ from the order used here and in the FORCE_CONSTANTS reader,
    ! so refuse them unless the atoms are already grouped
    if (iwf == isformat_w_vasp .or. iwf == isformat_w_abinit .or. iwf == isformat_w_elk .or.&
       iwf == isformat_w_siesta_struct .or. iwf == isformat_w_dftbp_hsd) then
       seen = .false.
       do i = 1, c%ncel
          if (i > 1) then
             if (c%atcel(i)%is /= c%atcel(i-1)%is .and. seen(c%atcel(i)%is)) then
                errmsg = "This file format writes the atoms grouped by species, which would change the &
                   &atom order of the displaced structures; reorder the atoms by species first or &
                   &use a template with a format that keeps the order (QE, cif, FHIaims, ...)"
                return
             end if
          end if
          seen(c%atcel(i)%is) = .true.
       end do
    end if

    ! the supercell and the displacement list
    call fc2_disp_setup(c,smat,dist,verbose,sc,seed,nlat,nsat,madj,lvec,lkey,&
       nindep,indep,ndisp,datom,ddir,errmsg,ti)
    if (len_trim(errmsg) > 0) return

    if (verbose .and. present(rklength)) then
       call sc%get_kpoints(rklength,nk)
       write (uout,'("  K-point grid for the supercell (RKLENGTH = ",A,"): ",3(A," "))') &
          string(rklength,'f',10,2), (string(nk(k)),k=1,3)
    end if

    ! write the undisplaced supercell and then the displaced ones,
    ! displacing one atom of the seed at a time and restoring it. QE
    ! inputs are written for a single-point calculation with forces.
    npad = max(3,len(string(ndisp)))
    fname = fc2_expand_star(template_,0,npad)
    call sc%write_any_file(fname,errmsg,iwformat=iwf,nosym=.true.,ti=ti,forces=.true.,&
       rklength=rklength)
    if (len_trim(errmsg) > 0) return
    if (verbose) then
       write (uout,'("+ Undisplaced supercell written to: ",A)') fname
       write (uout,'("+ List of displacements")')
       write (uout,'("# id  atom  --direction (frac)--  ------ displacement (Cartesian, bohr) ------  file")')
    end if
    seed%havesym = 0 ! the displaced structures do not have the symmetry of the supercell
    seed%findsym = 0
    allocate(fnames(ndisp))
    do i = 1, ndisp
       dc = fc2_dispvec(sc,ddir(:,i),dist)
       is = datom(i)
       y = seed%x(:,is)
       seed%x(:,is) = y + sc%c2x(dc)
       call scd%struct_new(seed,errmsg,ti=ti)
       seed%x(:,is) = y
       if (len_trim(errmsg) > 0) return
       fname = fc2_expand_star(template_,i,npad)
       fnames(i) = fname
       call scd%write_any_file(fname,errmsg,iwformat=iwf,nosym=.true.,ti=ti,forces=.true.,&
          rklength=rklength)
       if (len_trim(errmsg) > 0) return
       if (verbose) &
          write (uout,'(2X,A,X,A,2X,3(A,X),X,3(A,X),X,A)') string(i,4), string(is,5),&
             (string(ddir(k,i),4),k=1,3), (string(dc(k),'f',14,10,ioj_right),k=1,3), trim(fname)
    end do
    if (verbose) &
       write (uout,'("+ Written ",A," displaced supercells")') string(ndisp)

    ! The dataset file: everything READ_FORCES needs to interpret the
    ! forces of these structures. Without it they cannot be read back,
    ! so it is always written.
    call fc2_write_dataset(c,dataset,smat,dist,template_,ndisp,datom,ddir,fnames,errmsg,ti)
    if (len_trim(errmsg) > 0) return
    if (verbose) then
       write (uout,'("+ Displacement dataset written to: ",A)') trim(dataset)
       call fc2_output_template(template_,otemplate,errmsg2)
       if (len_trim(errmsg2) > 0) then
          write (uout,'("  Read the forces back with: VIBRATIONS READ_FORCES TEMPLATE <outputs>")')
       else
          write (uout,'("  Read the forces back with: VIBRATIONS READ_FORCES")')
          write (uout,'("  (the outputs are expected in: ",A,")")') trim(otemplate)
       end if
    end if

  end subroutine create_displacements

  !> Cartesian displacement (bohr) of length dist along the integer
  !> direction ddir, given in fractional coordinates of the supercell sc.
  function fc2_dispvec(sc,ddir,dist) result(dc)
    type(crystal), intent(in) :: sc
    integer, intent(in) :: ddir(3)
    real*8, intent(in) :: dist
    real*8 :: dc(3)

    dc = matmul(sc%m_x2c,real(ddir,8))
    dc = dc * dist / norm2(dc)

  end function fc2_dispvec

  !> Expand a file name template by replacing every * with the integer
  !> i, zero-padded to npad digits. Several stars are allowed (and all
  !> get the same index), so a template like "supercell-* /supercell-*.in"
  !> names one file per displacement inside its own directory.
  function fc2_expand_star(template,i,npad) result(fname)
    use tools_io, only: string
    character*(*), intent(in) :: template
    integer, intent(in) :: i, npad
    character(len=:), allocatable :: fname

    character(len=:), allocatable :: rest
    integer :: istar

    fname = ""
    rest = template
    do while (.true.)
       istar = index(rest,'*')
       if (istar == 0) exit
       fname = fname // rest(1:istar-1) // string(i,npad,pad0=.true.)
       rest = rest(istar+1:)
    end do
    fname = fname // rest

  end function fc2_expand_star

  !> Permutation of the complete atom list of c induced by the symmetry
  !> operation with rotation io and centering vector icv: perm(i) is the
  !> atom that the operation maps atom i onto. ok is false if some image
  !> could not be identified.
  subroutine fc2_atom_perm(c,io,icv,perm,ok)
    use global, only: symprec
    use param, only: icrd_crys
    type(crystal), intent(inout) :: c
    integer, intent(in) :: io, icv
    integer, intent(out) :: perm(:)
    logical, intent(out) :: ok

    integer :: i, j
    real*8 :: xx(3), dd

    ok = .false.
    do i = 1, c%ncel
       xx = matmul(c%rotm(1:3,1:3,io),c%atcel(i)%x) + c%rotm(:,4,io) + c%cen(:,icv)
       j = c%identify_atom(xx,icrd_crys,dist=dd,distmax=max(fc2_epspos,symprec))
       if (j == 0) return
       perm(i) = j
    end do
    ok = .true.

  end subroutine fc2_atom_perm

  !> Read the structure of a displaced supercell from file and compare
  !> it with the undisplaced supercell sc/seed: return in jmax the atom
  !> that moved the most and in dmax the difference (bohr) between that
  !> displacement and the expected one (dexp, Cartesian, on atom iat).
  !> If the file could not be read, return non-zero errmsg.
  subroutine fc2_check_disp(file,sc,seed,iat,dexp,jmax,dmax,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use tools_io, only: string
    character*(*), intent(in) :: file
    type(crystal), intent(in) :: sc
    type(crystalseed), intent(in) :: seed
    integer, intent(in) :: iat
    real*8, intent(in) :: dexp(3)
    integer, intent(out) :: jmax
    real*8, intent(out) :: dmax
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    type(crystalseed) :: dseed
    integer :: j
    real*8 :: dx(3), dd, dmx

    jmax = 0
    dmax = huge(1d0)
    call dseed%read_any_file(file,0,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) then
       errmsg = "Could not read the structure from " // trim(file) // ": " // errmsg
       return
    end if
    if (dseed%nat /= seed%nat) then
       errmsg = "The structure in " // trim(file) // " has " // string(dseed%nat) //&
          " atoms, but the supercell has " // string(seed%nat)
       return
    end if

    ! the atom that moved the most from the undisplaced supercell
    dmx = -1d0
    do j = 1, dseed%nat
       dx = dseed%x(:,j) - seed%x(:,j)
       dx = dx - nint(dx)
       dd = norm2(matmul(sc%m_x2c,dx))
       if (dd > dmx) then
          dmx = dd
          jmax = j
       end if
    end do

    ! how far that displacement is from the expected one
    dx = dseed%x(:,iat) - seed%x(:,iat)
    dx = dx - nint(dx)
    dmax = norm2(matmul(sc%m_x2c,dx) - dexp)

  end subroutine fc2_check_disp

  !> Write the displacement dataset file for a set of displaced
  !> structures just created: the reference structure c, the supercell
  !> matrix smat, the displacement length dist (bohr), the file name
  !> template, and the list of ndisp displacements (displaced supercell
  !> atom datom, direction ddir in supercell fractional coordinates,
  !> and the file written for each displacement). This is the file
  !> create_forces reads to know how the forces it is given were
  !> generated.
  subroutine fc2_write_dataset(c,file,smat,dist,template,ndisp,datom,ddir,fname,errmsg,ti)
    use tools_io, only: fopen_write, fclose, string, ioj_right
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: smat(3,3)
    real*8, intent(in) :: dist
    character*(*), intent(in) :: template
    integer, intent(in) :: ndisp
    integer, intent(in) :: datom(ndisp)
    integer, intent(in) :: ddir(3,ndisp)
    character(len=mlen), intent(in) :: fname(ndisp)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, k

    errmsg = ""
    lu = fopen_write(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Could not open the displacement dataset file for writing: " // trim(file)
       return
    end if

    write (lu,'("# critic2 VIBRATIONS displacement dataset")')
    write (lu,'("# the forces of these displaced structures are read back with VIBRATIONS READ_FORCES")')
    write (lu,'("VERSION 1")')
    write (lu,'("STRUCTURE ! the cell the displacements were generated from")')
    write (lu,'("  LATTICE ! lattice vectors, one per row, in bohr")')
    do i = 1, 3
       write (lu,'(4X,3(A," "))') (string(c%m_x2c(k,i),'f',22,14,ioj_right),k=1,3)
    end do
    write (lu,'("  NATOMS ",A)') string(c%ncel)
    do i = 1, c%ncel
       write (lu,'("  ATOM ",A,X,3(A," "))') string(c%spc(c%atcel(i)%is)%z,3,ioj_right),&
          (string(c%atcel(i)%x(k),'f',22,14,ioj_right),k=1,3)
    end do
    write (lu,'("ENDSTRUCTURE")')
    write (lu,'("SUPERCELL ! rows = supercell lattice vectors, in units of the cell vectors")')
    do i = 1, 3
       write (lu,'(4X,3(A," "))') (string(smat(i,k),4,ioj_right),k=1,3)
    end do
    write (lu,'("DISTANCE ",A," ! displacement length, always in bohr")') string(dist,'f',22,14)
    write (lu,'("TEMPLATE ",A)') trim(template)
    write (lu,'("NDISP ",A)') string(ndisp)
    write (lu,'("DISPLACEMENTS")')
    write (lu,'("# id  atom  --direction (frac)--  file")')
    do i = 1, ndisp
       write (lu,'(2X,A,X,A,2X,3(A,X),X,A)') string(i,4), string(datom(i),5),&
          (string(ddir(k,i),4),k=1,3), trim(fname(i))
    end do
    write (lu,'("ENDDISPLACEMENTS")')
    call fclose(lu)

  end subroutine fc2_write_dataset

  !> Read a displacement dataset file written by fc2_write_dataset into
  !> ds. Returns non-empty errmsg on error.
  subroutine fc2_read_dataset(file,ds,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline, lgetword, getword, isinteger,&
       isreal, equal, string
    character*(*), intent(in) :: file
    type(fc2_dataset), intent(out) :: ds
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word
    integer :: lu, lp, i, k, idx, idum
    logical :: ok

    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Could not open the displacement dataset file: " // trim(file)
       return
    end if
    errmsg = "Error reading the displacement dataset file: " // trim(file)

    do while (getline(lu,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"version")) then
          if (.not.isinteger(ds%version,line,lp)) goto 999
          if (ds%version /= 1) then
             errmsg = "Unknown version (" // string(ds%version) // ") of the displacement dataset file: " //&
                trim(file)
             goto 999
          end if
       elseif (equal(word,"structure")) then
          i = 0
          do while (getline(lu,line))
             lp = 1
             word = lgetword(line,lp)
             if (equal(word,"endstructure")) then
                exit
             elseif (equal(word,"lattice")) then
                do idx = 1, 3
                   if (.not.getline(lu,line)) goto 999
                   lp = 1
                   do k = 1, 3
                      if (.not.isreal(ds%m_x2c(k,idx),line,lp)) goto 999
                   end do
                end do
             elseif (equal(word,"natoms")) then
                if (allocated(ds%z)) then
                   errmsg = "Repeated NATOMS in the displacement dataset file: " // trim(file)
                   goto 999
                end if
                if (.not.isinteger(ds%nat,line,lp)) goto 999
                if (ds%nat < 1) goto 999
                allocate(ds%z(ds%nat),ds%x(3,ds%nat))
                ds%z = 0
                ds%x = 0d0
             elseif (equal(word,"atom")) then
                if (.not.allocated(ds%z)) then
                   errmsg = "ATOM before NATOMS in the displacement dataset file: " // trim(file)
                   goto 999
                end if
                i = i + 1
                if (i > ds%nat) then
                   errmsg = "Too many ATOM lines in the displacement dataset file: " // trim(file)
                   goto 999
                end if
                if (.not.isinteger(ds%z(i),line,lp)) goto 999
                do k = 1, 3
                   if (.not.isreal(ds%x(k,i),line,lp)) goto 999
                end do
             else
                errmsg = "Unknown keyword (" // trim(word) // ") in the STRUCTURE block of the &
                   &displacement dataset file: " // trim(file)
                goto 999
             end if
          end do
          if (.not.allocated(ds%z)) goto 999
          if (i /= ds%nat) then
             errmsg = "The STRUCTURE block of " // trim(file) // " announces " // string(ds%nat) //&
                " atoms but has " // string(i)
             goto 999
          end if
       elseif (equal(word,"supercell")) then
          do i = 1, 3
             if (.not.getline(lu,line)) goto 999
             lp = 1
             do k = 1, 3
                if (.not.isinteger(ds%smat(i,k),line,lp)) goto 999
             end do
          end do
       elseif (equal(word,"distance")) then
          if (.not.isreal(ds%dist,line,lp)) goto 999
       elseif (equal(word,"template")) then
          ds%template = getword(line,lp)
       elseif (equal(word,"ndisp")) then
          if (.not.isinteger(ds%ndisp,line,lp)) goto 999
          if (ds%ndisp < 1) goto 999
       elseif (equal(word,"displacements")) then
          if (ds%ndisp < 1) then
             errmsg = "DISPLACEMENTS before NDISP in the displacement dataset file: " // trim(file)
             goto 999
          end if
          if (allocated(ds%datom)) deallocate(ds%datom)
          if (allocated(ds%ddir)) deallocate(ds%ddir)
          if (allocated(ds%fname)) deallocate(ds%fname)
          allocate(ds%datom(ds%ndisp),ds%ddir(3,ds%ndisp),ds%fname(ds%ndisp))
          i = 0
          do while (getline(lu,line))
             lp = 1
             word = lgetword(line,lp)
             if (equal(word,"enddisplacements")) exit
             i = i + 1
             if (i > ds%ndisp) then
                errmsg = "Too many displacements in the displacement dataset file: " // trim(file)
                goto 999
             end if
             lp = 1
             ok = isinteger(idum,line,lp)
             ok = ok .and. isinteger(ds%datom(i),line,lp)
             do k = 1, 3
                ok = ok .and. isinteger(ds%ddir(k,i),line,lp)
             end do
             if (.not.ok) goto 999
             if (idum /= i) then
                errmsg = "The displacements in " // trim(file) // " are not in order (found index " //&
                   string(idum) // " on displacement " // string(i) // ")"
                goto 999
             end if
             ds%fname(i) = getword(line,lp)
          end do
          if (i /= ds%ndisp) then
             errmsg = "The displacement dataset file " // trim(file) // " announces " // string(ds%ndisp) //&
                " displacements but has " // string(i)
             goto 999
          end if
       else
          errmsg = "Unknown keyword (" // trim(word) // ") in the displacement dataset file: " // trim(file)
          goto 999
       end if
    end do
    call fclose(lu)

    ! everything must be there
    errmsg = ""
    if (ds%version == 0) then
       errmsg = "No VERSION in the displacement dataset file: " // trim(file)
    elseif (ds%nat == 0) then
       errmsg = "No STRUCTURE in the displacement dataset file: " // trim(file)
    elseif (all(ds%smat == 0)) then
       errmsg = "No SUPERCELL in the displacement dataset file: " // trim(file)
    elseif (ds%dist <= 0d0) then
       errmsg = "No (or non-positive) DISTANCE in the displacement dataset file: " // trim(file)
    elseif (.not.allocated(ds%datom)) then
       errmsg = "No DISPLACEMENTS in the displacement dataset file: " // trim(file)
    end if
    return

999 continue
    call fclose(lu)

  end subroutine fc2_read_dataset

  !> Check the displacement dataset ds against the current structure c
  !> and against the supercell matrix smat, displacement length dist,
  !> and displacement list (ndisp,datom,ddir) regenerated by
  !> fc2_disp_setup. Returns non-empty errmsg if any of them disagrees.
  subroutine fc2_check_dataset(c,ds,smat,dist,ndisp,datom,ddir,errmsg)
    use tools_io, only: string
    type(crystal), intent(in) :: c
    type(fc2_dataset), intent(in) :: ds
    integer, intent(in) :: smat(3,3)
    real*8, intent(in) :: dist
    integer, intent(in) :: ndisp
    integer, intent(in) :: datom(ndisp)
    integer, intent(in) :: ddir(3,ndisp)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i
    real*8 :: dev, x(3)

    errmsg = ""

    ! the reference structure
    if (ds%nat /= c%ncel) then
       errmsg = "The dataset was generated from a cell with " // string(ds%nat) // " atoms, but the &
          &current structure has " // string(c%ncel)
       return
    end if
    dev = maxval(abs(ds%m_x2c - c%m_x2c)) / max(maxval(abs(c%m_x2c)),1d-10)
    if (dev > fc2_epsdataset) then
       errmsg = "The lattice of the current structure differs from the one the displacements were &
          &generated from (relative deviation " // string(dev,'e',12,4) // ")"
       return
    end if
    do i = 1, c%ncel
       if (ds%z(i) /= c%spc(c%atcel(i)%is)%z) then
          errmsg = "Atom " // string(i) // " is " // string(c%spc(c%atcel(i)%is)%z) // " in the current &
             &structure but " // string(ds%z(i)) // " in the dataset; the atom order must be the same"
          return
       end if
       x = ds%x(:,i) - c%atcel(i)%x
       x = x - nint(x)
       if (any(abs(x) > fc2_epsdataset)) then
          errmsg = "Atom " // string(i) // " has moved since the displacements were generated (by " //&
             string(norm2(matmul(c%m_x2c,x)),'e',12,4) // " bohr)"
          return
       end if
    end do

    ! the supercell and the displacement length
    if (any(ds%smat /= smat)) then
       errmsg = "The supercell given (" // fc2_smatstr(smat) // ") is not the one in the dataset (" //&
          fc2_smatstr(ds%smat) // ")"
       return
    end if
    if (abs(ds%dist - dist) > fc2_epsdisp) then
       errmsg = "The displacement length given (" // string(dist,'f',12,8) // " bohr) is not the one in &
          &the dataset (" // string(ds%dist,'f',12,8) // " bohr)"
       return
    end if

    ! the displacement list
    if (ds%ndisp /= ndisp) then
       errmsg = "The dataset has " // string(ds%ndisp) // " displacements but " // string(ndisp) //&
          " were regenerated for this supercell; the symmetry may have changed between the two runs (SYMPREC, SYM, NOSYM?)"
       return
    end if
    do i = 1, ndisp
       if (ds%datom(i) /= datom(i) .or. any(ds%ddir(:,i) /= ddir(:,i))) then
          errmsg = "Displacement " // string(i) // " in the dataset (atom " // string(ds%datom(i)) //&
             ", direction " // string(ds%ddir(1,i)) // " " // string(ds%ddir(2,i)) // " " //&
             string(ds%ddir(3,i)) // ") is not the one regenerated here (atom " // string(datom(i)) //&
             ", direction " // string(ddir(1,i)) // " " // string(ddir(2,i)) // " " //&
             string(ddir(3,i)) // "); the symmetry may have changed between the two runs (SYMPREC, SYM, NOSYM?)"
          return
       end if
    end do

  end subroutine fc2_check_dataset

  !> A supercell matrix written on one line, for an error message.
  function fc2_smatstr(m) result(str)
    use tools_io, only: string
    integer, intent(in) :: m(3,3)
    character(len=:), allocatable :: str

    integer :: i, k

    str = ""
    do i = 1, 3
       if (i > 1) str = str // " /"
       do k = 1, 3
          str = str // " " // string(m(i,k))
       end do
    end do
    str = adjustl(str)

  end function fc2_smatstr

  !> Guess the file name template of the outputs of a set of displaced
  !> structure files whose template is template: the same template with
  !> the extension of the input replaced by the extension of the
  !> corresponding output. Only Quantum ESPRESSO is known here, which
  !> is also the only code whose forces can be read.
  subroutine fc2_output_template(template,otemplate,errmsg)
    use tools_io, only: lower
    character*(*), intent(in) :: template
    character(len=:), allocatable, intent(out) :: otemplate
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: idx
    character(len=:), allocatable :: aux

    errmsg = ""
    otemplate = ""
    aux = lower(trim(template))
    idx = len_trim(aux) - 2
    if (idx > 0) then
       if (aux(idx:) == ".in") then
          otemplate = template(1:idx-1) // ".out"
          return
       end if
    end if
    errmsg = "Could not guess the name of the output files from the template in the dataset (" //&
       trim(template) // "); give the outputs with TEMPLATE or LIST"

  end subroutine fc2_output_template

  !> Deduce the supercell matrix (rows = supercell lattice vectors, in
  !> units of the cell vectors, which is phonopy's convention) from the
  !> structure in file scfile, an integer supercell of c. Checks that
  !> the transformation is integer and that the two metrics are
  !> consistent. Returns non-empty errmsg on error.
  subroutine fc2_smat_from_cell(c,scfile,smat,errmsg,ti,sco,seedo)
    use crystalseedmod, only: crystalseed
    use tools_io, only: string
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: scfile
    integer, intent(out) :: smat(3,3)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti
    type(crystal), intent(inout), optional :: sco
    type(crystalseed), intent(inout), optional :: seedo

    real*8 :: rmat(3,3), guc(3,3), gsc(3,3), gchk(3,3), dev
    type(crystal) :: sc
    type(crystalseed) :: seed

    smat = 0
    call seed%read_any_file(scfile,-1,errmsg,ti=ti)
    if (len_trim(errmsg) > 0) return
    call sc%struct_new(seed,errmsg,noenv=.true.,ti=ti)
    if (len_trim(errmsg) > 0) return
    if (present(seedo)) seedo = seed
    if (present(sco)) sco = sc
    if (sc%ismolecule) then
       errmsg = "The supercell file " // trim(scfile) // " is a molecule"
       return
    end if

    ! columns of m_c2x*sc%m_x2c are the supercell vectors in cell coordinates
    rmat = transpose(matmul(c%m_c2x,sc%m_x2c))
    smat = nint(rmat)
    dev = maxval(abs(rmat - real(smat,8)))
    if (dev > fc2_epsint) then
       errmsg = "The cell in " // trim(scfile) // " is not an integer supercell of the current &
          &structure (deviation " // string(dev,'e',12,4) // ")"
       return
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
       return
    end if

  end subroutine fc2_smat_from_cell

  !> Read the forces on the nat atoms of a supercell from the output
  !> of a Quantum ESPRESSO calculation, in Hartree/bohr and in the
  !> atom order of the input. The last force block in the file is
  !> used. The mean force over the atoms (the drift) is subtracted. If
  !> error, return non-zero errmsg.
  subroutine fc2_read_forces(file,nat,f,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword, equal, isreal, string
    character*(*), intent(in) :: file
    integer, intent(in) :: nat
    real*8, intent(out) :: f(3,nat)
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, word
    integer :: lu, n, lp, i
    real*8 :: fx(3), fmean(3)
    logical :: found, ok

    ! Ry/bohr (QE) to Hartree/bohr
    real*8, parameter :: fac = 0.5d0

    errmsg = ""
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu <= 0) then
       errmsg = "Could not open the force file: " // trim(file)
       return
    end if

    found = .false.
    do while (getline_raw(lu,line))
       if (index(line,"Forces acting on atoms") == 0) cycle

       n = 0
       do while (getline_raw(lu,line))
          if (len_trim(line) == 0) then
             if (n > 0) exit
             cycle
          end if
          lp = 1
          word = lgetword(line,lp)
          if (.not.equal(word,"atom")) exit
          lp = index(line,"=")
          if (lp == 0) exit
          lp = lp + 1
          ok = isreal(fx(1),line,lp)
          ok = ok .and. isreal(fx(2),line,lp)
          ok = ok .and. isreal(fx(3),line,lp)
          if (.not.ok) then
             errmsg = "Error reading a force line in " // trim(file) // ": " // trim(line)
             goto 999
          end if
          n = n + 1
          if (n <= nat) f(:,n) = fx * fac
       end do

       if (n /= nat) then
          errmsg = "The force block in " // trim(file) // " has " // string(n) // " atoms, but " //&
             string(nat) // " were expected for this supercell"
          goto 999
       end if
       found = .true.
    end do
    call fclose(lu)
    lu = -1

    if (.not.found) then
       errmsg = "No force block found in " // trim(file) // ". Only Quantum ESPRESSO outputs &
          &(with tprnfor) can be read for now"
       return
    end if

    ! subtract the drift
    do i = 1, 3
       fmean(i) = sum(f(i,1:nat)) / real(nat,8)
       f(i,1:nat) = f(i,1:nat) - fmean(i)
    end do

    return
999 continue
    if (lu > 0) call fclose(lu)

  end subroutine fc2_read_forces

  !> Calculate the second-order force constants by finite differences
  !> from the forces of a set of already-run displaced supercells, and
  !> store them in c%vib. smat and dist must be the same supercell
  !> matrix and displacement length used to write the displaced
  !> structures with create_displacements. file locates the outputs:
  !> if it contains a *, it is a template with the * replaced by the
  !> zero-padded displacement index, and otherwise it is a text file
  !> with one file name per line, in displacement order.
  !> If error, return non-zero errmsg.
  module subroutine create_forces(c,file,dataset,verbose,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use tools_io, only: uout, string, fopen_read, fclose, getline_raw, getword
    use tools_math, only: matinv, eigsym
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    character*(*), intent(in) :: dataset
    logical, intent(in) :: verbose
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: line, file_
    character(len=mlen), allocatable :: fname(:)
    integer :: nlat, nsat, madj(3,3), nindep, ndisp, nop, smat(3,3)
    integer :: i, j, k, m, ia, is, js, id, ip, iq, io, icv, ier, npad, lu, nf, lp
    integer :: nd, ns, jmax
    real*8 :: gg(3,3), ggev(3,3), eval(3), atb(3,3), dc(3), rr(3,3), blk(3,3), dmax, dist
    logical :: ok
    integer, allocatable :: lvec(:,:), lkey(:,:), indep(:), datom(:), ddir(:,:)
    integer, allocatable :: perm(:,:), isite(:), ipsite(:,:)
    real*8, allocatable :: fall(:,:,:), rcart(:,:,:), amat(:,:), bmat(:,:), fcs(:,:,:,:)
    type(crystal) :: sc
    type(crystalseed) :: seed
    type(fc2_dataset) :: ds

    ! The displacement dataset written by create_displacements says how
    ! the displaced structures were generated: the supercell, the
    ! displacement length, the reference structure and the list of
    ! displacements. It is the only source of that information.
    call fc2_read_dataset(dataset,ds,errmsg,ti)
    if (len_trim(errmsg) > 0) return
    smat = ds%smat
    dist = ds%dist
    if (verbose) then
       write (uout,'("+ Displacement dataset read from: ",A)') trim(dataset)
       write (uout,'("  Displacement length: ",A," bohr (",A," ang)")') &
          string(dist,'f',10,6), string(dist*bohrtoa,'f',10,6)
    end if

    ! the supercell and the displacement list, exactly as create_displacements built them
    call fc2_disp_setup(c,smat,dist,verbose,sc,seed,nlat,nsat,madj,lvec,lkey,&
       nindep,indep,ndisp,datom,ddir,errmsg,ti)
    if (len_trim(errmsg) > 0) return

    ! the dataset knows the structure and the displacements the forces
    ! were calculated for; check they are the ones regenerated here
    call fc2_check_dataset(c,ds,smat,dist,ndisp,datom,ddir,errmsg)
    if (len_trim(errmsg) > 0) then
       errmsg = trim(errmsg) // " (displacement dataset " // trim(dataset) // ")"
       return
    end if

    ! where the outputs are: the argument if given, otherwise guessed
    ! from the template of the displaced structures in the dataset
    if (len_trim(file) > 0) then
       file_ = file
    else
       call fc2_output_template(ds%template,file_,errmsg)
       if (len_trim(errmsg) > 0) return
       if (verbose) &
          write (uout,'("  Output files guessed from the dataset template: ",A)') trim(file_)
    end if

    ! the files with the forces, one per displacement
    allocate(fname(ndisp))
    npad = max(3,len(string(ndisp)))
    if (index(file_,'*') > 0) then
       do i = 1, ndisp
          fname(i) = fc2_expand_star(file_,i,npad)
       end do
    else
       lu = fopen_read(file_,errstop=.false.,ti=ti)
       if (lu <= 0) then
          errmsg = "Could not open the list of force files: " // trim(file_)
          return
       end if
       nf = 0
       do while (getline_raw(lu,line))
          i = index(line,'#')
          if (i > 0) line = line(1:i-1)
          if (len_trim(line) == 0) cycle
          nf = nf + 1
          if (nf <= ndisp) then
             lp = 1
             fname(nf) = getword(line,lp)
          end if
       end do
       call fclose(lu)
       if (nf /= ndisp) then
          errmsg = "The list in " // trim(file_) // " has " // string(nf) // " file names, but " //&
             string(ndisp) // " displacements were generated for this supercell"
          return
       end if
    end if

    ! read the forces
    allocate(fall(3,nsat,ndisp))
    if (verbose) then
       write (uout,'("+ List of displacements and force files")')
       write (uout,'("# id  atom  --direction (frac)--  file")')
    end if
    do i = 1, ndisp
       call fc2_read_forces(fname(i),nsat,fall(:,:,i),errmsg,ti)
       if (len_trim(errmsg) > 0) return

       ! check the geometry in the file is the displacement we think it is;
       ! this catches a wrong supercell, a wrong displacement length, and
       ! files given in the wrong order
       call fc2_check_disp(fname(i),sc,seed,datom(i),fc2_dispvec(sc,ddir(:,i),dist),&
          jmax,dmax,errmsg,ti)
       if (len_trim(errmsg) > 0) return
       if (jmax /= datom(i) .or. dmax > fc2_epsdisp) then
          errmsg = "The structure in " // trim(fname(i)) // " is not the displacement expected for &
             &file " // string(i) // " (atom " // string(datom(i)) // " displaced by " //&
             string(norm2(fc2_dispvec(sc,ddir(:,i),dist)),'f',10,6) // " bohr): the largest &
             &displacement found is " // string(dmax,'e',12,4) // " bohr on atom " // string(jmax) //&
             ". Check the order of the files, and that they belong to the run described by the &
             &displacement dataset"
          return
       end if

       if (verbose) &
          write (uout,'(2X,A,X,A,2X,3(A,X),X,A)') string(i,4), string(datom(i),5),&
             (string(ddir(k,i),4),k=1,3), trim(fname(i))
    end do

    ! the atom permutations induced by the operations of the supercell
    nop = sc%neqv * sc%ncv
    allocate(perm(nsat,nop),rcart(3,3,nop),stat=ier)
    if (ier /= 0) then
       errmsg = "Could not allocate the symmetry permutation table (" //&
          string(nint(4d0*nsat*nop/1024d0**2)) // " Mb)"
       return
    end if
    k = 0
    do io = 1, sc%neqv
       do icv = 1, sc%ncv
          k = k + 1
          call fc2_atom_perm(sc,io,icv,perm(:,k),ok)
          if (.not.ok) then
             errmsg = "Could not build the atom permutation of a supercell symmetry operation"
             return
          end if
          rcart(:,:,k) = matmul(sc%m_x2c,matmul(sc%rotm(1:3,1:3,io),sc%m_c2x))
       end do
    end do

    ! solve for the rows of the independent atoms
    allocate(fcs(3,3,nindep,nsat),isite(nop),stat=ier)
    if (ier /= 0) then
       errmsg = "Could not allocate " // string(nint(72d0*nindep*nsat/1024d0**2)) //&
          " Mb for the force constants"
       return
    end if
    do i = 1, nindep
       is = indep(i)

       ! the operations that leave this atom in place (its site symmetry)
       ns = 0
       do k = 1, nop
          if (perm(is,k) == is) then
             ns = ns + 1
             isite(ns) = k
          end if
       end do

       ! the atom each site operation sends onto every other atom
       if (allocated(ipsite)) deallocate(ipsite)
       allocate(ipsite(nsat,ns))
       do k = 1, ns
          do j = 1, nsat
             ipsite(perm(j,isite(k)),k) = j
          end do
       end do

       ! the rotated displacements and the normal-equation matrix
       nd = count(datom(1:ndisp) == is)
       if (allocated(amat)) deallocate(amat)
       if (allocated(bmat)) deallocate(bmat)
       allocate(amat(3,nd*ns),bmat(3,nd*ns))
       m = 0
       do id = 1, ndisp
          if (datom(id) /= is) cycle
          dc = fc2_dispvec(sc,ddir(:,id),dist)
          do k = 1, ns
             m = m + 1
             amat(:,m) = matmul(rcart(:,:,isite(k)),dc)
          end do
       end do
       gg = matmul(amat,transpose(amat))

       ! guard the conditioning: matinv only fails on an exactly zero
       ! pivot, and gg is full rank by construction but can still be
       ! badly conditioned for an oblique supercell with a low-symmetry site
       ggev = gg
       call eigsym(ggev,3,eval,ier)
       if (ier /= 0 .or. eval(3) <= 0d0 .or. eval(1)/eval(3) < fc2_epscond) then
          errmsg = "Ill-conditioned least-squares matrix for the displacements of supercell atom " //&
             string(is) // " (eigenvalue ratio " // string(eval(1)/max(eval(3),tiny(1d0)),'e',12,4) //&
             "); the displacement directions do not span space well enough"
          return
       end if
       call matinv(gg,3,ier)
       if (ier /= 0) then
          errmsg = "Singular least-squares matrix for the displacements of supercell atom " // string(is)
          return
       end if

       ! solve for every atom of the supercell
       do js = 1, nsat
          m = 0
          do id = 1, ndisp
             if (datom(id) /= is) cycle
             do k = 1, ns
                m = m + 1
                bmat(:,m) = matmul(rcart(:,:,isite(k)),fall(:,ipsite(js,k),id))
             end do
          end do
          atb = matmul(amat,transpose(bmat))
          fcs(:,:,i,js) = -matmul(gg,atb)
       end do
    end do

    ! the rows of the remaining atoms of the cell, by symmetry
    call c%vib%end(keepvibs=.true.)
    allocate(c%vib%fc2(3,3,c%ncel,nsat),stat=ier)
    if (ier /= 0) then
       errmsg = "Could not allocate " // string(nint(72d0*c%ncel*nsat/1024d0**2)) //&
          " Mb for the force constants"
       return
    end if
    do ia = 1, c%ncel
       is = (ia-1)*nlat + 1

       ! this atom may be one of the independent ones
       ip = 0
       do i = 1, nindep
          if (indep(i) == is) then
             ip = i
             exit
          end if
       end do
       if (ip > 0) then
          c%vib%fc2(:,:,ia,:) = fcs(:,:,ip,:)
          cycle
       end if

       ! otherwise, an operation carrying it onto an independent atom
       iq = 0
       oper: do k = 1, nop
          do i = 1, nindep
             if (perm(is,k) == indep(i)) then
                iq = k
                ip = i
                exit oper
             end if
          end do
       end do oper
       if (iq == 0) then
          errmsg = "No symmetry operation carries atom " // string(is) // " of the supercell onto &
             &an atom whose force constants were calculated"
          return
       end if
       rr = rcart(:,:,iq)
       do js = 1, nsat
          blk = fcs(:,:,ip,perm(js,iq))
          c%vib%fc2(:,:,ia,js) = matmul(transpose(rr),matmul(blk,rr))
       end do
    end do

    ! wrap up
    c%vib%fc2_file = file_
    c%vib%hasfc2 = .true.
    c%vib%fc2_smat = smat
    c%vib%fc2_madj = madj
    c%vib%fc2_nlat = nlat
    c%vib%fc2_nsat = nsat
    c%vib%fc2_ncel = c%ncel
    c%vib%fc2_iscompact = (c%ncel < nsat)
    if (allocated(c%vib%fc2_lvec)) deallocate(c%vib%fc2_lvec)
    if (allocated(c%vib%fc2_lkey)) deallocate(c%vib%fc2_lkey)
    call move_alloc(lvec,c%vib%fc2_lvec)
    call move_alloc(lkey,c%vib%fc2_lkey)
    c%vib%fc2_acoustic = -1
    c%vib%fc2_vs_center = -1d0
    c%vib%fc2_vs_delta = -1d0
    call fc2_stamp_geometry(c%vib,c)

    if (verbose) &
       write (uout,'("+ Force constants calculated from ",A," displaced supercells")') string(ndisp)

  end subroutine create_forces

  !> What CREATE_DISPLACEMENTS would generate in the supercell smat
  !> (rows = supercell vectors in cell units), without building it: the
  !> number of displacements (ndisp) and of independent atoms
  !> (nindep).
  module subroutine disp_count(c,smat,ndisp,nindep,errmsg)
    use global, only: symprec
    use tools_math, only: idet3
    use param, only: icrd_crys
    class(crystal), intent(inout) :: c
    integer, intent(in) :: smat(3,3)
    integer, intent(out) :: ndisp, nindep
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: nlat, madj(3,3), nkeep, k, m, ia, id, nsym, nd
    integer :: ikeep(48), rp(3,3,48), irot(3,3,48), dsel(3,6)
    real*8 :: x0(3)
    logical :: issite
    logical, allocatable :: done(:)

    errmsg = ""
    ndisp = 0
    nindep = 0
    if (c%ismolecule) then
       errmsg = "the finite-difference force constants can only be used with crystals"
       return
    end if
    if (idet3(smat) <= 0) then
       errmsg = "The supercell matrix must have positive determinant"
       return
    end if
    call fc2_compatible_ops(c,smat,nlat,madj,nkeep,ikeep,rp)

    ! orbits, representatives and site symmetries (the same tolerance
    ! as reduceatoms uses to build the orbits of the supercell)
    allocate(done(c%ncel))
    done = .false.
    do ia = 1, c%ncel
       if (done(ia)) cycle
       nindep = nindep + 1
       nsym = 0
       do k = 1, nkeep
          issite = .false.
          do m = 1, c%ncv
             x0 = matmul(c%rotm(1:3,1:3,ikeep(k)),c%atcel(ia)%x) + c%rotm(1:3,4,ikeep(k)) + c%cen(:,m)
             id = c%identify_atom(x0,icrd_crys,distmax=symprec)
             if (id == 0) then
                errmsg = "a symmetry operation does not map the atoms onto themselves"
                return
             end if
             done(id) = .true.
             issite = issite .or. (id == ia)
          end do
          if (issite) then
             nsym = nsym + 1
             irot(:,:,nsym) = rp(:,:,k)
          end if
       end do
       call fc2_site_directions(nsym,irot(:,:,1:nsym),nd,dsel)
       ndisp = ndisp + nd
    end do

  end subroutine disp_count

  !> Default file name template and format for the displaced
  !> structures written by create_displacements. If template is not
  !> empty, it is used as is and the format is detected from its
  !> extension. Otherwise both are chosen from the source
  !> file of the structure c: the same or a related format when critic2
  !> can write it (for instance, a QE output gives QE inputs), and
  !> FHIaims geometry (geometry.in) otherwise. The root of the name is
  !> the source file name without its extension (also without a .scf or
  !> similar QE infix), followed by "-*".
  subroutine disp_default_template(c,template,template_,iwf)
    use tools_io, only: lower
    use param, only: dirsep, isformat_w_qein, isformat_w_vasp,&
       isformat_w_cif, isformat_w_abinit, isformat_w_elk, isformat_w_siesta_struct,&
       isformat_w_castepcell, isformat_w_crystal, isformat_w_dftbp_hsd, isformat_w_aimsin,&
       isformat_r_qein, isformat_r_qeout, isformat_r_vasp, isformat_r_cif, isformat_r_abinit,&
       isformat_r_elk, isformat_r_siesta, isformat_r_castepcell, isformat_r_castepgeom,&
       isformat_r_crystal, isformat_r_dmain
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: template
    character(len=:), allocatable, intent(out) :: template_
    integer, intent(out) :: iwf

    character(len=:), allocatable :: root, ext
    integer :: idx, istar

    if (len_trim(template) > 0) then
       ! detect the format from the template, with the * filled in (unknown
       ! extension -> isformat_w_unknown, the caller decides)
       template_ = trim(template)
       call struct_detect_write_format(fc2_expand_star(template_,0,3),iwf)
       return
    end if

    ! the root: source file name without directory and extension
    root = trim(c%file)
    idx = index(root,dirsep,back=.true.)
    if (idx > 0) root = root(idx+1:)
    if (len_trim(root) == 0) root = "supercell"
    idx = index(root,'.',back=.true.)
    if (idx > 1) root = root(1:idx-1)

    ! format and extension from the source format
    select case (c%isformat)
    case (isformat_r_qein, isformat_r_qeout)
       ! drop a .scf/.relax/... infix, so that x.scf.out gives x-*.scf.in
       idx = index(root,'.',back=.true.)
       if (idx > 1) then
          select case (lower(root(idx+1:)))
          case ("scf","relax","vc-relax","nscf","md","vc-md","pw")
             root = root(1:idx-1)
          end select
       end if
       ext = ".scf.in"
       iwf = isformat_w_qein
    case (isformat_r_vasp)
       ext = ""
       iwf = isformat_w_vasp
    case (isformat_r_cif)
       ext = ".cif"
       iwf = isformat_w_cif
    case (isformat_r_abinit)
       ext = ".abin"
       iwf = isformat_w_abinit
    case (isformat_r_elk)
       ext = ".elk"
       iwf = isformat_w_elk
    case (isformat_r_siesta)
       ext = ".STRUCT_IN"
       iwf = isformat_w_siesta_struct
    case (isformat_r_castepcell, isformat_r_castepgeom)
       ext = ".cell"
       iwf = isformat_w_castepcell
    case (isformat_r_crystal)
       ext = ".d12"
       iwf = isformat_w_crystal
    case (isformat_r_dmain)
       ext = ".hsd"
       iwf = isformat_w_dftbp_hsd
    case default
       ! FHIaims geometry (also for aims inputs and outputs)
       ext = ".in"
       iwf = isformat_w_aimsin
    end select
    template_ = root // "-*" // ext

  end subroutine disp_default_template

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
  !> Print the frequencies at the q-point with index id: the
  !> coordinates of the point, fractional and Cartesian, and the
  !> frequency of every mode in cm-1 and THz. A negative frequency is
  !> an imaginary (unstable) mode.
  module subroutine vibrations_print_freq(v,c,id)
    use global, only: iunit, iunitname0, dunit0
    use tools_io, only: uout, ferror, faterr, string, ioj_right
    use param, only: cm1tothz
    class(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    integer, intent(in) :: id

    integer :: i, j
    real*8 :: qc(3)

    ! check the q-point id
    if (id <= 0 .or. id > v%nqpt) &
       call ferror('vibrations_print_freq','q-point index out of range',faterr)

    ! the q-point, in fractional and Cartesian coordinates
    qc = c%rx2rc(v%qpt(:,id)) / dunit0(iunit)
    write (uout,'("# q-point ",A,": ",3(A," "),"(fractional)")') string(id),&
       (string(v%qpt(j,id),'f',12,7,ioj_right),j=1,3)
    write (uout,'("#",11X,3(A," "),"(Cartesian, ",A,"^-1)")') &
       (string(qc(j),'f',12,7,ioj_right),j=1,3), iunitname0(iunit)

    ! the frequencies
    write (uout,'("#Id   Frequency(cm^-1)   Frequency(THz)")')
    do i = 1, v%nfreq
       write (uout,'(3(A," "))') string(i,4,ioj_right),&
          string(v%freq(i,id),'f',16,6,ioj_right),&
          string(v%freq(i,id) * cm1tothz,'f',16,6,ioj_right)
    end do

  end subroutine vibrations_print_freq

  !> Write the q-points and frequencies currently stored to a file, as
  !> a table ready for plotting: the cumulative distance in reciprocal
  !> space (bohr^-1, the abscissa of a dispersion plot), the q-point in
  !> fractional coordinates, and the frequencies in cm^-1, one row per
  !> q-point. Returns non-empty errmsg on error.
  module subroutine vibrations_write_freq(v,c,file,verbose,errmsg)
    use tools_io, only: fopen_write, fclose, uout, string, ioj_right
    class(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    logical, intent(in) :: verbose
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: lu, i, j
    real*8 :: dq, qc(3), qcold(3)

    errmsg = ""
    if (.not.v%hasvibs .or. v%nqpt < 1) then
       errmsg = "No frequencies available to write"
       return
    end if

    lu = fopen_write(file,errstop=.false.)
    if (lu < 0) then
       errmsg = "Could not open the frequency file for writing: " // trim(file)
       return
    end if
    write (lu,'("# q-point frequencies calculated by critic2")')
    write (lu,'("# column 1: cumulative distance along the list of q-points (bohr^-1)")')
    write (lu,'("# columns 2-4: q-point (fractional coordinates of the reciprocal cell)")')
    write (lu,'("# columns 5-",A,": frequencies (cm^-1)")') string(4+v%nfreq)
    dq = 0d0
    do i = 1, v%nqpt
       qc = c%rx2rc(v%qpt(:,i))
       if (i > 1) dq = dq + norm2(qc - qcold)
       qcold = qc
       write (lu,'(*(A," "))') string(dq,'f',16,10,ioj_right),&
          (string(v%qpt(j,i),'f',14,8,ioj_right),j=1,3),&
          (string(v%freq(j,i),'f',16,6,ioj_right),j=1,v%nfreq)
    end do
    call fclose(lu)

    if (verbose) &
       write (uout,'("+ Frequencies of ",A," q-points written to: ",A)') string(v%nqpt), trim(file)

  end subroutine vibrations_write_freq

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
    ! into the self term.
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
  module subroutine vibrations_write_fc2(v,c,sline,verbose,errmsg)
    use tools_io, only: fopen_write, fclose, uout, string, ioj_right, lgetword, getword, equal
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: sline
    logical, intent(in) :: verbose
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: lu, ia, ja, il, jl, kl, is, js, ks, k, iunitfc, nrow, lp, lp0
    real*8 :: fac
    logical :: full_
    character(len=:), allocatable :: file_, word

    errmsg = ""
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) then
       errmsg = "No force constants available to write"
       return
    end if

    ! The file name: the first word, unless it is one of the options,
    ! in which case phonopy's default name is used.
    lp = 1
    if (fc2_unit_index(sline) /= 0) then
       file_ = "FORCE_CONSTANTS"
    else
       lp0 = lp
       word = lgetword(sline,lp0)
       if (len_trim(word) == 0 .or. equal(word,"full")) then
          file_ = "FORCE_CONSTANTS"
       else
          file_ = getword(sline,lp)
       end if
    end if

    ! options: the units of the file (Ry/bohr^2 by default, phonopy's
    ! units for Quantum ESPRESSO) and whether to write the full array
    iunitfc = 3
    full_ = .false.
    do while (.true.)
       lp0 = lp
       k = fc2_unit_index(sline,lp)
       if (k > 0) then
          iunitfc = k
          cycle
       elseif (k < 0) then
          errmsg = "Unknown force-constant units in UNITS; known units are " // fc2_unitlist
          return
       end if
       word = lgetword(sline,lp)
       if (len_trim(word) == 0) exit
       if (equal(word,"full")) then
          full_ = .true.
       else
          errmsg = "Unknown keyword in WRITE_FC2: " // trim(word)
          return
       end if
    end do
    fac = 1d0 / fc2_unit_factor(iunitfc)

    if (full_) then
       nrow = v%fc2_nsat
    else
       nrow = c%ncel
    end if

    lu = fopen_write(file_,errstop=.false.)
    if (lu < 0) then
       errmsg = "Could not open the force-constant file for writing: " // trim(file_)
       return
    end if
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

    if (verbose) then
       if (full_) then
          write (uout,'("+ Written FC2 file (",A,", full format): ",A)') &
             trim(fc2_unitname(iunitfc)), file_
       else
          write (uout,'("+ Written FC2 file (",A,", compact format): ",A)') &
             trim(fc2_unitname(iunitfc)), file_
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

  !> Record the geometry the force constants were obtained for, so that
  !> calculate_q can refuse to use them after the structure changes.
  subroutine fc2_stamp_geometry(v,c)
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c

    integer :: i

    if (allocated(v%fc2_x)) deallocate(v%fc2_x)
    allocate(v%fc2_x(3,c%ncel))
    do i = 1, c%ncel
       v%fc2_x(:,i) = c%atcel(i)%x
    end do
    v%fc2_m = c%m_x2c

  end subroutine fc2_stamp_geometry

  !> Build the table of shortest images of the vectors going from the
  !> atoms of the cell to the atoms of the supercell in which the force
  !> constants were calculated; the dynamical matrix needs them. For the
  !> pair (ia,js) it finds the images of the vector from cell atom ia to
  !> supercell atom js having the smallest length, over all translations
  !> by the supercell lattice, and stores them in cell fractional
  !> coordinates in v%fc2_svec (see the type for the indexing).
  !> If error, return non-zero errmsg.
  !>
  !> This routine was adapted from phonopy, by A. Togo.
  subroutine fc2_build_svec(v,c,errmsg)
    use tools, only: delaunay_reduction
    use types, only: realloc
    use tools_io, only: string
    class(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: ia, ja, jl, js, ip, k, n1, n2, n3, n4, nc, nsv, ier
    integer :: ncel, nlat, nsat
    real*8 :: scx2c(3,3), smatr(3,3), rmat(3,4), tt(3), y(3), yx(3), r0(3), dmin
    real*8 :: tcand(3,3**4), tcandx(3,3**4), dcand(3**4)
    real*8, allocatable :: sv(:,:)
    logical :: found

    errmsg = ""
    ncel = c%ncel
    nlat = v%fc2_nlat
    nsat = v%fc2_nsat
    smatr = real(v%fc2_smat,8)
    scx2c = matmul(c%m_x2c,transpose(smatr))

    ! the candidate translations: the four Delaunay vectors of the supercell
    ! lattice with coefficients -1, 0 and 1.
    call delaunay_reduction(scx2c,rmat)
    nc = 0
    do n1 = -1, 1
       do n2 = -1, 1
          do n3 = -1, 1
             do n4 = -1, 1
                tt = matmul(scx2c,n1*rmat(:,1) + n2*rmat(:,2) + n3*rmat(:,3) + n4*rmat(:,4))
                found = .false.
                do k = 1, nc
                   found = all(abs(tt - tcand(:,k)) < fc2_epssvec)
                   if (found) exit
                end do
                if (found) cycle
                nc = nc + 1
                tcand(:,nc) = tt
                tcandx(:,nc) = c%c2x(tt)
             end do
          end do
       end do
    end do

    ! allocate; the initial guess for the number of images is four per pair
    if (allocated(v%fc2_sptr)) deallocate(v%fc2_sptr)
    allocate(v%fc2_sptr(ncel*nsat+1),sv(3,4*ncel*nsat),stat=ier)
    if (ier /= 0) then
       errmsg = "Could not allocate the table of shortest images (" //&
          string(nint(52d0*ncel*nsat/1024d0**2)) // " Mb)"
       return
    end if

    nsv = 0
    do js = 1, nsat
       ja = (js-1)/nlat + 1
       jl = js - (ja-1)*nlat
       do ia = 1, ncel
          ip = (js-1)*ncel + ia
          v%fc2_sptr(ip) = nsv + 1

          ! the vector from atom ia to atom js, reduced into the supercell
          y = fc2_scpos(c%atcel(ja)%x - c%atcel(ia)%x,v%fc2_lvec(:,jl),v%fc2_madj,nlat)
          y = y - nint(y)
          r0 = matmul(scx2c,y)
          yx = matmul(y,smatr)

          ! keep the images at the smallest distance
          dmin = huge(1d0)
          do k = 1, nc
             dcand(k) = norm2(r0 + tcand(:,k))
             dmin = min(dmin,dcand(k))
          end do
          do k = 1, nc
             if (dcand(k) > dmin + fc2_epssvec) cycle
             nsv = nsv + 1
             if (nsv > size(sv,2)) call realloc(sv,3,2*nsv)
             sv(:,nsv) = yx + tcandx(:,k)
          end do
       end do
    end do
    v%fc2_sptr(ncel*nsat+1) = nsv + 1

    call realloc(sv,3,nsv)
    call move_alloc(sv,v%fc2_svec)

  end subroutine fc2_build_svec

  !> Check that the force constants stored in v are available and
  !> correspond to the structure c as it is now: same number of atoms,
  !> same lattice.
  subroutine fc2_check_current(v,c,errmsg)
    use tools_io, only: string
    type(vibrations), intent(in) :: v
    type(crystal), intent(in) :: c
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i
    real*8 :: dx(3)

    errmsg = ""
    if (.not.v%hasfc2 .or..not.allocated(v%fc2)) then
       errmsg = "No force constants available; use VIBRATIONS LOAD_FC2 or VIBRATIONS READ_FORCES first"
       return
    end if
    if (v%fc2_ncel /= c%ncel .or..not.allocated(v%fc2_x)) then
       errmsg = "The force constants correspond to a structure with " // string(v%fc2_ncel) //&
          " atoms in the cell but this one has " // string(c%ncel) // "; load them again"
       return
    end if
    if (any(abs(v%fc2_m - c%m_x2c) > fc2_epsgeo)) then
       errmsg = "The lattice has changed since the force constants were obtained; load them again"
       return
    end if
    do i = 1, c%ncel
       dx = v%fc2_x(:,i) - c%atcel(i)%x
       dx = dx - nint(dx)
       if (any(abs(dx) > fc2_epsgeo)) then
          errmsg = "Atom " // string(i) // " has moved since the force constants were obtained; &
             &load them again"
          return
       end if
    end do

  end subroutine fc2_check_current

  !> Calculate frequencies and eigenvectors from the FC2 at a single
  !> q-point (fractional coordinates in reciprocal space). If the
  !> optional arguments freqo or veco are present, return the
  !> frequencies (cm-1, ascending order) in freqo and the whole
  !> (3*ncel,3*ncel) eigenvector matrix in veco (columns = modes), and
  !> leave v alone. Otherwise, add the q-point to v, where the
  !> eigenvectors are reshaped to (3,ncel) per mode. If error, return
  !> non-zero errmsg.
  module subroutine vibrations_calculate_q(v,c,q,errmsg,freqo,veco)
    use tools_math, only: eigherm
    use tools_io, only: string
    use types, only: realloc
    use param, only: atmass, tpi, img
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: q(3)
    character(len=:), allocatable, intent(out) :: errmsg
    real*8, intent(inout), allocatable, optional :: freqo(:)
    complex*16, intent(inout), allocatable, optional :: veco(:,:)

    integer :: ia, ja, jl, js, ip, i, k, ier, ncel, nlat, nfreq
    complex*16 :: phase, dd(3,3)
    real*8, allocatable :: eval(:), sqrtm(:)
    complex*16, allocatable :: dm(:,:)

    errmsg = ""

    ! the force constants have to belong to the structure at hand
    call fc2_check_current(v,c,errmsg)
    if (len_trim(errmsg) > 0) return
    ncel = c%ncel
    nlat = v%fc2_nlat
    nfreq = 3 * ncel

    ! the shortest images, calculated once and kept with the force constants
    if (.not.allocated(v%fc2_svec)) then
       call fc2_build_svec(v,c,errmsg)
       if (len_trim(errmsg) > 0) return
    end if

    ! build the dynamical matrix
    allocate(dm(nfreq,nfreq),eval(nfreq),sqrtm(ncel),stat=ier)
    if (ier /= 0) then
       errmsg = "Could not allocate the dynamical matrix"
       return
    end if
    do ia = 1, ncel
       sqrtm(ia) = sqrt(atmass(c%spc(c%atcel(ia)%is)%z))
    end do
    do ia = 1, ncel
       do ja = 1, ncel
          dd = 0d0
          do jl = 1, nlat
             js = (ja-1)*nlat + jl
             ip = (js-1)*ncel + ia

             ! average the phases of the images at the same distance
             phase = 0d0
             do k = v%fc2_sptr(ip), v%fc2_sptr(ip+1)-1
                phase = phase + exp(img * tpi * dot_product(q,v%fc2_svec(:,k)))
             end do
             dd = dd + v%fc2(:,:,ia,js) * (phase / real(v%fc2_sptr(ip+1)-v%fc2_sptr(ip),8))
          end do
          dm(3*ia-2:3*ia,3*ja-2:3*ja) = dd / (sqrtm(ia) * sqrtm(ja))
       end do
    end do

    ! impose hermiticity and diagonalize; the eigenvectors are skipped
    ! when only the frequencies were asked for (the mesh sampling)
    dm = 0.5d0 * (dm + transpose(conjg(dm)))
    call eigherm(dm,nfreq,eval,ier,vectors=(present(veco).or..not.present(freqo)))
    if (ier /= 0) then
       errmsg = "Error diagonalizing the dynamical matrix at q = " //&
          string(q(1),'f',12,7) // string(q(2),'f',12,7) // string(q(3),'f',12,7)
       return
    end if
    eval = sign(sqrt(abs(eval)) * freqfactor,eval)

    ! return the result, or add the q-point to the vibrations object
    if (present(freqo).or.present(veco)) then
       if (present(freqo)) freqo = eval
       if (present(veco)) veco = dm
    else
       v%nqpt = v%nqpt + 1
       v%nfreq = nfreq
       if (.not.allocated(v%qpt)) then
          allocate(v%qpt(3,v%nqpt),v%freq(nfreq,v%nqpt),v%vec(3,ncel,nfreq,v%nqpt))
       elseif (v%nqpt > size(v%qpt,2)) then
          call realloc(v%qpt,3,2*v%nqpt)
          call realloc(v%freq,nfreq,2*v%nqpt)
          call realloc(v%vec,3,ncel,nfreq,2*v%nqpt)
       end if
       v%qpt(:,v%nqpt) = q
       v%freq(:,v%nqpt) = eval
       do i = 1, nfreq
          v%vec(:,:,i,v%nqpt) = reshape(dm(:,i),(/3,ncel/))
       end do
       v%hasvibs = .true.
    end if

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
  module subroutine vibrations_calculate_thermo(v,t,cutoff,zpe,fvib,svib,cv,nused,ntot,nimag,freqo)
    class(vibrations), intent(inout) :: v
    real*8, intent(in) :: t
    real*8, intent(in) :: cutoff
    real*8, intent(out) :: zpe, fvib, svib, cv
    integer, intent(out) :: nused, ntot, nimag
    real*8, intent(in), optional :: freqo(:,:)

    if (present(freqo)) then
       call thermo_sum(freqo,size(freqo,1),size(freqo,2),t,cutoff,zpe,fvib,svib,cv,nused,ntot,nimag)
    else
       call thermo_sum(v%freq,v%nfreq,v%nqpt,t,cutoff,zpe,fvib,svib,cv,nused,ntot,nimag)
    end if

  end subroutine vibrations_calculate_thermo

  !> Calculate the harmonic thermodynamic properties at temperature t
  !> (K) from the frequencies freq(1:nf,1:nq) (cm^-1), which are a
  !> sampling of the Brillouin zone with equal weights. Returns the
  !> zero-point energy and the vibrational Helmholtz free energy in
  !> kJ/mol, and the entropy and constant-volume heat capacity in
  !> J/K/mol, all per unit cell. Modes at or below the cutoff (cm^-1)
  !> are left out. The cutoff has a floor of thermo_epszero. nimag
  !> counts the modes below -max(cutoff,thermo_epsimag), which are
  !> genuinely imaginary rather than numerical zeros.
  !>
  !> This routine was adapted from phonopy, by A. Togo.
  subroutine thermo_sum(freq,nf,nq,t,cutoff,zpe,fvib,svib,cv,nused,ntot,nimag)
    use param, only: Rgas
    real*8, intent(in) :: freq(:,:)
    integer, intent(in) :: nf, nq
    real*8, intent(in) :: t, cutoff
    real*8, intent(out) :: zpe, fvib, svib, cv
    integer, intent(out) :: nused, ntot, nimag

    integer :: i, j
    real*8 :: nu, x, y, nut, nue, rt, l1mx, nutdiv, ym1, ff, cut, cutimag

    real*8, parameter :: small1 = 50000d0 * cminv_to_K / huge(1d0) ! protection against zerodiv in nu/(kB*T)
    real*8, parameter :: small2 = 0.5d0 * log(huge(1d0)) ! protection against overflow in exp(nu/kB*T)**2
    real*8, parameter :: small3 = 1d0 / sqrt(huge(1d0)) ! protection against 1/(e^(nu/kt)-1)^2 overflow

    zpe = 0d0
    fvib = 0d0
    svib = 0d0
    cv = 0d0
    nused = 0
    ntot = nf * nq
    nimag = 0
    rt = Rgas / 1000d0 * t ! RT in kJ/mol
    cut = max(cutoff,thermo_epszero)
    cutimag = max(cutoff,thermo_epsimag)

    do i = 1, nq
       do j = 1, nf
          ! leave out the modes at or below the cutoff (the acoustic
          ! branches at gamma, and any imaginary mode) without
          ! compensating for them
          nu = freq(j,i)
          if (nu < -cutimag) nimag = nimag + 1
          if (nu <= cut) cycle
          nused = nused + 1
          nut = nu * cminv_to_K      ! frequency in K
          nue = nu * cminv_to_kJmol  ! frequency in kJ/mol

          ! some quantities we will need, with protection
          if (t < small1) then
             x = 0d0
             nutdiv = huge(1d0)
          else
             nutdiv = nut / t
             x = exp(-nutdiv)
          end if
          l1mx = log(1d0 - x)

          ! zero-point energy and free energy
          zpe = zpe + 0.5d0 * nue
          fvib = fvib + 0.5d0 * nue + rt * l1mx

          ! entropy
          if (t > small1) &
             svib = svib + Rgas * (-l1mx + nutdiv * x / (1d0 - x))

          ! heat capacity
          if (t > small1 .and. nutdiv <= small2) then
             y = exp(nutdiv)
             ym1 = y - 1d0
             if (ym1 >= small3) &
                cv = cv + Rgas * nutdiv * nutdiv * y / (ym1 * ym1)
          end if
       end do
    end do

    ! per unit cell: the q-points all have the same weight
    ff = 1d0 / real(max(nq,1),8)
    zpe = zpe * ff
    fvib = fvib * ff
    svib = svib * ff
    cv = cv * ff

  end subroutine thermo_sum

  !> Calculate frequencies (cm^-1) on the uniform mesh nk(3), q = (i -
  !> 1 + qshift)/nk in fractional coordinates of the reciprocal cell.
  !> Returns freq(3*ncel,nk(1)*nk(2)*nk(3)). No frequencies are stored
  !> in v.
  module subroutine vibrations_mesh_freqs(v,c,nk,qshift,freq,errmsg)
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    integer, intent(in) :: nk(3)
    real*8, intent(in) :: qshift(3)
    real*8, allocatable, intent(inout) :: freq(:,:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i1, i2, i3, iq, nq, ierr
    real*8 :: q(3)
    real*8, allocatable :: f1(:)
    character(len=:), allocatable :: errmsg1, errmsg2

    errmsg = ""
    nq = nk(1) * nk(2) * nk(3)
    if (allocated(freq)) deallocate(freq)
    allocate(freq(3*c%ncel,nq),stat=ierr)
    if (ierr /= 0) then
       errmsg = "Could not allocate the frequencies of the mesh"
       return
    end if

    ! the force constants have to belong to the structure at hand; this
    ! is calculate_q's check, and it has to happen here too because
    ! fc2_build_svec below indexes the cell atoms with the force-constant
    ! dimensions
    call fc2_check_current(v,c,errmsg)
    if (len_trim(errmsg) > 0) return

    ! calculate_q builds the shortest-image table the first time it is
    ! called; do it here, or the threads below race to build it
    if (.not.allocated(v%fc2_svec)) then
       call fc2_build_svec(v,c,errmsg)
       if (len_trim(errmsg) > 0) return
    end if

    ! every (i1,i2,i3) writes its own column of freq, so only the error
    ! report needs serializing
    errmsg1 = ""
    !$omp parallel do collapse(3) private(i1,i2,i3,iq,q,f1,errmsg2) schedule(dynamic)
    do i1 = 1, nk(1)
       do i2 = 1, nk(2)
          do i3 = 1, nk(3)
             iq = ((i1-1)*nk(2) + (i2-1))*nk(3) + i3
             q = (real((/i1,i2,i3/)-1,8) + qshift) / real(nk,8)
             call v%calculate_q(c,q,errmsg2,freqo=f1)
             if (len_trim(errmsg2) > 0) then
                !$omp critical (meshfreq)
                errmsg1 = errmsg2
                !$omp end critical (meshfreq)
             else
                freq(:,iq) = f1
             end if
          end do
       end do
    end do
    !$omp end parallel do
    errmsg = errmsg1

  end subroutine vibrations_mesh_freqs

  !> Write the phonon density of states obtained from the frequencies
  !> freq(:,1:nq) (cm^-1) by Gaussian smearing of width sigma (cm^-1),
  !> on npts points spanning the sampled range with a margin. The
  !> integral of the DOS is the number of modes of the cell.
  module subroutine vibrations_write_dos(v,c,file,nk,qshift,sigma,npts,verbose,errmsg,freqo)
    class(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: nk(3)
    real*8, intent(in) :: qshift(3)
    real*8, intent(in) :: sigma
    integer, intent(in) :: npts
    logical, intent(in) :: verbose
    character(len=:), allocatable, intent(out) :: errmsg
    real*8, intent(in), optional :: freqo(:,:)

    if (present(freqo)) then
       call dos_run(c,freqo,size(freqo,1),size(freqo,2),file,nk,qshift,sigma,npts,verbose,errmsg)
    else
       call dos_run(c,v%freq,v%nfreq,v%nqpt,file,nk,qshift,sigma,npts,verbose,errmsg)
    end if

  end subroutine vibrations_write_dos

  !> Write the phonon density of states of the frequencies
  !> freq(1:nf,1:nq) (cm^-1) to a file, on npts points. A positive
  !> sigma smears every mode with a Gaussian of that width (cm^-1);
  !> otherwise the linear tetrahedron method is used, which needs the
  !> frequencies to be those of the uniform mesh nk(3), in the order
  !> vibrations_mesh_freqs generates them. The density is normalized so
  !> that its integral is the number of modes of the cell.
  subroutine dos_run(c,freq,nf,nq,file,nk,qshift,sigma,npts,verbose,errmsg)
    use tools_io, only: fopen_write, fclose, uout, string, ioj_right
    type(crystal), intent(in) :: c
    real*8, intent(in) :: freq(:,:)
    integer, intent(in) :: nf, nq
    character*(*), intent(in) :: file
    integer, intent(in) :: nk(3)
    real*8, intent(in) :: qshift(3)
    real*8, intent(in) :: sigma
    integer, intent(in) :: npts
    logical, intent(in) :: verbose
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: k, lu, ierr
    real*8 :: fmin, fmax, step, sums
    logical :: dotetra
    real*8, allocatable :: dos(:)

    errmsg = ""
    if (nq < 1 .or. nf < 1) then
       errmsg = "No frequencies available for the phonon DOS"
       return
    end if
    if (npts < 2) then
       errmsg = "The DOS needs at least two points"
       return
    end if

    ! the tetrahedron method needs the mesh and its connectivity; a
    ! positive smearing width asks for Gaussians instead
    dotetra = (sigma <= 0d0)
    if (dotetra) then
       if (any(nk <= 0) .or. nk(1)*nk(2)*nk(3) /= nq) then
          errmsg = "The tetrahedron method needs the frequencies of a uniform mesh"
          return
       end if
    end if

    ! the range of the output grid
    fmin = minval(freq(1:nf,1:nq))
    fmax = maxval(freq(1:nf,1:nq))
    if (dotetra) then
       step = 0.02d0 * (fmax - fmin) + 1d-3
    else
       step = 5d0 * sigma
       fmin = min(fmin,0d0)
    end if
    fmin = fmin - step
    fmax = fmax + step
    step = (fmax - fmin) / real(npts-1,8)

    allocate(dos(npts),stat=ierr)
    if (ierr /= 0) then
       errmsg = "Could not allocate the phonon DOS"
       return
    end if
    dos = 0d0

    if (dotetra) then
       call dos_tetrahedra(c,freq,nf,nq,nk,fmin,step,npts,dos)
    else
       call dos_gaussian(freq,nf,nq,sigma,fmin,step,npts,dos)
    end if

    ! write it out
    lu = fopen_write(file,errstop=.false.)
    if (lu < 0) then
       errmsg = "Could not open the phonon DOS file for writing: " // trim(file)
       return
    end if
    sums = sum(dos) * step
    write (lu,'("# Phonon density of states calculated by critic2")')
    if (all(abs(qshift) < 1d-10)) then
       write (lu,'("# ",A," q-points (",A,"x",A,"x",A," gamma-centered mesh), ",A," modes per cell")') &
          string(nq), (string(nk(k)),k=1,3), string(nf)
    else
       write (lu,'("# ",A," q-points (",A,"x",A,"x",A," mesh shifted by ",3(A," "),&
          &"of a mesh step), ",A," modes per cell")') string(nq), (string(nk(k)),k=1,3),&
          (string(qshift(k),'f',8,4),k=1,3), string(nf)
    end if
    if (dotetra) then
       write (lu,'("# linear tetrahedron method")')
    else
       write (lu,'("# Gaussian smearing of ",A," cm^-1")') string(sigma,'f',10,4)
    end if
    write (lu,'("# integral of the density = ",A," (should be the ",A," modes of the cell)")') &
       string(sums,'f',12,6), string(nf)
    write (lu,'("# column 1: frequency (cm^-1)   column 2: density of states (1/cm^-1, per cell)")')
    do k = 1, npts
       write (lu,'(2(A," "))') string(fmin + real(k-1,8)*step,'f',16,6,ioj_right),&
          string(dos(k),'f',18,10,ioj_right)
    end do
    call fclose(lu)

    if (verbose) then
       if (dotetra) then
          write (uout,'("+ Phonon DOS (linear tetrahedron method) written to: ",A)') trim(file)
       else
          write (uout,'("+ Phonon DOS (Gaussian smearing, sigma = ",A," cm^-1) written to: ",A)') &
             string(sigma,'f',10,4), trim(file)
       end if
       write (uout,'("  Integral of the density: ",A," (",A," modes per cell)")') &
          string(sums,'f',12,6), string(nf)
    end if

  end subroutine dos_run

  !> Accumulate the density of states of the frequencies
  !> freq(1:nf,1:nq) (cm^-1) on the npts-point grid that starts at
  !> fmin with spacing step, by Gaussian smearing of width sigma. Bin
  !> k spans half a step on either side of its frequency. The
  !> normalization is one state per mode and per cell.
  subroutine dos_gaussian(freq,nf,nq,sigma,fmin,step,npts,dos)
    real*8, intent(in) :: freq(:,:)
    integer, intent(in) :: nf, nq
    real*8, intent(in) :: sigma, fmin, step
    integer, intent(in) :: npts
    real*8, intent(inout) :: dos(npts)

    integer :: i, j, k, k0, k1
    real*8 :: fac, s2, elo, ehi

    fac = 0.5d0 / (real(nq,8) * step)
    s2 = sigma * sqrt(2d0)
    do i = 1, nq
       do j = 1, nf
          k0 = max(int((freq(j,i) - 5d0*sigma - fmin) / step),1)
          k1 = min(int((freq(j,i) + 5d0*sigma - fmin) / step) + 2,npts)
          elo = erf((fmin + (real(k0,8)-1.5d0)*step - freq(j,i)) / s2)
          do k = k0, k1
             ehi = erf((fmin + (real(k,8)-0.5d0)*step - freq(j,i)) / s2)
             dos(k) = dos(k) + fac * (ehi - elo)
             elo = ehi
          end do
       end do
    end do

  end subroutine dos_gaussian

  !> Accumulate the density of states of the frequencies freq(1:nf,1:nq)
  !> (cm^-1) of the uniform mesh nk(3) with the linear tetrahedron
  !> method (Lehmann and Taut; Jepsen and Andersen; Bloechl, Jepsen and
  !> Andersen, Phys. Rev. B 49, 16223 (1994), appendix A). Each cell of
  !> the mesh is cut into six tetrahedra sharing the shortest of its
  !> four body diagonals, the frequency is taken to be linear inside
  !> each of them.
  subroutine dos_tetrahedra(c,freq,nf,nq,nk,fmin,step,npts,dos)
    type(crystal), intent(in) :: c
    real*8, intent(in) :: freq(:,:)
    integer, intent(in) :: nf, nq, nk(3)
    real*8, intent(in) :: fmin, step
    integer, intent(in) :: npts
    real*8, intent(inout) :: dos(npts)

    ! the six tetrahedra of a cell are the six ways of going from one
    ! end of the main diagonal to the other along the cell edges
    integer, parameter :: ntetra = 6
    integer, parameter :: perm(3,ntetra) = reshape((/&
       1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1/),(/3,ntetra/))
    real*8, parameter :: epsdeg = 1d-8 ! degenerate corner frequencies (cm^-1)

    integer :: i, j, k, ic, it, ib, j1, j2, j3, id, k0, k1
    integer :: off(3,4,ntetra), flip(3), iq(4), ix(3)
    real*8 :: am(3,3), dg(3), dlen, dmin, x(3)
    real*8 :: e(4), om, nlo, nhi, fac

    ! the lattice vectors of one cell of the mesh, in Cartesian
    ! reciprocal coordinates
    do k = 1, 3
       x = 0d0
       x(k) = 1d0 / real(nk(k),8)
       am(:,k) = c%rx2rc(x)
    end do

    ! the main diagonal: the shortest of the four, so that the
    ! tetrahedra are as regular as the mesh allows
    id = 0
    dmin = huge(1d0)
    do i = 0, 3
       dg = am(:,1) + am(:,2) + am(:,3)
       if (i > 0) dg = dg - 2d0 * am(:,i)
       dlen = dot_product(dg,dg)
       if (dlen < dmin) then
          dmin = dlen
          id = i
       end if
    end do
    flip = 0
    if (id > 0) flip(id) = 1

    ! the corner offsets of the six tetrahedra, in the frame where the
    ! main diagonal runs from (0,0,0) to (1,1,1), then reflected onto
    ! the diagonal actually chosen
    do it = 1, ntetra
       off(:,1,it) = 0
       off(:,2,it) = 0
       off(perm(1,it),2,it) = 1
       off(:,3,it) = off(:,2,it)
       off(perm(2,it),3,it) = 1
       off(:,4,it) = 1
       do ic = 1, 4
          do k = 1, 3
             if (flip(k) == 1) off(k,ic,it) = 1 - off(k,ic,it)
          end do
       end do
    end do

    ! every tetrahedron is one sixth of one cell of the mesh
    fac = 1d0 / (6d0 * real(nq,8))

    do j1 = 0, nk(1)-1
       do j2 = 0, nk(2)-1
          do j3 = 0, nk(3)-1
             do it = 1, ntetra
                ! the four corners, wrapped around the mesh
                do ic = 1, 4
                   ix(1) = modulo(j1 + off(1,ic,it),nk(1))
                   ix(2) = modulo(j2 + off(2,ic,it),nk(2))
                   ix(3) = modulo(j3 + off(3,ic,it),nk(3))
                   iq(ic) = (ix(1)*nk(2) + ix(2))*nk(3) + ix(3) + 1
                end do

                do ib = 1, nf
                   ! the corner frequencies, sorted
                   do ic = 1, 4
                      e(ic) = freq(ib,iq(ic))
                   end do
                   do i = 2, 4
                      om = e(i)
                      j = i - 1
                      do while (j >= 1)
                         if (e(j) <= om) exit
                         e(j+1) = e(j)
                         j = j - 1
                      end do
                      e(j+1) = om
                   end do

                   ! a completely flat tetrahedron is a delta function
                   if (e(4) - e(1) < epsdeg) then
                      k = nint((e(1) - fmin) / step) + 1
                      if (k >= 1 .and. k <= npts) dos(k) = dos(k) + fac / step
                      cycle
                   end if

                   ! Add the exact number of states the tetrahedron puts
                   ! in each bin, so that the sum rule holds whatever the
                   ! spacing of the output grid: bin k spans half a step
                   ! on either side of its frequency.
                   k0 = max(int((e(1) - fmin) / step),1)
                   k1 = min(int((e(4) - fmin) / step) + 2,npts)
                   nlo = tetra_nstates(fmin + (real(k0,8)-1.5d0)*step,e)
                   do k = k0, k1
                      nhi = tetra_nstates(fmin + (real(k,8)-0.5d0)*step,e)
                      dos(k) = dos(k) + fac * (nhi - nlo) / step
                      nlo = nhi
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine dos_tetrahedra

  !> Fraction of a tetrahedron whose frequency is below om, given its
  !> four corner frequencies e(4) sorted in ascending order and not all
  !> equal.
  pure function tetra_nstates(om,e) result(n)
    real*8, intent(in) :: om, e(4)
    real*8 :: n

    if (om <= e(1)) then
       n = 0d0
    elseif (om < e(2)) then
       ! e1 < om < e2, so e2-e1 > 0 and the other differences are larger
       n = (om - e(1))**3 / ((e(2)-e(1)) * (e(3)-e(1)) * (e(4)-e(1)))
    elseif (om < e(3)) then
       ! e2 <= om < e3, so e3-e2 > 0 and e3-e1, e4-e2, e4-e1 are larger
       n = (om - e(2))**2 / ((e(4)-e(2)) * (e(3)-e(2))) +&
          (om - e(1)) * (e(4) - om) * (om - e(2)) / ((e(4)-e(1)) * (e(4)-e(2)) * (e(3)-e(2))) +&
          (om - e(1))**2 * (e(3) - om) / ((e(4)-e(1)) * (e(3)-e(1)) * (e(3)-e(2)))
    elseif (om < e(4)) then
       ! e3 <= om < e4, so e4-e3 > 0 and e4-e2, e4-e1 are larger
       n = 1d0 - (e(4) - om)**3 / ((e(4)-e(1)) * (e(4)-e(2)) * (e(4)-e(3)))
    else
       n = 1d0
    end if

  end function tetra_nstates

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

    character(len=:), allocatable :: errmsg
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
          call v%calculate_q(c,(/0d0,0d0,0d0/),errmsg)
          if (len_trim(errmsg) > 0) return
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
  !> optional supercell specification (one integer for an n x n x n
  !> supercell, three for a diagonal one, nine for a general one in the
  !> NEWCELL order, or in phonopy's DIM order if the keyword PHONOPY is
  !> also given; or the name of a file containing the supercell;
  !> default is no supercell), and the optional keyword ASR to enforce
  !> the acoustic sum rule after reading. If error, return non-zero
  !> errmsg.
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
    type(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file, sline
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: word, line, scfile, file_
    integer :: lu, lp, lp0, i, j, ia, il, ip, jp, is, js, iap, ilp, irow, idum
    integer :: iunitfc, iun, ndim, nlat, nsat, nrow, n1, n2, ierr
    integer :: smat(3,3), madj(3,3), idim(9)
    real*8 :: fc2factor, rmat(3,3), dev
    real*8 :: t(3), xsc(3), dx(3)
    logical :: iscompact, haveseed, doasr, found, phonopydim, flipped
    integer, allocatable :: lvec(:,:), lkey(:,:), p2s(:), s2row(:), perm(:)
    logical, allocatable :: coldone(:)
    real*8, allocatable :: fcrow(:,:,:,:)
    type(crystal) :: sc
    type(crystalseed) :: seed

    ! initialize
    lu = -1
    scfile = ""
    if (c%ismolecule) then
       errmsg = "FORCE_CONSTANTS (phonopy) files can only be read for crystals"
       return
    end if

    ! The file name. It is empty when the FC2 keyword was used: then
    ! the first word of the option string is the file name, unless it
    ! is one of the options, in which case phonopy's default name is
    ! taken.
    lp = 1
    file_ = file
    if (len_trim(file_) == 0) then
       if (fc2_unit_index(sline) /= 0) then
          file_ = "FORCE_CONSTANTS"
       else
          lp0 = lp
          word = lgetword(sline,lp0)
          if (len_trim(word) == 0 .or. equal(word,"phonopy") .or. equal(word,"asr") .or.&
             equal(word,"acoustic") .or. equal(word,"acoustic_sum_rules")) then
             file_ = "FORCE_CONSTANTS"
          else
             lp0 = lp
             if (isinteger(idum,sline,lp0)) then
                file_ = "FORCE_CONSTANTS"
             else
                file_ = getword(sline,lp)
             end if
          end if
       end if
    end if
    errmsg = "Error reading FORCE_CONSTANTS file: " // trim(file_)

    ! supercell specification and options, in any order
    ndim = 0
    iunitfc = 0
    haveseed = .false.
    doasr = .false.
    phonopydim = .false.
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

       ! the units of the file: mandatory, the file does not carry them
       iun = fc2_unit_index(sline,lp)
       if (iun > 0) then
          iunitfc = iun
          cycle
       elseif (iun < 0) then
          errmsg = "Unknown force-constant units in UNITS; known units are " // fc2_unitlist
          goto 999
       end if

       word = lgetword(sline,lp)
       if (len_trim(word) == 0) exit
       if (equal(word,"asr").or.equal(word,"acoustic").or.equal(word,"acoustic_sum_rules")) then
          doasr = .true.
       elseif (equal(word,"phonopy")) then
          phonopydim = .true.
       elseif (haveseed) then
          errmsg = "Unknown keyword in VIBRATIONS LOAD (FORCE_CONSTANTS): " // trim(word)
          goto 999
       else
          ! the supercell, read from a structure file
          lp = lp0
          scfile = getword(sline,lp)
          haveseed = .true.
       end if
    end do
    if (iunitfc == 0) then
       errmsg = "The units of the force constants are required: give the code that wrote the file &
          &(vasp, qe, aims, alamode, ...) or UNITS followed by the units themselves (" //&
          fc2_unitlist // ")"
       goto 999
    end if
    fc2factor = fc2_unit_factor(iunitfc)

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
       call fc2_smat_from_cell(c,scfile,smat,errmsg,ti,sco=sc,seedo=seed)
       if (len_trim(errmsg) > 0) goto 999
    elseif (ndim > 0) then
       call supercell_matrix_from_ints(ndim,idim,phonopydim,smat,flipped,errmsg)
       if (len_trim(errmsg) > 0) goto 999
       if (flipped) then
          errmsg = "The supercell matrix has negative determinant; phonopy requires a positive one"
          goto 999
       end if
    end if

    ! lattice points of the supercell, in phonopy order
    call fc2_lattice_points(smat,nlat,madj,lvec,lkey,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    nsat = c%ncel * nlat

    ! open the file and read the header
    lu = fopen_read(file_,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file_)
       goto 999
    end if
    errmsg = "Error reading the header of FORCE_CONSTANTS file: " // trim(file_)
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
       errmsg = "Inconsistent number of atoms in " // trim(file_) // ": the file has (" // string(n1) //&
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
    errmsg = "Error reading the force constants from file: " // trim(file_)
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
       do ia = 1, c%ncel
          do il = 1, nlat
             is = (ia-1)*nlat + il
             xsc = fc2_scpos(c%atcel(ia)%x,lvec(:,il),madj,nlat)
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
    v%fc2_file = file_
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
    call fc2_stamp_geometry(v,c)

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
  !> matrix smat. Returns the number of lattice points (nlat = det(smat)), the
  !> integer adjugate of smat (madj, with smat*madj = nlat*I), the
  !> lattice points in cell fractional coordinates (lvec) and their
  !> residue keys modulo the supercell lattice (lkey). Two lattice
  !> points are equivalent modulo the supercell lattice if and only if
  !> they have the same key.
  subroutine fc2_lattice_points(smat,nlat,madj,lvec,lkey,errmsg)
    use tools_io, only: string
    use tools_math, only: idet3, iadj3
    integer, intent(in) :: smat(3,3)
    integer, intent(out) :: nlat
    integer, intent(out) :: madj(3,3)
    integer, allocatable, intent(inout) :: lvec(:,:), lkey(:,:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, k, ia, ib, ic, nl
    integer :: cmin(3), cmax(3), frame(3), vv(3), ll(3), kk(3)
    logical :: found

    errmsg = ""

    ! determinant and adjugate
    nlat = idet3(smat)
    if (nlat <= 0) then
       errmsg = "The supercell matrix must have positive determinant (found " // string(nlat) // ")"
       return
    end if
    madj = iadj3(smat)

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

  !> Index (in fc2_unitname) of the units of a force-constant file,
  !> named at position lp of sline: either the code that wrote the file
  !> (vasp, qe, aims, ...), the units themselves (eV/ang^2,
  !> Ry/bohr^2, ...), or the keyword UNITS followed by either of the
  !> two. Returns zero if the word is neither, and -1 if UNITS is
  !> followed by something unknown. If lp is present, the search starts
  !> there and lp is advanced past the words consumed (left untouched
  !> if the result is zero).
  function fc2_unit_index(sline,lp) result(iunit)
    use tools_io, only: lgetword, equal, lower
    character*(*), intent(in) :: sline
    integer, intent(inout), optional :: lp
    integer :: iunit

    integer :: i, lp_
    logical :: haveunits
    character(len=:), allocatable :: word

    iunit = 0
    lp_ = 1
    if (present(lp)) lp_ = lp
    word = lgetword(sline,lp_)
    haveunits = equal(word,"units")
    if (haveunits) word = lgetword(sline,lp_)

    do i = 1, fc2_ngen
       if (equal(word,trim(fc2_genname(i)))) then
          iunit = fc2_genunit(i)
          if (present(lp)) lp = lp_
          return
       end if
    end do
    do i = 1, fc2_nunit
       if (equal(word,lower(trim(fc2_unitname(i))))) then
          iunit = i
          if (present(lp)) lp = lp_
          return
       end if
    end do
    if (haveunits) then
       iunit = -1
       if (present(lp)) lp = lp_
    end if

  end function fc2_unit_index

  !> Factor that converts force constants in the units iunit (an index
  !> into fc2_unitname) to the internal units (Hartree/bohr^2).
  function fc2_unit_factor(iunit) result(fac)
    use param, only: hartoev, bohrtoa
    integer, intent(in) :: iunit
    real*8 :: fac

    select case (iunit)
    case (1)
       fac = bohrtoa**2 / hartoev ! eV/ang^2
    case (2)
       fac = bohrtoa / hartoev ! eV/(ang*bohr)
    case (3)
       fac = 0.5d0 ! Ry/bohr^2
    case (4)
       fac = 0.5d-3 ! mRy/bohr^2
    case (5)
       fac = 1d0 ! Hartree/bohr^2
    case (6)
       fac = bohrtoa ! Hartree/(ang*bohr)
    case default
       fac = 1d0
    end select

  end function fc2_unit_factor

  !> Position, in fractional coordinates of the supercell, of the image
  !> of the cell atom at x (cell fractional coordinates) at the lattice
  !> point lv. madj and nlat are the adjugate and determinant of the
  !> supercell matrix (see fc2_lattice_points), so madj/nlat is its
  !> exact inverse. This is the position convention shared by
  !> create_displacements and the FORCE_CONSTANTS reader.
  pure function fc2_scpos(x,lv,madj,nlat) result(xs)
    real*8, intent(in) :: x(3)
    integer, intent(in) :: lv(3)
    integer, intent(in) :: madj(3,3)
    integer, intent(in) :: nlat
    real*8 :: xs(3)

    xs = matmul(x + real(lv,8),real(madj,8)) / real(nlat,8)

  end function fc2_scpos

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
