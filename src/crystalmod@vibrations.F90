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

  !xx! private procedures
  ! subroutine vibrations_detect_format(file,ivformat)
  ! subroutine read_matdyn_modes(v,c,file,ivformat,errmsg,ti)
  ! subroutine read_qe_dyn(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_ascii(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_yaml(v,c,file,errmsg,ti)
  ! subroutine read_phonopy_hdf5(v,c,file,errmsg)
  ! subroutine read_phonopy_fc2(v,c,file,sline,errmsg,ti)
  ! subroutine read_crystal_out(v,c,file,errmsg,ti)
  ! subroutine read_gaussian_log(v,c,file,errmsg,ti)
  ! subroutine read_gaussian_fchk(v,c,file,errmsg,ti)

contains

  !> Terminate a vibrations object
  module subroutine vibrations_end(v)
    use param, only: ivformat_unknown
    class(vibrations), intent(inout) :: v

    v%isinit = .false.
    v%hasfc2 = .false.
    v%file = ""
    v%ivformat = ivformat_unknown
    v%nqpt = 0
    if (allocated(v%fc2)) deallocate(v%fc2)
    if (allocated(v%qpt)) deallocate(v%qpt)
    if (allocated(v%freq)) deallocate(v%freq)
    if (allocated(v%vec)) deallocate(v%vec)
    v%nfreq = 0

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
       ivformat_phonopy_fc2
    use types, only: realloc
    class(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file, sline
    integer, intent(in) :: ivformat
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: ivf

    ! detect the format
    if (ivformat == ivformat_unknown) then
       call vibrations_detect_format(file,ivf)
       if (ivf == ivformat_unknown) then
          errmsg = "Unknown vibration file format: " // trim(file)
          return
       end if
    else
       ivf = ivformat
    end if

    ! initialize
    call v%end()

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
    else
       errmsg = "Unknown vibration file format: " // trim(file)
       return
    end if
    if (len(errmsg) > 0) return
    if (v%isinit) then
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

  !> Print information about the stored FC2 to the standard output.
  !> The structural info comes from crystal structures c.
  module subroutine vibrations_print_fc2(v,c,disteps,fc2eps,environ)
    use tools_io, only: uout, ioj_center, ioj_right, string
    use global, only: iunit, iunitname0, dunit0
    use tools, only: mergesort
    use param, only: icrd_crys
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in), optional :: disteps, fc2eps
    logical, intent(in), optional :: environ

    integer :: i, j, k, idum, jdum, id
    integer :: npair, cidx, nidx
    real*8 :: maxdisteps, dist, disteps_, fc2eps_, up2d, xx(3), norm
    logical :: isintra, environ_
    real*8, allocatable :: dista(:)
    integer, allocatable :: idx(:), i2(:,:), pair(:,:)
    logical, allocatable :: isdone(:), isdonej(:)
    integer :: nat
    integer, allocatable :: eid(:)
    real*8, allocatable :: distr(:)
    integer, allocatable :: lvec(:,:)

    ! header and checks
    write (uout,'("+ PRINT 2nd-order force constant information (FC2)")')
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) then
       write (uout,'("!! No FC2 info has been loaded, exiting.")')
       return
    end if

    ! process options
    disteps_ = huge(1d0)
    fc2eps_ = 0d0
    environ_ = .false.
    if (present(disteps)) then
       disteps_ = disteps
       if (disteps_ < huge(1d0)) &
          write (uout,'("# Cutoff distance for pairs = ",A," ",A)') &
             string(disteps_,'f',decimal=5), iunitname0(iunit)
    end if
    if (present(fc2eps)) then
       fc2eps_ = fc2eps
       if (fc2eps_ > 0d0) &
          write (uout,'("# Cutoff for FC2 values = ",A)') string(fc2eps_,'e',decimal=5)
    end if
    if (present(environ)) environ_ = environ

    ! allocate information about atom pairs, and sort by distance
    npair = c%ncel*(c%ncel + 1) / 2
    allocate(dista(npair),idx(npair),i2(2,npair))
    allocate(pair(2,npair))
    k = 0
    do i = 1, c%ncel
       do j = i, c%ncel
          k = k + 1
          dista(k) = c%eql_distance(c%atcel(i)%x, c%atcel(j)%x)
          idx(k) = k
          i2(1,k) = i
          i2(2,k) = j
       end do
    end do
    maxdisteps = maxval(dista)
    call mergesort(dista,idx,1,npair)

    if (.not.environ) then
       ! print the whole fc2 matrix
       write (uout,'("# Id1 At1  Id2 At2  dist(bohr) isintra norm          fc2xx           fc2xy           fc2xz&
          &           fc2yx           fc2yy           fc2yz           fc2zx           fc2zy           fc2zz")')
       k = 0
       do idum = 1, c%ncel
          do jdum = idum, c%ncel
             k = k + 1
             i = i2(1,idx(k))
             j = i2(2,idx(k))
             dist = dista(idx(k))
             isintra = (c%idatcelmol(1, i) == c%idatcelmol(1, j))
             if (dist <= disteps_ .and. any(abs(v%fc2(:,:,i,j)) >= fc2eps_)) then
                norm = sqrt(sum(v%fc2(:,:,i,j)**2))
                write (uout,'(99(A," "))') string(i, 5, ioj_center), string(c%spc(c%atcel(i)%is)%name,2),&
                   string(j, 5, ioj_center), string(c%spc(c%atcel(j)%is)%name,2),&
                   string(dist*dunit0(iunit),'f',12,6,ioj_right), string(isintra), string(norm,'e',15,8,ioj_right),&
                   string(v%fc2(1,1,i,j),'e',15,8,ioj_right), string(v%fc2(1,2,i,j),'e',15,8,ioj_right),&
                   string(v%fc2(1,3,i,j),'e',15,8,ioj_right), string(v%fc2(2,1,i,j),'e',15,8,ioj_right),&
                   string(v%fc2(2,2,i,j),'e',15,8,ioj_right), string(v%fc2(2,3,i,j),'e',15,8,ioj_right),&
                   string(v%fc2(3,1,i,j),'e',15,8,ioj_right), string(v%fc2(3,2,i,j),'e',15,8,ioj_right),&
                   string(v%fc2(3,3,i,j),'e',15,8,ioj_right)
             end if
          end do
       end do
    else
       ! print the FC2 in environments of nneq atoms
       allocate(isdone(c%nneq),isdonej(c%ncel))
       isdone = .false.
       up2d = min(disteps_,maxdisteps+1d-4)

       ! run over non-equivalent atoms, i is the complete list index
       write (uout,*)
       do i = 1, c%ncel
          id = c%atcel(i)%idx
          if (isdone(id)) cycle
          isdone(id) = .true.

          ! find the atomic environment
          call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.true.,nat,eid,distr,lvec,up2d=up2d)

          write (uout,'("+ Environment of atom ",A," (spc=",A,", nid=",A,") at ",3(A," "))') &
             string(i), string(c%spc(c%at(id)%is)%name), string(id), &
             (string(c%atcel(i)%x(j),'f',length=10,decimal=6),j=1,3)
          write (uout,'("# Up to distance (",A,"): ",A)') iunitname0(iunit), string(up2d*dunit0(iunit),'f',length=10,decimal=6)
          write (uout,'("# Number of atoms in the environment: ",A)') string(nat)
          write (uout,'("# nid = non-equivalent list atomic ID. id = complete list ID plus lattice vector (lvec).")')
          write (uout,'("# name = atomic name. dist = distance. norm = Frobenius norm FC2. xx,... = FC2 components")')
          write (uout,'("#nid   id      lvec     name  dist(",A,")  isintra FC2: norm        xx              &
             &xy              xz              yx              yy              yz              &
             &zx              zy              zz")') iunitname0(iunit)
          isdonej = .false.
          do j = 1, nat
             cidx = eid(j)
             if (isdonej(cidx)) cycle
             isdonej(cidx) = .true.

             nidx = c%atcel(cidx)%idx
             if (c%ismolecule) then
                xx = (c%atcel(cidx)%r + c%molx0) * dunit0(iunit)
             else
                xx = c%atcel(cidx)%x + lvec(:,j)
             end if

             isintra = (c%idatcelmol(1,i) == c%idatcelmol(1,cidx))
             norm = sqrt(sum(v%fc2(:,:,i,cidx)**2))
             write (uout,'("  ",2(A," "),"(",A," ",A," ",A,")",99(" ",A))') string(nidx,4,ioj_center), string(cidx,4,ioj_center),&
                (string(lvec(k,j),2,ioj_right),k=1,3), string(c%spc(c%atcel(cidx)%is)%name,7,ioj_center),&
                string(distr(j)*dunit0(iunit),'f',12,6,4), string(isintra), string(norm,'e',15,8,ioj_right),&
                string(v%fc2(1,1,i,cidx),'e',15,8,ioj_right), string(v%fc2(1,2,i,cidx),'e',15,8,ioj_right),&
                string(v%fc2(1,3,i,cidx),'e',15,8,ioj_right), string(v%fc2(2,1,i,cidx),'e',15,8,ioj_right),&
                string(v%fc2(2,2,i,cidx),'e',15,8,ioj_right), string(v%fc2(2,3,i,cidx),'e',15,8,ioj_right),&
                string(v%fc2(3,1,i,cidx),'e',15,8,ioj_right), string(v%fc2(3,2,i,cidx),'e',15,8,ioj_right),&
                string(v%fc2(3,3,i,cidx),'e',15,8,ioj_right)
          end do
          write (uout,*)
       end do
    end if

  end subroutine vibrations_print_fc2

  !xx! private procedures

  !> Detect the format for a file containing molecular or
  !> crystal vibrations.
  subroutine vibrations_detect_format(file,ivformat)
    use param, only: ivformat_unknown, ivformat_matdynmodes, ivformat_matdyneig,&
       ivformat_qedyn, ivformat_phonopy_ascii, ivformat_phonopy_yaml,&
       ivformat_phonopy_hdf5, ivformat_crystal_out, ivformat_gaussian_log,&
       ivformat_gaussian_fchk, ivformat_phonopy_fc2, dirsep
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
    use param, only: atmass, isformat_qeout, ivformat_matdynmodes
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: ivformat
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: ivf, lu, nline, nline1, nat, idx
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
    if (c%isformat /= isformat_qeout) then
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
    v%ivformat = ivf

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

    ! wrap up
    errmsg = ""
    v%isinit = .true.
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
    use param, only: atmass, isformat_qeout, ivformat_qedyn
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
    if (c%isformat /= isformat_qeout) then
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
    v%isinit = .true.
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
    v%isinit = .true.
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
    use tools_io, only: fopen_read, fclose, getline_raw
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
             if (.not.ok) goto 999
             do j = 1, 3
                ok = getline_raw(lu,line,.false.)
                if (.not.ok) goto 999
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
    v%isinit = .true.
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
    v%isinit = .true.
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
  subroutine read_phonopy_fc2(v,c,file,sline,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword, equal
    use param, only: ivformat_phonopy_fc2
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file, sline
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: word
    integer :: lp, lu, nat, mat, i, j, idum, jdum
    real*8 :: fc2factor

    ! initialize
    errmsg = "Error reading FORCE_CONSTANTS file: " // trim(file)
    lu = -1

    ! interpret format and set conversion factor
    lp = 1
    word = lgetword(sline,lp)
    if (equal(word,"qe")) then
       fc2factor = 0.5d0
    elseif (equal(word,"vasp")) then
       write (*,*) "fixme: vasp in FC2 reader"
       stop 1
    elseif (equal(word,"fhiaims")) then
       write (*,*) "fixme: fhiaims in FC2 reader"
       stop 1
    else
       errmsg = "A generator keyword (qe,vasp,fhiaims,...) is required to read FORCE_CONSTANTS"
       goto 999
    end if

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! read the number of atoms from the fc2 file
    read(lu,*,err=999,end=999) nat, mat
    if (c%ncel /= nat .or. c%ncel /= mat ) then
       errmsg = 'Inconsistent atom number in FORCE_CONSTANTS file'
       goto 999
    end if

    ! prepare the fc2 matrix
    errmsg = "Error reading force constants from FORCE_CONSTANTS file"
    if (allocated(v%fc2)) deallocate(v%fc2)
    allocate(v%fc2(3,3,nat,nat))
    do idum = 1, c%ncel
       do jdum = 1, c%ncel
          read (lu,*,err=999,end=999) i, j
          read (lu,*,err=999,end=999) v%fc2(1,1,i,j), v%fc2(2,1,i,j), v%fc2(3,1,i,j)
          read (lu,*,err=999,end=999) v%fc2(1,2,i,j), v%fc2(2,2,i,j), v%fc2(3,2,i,j)
          read (lu,*,err=999,end=999) v%fc2(1,3,i,j), v%fc2(2,3,i,j), v%fc2(3,3,i,j)
       end do
    end do
    v%fc2 = v%fc2 * fc2factor

    ! wrap up
    v%file = file
    v%ivformat = ivformat_phonopy_fc2
    v%hasfc2 = .true.
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_phonopy_fc2

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
    ! integer :: jfreq, iqpt ! checking normalization
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

    ! wrap up
    v%isinit = .true.
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
    use param, only: ivformat_gaussian_log, atmass
    type(vibrations), intent(inout) :: v
    type(crystal), intent(in) :: c
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok, okl(3)
    character(len=:), allocatable :: line
    integer :: lu, i, j, nread, iz, nn, lp, idum(2)
    real*8 :: xdum, xread(9)
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

    nread = 0
    do while (.true.)
       ! advance to the next frequencies header
       ok = .false.
       do while (getline_raw(lu,line,.false.))
          if (len(line) >= 15) then
             if (index(line," Frequencies --") > 0) then
                ok = .true.
                exit
             end if
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
       do i = 1, 4
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
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
    v%isinit = .true.
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
    v%isinit = .true.
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_gaussian_fchk

end submodule vibrationsmod
