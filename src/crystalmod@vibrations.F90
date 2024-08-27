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
submodule (crystalmod) vibrations
  implicit none

  !xx! private procedures
  ! subroutine read_matdyn_modes(c,file,ivformat,errmsg,ti)
  ! subroutine read_qe_dyn(c,file,errmsg,ti)
  ! subroutine read_phonopy_ascii(c,vib,file,errmsg,ti)
  ! subroutine read_phonopy_yaml(c,vib,file,errmsg,ti)
  ! subroutine read_phonopy_hdf5(c,vib,file,errmsg)
  ! subroutine read_crystal_out(c,vib,file,errmsg,ti)
  ! subroutine read_gaussian_fchk(c,vib,file,errmsg,ti)

contains

  !> Remove the vibration data from the crystal object.
  module subroutine clear_vibrations(c)
    class(crystal), intent(inout) :: c
    if (allocated(c%vib)) deallocate(c%vib)
  end subroutine clear_vibrations

  !> Read a file (file) containing the vibrational information for
  !> this system and add the q-points to c%vib. The vibration data
  !> format is detected from the file extension if ivformat is
  !> isformat_unknown or format ivformat is used otherwise. Returns
  !> non-zero errmsg if error.
  module subroutine read_vibrations_file(c,file,ivformat,errmsg,ti)
    use param, only: ivformat_unknown, ivformat_matdynmodes, ivformat_matdyneig,&
       ivformat_qedyn, ivformat_phonopy_ascii, ivformat_phonopy_yaml, ivformat_phonopy_hdf5,&
       ivformat_crystal_out, ivformat_gaussian_log, ivformat_gaussian_fchk
    use crystalmod, only: vibrations
    use crystalseedmod, only: vibrations_detect_format
    use types, only: realloc
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: ivformat
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: ivf, nqpt_
    type(vibrations) :: vib

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

    ! read the vibrations
    if (ivf == ivformat_matdynmodes .or. ivf == ivformat_matdyneig) then
       call read_matdyn_modes(c,vib,file,ivf,errmsg,ti)
    elseif (ivf == ivformat_qedyn) then
       call read_qe_dyn(c,vib,file,errmsg,ti)
    elseif (ivf == ivformat_phonopy_ascii) then
       call read_phonopy_ascii(c,vib,file,errmsg,ti)
    elseif (ivf == ivformat_phonopy_yaml) then
       call read_phonopy_yaml(c,vib,file,errmsg,ti)
    elseif (ivf == ivformat_phonopy_hdf5) then
       call read_phonopy_hdf5(c,vib,file,errmsg)
    elseif (ivf == ivformat_crystal_out) then
       call read_crystal_out(c,vib,file,errmsg)
    elseif (ivf == ivformat_gaussian_log) then
       call read_gaussian_log(c,vib,file,errmsg)
    elseif (ivf == ivformat_gaussian_fchk) then
       call read_gaussian_fchk(c,vib,file,errmsg)
    else
       errmsg = "Unknown vibration file format: " // trim(file)
       return
    end if
    if (len(errmsg) > 0) return
    if (vib%nfreq == 0) then
       errmsg = "No frequencies found in file: " // trim(file)
       return
    end if
    if (vib%nqpt == 0) then
       errmsg = "No q-points found in file: " // trim(file)
       return
    end if

    ! add to c%vib
    if (.not.allocated(c%vib)) then
       allocate(c%vib)
       c%vib%nqpt = 0
       c%vib%nfreq = vib%nfreq
    end if
    if (c%vib%nfreq /= vib%nfreq) then
       errmsg = "Inconsistent number of frequencies: " // trim(file)
       return
    end if
    nqpt_ = c%vib%nqpt
    c%vib%file = vib%file
    c%vib%ivformat = vib%ivformat
    c%vib%nqpt = nqpt_ + vib%nqpt
    call realloc(c%vib%qpt,3,c%vib%nqpt)
    c%vib%qpt(:,nqpt_+1:c%vib%nqpt) = vib%qpt
    call realloc(c%vib%freq,c%vib%nfreq,c%vib%nqpt)
    c%vib%freq(:,nqpt_+1:c%vib%nqpt) = vib%freq
    call realloc(c%vib%vec,3,c%ncel,c%vib%nfreq,c%vib%nqpt)
    c%vib%vec(:,:,:,nqpt_+1:c%vib%nqpt) = vib%vec

  end subroutine read_vibrations_file

  !xx! private procedures

  !> Read vibration data with Quantum ESPRESSO matdyn.modes and
  !> matdyn.eig format from a file, and return it in vib. If error,
  !> return non-zero errmsg.
  subroutine read_matdyn_modes(c,vib,file,ivformat,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalmod, only: vibrations
    use crystalseedmod, only: read_alat_from_qeout
    use param, only: atmass, isformat_qeout, ivformat_matdynmodes
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
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
    vib%file = file
    vib%ivformat = ivf

    ! first pass: determine nqpt and nfreq
    nline = 0
    vib%nqpt = 0
    vib%nfreq = 0
    do while (getline_raw(lu,line,.false.))
       nline = nline + 1
       if (len(line) >= 4) then
          if (line(1:4) == " q =") vib%nqpt = vib%nqpt + 1
       end if
       if (vib%nqpt == 1) then
          if (len(line) >= 9) then
             if (line(1:9) == "     freq") then
                vib%nfreq = vib%nfreq + 1
                if (vib%nfreq == 1) then
                   nline1 = nline
                elseif (vib%nfreq == 2) then
                   nat = nline - nline1 - 1
                end if
             end if
          end if
       end if
    end do
    if (vib%nfreq == 0) then
       errmsg = "Found no frequencies in file: " // trim(file)
       goto 999
    end if
    if (vib%nqpt == 0) then
       errmsg = "Found no q-points in file: "  // trim(file)
       goto 999
    end if
    if (nat /= c%ncel) then
       errmsg = "Number of atomic displacements inconsistent with number of atoms in crystal structure"
       goto 999
    end if

    ! allocate arrays
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(3,vib%nqpt))
    if (allocated(vib%freq)) deallocate(vib%freq)
    allocate(vib%freq(vib%nfreq,vib%nqpt))
    if (allocated(vib%vec)) deallocate(vib%vec)
    allocate(vib%vec(3,c%ncel,vib%nfreq,vib%nqpt))

    ! second pass: read the information
    errmsg = "Error reading modes file: " // trim(file)
    rewind(lu)
    vib%nqpt = 0
    do while (getline_raw(lu,line,.false.))
       ok = .false.
       if (len(line) >= 4) ok = (line(1:4) == " q =")
       if (ok) then
          vib%nqpt = vib%nqpt + 1
          read (line(5:),*,end=999,err=999) vib%qpt(:,vib%nqpt)
          vib%qpt(:,vib%nqpt) = c%rc2rx(vib%qpt(:,vib%nqpt)) / alat

          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          do ifreq = 1, vib%nfreq
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) goto 999
             idx = index(line,'=',back=.true.)
             read(line(idx+1:),*) vib%freq(ifreq,vib%nqpt)
             do iat = 1, c%ncel
                ok = getline_raw(lu,line,.false.)
                if (.not.ok) goto 999
                idx = index(line,'(')
                read(line(idx+1:),*,end=999,err=999) xdum
                vib%vec(1,iat,ifreq,vib%nqpt) = cmplx(xdum(1),xdum(2),8)
                vib%vec(2,iat,ifreq,vib%nqpt) = cmplx(xdum(3),xdum(4),8)
                vib%vec(3,iat,ifreq,vib%nqpt) = cmplx(xdum(5),xdum(6),8)
             end do
          end do
       end if
    end do

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    if (ivformat == ivformat_matdynmodes) then
       do iat = 1, c%ncel
          iz = c%spc(c%atcel(iat)%is)%z
          vib%vec(:,iat,:,:) = vib%vec(:,iat,:,:) * sqrt(atmass(iz))
       end do
    end if

    ! normalize
    do iqpt = 1, vib%nqpt
       do ifreq = 1, vib%nfreq
          vib%vec(:,:,ifreq,iqpt) = vib%vec(:,:,ifreq,iqpt) / &
             sqrt(sum(vib%vec(:,:,ifreq,iqpt)*conjg(vib%vec(:,:,ifreq,iqpt))))
       end do
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, vib%nqpt
    !    do ifreq = 1, vib%nfreq
    !       do jfreq = 1, vib%nfreq
    !          summ = sum(vib%vec(:,:,ifreq,iqpt)*conjg(vib%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_matdyn_modes

  !> Read vibration data with Quantum ESPRESSO *.dyn* file format, and
  !> return it in vib. If error, return non-zero errmsg.
  subroutine read_qe_dyn(c,vib,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalmod, only: vibrations
    use crystalseedmod, only: read_alat_from_qeout
    use param, only: atmass, isformat_qeout, ivformat_qedyn
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
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
    vib%file = file
    vib%ivformat = ivformat_qedyn

    ! advance to the header
    vib%nqpt = 1
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(3,vib%nqpt))
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
    read (line(11:),*,end=999,err=999) vib%qpt(:,1)
    vib%qpt(:,1) = c%rc2rx(vib%qpt(:,1)) / alat

    ! allocate frequencies
    vib%nfreq = 3 * c%ncel
    if (allocated(vib%freq)) deallocate(vib%freq)
    allocate(vib%freq(vib%nfreq,vib%nqpt))
    if (allocated(vib%vec)) deallocate(vib%vec)
    allocate(vib%vec(3,c%ncel,vib%nfreq,vib%nqpt))

    ! read the frequencies
    ok = getline_raw(lu,line,.false.)
    ok = ok .and. getline_raw(lu,line,.false.)
    if (.not.ok) goto 999
    do i = 1, vib%nfreq
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999
       if (len(line) < 42) goto 999
       read(line(43:),*,end=999,err=999) vib%freq(i,1)

       do j = 1, c%ncel
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          if (len(line) < 2) goto 999
          read(line(3:),*,end=999,err=999) xdum
          vib%vec(1,j,i,1) = cmplx(xdum(1),xdum(2),8)
          vib%vec(2,j,i,1) = cmplx(xdum(3),xdum(4),8)
          vib%vec(3,j,i,1) = cmplx(xdum(5),xdum(6),8)
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
       vib%vec(:,i,:,:) = vib%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, vib%nfreq
       vib%vec(:,:,i,1) = vib%vec(:,:,i,1) / &
          sqrt(sum(vib%vec(:,:,i,1)*conjg(vib%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, vib%nqpt
    !    do ifreq = 1, vib%nfreq
    !       do jfreq = 1, vib%nfreq
    !          summ = sum(vib%vec(:,:,ifreq,iqpt)*conjg(vib%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_qe_dyn

  !> Read vibration data from a phonopy ascii file (ANIME keyword),
  !> and return it in vib. If error, return non-zero errmsg.
  subroutine read_phonopy_ascii(c,vib,file,errmsg,ti)
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalmod, only: vibrations
    use param, only: atmass, ivformat_phonopy_ascii, cm1tothz
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
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
    vib%file = file
    vib%ivformat = ivformat_phonopy_ascii

    ! allocate and initialize
    vib%nqpt = 1
    vib%nfreq = 3 * c%ncel
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(3,vib%nqpt))
    if (allocated(vib%freq)) deallocate(vib%freq)
    allocate(vib%freq(vib%nfreq,vib%nqpt))
    if (allocated(vib%vec)) deallocate(vib%vec)
    allocate(vib%vec(3,c%ncel,vib%nfreq,vib%nqpt))

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
    read (line,*,err=999,end=999) vib%qpt(:,1), vib%freq(1,1)

    ! get the rest of the info
    do ifreq = 1, vib%nfreq
       ! read the eigenvector
       do i = 1, c%ncel
          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          call strip(line)
          read (line,*,err=999,end=999) xdum
          vib%vec(1,i,ifreq,1) = cmplx(xdum(1),xdum(4),8)
          vib%vec(2,i,ifreq,1) = cmplx(xdum(2),xdum(5),8)
          vib%vec(3,i,ifreq,1) = cmplx(xdum(3),xdum(6),8)
       end do

       ! advance to the next metadata
       if (ifreq < vib%nfreq) then
          ok = getline_raw(lu,line,.false.)
          ok = ok .and. getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          idx = index(line,"[")
          line = line(idx+1:)
          call strip(line)
          read (line,*,err=999,end=999) xdum(1:3), vib%freq(ifreq+1,1)
       end if
    end do

    ! THz to cm-1
    vib%freq = vib%freq / cm1tothz

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       vib%vec(:,i,:,:) = vib%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, vib%nfreq
       vib%vec(:,:,i,1) = vib%vec(:,:,i,1) / &
          sqrt(sum(vib%vec(:,:,i,1)*conjg(vib%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, vib%nqpt
    !    do ifreq = 1, vib%nfreq
    !       do jfreq = 1, vib%nfreq
    !          summ = sum(vib%vec(:,:,ifreq,iqpt)*conjg(vib%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
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
  subroutine read_phonopy_yaml(c,vib,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalmod, only: vibrations
    use param, only: ivformat_phonopy_yaml, cm1tothz
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
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
    vib%file = file
    vib%ivformat = ivformat_phonopy_yaml

    ! allocate and initialize
    vib%nqpt = 0
    vib%nfreq = 3 * c%ncel
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(3,10))
    if (allocated(vib%freq)) deallocate(vib%freq)
    allocate(vib%freq(vib%nfreq,10))
    if (allocated(vib%vec)) deallocate(vib%vec)
    allocate(vib%vec(3,c%ncel,vib%nfreq,10))

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
          vib%nqpt = vib%nqpt + 1
          ifreq = 0
          if (vib%nqpt > size(vib%qpt,2)) then
             call realloc(vib%qpt,3,2*vib%nqpt)
             call realloc(vib%freq,vib%nfreq,2*vib%nqpt)
             call realloc(vib%vec,3,c%ncel,vib%nfreq,2*vib%nqpt)
          end if

          ! a new q-point
          idx = index(line,":")
          line = line(idx+1:)
          call strip(line)
          read (line,*,err=999,end=999) vib%qpt(:,vib%nqpt)

       elseif (index(line,"frequency") > 0) then
          ! a new frequency (band)
          ifreq = ifreq + 1
          if (ifreq > vib%nfreq) then
             errmsg = "Inconsistent number of frequencies"
             goto 999
          end if
          idx = index(line,":")
          line = line(idx+1:)
          read (line,*,err=999,end=999) vib%freq(ifreq,vib%nqpt)

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
                vib%vec(j,i,ifreq,vib%nqpt) = cmplx(xdum(1),xdum(2),8)
             end do
          end do
       end if
    end do
    call realloc(vib%qpt,3,vib%nqpt)
    call realloc(vib%freq,vib%nfreq,vib%nqpt)
    call realloc(vib%vec,3,c%ncel,vib%nfreq,vib%nqpt)
    if (.not.readvec) then
       errmsg = "the yaml file contains no eigenvector information"
       goto 999
    end if

    ! THz to cm-1
    vib%freq = vib%freq / cm1tothz

    ! normalize
    do i = 1, vib%nfreq
       norm = sqrt(sum(vib%vec(:,:,i,1)*conjg(vib%vec(:,:,i,1))))
       if (norm < eps_error) then
          errmsg = "zero-length eigenvector found (the yaml file has eigenvector?)"
          goto 999
       end if
       vib%vec(:,:,i,1) = vib%vec(:,:,i,1) / norm
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, vib%nqpt
    !    do ifreq = 1, vib%nfreq
    !       do jfreq = 1, vib%nfreq
    !          summ = sum(vib%vec(:,:,ifreq,iqpt)*conjg(vib%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
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
  !> QPOINTS, plus EIGENVECTORS keywords), and return it in vib. If
  !> error, return non-zero errmsg.
  subroutine read_phonopy_hdf5(c,vib,file,errmsg)
#ifdef HAVE_HDF5
    use iso_c_binding, only: c_loc, c_ptr
    use hdf5
#endif
    use param, only: ivformat_phonopy_hdf5, cm1tothz
    use crystalmod, only: vibrations
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
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
    vib%file = file
    vib%ivformat = ivformat_phonopy_hdf5

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
    vib%nqpt = int(dim(2))
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(dim(1),dim(2)))
    call h5dread_f(dsetid,H5T_NATIVE_DOUBLE,vib%qpt,dim,ier)
    if (ier < 0) goto 999

    ! read frequencies
    call h5dopen_f(fid,"frequency",dsetid,ier)
    if (ier < 0) goto 999
    call h5dget_space_f(dsetid,spaceid,ier)
    if (ier < 0) goto 999
    CALL h5sget_simple_extent_dims_f(spaceid,dim,maxdim,ier)
    if (ier < 0) goto 999
    vib%nfreq = int(dim(1))
    if (3*c%ncel /= vib%nfreq) then
       errmsg = "Inconsistent number of frequencies: " // trim(file)
       goto 999
    end if
    if (allocated(vib%freq)) deallocate(vib%freq)
    allocate(vib%freq(dim(1),dim(2)))
    call h5dread_f(dsetid,H5T_NATIVE_DOUBLE,vib%freq,dim,ier)
    if (ier < 0) goto 999

    ! THz to cm-1
    vib%freq = vib%freq / cm1tothz

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
    allocate(vaux(vib%nfreq,3,c%ncel,vib%nqpt))
    cptr = c_loc(vaux)
    call h5dread_f(dsetid,ctypeid,cptr,ier)
    if (ier < 0) goto 999
    allocate(vib%vec(3,c%ncel,vib%nfreq,vib%nqpt))
    do iqpt = 1, vib%nqpt
       do ifreq = 1, vib%nfreq
          do i = 1, c%ncel
             vib%vec(:,i,ifreq,iqpt) = vaux(ifreq,:,i,iqpt)
          end do
       end do
    end do
    deallocate(vaux)

    ! wrap up
    errmsg = ""
999 continue
    call h5close_f(ier)

#else
    errmsg = "Critic2 must be compiled against the HDF5 library to read phonopy hdf5 files"
#endif

  end subroutine read_phonopy_hdf5

  !> Read vibration data from a crystal output file, and return it in
  !> vib. If error, return non-zero errmsg.
  subroutine read_crystal_out(c,vib,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalmod, only: vibrations
    use param, only: ivformat_crystal_out, atmass
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
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
    vib%file = file
    vib%ivformat = ivformat_crystal_out

    ! allocate and initialize -> only gamma point
    vib%nqpt = 1
    vib%nfreq = 3 * c%ncel
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(3,1))
    vib%qpt = 0d0
    if (allocated(vib%freq)) deallocate(vib%freq)
    allocate(vib%freq(vib%nfreq,1))
    if (allocated(vib%vec)) deallocate(vib%vec)
    allocate(vib%vec(3,c%ncel,vib%nfreq,1))

    ! advance to the header
    do while (getline_raw(lu,line,.false.))
       if (index(line,"NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES (IN BOHR)") > 0) exit
    end do
    ok = getline_raw(lu,line,.false.)
    if (.not.ok) goto 999

    ! read the info
    do i = 1, vib%nfreq, 6
       n0 = i
       n1 = min(i+5,vib%nfreq)

       ! read the frequency line
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999
       line = line(14:)
       read (line,*,end=999,err=999) vib%freq(n0:n1,1)
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999

       ! read the eigenvector
       do j = 1, c%ncel
          do k = 1, 3
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) goto 999
             line = line(14:)
             read (line,*,end=999,err=999) xdum
             vib%vec(k,j,n0:n1,1) = xdum
          end do
       end do
       ok = getline_raw(lu,line,.false.)
       if (.not.ok) goto 999
    end do

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       vib%vec(:,i,:,:) = vib%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, vib%nfreq
       vib%vec(:,:,i,1) = vib%vec(:,:,i,1) / &
          sqrt(sum(vib%vec(:,:,i,1)*conjg(vib%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do iqpt = 1, vib%nqpt
    !    do ifreq = 1, vib%nfreq
    !       do jfreq = 1, vib%nfreq
    !          summ = sum(vib%vec(:,:,ifreq,iqpt)*conjg(vib%vec(:,:,jfreq,iqpt)))
    !          write (*,*) ifreq, jfreq, summ
    !       end do
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_crystal_out

  !> Read vibration data from a Gaussian output file, and return it in
  !> vib. If error, return non-zero errmsg.
  subroutine read_gaussian_log(c,vib,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw, isreal
    use crystalmod, only: vibrations
    use param, only: ivformat_gaussian_log, atmass
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
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
    vib%file = file
    vib%ivformat = ivformat_gaussian_log

    ! allocate and initialize -> only gamma point
    vib%nqpt = 1
    vib%nfreq = 3 * c%ncel
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(3,1))
    vib%qpt = 0d0
    if (allocated(vib%freq)) deallocate(vib%freq)
    allocate(vib%freq(vib%nfreq,1))
    if (allocated(vib%vec)) deallocate(vib%vec)
    allocate(vib%vec(3,c%ncel,vib%nfreq,1))

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
          if (okl(i)) vib%freq(nread+i,1) = xdum
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
             vib%vec(:,i,nread+j,1) = xread(3*(j-1)+1:3*j)
          end do
       end do
       nread = nread + nn
    end do
    vib%nfreq = nread
    call realloc(vib%freq,vib%nfreq,vib%nqpt)
    call realloc(vib%vec,3,c%ncel,vib%nfreq,vib%nqpt)

    ! convert to mass-weighed coordinates (orthonormal eigenvectors)
    ! (note: units are incorrect, but we are normalizing afterwards)
    do i = 1, c%ncel
       iz = c%spc(c%atcel(i)%is)%z
       vib%vec(:,i,:,:) = vib%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, vib%nfreq
       vib%vec(:,:,i,1) = vib%vec(:,:,i,1) / &
          sqrt(sum(vib%vec(:,:,i,1)*conjg(vib%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do ifreq = 1, vib%nfreq
    !    do jfreq = 1, vib%nfreq
    !       summ = sum(vib%vec(:,:,ifreq,1)*conjg(vib%vec(:,:,jfreq,1)))
    !       write (*,*) ifreq, jfreq, summ
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_gaussian_log

  !> Read vibration data from a Gaussian formatted checkpoint file,
  !> and return it in vib. If error, return non-zero errmsg.
  subroutine read_gaussian_fchk(c,vib,file,errmsg,ti)
    use types, only: realloc
    use tools_io, only: fopen_read, fclose, getline_raw, isreal
    use crystalmod, only: vibrations
    use param, only: ivformat_gaussian_fchk, atmass
    class(crystal), intent(inout) :: c
    type(vibrations), intent(inout) :: vib
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    logical :: ok
    character(len=:), allocatable :: line
    integer :: lu, i, j, k, nread, iz, nn, ll, n0, n1, nread, nnm
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
    vib%file = file
    vib%ivformat = ivformat_gaussian_fchk

    ! allocate and initialize -> only gamma point
    vib%nqpt = 1
    vib%nfreq = 0
    if (allocated(vib%qpt)) deallocate(vib%qpt)
    allocate(vib%qpt(3,1))
    vib%qpt = 0d0
    if (allocated(vib%freq)) deallocate(vib%freq)
    if (allocated(vib%vec)) deallocate(vib%vec)

    do while (getline_raw(lu,line,.false.))
       ll = len(line)

       ! read the number of frequencies
       if (ll >= 22) then
          if (line(1:22) == "Number of Normal Modes") then
             line = line(23:)
             read (line,*,end=999,err=999) ch, vib%nfreq
             allocate(vib%freq(vib%nfreq,1))
             allocate(vib%vec(3,c%ncel,vib%nfreq,1))
          end if
       end if

       ! read frequencies
       if (ll >= 6) then
          if (line(1:6) == "Vib-E2") then
             do i = 1, vib%nfreq, 5
                n0 = i
                n1 = min(i+4,vib%nfreq)
                ok = getline_raw(lu,line,.false.)
                if (.not.ok) goto 999
                read (line,*,end=999,err=999) vib%freq(n0:n1,1)
             end do
          end if
       end if

       ! read modes
       if (ll >= 9) then
          if (line(1:9) == "Vib-Modes") then

             nn = -1
             do i = 1, vib%nfreq
                do j = 1, c%ncel
                   do k = 1, 3
                      nn = nn + 1
                      nnm = mod(nn,5) + 1
                      if (mod(nn,5) == 0) then
                         nread = min(vib%nfreq * c%ncel * 3 - nn,5)
                         ok = getline_raw(lu,line,.false.)
                         if (.not.ok) goto 999
                         read (line,*,end=999,err=999) xdum(1:nread)
                      end if
                      vib%vec(k,j,i,1) = xdum(nnm)
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
       vib%vec(:,i,:,:) = vib%vec(:,i,:,:) * sqrt(atmass(iz))
    end do

    ! normalize
    do i = 1, vib%nfreq
       vib%vec(:,:,i,1) = vib%vec(:,:,i,1) / &
          sqrt(sum(vib%vec(:,:,i,1)*conjg(vib%vec(:,:,i,1))))
    end do

    ! ! checking normalization
    ! write (*,*) "checking normalization..."
    ! do ifreq = 1, vib%nfreq
    !    do jfreq = 1, vib%nfreq
    !       summ = sum(vib%vec(:,:,ifreq,1)*conjg(vib%vec(:,:,jfreq,1)))
    !       write (*,*) ifreq, jfreq, summ
    !    end do
    ! end do

    ! wrap up
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (lu >= 0) call fclose(lu)

  end subroutine read_gaussian_fchk

end submodule vibrations
