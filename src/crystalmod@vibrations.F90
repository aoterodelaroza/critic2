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
  module subroutine vibrations_end(v,keepfc2)
    use param, only: ivformat_unknown
    class(vibrations), intent(inout) :: v
    logical, intent(in), optional :: keepfc2

    logical :: keepfc2_

    keepfc2_ = .false.
    if (present(keepfc2)) keepfc2_ = keepfc2

    v%hasvibs = .false.
    v%file = ""
    v%ivformat = ivformat_unknown
    v%nqpt = 0
    if (allocated(v%qpt)) deallocate(v%qpt)
    if (allocated(v%freq)) deallocate(v%freq)
    if (allocated(v%vec)) deallocate(v%vec)
    v%nfreq = 0
    if (.not.keepfc2_) then
       v%hasfc2 = .false.
       v%fc2_vs_delta = -1d0
       v%fc2_gamma_ac = -1d0
       if (allocated(v%fc2)) deallocate(v%fc2)
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
       ivformat_phonopy_fc2
    use types, only: realloc
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
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
       ivformat_crystal_out, ivformat_gaussian_log,&
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
    else
       write (uout,*)
       write (uout,'("# Second-order force constants NOT available")')
    end if

  end subroutine vibrations_print_summary

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
       write (uout,'("# Id1,Id2 = complete cell IDs. At1,At2 = atomic symbols. isintra = is intramolecular?")')
       write (uout,'("# norm = Frobenius norm FC2. xx,... = FC2 components")')
       write (uout,'("# All FC2 in Hartree/bohr^2.")')
       write (uout,'("# Id1 At1  Id2 At2  dist(bohr) isintra FC2: norm      xx              xy              xz   &
          &           yx              yy              yz              zx              zy              zz   ")')
       k = 0
       do idum = 1, c%ncel
          do jdum = idum, c%ncel
             k = k + 1
             i = i2(1,idx(k))
             j = i2(2,idx(k))
             dist = dista(idx(k))
             isintra = (c%idatcelmol(1, i) == c%idatcelmol(1, j))
             norm = sqrt(sum(v%fc2(:,:,i,j)**2))
             if (dist <= disteps_ .and. norm >= fc2eps_) then
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
          write (uout,'("# name = atomic name. dist = distance. isintra = is intramolecular?")')
          write (uout,'("# norm = Frobenius norm FC2. xx,... = FC2 components")')
          write (uout,'("# All FC2 in Hartree/bohr^2.")')
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
             if (distr(j) <= disteps_ .and. norm >= fc2eps_) then
                write (uout,'("  ",2(A," "),"(",A," ",A," ",A,")",99(" ",A))') string(nidx,4,ioj_center), string(cidx,4,ioj_center),&
                   (string(lvec(k,j),2,ioj_right),k=1,3), string(c%spc(c%atcel(cidx)%is)%name,7,ioj_center),&
                   string(distr(j)*dunit0(iunit),'f',12,6,4), string(isintra), string(norm,'e',15,8,ioj_right),&
                   string(v%fc2(1,1,i,cidx),'e',15,8,ioj_right), string(v%fc2(1,2,i,cidx),'e',15,8,ioj_right),&
                   string(v%fc2(1,3,i,cidx),'e',15,8,ioj_right), string(v%fc2(2,1,i,cidx),'e',15,8,ioj_right),&
                   string(v%fc2(2,2,i,cidx),'e',15,8,ioj_right), string(v%fc2(2,3,i,cidx),'e',15,8,ioj_right),&
                   string(v%fc2(3,1,i,cidx),'e',15,8,ioj_right), string(v%fc2(3,2,i,cidx),'e',15,8,ioj_right),&
                   string(v%fc2(3,3,i,cidx),'e',15,8,ioj_right)
             end if
          end do
          write (uout,*)
       end do
    end if

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

    integer :: i, j, iat, jat
    real*8 :: summ
    logical :: verbose_

    ! return if no FC2 is available
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    ! optional parameters
    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose

    ! apply acoustic sum rules
    do iat = 1, c%ncel
       do i = 1, 3
          do j = 1, 3
             summ = sum(v%fc2(i,j,iat,:)) / c%ncel
             v%fc2(i,j,iat,:) = v%fc2(i,j,iat,:) - summ
          end do
       end do
    end do
    do iat = 1, c%ncel
       do i = 1, 3
          do j = 1, 3
             summ = sum(v%fc2(i,j,:,iat)) / c%ncel
             v%fc2(i,j,:,iat) = v%fc2(i,j,:,iat) - summ
          end do
       end do
    end do

    do iat = 1, c%ncel
       do jat = iat+1, c%ncel
          do i = 1, 3
             do j = 1, 3
                v%fc2(i,j,iat,jat) = 0.5d0 * (v%fc2(i,j,iat,jat) + v%fc2(j,i,jat,iat))
                v%fc2(j,i,jat,iat) = v%fc2(i,j,iat,jat)
             end do
          end do
       end do
    end do
    if (verbose_) then
       write (uout,'("+ FC2 ACOUSTIC_SUM_RULES: apply acoustic sum rules to FC2")')
       write (uout,'("  Resulting acoustic sum = ",A)') string(sum(v%fc2),'e',10,5)
    end if

  end subroutine vibrations_apply_acoustic

  !> Write the FC2 to a phonopy-style file. If no file is given or if
  !> file is empty, use FORCE_CONSTANTS.
  module subroutine vibrations_write_fc2(v,c,file,verbose)
    use tools_io, only: fopen_write, fclose, uout, string, ioj_right
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(in), optional :: file
    logical, intent(in), optional :: verbose

    integer :: lu, i, j, k
    logical :: verbose_
    character(len=:), allocatable :: file_

    ! return if no FC2 is available
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    ! optional parameters
    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose
    file_ = "FORCE_CONSTANTS"
    if (present(file)) then
       if (len(file) > 0) file_ = file
    end if

    ! write the file
    lu = fopen_write(file_)
    write (lu,'(A,X,A)') string(c%ncel), string(c%ncel)
    do i = 1, c%ncel
       do j = 1, c%ncel
          write (lu,'(A,X,A)') string(i), string(j)
          write (lu,'(3(A,X))') (string(2*v%fc2(1,k,i,j),'f',22,15,ioj_right),k=1,3)
          write (lu,'(3(A,X))') (string(2*v%fc2(2,k,i,j),'f',22,15,ioj_right),k=1,3)
          write (lu,'(3(A,X))') (string(2*v%fc2(3,k,i,j),'f',22,15,ioj_right),k=1,3)
       end do
    end do
    call fclose(lu)

    write (uout,'("+ Written FC2 file (QE format, Ry/bohr^2): ",A)') file_

  end subroutine vibrations_write_fc2

  !> Calculate frequencies and eigenvectors from the FC2 for a single
  !> q (fractional coordiantes in reciprocal space); adds the
  !> resulting frequencies and eigenvectors to the v type.  If the
  !> optional arguments freqo and veco are present, return the
  !> frequencies and eigenvectors in these variables instead of
  !> overwriting the vibrational info in v. Frequencies in output are
  !> in ascending order.
  module subroutine vibrations_calculate_q(v,c,q,freqo,veco)
    use param, only: icrd_crys, atmass, tpi, img
    use types, only: realloc
    use tools_math, only: eigherm
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: q(3)
    real*8, intent(inout), allocatable, optional :: freqo(:)
    complex*16, intent(inout), allocatable, optional :: veco(:,:)

    complex*16 :: dm_local(3,3), phase
    real*8 :: sqrt_ij
    integer :: nat
    real*8 :: x1(3), x2(3), maxdisteps
    integer :: i, jj, j
    integer, allocatable :: nida(:), lvec(:,:), mult(:)
    real*8, allocatable:: eval(:), dist(:), freq(:), rused(:)
    complex*16, allocatable :: dm(:,:)
    logical :: varoutput

    ! conversion factor
    ! sqrt(Hartree / bohr^2 / amu) / 1e12 -> THz
    ! sqrt(4.3597447222071e-18) / 5.29177210903e-11 / sqrt(1.66053906660e-27) /  2 / pi / c / 100
    real*8, parameter :: factor = 5140.487143715827d0
    real*8, parameter :: epsgen = 1d-4
    real*8, parameter :: epsneigh = 1d-5

    ! return if no FC2 is available
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    ! process input arguments
    varoutput = present(freqo) .and. present(veco)

    ! maximum distance of all pairs of atoms in the unit cell
    maxdisteps = 0d0
    do i = 1, c%ncel
       do j = i+1, c%ncel
          maxdisteps = max(maxdisteps,c%eql_distance(c%atcel(i)%x, c%atcel(j)%x))
       end do
    end do

    ! build the dynamical matrix for this q-point (Parlinski's recipe)
    allocate(dm(c%ncel*3,c%ncel*3),rused(c%ncel),mult(c%ncel))
    dm = 0d0
    do i = 1, c%ncel
       rused = huge(1d0)
       mult = 0
       call c%list_near_atoms(c%atcel(i)%x,icrd_crys,.true.,nat,nida,dist=dist,lvec=lvec,up2d=maxdisteps+epsgen)
       do jj = 1, nat
          j = nida(jj)
          if (dist(jj) <= rused(j) + epsneigh) then
             x1 = c%atcel(i)%x
             x2 = c%atcel(j)%x + lvec(:,jj)

             sqrt_ij = sqrt(atmass(c%spc(c%atcel(i)%is)%z) * atmass(c%spc(c%atcel(j)%is)%z))
             phase = exp(-tpi * img * dot_product(q,x1 - x2))

             dm_local = v%fc2(:,:,i,j) * phase / sqrt_ij
             dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) = dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) + dm_local

             rused(j) = min(rused(j),dist(jj))
             mult(j) = mult(j) + 1
          end if
       end do

       ! apply the multiplicities
       do j = 1, c%ncel
          dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) = dm(3*(i-1)+1:3*i,3*(j-1)+1:3*j) / mult(j)
       end do
    end do

    ! impose hermiticity
    dm = (dm + transpose(conjg(dm)))/2

    ! diagonalize
    allocate(eval(c%ncel*3),freq(c%ncel*3))
    call eigherm(dm,c%ncel*3,eval)

    ! calculate the frequencies
    freq = sqrt(abs(eval)) * factor

    ! save output
    if (varoutput) then
       freqo = freq
       veco = dm
    else
       ! save the info
       v%nqpt = v%nqpt + 1
       if (.not.allocated(v%qpt)) allocate(v%qpt(3,v%nqpt))
       if (.not.allocated(v%freq)) then
          v%nfreq = 3 * c%ncel
          allocate(v%freq(3*c%ncel,v%nqpt))
       end if
       if (.not.allocated(v%vec)) allocate(v%vec(3,c%ncel,3*c%ncel,v%nqpt))
       if (v%nqpt > size(v%qpt,2)) call realloc(v%qpt,3,2*v%nqpt)
       if (v%nqpt > size(v%freq,2)) call realloc(v%freq,size(v%freq,1),2*v%nqpt)
       if (v%nqpt > size(v%vec,4)) call realloc(v%vec,size(v%vec,1),size(v%vec,2),size(v%vec,3),2*v%nqpt)
       v%qpt(:,v%nqpt) = q
       v%freq(:,v%nqpt) = freq
       do i = 1, 3*c%ncel
          v%vec(:,:,i,v%nqpt) = reshape(dm(:,i),(/3,c%ncel/))
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
    use param, only: fact, cm1tohz, bohrtom
    class(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: q(3)
    real*8, intent(out) :: vs(3)

    real*8 :: normq, qn(3), qthis(3), h, fd0(3), fdn, ratio, fdiff(3)
    integer :: i, j
    real*8, allocatable :: freq(:), ff(:,:)
    complex*16, allocatable :: vec(:,:)

    real*8, parameter :: fdiff_thr = 0.05d0 ! threshold for bracketing, cm-1
    real*8, parameter :: maxlen_ini = 1d-8 ! initial bracketing value
    real*8, parameter :: maxlen_step = 2d0 ! bracketing step
    real*8, parameter :: epsnq = 1d-10 ! cutoff for zero q-point
    integer, parameter :: npts = 3 ! can be up to 6
    real*8, parameter :: epsfd = 1d-3 ! error if finite differences disagree more than this

    ! abramowitz, stegun, numerical differentiation formulas, table 25.2 (first derivatives)
    real*8, parameter :: coefd(6,2:5) = reshape((/&
       -3d0, 4d0, -1d0, 0d0, 0d0, 0d0,&
       -11d0, 18d0, -9d0, 2d0, 0d0, 0d0,&
       -50d0, 96d0, -72d0, 32d0, -6d0, 0d0,&
       -274d0, 600d0, -600d0, 400d0, -150d0, 24d0/),(/6,4/))

    ! return if no FC2 is available
    if (.not.v%hasfc2.or..not.allocated(v%fc2)) return

    ! normalize q
    normq = norm2(q)
    if (normq < epsnq) &
       call ferror('vibrations_calculate_vs','zero-length q-point',faterr)
    qn = q / normq

    ! apply acoustic sum rules
    call v%apply_acoustic(c)

    ! find frequencies at gamma and maxlen if they are not available
    if (v%fc2_vs_delta < 0d0) then
       ! calculate the frequencies at gamma
       qthis = 0d0
       call v%calculate_q(c,qthis,freqo=freq,veco=vec)
       v%fc2_gamma_ac = freq(1:3)

       ! bracketing for maxlen
       fdiff = 0d0
       v%fc2_vs_delta = maxlen_ini
       do while (any(fdiff < fdiff_thr))
          v%fc2_vs_delta = v%fc2_vs_delta * maxlen_step
          h = v%fc2_vs_delta / real(npts-1,8)
          qthis = h * qn
          qthis = c%rc2rx(qthis)
          call v%calculate_q(c,qthis,freqo=freq,veco=vec)
          fdiff = abs(freq(1:3) - v%fc2_gamma_ac)
       end do
    end if

    ! build the list of q and run the q-point calculation
    allocate(ff(npts,3))
    ff(1,:) = v%fc2_gamma_ac
    do i = 2, npts
       h = v%fc2_vs_delta / real(npts-1,8)
       qthis = real(i-1,8) * h * qn
       qthis = c%rc2rx(qthis)
       call v%calculate_q(c,qthis,freqo=freq,veco=vec)
       ff(i,:) = freq(1:3)
    end do
    if (allocated(freq)) deallocate(freq)
    if (allocated(vec)) deallocate(vec)

    ! calculate finite differences, verify the error is low
    do i = 1, 3
       fd0(i) = (ff(2,i) - ff(1,i)) / h
       do j = 2, npts-1
          fdn = sum(ff(:,i) * coefd(1:npts,j)) / (h * fact(j))
          ratio = abs(fdn-fd0(i)) / abs(fdn+fd0(i))
          if (ratio > epsfd) &
             call ferror('vibrations_calculate_vs','error calculating finite differences',faterr)
       end do
    end do

    ! convert from cm-1*bohr to m/s
    vs = fd0 * cm1tohz * bohrtom

  end subroutine vibrations_calculate_vs

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

    real*8, parameter :: cminv_to_kJmol = 1.196265656955833d-02 ! c * h * NA / 10
    real*8, parameter :: cminv_to_K = 1.438776959983815d0 ! c * h * 100 / kb
    real*8, parameter :: R = 8.314462618d0 ! molar gas constant (kB/NA, J/K/mol)
    real*8, parameter :: small1 = 50000d0 * cminv_to_K / huge(1d0) ! protection against zerodiv in nu/(kB*T)
    real*8, parameter :: small2 = 0.5d0 * log(huge(1d0)) ! protection against overflow in exp(nu/kB*T)**2
    real*8, parameter :: small3 = sqrt(1d0 / huge(1d0)) ! protection against 1/(e^(nu/kt)-1)^2 overflow

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
  subroutine read_phonopy_fc2(v,c,file,sline,errmsg,ti)
    use crystalseedmod, only: crystalseed
    use global, only: eval_next
    use types, only: realloc
    use tools_math, only: matinv, det3
    use tools_io, only: fopen_read, fclose, getline_raw, lgetword, equal, getword
    use param, only: ivformat_phonopy_fc2, hartoev, bohrtoa
    type(vibrations), intent(inout) :: v
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: file, sline
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    character(len=:), allocatable :: word
    integer :: lp, lp2, lu, nat, mat, i, j, idum, jdum
    real*8 :: fc2factor, rdum(4), x0(3,3), rmat(3,3)
    logical :: ok
    integer :: isin ! 0 = full; 1 = left is cell; 2 = right is cell
    integer, allocatable :: atop(:,:,:), ineq(:), ineqrot(:), ineqcen(:)
    type(crystal) :: sc
    type(crystalseed) :: seed
    logical, allocatable :: done(:,:), dneq(:)
    real*8, allocatable :: rotc(:,:,:)
    integer :: irot, icv
    real*8 :: fc2_(3,3), rot(3,3)
    real*8 :: xi_(3), xj_(3), xdif(3)

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
    elseif (equal(word,"aims").or.equal(word,"fhiaims")) then
       fc2factor = bohrtoa**2 / hartoev
    else
       errmsg = "A generator keyword (qe,vasp,fhiaims,...) is required to read FORCE_CONSTANTS"
       goto 999
    end if

    ! read the supercell transformation or the supercell file
    lp2 = lp
    x0 = 0d0
    ok = eval_next(rdum(1),sline,lp)
    ok = ok .and. eval_next(rdum(2),sline,lp)
    ok = ok .and. eval_next(rdum(3),sline,lp)
    if (ok) then
       ok = eval_next(rdum(4),sline,lp)
       if (ok) then
          x0(:,1) = rdum(1:3)
          x0(1,2) = rdum(4)
          ok = eval_next(x0(2,2),sline,lp)
          ok = ok .and. eval_next(x0(3,2),sline,lp)
          ok = ok .and. eval_next(x0(1,3),sline,lp)
          ok = ok .and. eval_next(x0(2,3),sline,lp)
          ok = ok .and. eval_next(x0(3,3),sline,lp)
          if (.not.ok) then
             errmsg = "Wrong syntax for supercell transformation"
             goto 999
          end if
       else
          do i = 1, 3
             x0(i,i) = rdum(i)
          end do
       end if
       ! create the supercell
       sc = c
       call sc%newcell(x0)
    else
       lp = lp2
       word = getword(sline,lp)
       call seed%read_any_file(word,-1,errmsg)
       if (len_trim(errmsg) > 0) goto 999
       call sc%struct_new(seed,.true.)
       if (sc%ismolecule) then
          errmsg = "Supercell file in VIBRATIONS LOAD is a molecule"
          goto 999
       end if
    end if

    ! open file
    lu = fopen_read(file,ti=ti)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! read the number of atoms from the fc2 file
    read(lu,*,err=999,end=999) nat, mat
    isin = -1
    if (nat == mat .and. nat == c%ncel) then
       isin = 0
    elseif (nat /= mat .and. mat == c%ncel) then
       isin = 1
    elseif (nat /= mat .and. nat == c%ncel) then
       isin = 2
       write (*,*) "fixme!! isin = 2"
       stop 1
    else
       errmsg = 'Inconsistent atom number in FORCE_CONSTANTS file'
       goto 999
    end if

    ! read the fc2 matrix
    if (isin > 0) then
       allocate(done(mat,mat))
       done = .false.
    end if
    errmsg = "Error reading force constants from FORCE_CONSTANTS file"
    if (allocated(v%fc2)) deallocate(v%fc2)
    allocate(v%fc2(3,3,mat,mat))
    v%fc2 = 0d0
    do idum = 1, nat
       do jdum = 1, mat
          read (lu,*,err=999,end=999) i, j
          read (lu,*,err=999,end=999) rmat(1,:)
          read (lu,*,err=999,end=999) rmat(2,:)
          read (lu,*,err=999,end=999) rmat(3,:)

          v%fc2(:,:,i,j) = rmat
          if (isin > 0) done(i,j) = .true.
       end do
    end do
    v%fc2 = v%fc2 * fc2factor

    ! for which atoms do we have full FC2 info?
    allocate(dneq(c%ncel))
    dneq = .false.
    do i = 1, c%ncel
       dneq(i) = all(done(i,:))
    end do
    deallocate(done)

    ! isin == 1 corresponds to the case when the crystal structure is
    ! a supercell and we are only given the force constants for the
    ! non-equivalent atoms on the left atom.
    if (isin == 1) then
       ! calculate the effect of every symmetry operation on every atom
       ! xxxx FIXME?
       allocate(atop(c%ncel,c%neqv,c%ncv))
       atop = 0
       do i = 1, c%ncel
          do irot = 1, c%neqv
             xj_ = matmul(c%rotm(1:3,1:3,irot),c%atcel(i)%x) + c%rotm(:,4,irot)
             do icv = 1, c%ncv
                xi_ = xj_ + c%cen(:,icv)
                do j = 1, c%ncel
                   xdif = c%atcel(j)%x - xi_
                   if (all(abs(xdif - nint(xdif)) < 1d-2)) then
                      atop(i,irot,icv) = j
                      exit
                   end if
                end do
             end do
          end do
       end do
       if (any(atop == 0)) then
          errmsg = "Error calculating atom mapping from symmetry operations; check structure"
          return
       end if

       ! calculate the representative atoms
       ! calculate the symmetry operations that lead to each atom from the representatives
       !   ineqrot(i) and ineqcen(i) are the symmetry operation IDs such that:
       !   rot * i + cen = j
       allocate(ineq(c%ncel),ineqrot(c%ncel),ineqcen(c%ncel))
       ineq = c%ncel + 1
       do icv = 1, c%ncv
          do irot = 1, c%neqv
             do i = 1, c%ncel
                j = atop(i,irot,icv)
                if (dneq(j) .and. j < ineq(i)) then
                   ineq(i) = j
                   ineqrot(i) = irot
                   ineqcen(i) = icv
                end if
             end do
          end do
       end do
       if (any(ineq == c%ncel + 1)) then
          errmsg = "Error calculating non-equivalent atom mapping; check structure"
          return
       end if

       ! calculate the rotation matrices in cartesian coordinates
       allocate(rotc(3,3,c%neqv))
       do irot = 1, c%neqv
          rotc(:,:,irot) = matmul(c%m_x2c,matmul(c%rotm(1:3,1:3,irot),c%m_c2x))
       end do

       ! calculate all remaining fc2
       do i = 1, c%ncel
          if (dneq(i)) cycle
          rot = rotc(:,:,ineqrot(i))
          do j = 1, c%ncel
             fc2_ = matmul(transpose(rot),matmul(v%fc2(:,:,ineq(i),atop(j,ineqrot(i),ineqcen(i))),rot))
             v%fc2(:,:,i,j) = fc2_
          end do
       end do

       ! clean up
       deallocate(atop,ineq,ineqrot,ineqcen,rotc,dneq)
    end if

    ! wrap up
    v%file = file
    v%ivformat = ivformat_phonopy_fc2
    v%hasfc2 = .true.
    v%fc2_vs_delta = -1d0
    v%fc2_gamma_ac = -1d0
    v%hasvibs = .false.
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
          if (len(line) >= 30) then
             if (index(line," Harmonic frequencies (cm**-1)") > 0) nread = 0
          end if
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

end submodule vibrationsmod
