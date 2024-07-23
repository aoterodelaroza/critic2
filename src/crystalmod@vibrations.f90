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
submodule (crystalmod) edit
  implicit none

contains

  !> Remove the vibration data from the crystal object.
  module subroutine clear_vibrations(c)
    class(crystal), intent(inout) :: c
    if (allocated(c%vib)) deallocate(c%vib)
  end subroutine clear_vibrations

  !> Read a file containing the vibrational information for this
  !> system. Populates c%vib.
  module subroutine read_vibrations_file(c,file,ivformat,errmsg,ti)
    use tools_math, only: matinv
    use tools_io, only: fopen_read, fclose, getline_raw
    use crystalseedmod, only: vibrations_detect_format
    use param, only: isformat_unknown, pi
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    integer, intent(in) :: ivformat
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: ivf, lu, nline, nline1, nat, idx
    logical :: ok
    character(len=:), allocatable :: line
    integer :: ifreq, iat, ier
    real*8 :: xdum(6)

    ! initialize
    if (allocated(c%vib)) deallocate(c%vib)
    lu = -1

    ! detect the format
    if (ivformat == isformat_unknown) then
       call vibrations_detect_format(file,ivf,ti)
       if (ivf == isformat_unknown) then
          errmsg = "Unknown vibration file format: " // trim(file)
          return
       end if
    else
       ivf = ivformat
    end if

    !!xxxx!! for now, only matdyn.modes is understood !!xxxx!!

    ! open file
    lu = fopen_read(file)
    if (lu <= 0) then
       errmsg = "File not found: " // trim(file)
       goto 999
    end if

    ! prepare container for data
    allocate(c%vib)
    c%vib%file = file

    ! first pass: determine nqpt and nfreq
    nline = 0
    c%vib%qpt_digits = 3
    c%vib%nqpt = 0
    c%vib%nfreq = 0
    do while (getline_raw(lu,line,.false.))
       nline = nline + 1
       if (len(line) >= 4) then
          if (line(1:4) == " q =") c%vib%nqpt = c%vib%nqpt + 1
       end if
       if (c%vib%nqpt == 1) then
          if (len(line) >= 9) then
             if (line(1:9) == "     freq") then
                c%vib%nfreq = c%vib%nfreq + 1
                if (c%vib%nfreq == 1) then
                   nline1 = nline
                elseif (c%vib%nfreq == 2) then
                   nat = nline - nline1 - 1
                end if
             end if
          end if
       end if
    end do
    if (c%vib%nfreq == 0) then
       errmsg = "Found no frequencies in file: " // trim(file)
       goto 999
    end if
    if (c%vib%nqpt == 0) then
       errmsg = "Found no q-points in file: "  // trim(file)
       goto 999
    end if
    if (nat /= c%ncel) then
       errmsg = "Number of atomic displacements inconsistent with number of atoms in crystal structure"
       goto 999
    end if

    ! allocate arrays
    allocate(c%vib%qpt(3,c%vib%nqpt))
    allocate(c%vib%freq(c%vib%nfreq,c%vib%nqpt))
    allocate(c%vib%vec(3,c%ncel,c%vib%nfreq,c%vib%nqpt))

    ! second pass: read the information
    errmsg = "Error reading modes file: " // trim(file)
    rewind(lu)
    c%vib%nqpt = 0
    do while (getline_raw(lu,line,.false.))
       ok = .false.
       if (len(line) >= 4) ok = (line(1:4) == " q =")
       if (ok) then
          c%vib%nqpt = c%vib%nqpt + 1
          read (line(5:),*,end=999,err=999) c%vib%qpt(:,c%vib%nqpt)
          c%vib%qpt(:,c%vib%nqpt) = c%rc2rx(c%vib%qpt(:,c%vib%nqpt))

          ok = getline_raw(lu,line,.false.)
          if (.not.ok) goto 999
          do ifreq = 1, c%vib%nfreq
             ok = getline_raw(lu,line,.false.)
             if (.not.ok) goto 999
             idx = index(line,'=',back=.true.)
             read(line(idx+1:),*) c%vib%freq(ifreq,c%vib%nqpt)
             do iat = 1, c%ncel
                ok = getline_raw(lu,line,.false.)
                if (.not.ok) goto 999
                idx = index(line,'(')
                read(line(idx+1:),*,end=999,err=999) xdum
                c%vib%vec(1,iat,ifreq,c%vib%nqpt) = cmplx(xdum(1),xdum(2),16)
                c%vib%vec(2,iat,ifreq,c%vib%nqpt) = cmplx(xdum(3),xdum(4),16)
                c%vib%vec(3,iat,ifreq,c%vib%nqpt) = cmplx(xdum(5),xdum(6),16)
             end do
          end do
       end if
    end do

    ! wrap up
    errmsg = ""
    call fclose(lu)

    return
999 continue
    if (allocated(c%vib)) deallocate(c%vib)
    if (lu >= 0) call fclose(lu)

  end subroutine read_vibrations_file

end submodule edit
