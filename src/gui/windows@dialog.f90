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

! Windows, dialog.
submodule (windows) dialog
  use interfaces_cimgui
  implicit none

contains

  !> Draw the open files dialog.
  module subroutine draw_dialog(w)
    use commands, only: com, idcom, icom
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use windows, only: win
    use systems, only: add_systems_from_name, launch_initialization_thread,&
       system_shorten_names, sys_init, sys, ok_system
    use c_interface_module, only: C_F_string_alloc, c_free
    use tools_io, only: ferror, faterr, fopen_write, fclose, uout
    use param, only: dirsep, bohrtoa
    class(window), intent(inout), target :: w

    type(ImVec2) :: minsize, maxsize
    type(IGFD_Selection_Pair), pointer :: s(:)
    type(IGFD_Selection) :: sel
    type(c_ptr) :: cstr
    integer(c_size_t) :: i
    character(len=:), allocatable :: name, path, str, errmsg
    logical :: readlastonly, doquit
    integer :: lu, ios, idx

    ! set initial, minimum, and maximum sizes
    minsize%x = 0._c_float
    minsize%y = 0._c_float
    maxsize%x = 1e10_c_float
    maxsize%y = 1e10_c_float

    ! handle cases when the dialog must quit
    doquit = .false.
    if (w%dialog_purpose == wpurp_dialog_openvibfile) then
       ! open vibrations file dialog => quit if the system is gone
       doquit = (.not.ok_system(w%isys,sys_init))
    elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile.or.w%dialog_purpose == wpurp_dialog_saveimagefile.or.&
       w%dialog_purpose == wpurp_dialog_openfieldfile.or.w%dialog_purpose == wpurp_dialog_openonefilemodal) then
       ! open library file, save image file, open field file, open one file modal
       ! => quit if the caller window is gone
       doquit = .not.win(w%idparent)%isinit
    end if

    ! process the dialog
    if (IGFD_DisplayDialog(w%dptr,c_loc(w%name),w%flags,minsize,maxsize)) then
       ! the dialog has been closed
       if (IGFD_IsOk(w%dptr).and..not.doquit) then
          ! with an OK, gather information
          if (w%dialog_purpose == wpurp_dialog_openfiles) then
             !! open files dialog !!
             ! open all files selected and add the new systems
             sel = IGFD_GetSelection(w%dptr)
             call c_f_pointer(sel%table,s,(/sel%count/))
             cstr = IGFD_GetCurrentPath(w%dptr)
             call C_F_string_alloc(cstr,path)
             call c_free(cstr)
             do i = 1, sel%count
                call C_F_string_alloc(s(i)%fileName,name)
                name = trim(path) // dirsep // trim(name)
                readlastonly = w%dialog_data%readlastonly
                call add_systems_from_name(name,w%dialog_data%mol,w%dialog_data%isformat,&
                   readlastonly,w%dialog_data%rborder/bohrtoa,logical(w%dialog_data%molcubic))
             end do

             ! initialize
             call launch_initialization_thread()

             ! shorten the names
             call system_shorten_names()

          elseif (w%dialog_purpose == wpurp_dialog_savelogfile) then
             !! save log file dialog !!
             cstr = IGFD_GetFilePathName(w%dptr)
             call C_F_string_alloc(cstr,name)
             call c_free(cstr)

             lu = fopen_write(name,errstop=.false.)
             if (lu < 0) then
                write (uout,'("Could not open file for writing: ",A)') trim(name)
             else
                if (idcom == 0 .and. allocated(outputb)) then
                   write(lu,'(A)',iostat=ios) outputb(1:lob)
                elseif (idcom > 0) then
                   write(lu,'(A)',iostat=ios) com(icom(idcom))%output
                end if
                if (ios /= 0) then
                   write (uout,'("Error writing to file: ",A)') trim(name)
                end if
                call fclose(lu)
             end if
          elseif (w%dialog_purpose == wpurp_dialog_openlibraryfile.or.&
             w%dialog_purpose == wpurp_dialog_openfieldfile.or.w%dialog_purpose == wpurp_dialog_openonefilemodal) then
             ! !! new structure file dialog !!
             ! cstr = IGFD_GetFilePathName(w%dptr)
             ! call C_F_string_alloc(cstr,name)
             ! call c_free(cstr)
             ! w%okfile = trim(name)
             win(w%idparent)%okfile_set = .true.
             win(w%idparent)%okfile_read = .true.
             win(w%idparent)%okfile = ""
             win(w%idparent)%itoken = w%itoken

             ! open all files selected and add the new systems
             sel = IGFD_GetSelection(w%dptr)
             call c_f_pointer(sel%table,s,(/sel%count/))
             cstr = IGFD_GetCurrentPath(w%dptr)
             call C_F_string_alloc(cstr,path)
             call c_free(cstr)
             do i = 1, sel%count
                call C_F_string_alloc(s(i)%fileName,name)
                name = trim(path) // dirsep // trim(name)
                win(w%idparent)%okfile = win(w%idparent)%okfile // name
                if (i < sel%count) &
                   win(w%idparent)%okfile = win(w%idparent)%okfile // c_null_char
             end do
             win(w%idparent)%okfile_format = w%dialog_data%isformat

          elseif (w%dialog_purpose == wpurp_dialog_openvibfile) then
             w%okfile = ""

             ! open all files selected and add the new systems
             sel = IGFD_GetSelection(w%dptr)
             call c_f_pointer(sel%table,s,(/sel%count/))
             cstr = IGFD_GetCurrentPath(w%dptr)
             call C_F_string_alloc(cstr,path)
             call c_free(cstr)
             do i = 1, sel%count
                call C_F_string_alloc(s(i)%fileName,name)
                name = trim(path) // dirsep // trim(name)
                w%okfile = w%okfile // name // c_null_char
             end do

             ! open the vibrations file
             if (ok_system(w%isys,sys_init)) then
                call sys(w%isys)%c%vib%end()
                str = w%okfile
                do while (.true.)
                   idx = index(str,c_null_char)
                   if (idx == 0) exit
                   call sys(w%isys)%c%vib%read_file(sys(w%isys)%c,str(1:idx-1),"",w%dialog_data%isformat,errmsg)
                   if (len_trim(w%errmsg) > 0) &
                      write (uout,'(A)') errmsg
                   str = str(idx+1:)
                end do
             end if

          elseif (w%dialog_purpose == wpurp_dialog_saveimagefile) then
             !! save image file dialog !!
             win(w%idparent)%okfile_set = .true.

             cstr = IGFD_GetFilePathName(w%dptr)
             call C_F_string_alloc(cstr,name)
             call c_free(cstr)
             win(w%idparent)%okfile = trim(name)

             cstr = IGFD_GetCurrentFilter(w%dptr)
             call C_F_string_alloc(cstr,name)
             call c_free(cstr)
             win(w%idparent)%okfilter = trim(name)
          else
             call ferror('draw_dialog','unknown dialog purpose',faterr)
          end if
       end if

       ! close the dialog and terminate the window
       call IGFD_CloseDialog(w%dptr)
       call w%end()
    end if

    ! exit if focused and received the close keybinding, or if forced by some other window
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) .or. &
       is_bind_event(BIND_CLOSE_ALL_DIALOGS) .or. doquit) &
       call IGFD_ForceQuit(w%dptr)

    ! exit if focused and received the OK keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) &
       call IGFD_ForceOK(w%dptr)

  end subroutine draw_dialog

end submodule dialog
