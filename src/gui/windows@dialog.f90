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
    use tools_io, only: ferror, faterr, string, fopen_write, fclose, uout
    use param, only: dirsep, bohrtoa
    class(window), intent(inout), target :: w

    type(ImVec2) :: minsize, maxsize
    type(IGFD_Selection_Pair), pointer :: s(:)
    type(IGFD_Selection) :: sel
    type(c_ptr) :: cstr
    integer(c_size_t) :: i
    character(len=:), allocatable :: name, path, str, errmsg
    logical :: readlastonly, doquit
    integer :: lu, ios, idx, idp

    ! set initial, minimum, and maximum sizes
    minsize%x = 0._c_float
    minsize%y = 0._c_float
    maxsize%x = 1e10_c_float
    maxsize%y = 1e10_c_float

    ! the window that opened this dialog, if it is still there
    idp = w%parent()

    ! handle cases when the dialog must quit
    doquit = .false.
    if (w%purpose == wpurp_dialog_openvibfile) then
       ! open vibrations file dialog => quit if the system is gone
       doquit = (.not.ok_system(w%isys,sys_init))
    elseif (w%purpose == wpurp_dialog_openlibraryfile.or.w%purpose == wpurp_dialog_saveimagefile.or.&
       w%purpose == wpurp_dialog_openfieldfile.or.w%purpose == wpurp_dialog_openonefilemodal.or.&
       w%purpose == wpurp_dialog_savefile.or.w%purpose == wpurp_dialog_selectdir) then
       ! open library file, save image file, open field file, open one file modal,
       ! save structure file, select directory => quit if the caller window is gone
       doquit = .true.
       doquit = (idp == 0)
    end if

    ! process the dialog
    if (IGFD_DisplayDialog(w%dptr,c_loc(w%name),w%flags,minsize,maxsize)) then
       ! the dialog has been closed
       if (IGFD_IsOk(w%dptr).and..not.doquit) then
          ! with an OK, gather information
          if (w%purpose == wpurp_dialog_openfiles) then
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
             if (sel%count == 0) then
                ! nothing clicked in the list: use the typed file name
                name = typed_filename()
                if (len(name) > 0) then
                   readlastonly = w%dialog_data%readlastonly
                   call add_systems_from_name(name,w%dialog_data%mol,w%dialog_data%isformat,&
                      readlastonly,w%dialog_data%rborder/bohrtoa,logical(w%dialog_data%molcubic))
                end if
             end if

             ! initialize
             call launch_initialization_thread()

             ! shorten the names
             call system_shorten_names()

          elseif (w%purpose == wpurp_dialog_savelogfile) then
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
          elseif (w%purpose == wpurp_dialog_openlibraryfile.or.&
             w%purpose == wpurp_dialog_openfieldfile.or.w%purpose == wpurp_dialog_openonefilemodal) then
             ! !! new structure file dialog !!
             ! cstr = IGFD_GetFilePathName(w%dptr)
             ! call C_F_string_alloc(cstr,name)
             ! call c_free(cstr)
             ! w%okfile = trim(name)
             win(idp)%okfile_set = .true.
             win(idp)%okfile_read = .true.
             win(idp)%okfile = ""
             win(idp)%itoken = w%itoken

             ! open all files selected and add the new systems
             sel = IGFD_GetSelection(w%dptr)
             call c_f_pointer(sel%table,s,(/sel%count/))
             cstr = IGFD_GetCurrentPath(w%dptr)
             call C_F_string_alloc(cstr,path)
             call c_free(cstr)
             do i = 1, sel%count
                call C_F_string_alloc(s(i)%fileName,name)
                name = trim(path) // dirsep // trim(name)
                win(idp)%okfile = win(idp)%okfile // name
                if (i < sel%count) &
                   win(idp)%okfile = win(idp)%okfile // c_null_char
             end do
             if (sel%count == 0) then
                ! nothing clicked in the list: use the typed file name
                win(idp)%okfile = typed_filename()
             end if
             win(idp)%okfile_format = w%dialog_data%isformat

          elseif (w%purpose == wpurp_dialog_openvibfile) then
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
             if (sel%count == 0) then
                ! nothing clicked in the list: use the typed file name
                name = typed_filename()
                if (len(name) > 0) w%okfile = name // c_null_char
             end if

             ! open the vibrations file
             if (ok_system(w%isys,sys_init)) then
                call sys(w%isys)%c%vib%end()
                str = w%okfile
                do while (.true.)
                   idx = index(str,c_null_char)
                   if (idx == 0) exit
                   call sys(w%isys)%c%vib%read_file(sys(w%isys)%c,str(1:idx-1),"",w%dialog_data%isformat,errmsg)
                   if (len_trim(errmsg) > 0) &
                      write (uout,'(A)') "WARNING : " // trim(errmsg)
                   str = str(idx+1:)
                end do
             end if

          elseif (w%purpose == wpurp_dialog_saveimagefile.or.w%purpose == wpurp_dialog_savefile) then
             !! save image file or save structure file dialog !!
             win(idp)%okfile_set = .true.

             cstr = IGFD_GetFilePathName(w%dptr)
             call C_F_string_alloc(cstr,name)
             call c_free(cstr)
             win(idp)%okfile = trim(name)

             ! the filter carries the image format (save image file only)
             if (w%purpose == wpurp_dialog_saveimagefile) then
                cstr = IGFD_GetCurrentFilter(w%dptr)
                call C_F_string_alloc(cstr,name)
                call c_free(cstr)
                win(idp)%okfilter = trim(name)
             end if
          elseif (w%purpose == wpurp_dialog_selectdir) then
             !! select directory dialog !!
             ! in directory mode the dialog has no file name, so the
             ! selected directory is the current path
             win(idp)%okfile_set = .true.

             cstr = IGFD_GetCurrentPath(w%dptr)
             call C_F_string_alloc(cstr,path)
             call c_free(cstr)
             win(idp)%okfile = trim(path)
          else
             call ferror('draw_dialog','unknown dialog purpose: ' // string(w%purpose),faterr)
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

  contains
    !> The file name typed in the dialog's File Name box, for when the
    !> selection is empty (no row was clicked in the list). An absolute
    !> typed name is used as is; a relative one is composed with the
    !> current path by the dialog. Empty if nothing was typed.
    function typed_filename() result(fname)
      character(len=:), allocatable :: fname

      type(c_ptr) :: cstr_
      logical :: isabs

      cstr_ = IGFD_GetCurrentFileName(w%dptr)
      call C_F_string_alloc(cstr_,fname)
      call c_free(cstr_)
      fname = trim(fname)
      if (len(fname) == 0) return

      isabs = (fname(1:1) == "/" .or. fname(1:1) == dirsep)
      if (.not.isabs .and. len(fname) >= 2) isabs = (fname(2:2) == ":")
      if (.not.isabs) then
         cstr_ = IGFD_GetFilePathName(w%dptr)
         call C_F_string_alloc(cstr_,fname)
         call c_free(cstr_)
         fname = trim(fname)
      end if

    end function typed_filename
  end subroutine draw_dialog

end submodule dialog
