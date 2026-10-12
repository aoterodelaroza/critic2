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

! Persistent user preferences and window layout, in the standard
! per-user directories of each operating system (user_file). The
! settings of the preferences window are kept in a JSON file, saved
! automatically. The window layout (the ImGui ini plus which permanent
! windows are open) is read at startup and written only when the user
! asks for it.
submodule (gui_main) settings
  implicit none

  ! version of the settings file format written
  integer, parameter :: settings_file_version = 1

  ! When the settings file was rejected (unreadable, or a newer
  ! format), the preferences in memory at that point, serialized. The
  ! file is not overwritten while they stay unchanged.
  logical :: rejected = .false.
  character(len=:), allocatable :: rejected_json

  ! section of the layout file with the open state of the permanent
  ! windows (ImGui skips sections it does not know), and the keys in it
  character(len=*), parameter :: layout_section = "[critic2][Windows]"
  character(len=*), parameter :: layout_keys(nlayout) = (/character(len=6) ::&
     "Tree","View","Input","Output"/)

  !xx! private procedures
  ! function user_file(state,basename)
  ! function xdg_dir(var,suffix)
  ! function settings_string()
  ! subroutine write_file(file,str,what)
  ! subroutine settings_xfer(x,p)

contains

  !> Read the user preferences from the settings file, if there is
  !> one. A missing or unreadable key keeps its current value (the
  !> defaults, at startup). A file that cannot be parsed or has a
  !> newer format is reported and ignored, and it is not overwritten
  !> unless the preferences are changed.
  module subroutine settings_load()
    use json_module, only: json_core, json_value
    use jsonxfer, only: json_xfer
    use tools_io, only: ferror, warning
    type(json_core), target :: json
    type(json_value), pointer :: p
    type(json_xfer) :: x
    character(len=:), allocatable :: file, errmsg
    integer :: iver
    logical :: found, ok

    if (.not.persist_settings) return
    file = user_file(.false.,"settings.json")
    if (len(file) == 0) return
    inquire(file=file,exist=found)
    if (.not.found) return

    nullify(p)
    call json%initialize()
    call json%load(file=file,p=p)
    if (json%failed()) then
       call json%check_for_errors(ok,errmsg)
       rejected = .true.
    else
       call json%get(p,'critic2_settings',iver,found)
       if (.not.found .or. iver > settings_file_version) then
          errmsg = "unknown format"
          rejected = .true.
       end if
    end if

    if (rejected) then
       call ferror('settings_load','Settings file ' // file // ' not read (' // errmsg //&
          '); it will not be overwritten unless the preferences are changed',warning)
       rejected_json = settings_string()
    else
       x%json => json
       x%reading = .true.
       call settings_xfer(x,p)
    end if

    if (json%failed()) call json%clear_exceptions()
    if (associated(p)) call json%destroy(p)

  end subroutine settings_load

  !> Write the current user preferences to the settings file,
  !> creating its directory if necessary. If the file was rejected at
  !> startup, it is written only once the preferences change.
  module subroutine settings_save()
    character(len=:), allocatable :: file, str

    if (.not.persist_settings) return
    file = user_file(.false.,"settings.json")
    if (len(file) == 0) return

    str = settings_string()
    if (rejected) then
       if (str == rejected_json) return
       rejected = .false.
    end if
    call write_file(file,str,"Settings")

  end subroutine settings_save

  !> Reset all preferences to their default values (as if there were
  !> no settings file), apply them to the open scenes, and save them,
  !> also over a rejected file.
  module subroutine settings_reset()

    call set_default_ui_settings()
    call apply_color_settings()
    rejected = .false.
    call settings_save()

  end subroutine settings_reset

  !> Apply the element, critical point, and gradient path colors to
  !> the representations in the open systems and alternate views.
  module subroutine apply_color_settings()
    use systems, only: nsys, sysc
    use windows, only: nwin, win, wintype_view
    integer :: i

    do i = 1, nsys
       call sysc(i)%sc%reset_atom_colors()
    end do
    do i = 1, nwin
       if (win(i)%type == wintype_view .and. .not.win(i)%ismain .and. associated(win(i)%sc)) &
          call win(i)%sc%reset_atom_colors()
    end do

  end subroutine apply_color_settings

  !> The window indices of the permanent windows whose open state is
  !> kept in the layout, in the order of layout_keys.
  module function layout_windows() result(iw)
    use windows, only: iwin_tree, iwin_view, iwin_console_input, iwin_console_output
    integer :: iw(nlayout)

    iw = (/iwin_tree,iwin_view,iwin_console_input,iwin_console_output/)

  end function layout_windows

  !> Read the window layout saved by the user, if there is one: give
  !> ImGui the window positions and docking, and open or close the
  !> permanent windows as they were. Call before the first frame.
  module subroutine layout_load()
    use iso_c_binding, only: c_char, c_null_char, c_loc
    use interfaces_cimgui, only: igLoadIniSettingsFromDisk
    use windows, only: win
    use tools_io, only: fopen_read, getline_raw, fclose
    character(kind=c_char,len=:), allocatable, target :: cfile
    character(len=:), allocatable :: file, line
    integer :: lu, k, ieq, iw(nlayout)
    logical :: found, insection

    if (.not.persist_settings) return
    file = user_file(.true.,"imgui.ini")
    if (len(file) == 0) return
    inquire(file=file,exist=found)
    if (.not.found) return

    ! the ImGui settings
    cfile = file // c_null_char
    call igLoadIniSettingsFromDisk(c_loc(cfile))

    ! open or close the permanent windows; blanks and carriage
    ! returns are ignored
    lu = fopen_read(file,errstop=.false.)
    if (lu < 0) return
    iw = layout_windows()
    insection = .false.
    do while (getline_raw(lu,line))
       line = strip_blanks(line)
       if (len(line) == 0) cycle
       if (line(1:1) == "[") then
          insection = (line == layout_section)
       elseif (insection) then
          ieq = index(line,"=")
          if (ieq == 0) cycle
          do k = 1, nlayout
             if (line(1:ieq-1) == trim(layout_keys(k))) win(iw(k))%isopen = (line(ieq+1:) == "1")
          end do
       end if
    end do
    call fclose(lu)

  contains
    ! str without spaces, tabs, and carriage returns
    function strip_blanks(str) result(res)
      character*(*), intent(in) :: str
      character(len=:), allocatable :: res

      integer :: i

      res = ""
      do i = 1, len(str)
         if (index(" " // char(9) // char(13),str(i:i)) == 0) res = res // str(i:i)
      end do

    end function strip_blanks
  end subroutine layout_load

  !> Write the current window layout to the layout file, creating
  !> its directory if necessary.
  module subroutine layout_save()
    use iso_c_binding, only: c_size_t, c_ptr
    use c_interface_module, only: c_f_string_alloc
    use interfaces_cimgui, only: igSaveIniSettingsToMemory
    use windows, only: win
    use param, only: newline
    character(len=:), allocatable :: file, str
    integer(c_size_t) :: n
    integer :: k, iw(nlayout)
    type(c_ptr) :: cstr

    if (.not.persist_settings) return
    file = user_file(.true.,"imgui.ini")
    if (len(file) == 0) return

    ! the ImGui settings, followed by the open state of the permanent windows
    cstr = igSaveIniSettingsToMemory(n)
    call c_f_string_alloc(cstr,str)
    str = str // layout_section // newline
    iw = layout_windows()
    do k = 1, nlayout
       str = str // trim(layout_keys(k)) // "=" // merge("1","0",logical(win(iw(k))%isopen)) // newline
    end do
    call write_file(file,str,"Window layout")

  end subroutine layout_save

  !xx! private procedures

  !> Path of the file basename in the critic2 directory under the
  !> standard per-user configuration directory (state = .false.) or
  !> state directory (.true.) of each operating system:
  !>   Linux: $XDG_CONFIG_HOME (~/.config) or $XDG_STATE_HOME (~/.local/state)
  !>   Windows: %APPDATA% or %LOCALAPPDATA%
  !>   macOS: ~/Library/Application Support (both)
  !> Empty if the directory cannot be found. It is not created.
  function user_file(state,basename) result(file)
    use tools_io, only: getenv_alloc
    use param, only: dirsep
    logical, intent(in) :: state
    character*(*), intent(in) :: basename
    character(len=:), allocatable :: file

    character(len=:), allocatable :: dir

#if defined(_WIN32)
    dir = ""
    if (state) dir = getenv_alloc("LOCALAPPDATA")
    if (len(dir) == 0) dir = getenv_alloc("APPDATA")
#elif defined(__APPLE__)
    dir = getenv_alloc("HOME")
    if (len(dir) > 0) dir = dir // "/Library/Application Support"
#else
    if (state) then
       dir = xdg_dir("XDG_STATE_HOME","/.local/state")
    else
       dir = xdg_dir("XDG_CONFIG_HOME","/.config")
    end if
#endif
    file = ""
    if (len(dir) > 0) file = dir // dirsep // "critic2" // dirsep // basename

  end function user_file

#if !defined(_WIN32) && !defined(__APPLE__)
  !> The XDG base directory in environment variable var or, if it is
  !> not set or not an absolute path (XDG spec), its default: the home
  !> directory followed by suffix. Empty if neither is available.
  function xdg_dir(var,suffix) result(dir)
    use tools_io, only: getenv_alloc
    character*(*), intent(in) :: var, suffix
    character(len=:), allocatable :: dir

    dir = getenv_alloc(var)
    if (len(dir) > 0) then
       if (dir(1:1) /= "/") dir = ""
    end if
    if (len(dir) == 0) then
       dir = getenv_alloc("HOME")
       if (len(dir) > 0) dir = dir // suffix
    end if

  end function xdg_dir
#endif

  !> The current user preferences, as the contents of a settings file.
  function settings_string() result(str)
    use json_module, only: json_core, json_value
    use jsonxfer, only: json_xfer
    character(len=:), allocatable :: str

    type(json_core), target :: json
    type(json_value), pointer :: p
    type(json_xfer) :: x

    call json%initialize(compress_vectors=.true.)
    call json%create_object(p,'')
    call json%add(p,'critic2_settings',settings_file_version)
    x%json => json
    x%reading = .false.
    call settings_xfer(x,p)
    call json%serialize(p,str)
    call json%destroy(p)

  end function settings_string

  !> Write the string str to file, creating its directory if
  !> necessary. Warn, naming the contents (what), if it fails.
  subroutine write_file(file,str,what)
    use tools_io, only: falloc, fclose, ferror, warning, mkpath
    use param, only: dirsep
    character*(*), intent(in) :: file, str, what

    character(len=:), allocatable :: errmsg
    integer :: lu, ios

    call mkpath(file(1:index(file,dirsep,.true.)-1),errmsg)
    if (len(errmsg) > 0) then
       call ferror('write_file',what // ' not saved: ' // errmsg,warning)
       return
    end if
    lu = falloc()
    open(unit=lu,file=file,access='stream',form='unformatted',action='write',status='replace',iostat=ios)
    if (ios == 0) write(lu,iostat=ios) str
    call fclose(lu)
    if (ios /= 0) call ferror('write_file',what // ' not saved to ' // file,warning)

  end subroutine write_file

  !> Transfer the user preferences between the variables and the json
  !> object p, in the direction given by x%reading. Grouped as in the
  !> categories of the preferences window.
  subroutine settings_xfer(x,p)
    use json_module, only: json_value
    use jsonxfer, only: json_xfer
    use keybindings, only: BIND_NUM, keybind, modbind, bindnames, groupbind, groupnames,&
       set_bind, is_valid_bind
    use representations, only: iso_defaultlevel, iso_nlevel, cps_abbrev
    use tools_io, only: string, nameguess, lower
    use systems, only: always_read_virtuals
    use interfaces_cimgui, only: ImGuiKey_NamedKey_BEGIN, ImGuiKey_NamedKey_COUNT
    type(json_xfer), intent(inout) :: x
    type(json_value), pointer, intent(in) :: p

    type(json_value), pointer :: q, c
    logical :: found
    integer :: i, iv(2), ikeys(2)
    integer, allocatable :: ivv(:)
    character(len=:), allocatable :: name

    ! names of the tree background colors, by system kind (ColorTableCellBg)
    character(len=*), parameter :: treebg_name(0:8) = (/character(len=24) ::&
       "crystal","layered_crystal","chain_crystal","molecular_crystal","surface",&
       "chain","molecule_in_box","molecule","molecular_cluster"/)

    ! interface
    call x%obj(p,'interface',q,found)
    if (found) then
       call x%x(q,'font_scale',io%FontGlobalScale,0.3d0,3d0)
       call x%x(q,'tooltips',tooltip_enabled)
       call x%x(q,'toolbar_labels',toolbar_labels)
       call x%x(q,'tooltip_delay',tooltip_delay,0d0,5d0)
       call x%x(q,'tooltip_width',tooltip_wrap_factor,0d0,1000d0)
       call x%x(q,'tree_select_updates_input_console',tree_select_updates_inpcon)
       call x%x(q,'tree_select_updates_view',tree_select_updates_view)
       call x%x(q,'isosurface_level',iso_defaultlevel,1,iso_nlevel)
    end if

    ! key bindings, each under the name of its group and its
    ! description, as [key, modifiers]. The keys are ImGui key codes,
    ! so the block is skipped if the ImGui key enumeration changed.
    ! The names have characters with meaning in json-fortran paths,
    ! so the entries are looked up by child name, not by path.
    call x%obj(p,'keybindings',q,found)
    if (found) then
       ikeys = (/ImGuiKey_NamedKey_BEGIN,ImGuiKey_NamedKey_COUNT/)
       ! the stamp must be in the file to trust its key codes
       if (x%reading) then
          iv = -1
       else
          iv = ikeys
       end if
       call x%x(q,'imgui_keys',iv)
       if (all(iv == ikeys)) then
          do i = 1, BIND_NUM
             name = trim(groupnames(groupbind(i))) // " / " // trim(bindnames(i))
             if (.not.x%reading) then
                call x%json%add(q,name,(/keybind(i),modbind(i)/))
                cycle
             end if
             call x%json%get_child(q,name,c,found)
             if (.not.found) cycle
             call x%json%get(c,ivv)
             if (x%json%failed()) then
                call x%json%clear_exceptions()
                cycle
             end if
             if (size(ivv) /= 2) cycle
             if (is_valid_bind(i,ivv(1),ivv(2))) call set_bind(i,ivv(1),ivv(2))
          end do
       end if
    end if

    ! colors, one hex string per name
    call x%obj(p,'colors',q,found)
    if (found) then
       do i = 0, 8
          call xcolor(q,'tree_' // trim(treebg_name(i)),ColorTableCellBg(:,i))
       end do
       call xcolor(q,'tree_selected_border',ColorTableSelectedBorder)
       call xcolor(q,'tree_multiselected_rows',ColorTableHighlightRow)
       call xcolor(q,'tree_reference_field_border',ColorFieldSelected)
       call xcolor(q,'hovered_atom',ColorHighlightScene)
       call xcolor(q,'selected_atom',ColorHighlightSelectScene)
       do i = 1, 4
          call xcolor(q,'measure_atom_' // string(i),ColorMeasureSelect(:,i))
       end do
       call xcolor(q,'gradient_paths',ColorGpath)
    end if
    call x%obj(p,'element_colors',q,found)
    if (found) then
       do i = 0, maxzat0
          call xcolor(q,trim(nameguess(i,.true.)),ColorElement(:,i))
       end do
    end if
    call x%obj(p,'critical_point_colors',q,found)
    if (found) then
       do i = 0, 3
          call xcolor(q,trim(cps_abbrev(i)),ColorCP(:,i))
       end do
    end if

    ! readers
    call x%obj(p,'readers',q,found)
    if (found) then
       call x%x(q,'always_read_virtuals',always_read_virtuals)
    end if

  contains
    ! transfer the color v, as a hex string in the file (#rrggbb or #rrggbbaa)
    subroutine xcolor(q,key,v)
      type(json_value), pointer, intent(in) :: q
      character*(*), intent(in) :: key
      real(c_float), intent(inout) :: v(:)

      character(len=:), allocatable :: str
      integer :: k, ic(size(v)), ios

      if (x%reading) then
         str = ""
         call x%xs(q,key,str)
         if (len(str) /= 1 + 2 * size(v)) return
         if (str(1:1) /= "#") return
         if (verify(str(2:),"0123456789abcdefABCDEF") /= 0) return
         do k = 1, size(v)
            read(str(2*k:2*k+1),'(Z2)',iostat=ios) ic(k)
            if (ios /= 0) return
         end do
         v = real(ic,c_float) / 255._c_float
      else
         allocate(character(len=1+2*size(v)) :: str)
         str(1:1) = "#"
         do k = 1, size(v)
            write(str(2*k:2*k+1),'(Z2.2)') nint(min(max(v(k),0._c_float),1._c_float) * 255._c_float)
         end do
         str = lower(str)
         call x%xs(q,key,str)
      end if

    end subroutine xcolor
  end subroutine settings_xfer

end submodule settings
