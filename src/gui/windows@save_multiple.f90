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

! Routines for the save multiple structures window
submodule (windows) save_multiple
  use interfaces_cimgui
  use param, only: mlen
  implicit none

  ! columns of the file table; the ids double as the column indices
  integer(c_int), parameter :: ic_savemult_id = 0
  integer(c_int), parameter :: ic_savemult_type = 1
  integer(c_int), parameter :: ic_savemult_name = 2
  integer(c_int), parameter :: ic_savemult_file = 3
  integer(c_int), parameter :: ic_savemult_NUMCOLUMNS = 4

  ! characters replaced by an underscore in the name of a system (%n)
  character(len=*), parameter :: badchars = " |/:*?<>" //&
     achar(92) // achar(34) // achar(39) // char(9)

  ! modulus for the hash of the system list (keeps it well inside int8)
  integer*8, parameter :: sighashmod = 2147483647_8

  ! settings remembered from the last successful save (this session).
  ! The 50 rklength default is a deliberate GUI choice (the CLI writers
  ! default to 40 when no rklength is given).
  integer(c_int) :: lastfmt = 1
  character(len=:), allocatable :: lastpattern
  real(c_float) :: lastrk = 50._c_float
  logical :: lastnosym = .false.
  logical :: lastcartesian = .false.
  logical :: lastdocell = .true.

contains

  !> Draw the save multiple structures window
  module subroutine draw_save_multiple(w)
    use windows, only: wintype_dialog, wpurp_dialog_selectdir
    use systems, only: sys, sysc, sys_init, ok_system, are_threads_running,&
       launch_initialization_thread
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_calcheight, iw_tooltip,&
       iw_checkbox, iw_combo_simple, iw_dragfloat_realc, iw_close_event,&
       iw_setpos_bottomright, iw_inputtext, iw_table_column, iw_radiobutton,&
       get_current_working_dir
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG
    use icons, only: icon_tex, icon_ui_cell, icon_vm_bond, rgba_fam_crys, rgba_fam_molgeom
    use gui_main, only: fontsize
    use tools_io, only: string, uout
    use param, only: dirsep
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1
    character(len=:), allocatable :: file, firsterr, warnmsg, savemsg
    integer :: i, isys, ifmt, iaux, nfail, nwritten, nskip, nleft, nbot
    integer(c_int) :: flags
    type(ImVec2) :: uv0, uv1
    type(ImVec4) :: nobord
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    logical :: dooverwrite
    logical :: doquit, ok, okvalid, userk, usecell, changed, dirok
    logical :: anycrystal, thisrk
    logical(c_bool) :: ldum
    type(ImVec2) :: sz

    logical, save :: ttshown = .false. ! tooltip flag

    ! whole-texture image, no border (for the kind icons in the table)
    uv0 = ImVec2(0._c_float,0._c_float)
    uv1 = ImVec2(1._c_float,1._c_float)
    nobord = ImVec4(0._c_float,0._c_float,0._c_float,0._c_float)

    ! build the format combo options on first use; unlike the save-as
    ! window there is no auto-detect entry, as the pattern carries no
    ! extension, so skip the first entry of the shared combo string
    call build_write_format_combo()

    ! this window works on the tree, so it closes when the tree is gone
    doquit = (w%parent() == 0)

    ! initialize state; the cache is stale on reopen (the files on disk
    ! may have changed in the meantime), so drop it
    if (w%firstpass) then
       if (allocated(w%sm%lastsig)) deallocate(w%sm%lastsig)
       w%errmsg = ""
       w%okmsg = ""
       w%savemult_format = lastfmt
       w%savemult_rk = lastrk
       w%savemult_nosym = lastnosym
       w%savemult_cartesian = lastcartesian
       w%savemult_docell = lastdocell
       w%savemult_exist = 0
       w%sm%novr = 0
       w%sm%editsys = 0
       if (allocated(lastpattern)) then
          w%savemult_pattern = lastpattern
       else
          w%savemult_pattern = "%f-*"
       end if
       ! write to the last directory used by a save window, or to the cwd
       w%okfile = ""
       if (allocated(okfile_lastdir)) w%okfile = okfile_lastdir
       if (len_trim(w%okfile) == 0) w%okfile = get_current_working_dir()
       ! use the tree selection if the user made one; a single selected
       ! system is just the current one, and a poor default for this window
       w%savemult_scope = 0
       if (tree_nselected() <= 1) w%savemult_scope = 1
    end if

    ! the systems being written
    call iw_text("Systems",highlight=.true.)
    str1 = "Selected in the tree" // c_null_char // "All systems in the tree" // c_null_char
    call igPushItemWidth(iw_calcwidth(24,1))
    call iw_combo_simple("##savemultscope",str1,w%savemult_scope)
    call igPopItemWidth()
    call iw_tooltip("Systems whose structures will be written. Select systems in the&
       & tree with control-click, shift-click, or the buttons",ttshown)

    ! the output directory
    call iw_text("Directory",highlight=.true.)
    ldum = iw_inputtext("##savemultdir",bufsize=1023,texta=w%okfile,width=36,notlive=.true.)
    call iw_tooltip("Directory where the structure files will be written",ttshown)
    if (iw_button("Browse...",sameline=.true.)) &
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_selectdir,idparent=w%id,orraise=-1)
    call iw_tooltip("Choose the directory in the file browser",ttshown)
    dirok = (len_trim(w%okfile) > 0)
    if (dirok) then
       ! a missing directory is only a warning: inquire() on a directory is
       ! not portable, and the write itself reports the real error. Cached
       ! on okfile, like the other existence checks: this runs every frame,
       ! and a stat() per frame stalls the render thread on a slow mount
       if (.not.allocated(w%savemult_lastcheck)) w%savemult_lastcheck = ""
       if (w%okfile /= w%savemult_lastcheck) then
          w%savemult_lastcheck = w%okfile
          inquire(file=savedir(w),exist=w%savemult_direxists)
       end if
       if (.not.w%savemult_direxists) &
          call iw_text("This directory does not exist",danger=.true.,wrap=.true.)
    else
       call iw_text("Choose a directory to write the files to",danger=.true.,wrap=.true.)
    end if

    ! format combo
    call iw_text("Format",highlight=.true.)
    call igPushItemWidth(iw_calcwidth(45,1))
    call iw_combo_simple("##savemultformat",write_format_combostr(icombo_fmt1:),w%savemult_format,startsatone=.true.)
    call igPopItemWidth()
    call iw_tooltip("File format for the written structures",ttshown)
    ifmt = fmtperm(max(min(w%savemult_format,isformat_w_max),1))

    ! the file name pattern
    call iw_text("File Names",highlight=.true.)
    call draw_pattern_help()
    ldum = iw_inputtext("##savemultpattern",bufsize=1023,texta=w%savemult_pattern,width=36,&
       notlive=.true.)
    call iw_tooltip("Pattern for the file w%sm%names, without extension&
       & (see the ? for the substitutions)",ttshown)

    ! rebuild the list of target file w%sm%names, if the settings changed
    call rebuild_names(w,ifmt,changed)
    if (changed) then
       w%okmsg = ""
       w%errmsg = ""
       ! the user asked for these systems, so get the ones that are only
       ! waiting for a worker moving
       if (w%sm%nloading > 0 .and. .not.are_threads_running()) call launch_initialization_thread()
    end if

    ! the list can mix molecules and crystals, so the options that only
    ! apply to a crystal are offered when at least one of the systems is
    ! one; they are then used only for those systems
    anycrystal = .false.
    do i = 1, w%sm%nnames
       if (.not.sys(w%sm%nameidx(i))%c%ismolecule) then
          anycrystal = .true.
          exit
       end if
    end do

    ! the k-point grid (rklength) applies to espresso, CASTEP, and to
    ! FHIaims crystals (for a molecule, FHIaims would write a
    ! meaningless k-point grid to the companion control file)
    userk = (ifmt == isformat_w_qein .or. ifmt == isformat_w_castepcell .or.&
       (ifmt == isformat_w_aimsin .and. anycrystal))

    ! the unit cell can be shown in 3D model files (crystals only)
    usecell = (ifmt == isformat_w_obj .or. ifmt == isformat_w_ply .or. ifmt == isformat_w_off)&
       .and. anycrystal

    ! options, only for the formats that support them
    if (userk .or. usecell .or. ifmt == isformat_w_cif .or. ifmt == isformat_w_crystal .or.&
       ifmt == isformat_w_shelx) &
       call iw_text("Options",highlight=.true.)

    ! k-point grid (rklength)
    if (userk) then
       ldum = iw_dragfloat_realc("RKlength##savemultrk",x1=w%savemult_rk,speed=1._c_float,&
          min=1._c_float,max=200._c_float,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Length parameter for calculating the k-point grid",ttshown)
    end if

    ! do not use symmetry
    if (ifmt == isformat_w_cif .or. ifmt == isformat_w_crystal .or.&
       ifmt == isformat_w_shelx) then
       ldum = iw_checkbox("Do not use symmetry##savemultnosym",w%savemult_nosym)
       call iw_tooltip("Write the structures without using the space group symmetry (P1)",ttshown)
    end if

    ! Cartesian coordinates (crystals only)
    if (ifmt == isformat_w_aimsin .and. anycrystal) then
       ldum = iw_checkbox("Cartesian coordinates##savemultcartesian",w%savemult_cartesian)
       call iw_tooltip("Write Cartesian atomic coordinates (atom) instead of&
          & fractional coordinates (atom_frac)",ttshown)
    end if

    ! show the unit cell in the 3D model
    if (usecell) then
       ldum = iw_checkbox("Show unit cell##savemultdocell",w%savemult_docell)
       call iw_tooltip("Draw the unit cell edges in the written 3D model files",ttshown)
    end if

    ! the single warning that goes under the table. It is built before the
    ! table because the table is sized to fill whatever the block below it
    ! does not use, so its height has to be known first
    warnmsg = ""
    dooverwrite = .false.
    if (w%sm%nloading > 0) then
       warnmsg = string(w%sm%nloading) // " more systems are still loading"
    elseif (w%sm%nhidden > 0) then
       warnmsg = string(w%sm%nhidden) // " systems in collapsed groups are skipped"
    elseif (w%sm%ndup > 0) then
       warnmsg = string(w%sm%ndup) // " file w%sm%names are repeated: add %i to the pattern"
    elseif (w%sm%nexist > 0) then
       warnmsg = string(w%sm%nexist) // " of these files exist already"
    end if
    ! the choice is asked for whenever files exist, whichever warning won
    ! the chain above: Save gates on it and would otherwise be unanswerable
    dooverwrite = (w%sm%nexist > 0)
    ! whether the structures can be written, and if not, the
    ! reason. Costs one frame of lag on the checkbox.
    okvalid = (w%sm%nnames > 0) .and. (w%sm%ndup == 0) .and. dirok .and. .not.doquit .and.&
       .not.are_threads_running()
    if (okvalid) okvalid = (w%sm%nexist == 0) .or. (w%savemult_exist /= 0)
    savemsg = ""
    if (.not.okvalid .and. .not.doquit) then
       if (are_threads_running()) then
          savemsg = "Cannot write: the systems are still being loaded"
       elseif (.not.dirok) then
          savemsg = "Cannot write: choose a directory first"
       elseif (w%sm%nnames == 0) then
          savemsg = "Cannot write: there are no systems to write"
       elseif (w%sm%ndup > 0) then
          savemsg = "Cannot write: the file w%sm%names are not unique"
       else
          savemsg = "Cannot write: say what to do with the files that exist"
       end if
    end if

    nbot = 1 ! the button row
    if (len_trim(warnmsg) > 0) nbot = nbot + 1
    if (dooverwrite) nbot = nbot + 1
    if (len_trim(w%errmsg) > 0 .or. len_trim(w%okmsg) > 0) nbot = nbot + 1
    if (len_trim(savemsg) > 0) nbot = nbot + 1

    ! the list of files that will be written
    call iw_text("Files (" // string(w%sm%nnames) // ")",highlight=.true.)
    str1 = "##savemulttable" // c_null_char
    flags = ImGuiTableFlags_Borders
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_RowBg)
    sz%x = 0._c_float
    sz%y = -iw_calcheight(nbot,0,.false.)
    if (igBeginTable(c_loc(str1),ic_savemult_NUMCOLUMNS,flags,sz,0._c_float)) then
       call iw_table_column("Id",id=ic_savemult_id,flags=ImGuiTableColumnFlags_WidthFixed)
       flags = ior(ImGuiTableColumnFlags_WidthFixed,ImGuiTableColumnFlags_NoHeaderLabel)
       call iw_table_column("(kind)##0savemultkind",id=ic_savemult_type,flags=flags,&
          width=max(4._c_float,fontsize%y + 2._c_float))
       call iw_table_column("System",id=ic_savemult_name,flags=ImGuiTableColumnFlags_WidthStretch)
       call iw_table_column("File",id=ic_savemult_file,flags=ImGuiTableColumnFlags_WidthStretch)
       call igTableSetupScrollFreeze(0,1) ! the header row is always visible
       call igTableHeadersRow()

       if (w%sm%nnames > 0) then
          ! the rows go through a clipper: the list is as long as the tree,
          ! and only the visible rows need to be emitted
          clipper = ImGuiListClipper_ImGuiListClipper()
          call ImGuiListClipper_Begin(clipper,w%sm%nnames,-1._c_float)
          do while (ImGuiListClipper_Step(clipper))
             call c_f_pointer(clipper,clipper_f)
             do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                isys = w%sm%nameidx(i)
                call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                if (igTableSetColumnIndex(ic_savemult_id)) &
                   call iw_text(string(isys))
                if (igTableSetColumnIndex(ic_savemult_type)) &
                   call draw_kind_icon(ifmt,isys)
                if (igTableSetColumnIndex(ic_savemult_name)) &
                   call iw_text(trim(sysc(isys)%seed%name))
                if (igTableSetColumnIndex(ic_savemult_file)) &
                   call draw_file_cell(i,isys,ifmt)
             end do
          end do
          call ImGuiListClipper_End(clipper)
          call ImGuiListClipper_destroy(clipper)
       else
          call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
          if (igTableSetColumnIndex(ic_savemult_name)) &
             call iw_text("No systems to write")
       end if
       call igEndTable()
    end if

    ! the warning, and the overwrite acknowledgement it may need
    if (len_trim(warnmsg) > 0) then
       call iw_text(warnmsg,danger=.true.,wrap=.true.)
    end if
    if (dooverwrite) then
       call iw_text("Write them?",highlight=.true.,alignframe=.true.)
       ldum = iw_radiobutton("Overwrite##savemultexist",int=w%savemult_exist,intval=1_c_int,&
          sameline=.true.)
       call iw_tooltip("Write over the files that already exist in this directory",ttshown)
       ldum = iw_radiobutton("Skip##savemultexist",int=w%savemult_exist,intval=2_c_int,&
          sameline=.true.)
       call iw_tooltip("Leave the files that already exist alone and write only the others",ttshown)
    end if

    ! maybe the error or result message
    if (len_trim(w%errmsg) > 0) then
       call iw_text(w%errmsg,danger=.true.,wrap=.true.)
    elseif (len_trim(w%okmsg) > 0) then
       call iw_text(w%okmsg,highlight=.true.,wrap=.true.)
    end if

    ! why the Save button is disabled, on the line just above it
    if (len_trim(savemsg) > 0) call iw_text(savemsg,danger=.true.,wrap=.true.)

    ! right-align and bottom-align for the rest of the contents: the two
    ! buttons of the row below carry 9 characters between them
    call iw_setpos_bottomright(9,2)

    ! final buttons: save
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) .and. okvalid
    ok = ok .or. iw_button("Save",disabled=.not.okvalid)
    if (ok) then
       w%errmsg = ""
       w%okmsg = ""
       ! re-check the directory warning next frame: the write may have
       ! materialized the directory (or the user created it meanwhile)
       w%savemult_lastcheck = ""
       nfail = 0
       nwritten = 0
       nskip = 0
       nleft = 0
       firsterr = ""
       do i = 1, w%sm%nnames
          isys = w%sm%nameidx(i)
          ! a system may have been closed since the list was built
          if (.not.ok_system(isys,sys_init)) then
             nskip = nskip + 1
             cycle
          end if
          ! leave the files that are there already, if that was the choice
          if (w%savemult_exist == 2 .and. w%sm%nameexists(i)) then
             nleft = nleft + 1
             cycle
          end if
          file = savedir(w) // dirsep // trim(w%sm%names(i))
          ! FHIaims writes the k-point grid to a companion control file,
          ! which is meaningless for a molecule: decide per system
          thisrk = userk .and. .not.writes_molecule(ifmt,isys)
          if (thisrk) then
             call sys(isys)%c%write_any_file(file,w%errmsg,iwformat=ifmt,&
                rklength=real(w%savemult_rk,8),nosym=w%savemult_nosym,&
                cartesian=w%savemult_cartesian,docell=w%savemult_docell)
          else
             call sys(isys)%c%write_any_file(file,w%errmsg,iwformat=ifmt,&
                nosym=w%savemult_nosym,cartesian=w%savemult_cartesian,&
                docell=w%savemult_docell)
          end if
          if (len_trim(w%errmsg) > 0) then
             nfail = nfail + 1
             if (nfail == 1) firsterr = trim(w%errmsg)
             w%errmsg = ""
          else
             nwritten = nwritten + 1
             write (uout,'("Saved structure file: ",A)') trim(file)
          end if
       end do

       ! remember the settings and report; a failure keeps the window open
       lastfmt = w%savemult_format
       lastpattern = w%savemult_pattern
       lastrk = w%savemult_rk
       lastnosym = w%savemult_nosym
       lastcartesian = w%savemult_cartesian
       lastdocell = w%savemult_docell
       if (nwritten > 0) okfile_lastdir = savedir(w)
       if (nfail == 0) then
          w%okmsg = string(nwritten) // " structure files written"
          if (nleft > 0) &
             w%okmsg = w%okmsg // ", " // string(nleft) // " already there"
          if (nskip > 0) &
             w%okmsg = w%okmsg // ", " // string(nskip) // " not loaded"
          write (uout,'("Saved ",A," structure files to: ",A)') string(nwritten), savedir(w)
       else
          file = string(nfail) // " of " // string(w%sm%nnames) // " files could not be written"
          if (len_trim(firsterr) > 0) file = file // ": " // trim(firsterr)
          w%errmsg = file
       end if
       ! the files on disk changed: refresh the existence check in place
       call update_existence(w)
       if (w%savemult_exist == 0) w%savemult_exist = 1
    end if

    ! final buttons: close
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  contains
    !> Draw the icon that says whether system isys comes out of format
    !> ifmt as a molecule or as a crystal, in the current table cell.
    subroutine draw_kind_icon(ifmt,isys)
      integer, intent(in) :: ifmt, isys

      integer(c_int) :: itex
      real(c_float) :: rgba(4)
      character(len=:), allocatable :: str2

      if (writes_molecule(ifmt,isys)) then
         itex = icon_tex(icon_vm_bond)
         rgba = rgba_fam_molgeom
         str2 = "Written as a molecule"
      else
         itex = icon_tex(icon_ui_cell)
         rgba = rgba_fam_crys
         str2 = "Written as a periodic crystal"
      end if
      if (itex == 0) return
      call igImage(int(itex,c_intptr_t),ImVec2(fontsize%y,fontsize%y),uv0,uv1,&
         ImVec4(rgba(1),rgba(2),rgba(3),rgba(4)),nobord)
      ! no ttshown: an image carries no ImGui ID, so HoveredIdTimer never
      ! runs on it and the delayed form of the tooltip would never fire
      call iw_tooltip(str2)

    end subroutine draw_kind_icon

    !> Draw the File cell of row i (system isys): the file name, or an
    !> input box when this is the row being edited. A double-click starts
    !> the edit, Enter or moving away commits it, and an empty name gives
    !> the row back to the pattern.
    subroutine draw_file_cell(i,isys,ifmt)
      integer, intent(in) :: i, isys, ifmt

      logical :: committed
      real(c_float) :: pos
      type(ImVec2) :: szero
      character(kind=c_char,len=:), allocatable, target :: strl

      szero = ImVec2(0._c_float,0._c_float)

      if (w%sm%editsys == isys) then
         call igPushItemWidth(-1._c_float)
         committed = iw_inputtext("##savemultfile",bufsize=mlen-1,texta=w%sm%editbuf,&
            grabfocus=w%sm%editfocus,notlive=.true.,flags=ImGuiInputTextFlags_AutoSelectAll)
         if (committed) call override_set(w%sm,isys,typed_root(w%sm%editbuf,ifmt))
         if (committed .or. igIsItemDeactivated()) w%sm%editsys = 0
         call igPopItemWidth()
         w%sm%editfocus = .false.
      else
         ! a zero-width selectable under the text
         pos = igGetCursorPosX()
         strl = "##savemultfile" // string(isys) // c_null_char
         if (igSelectable_Bool(c_loc(strl),.false._c_bool,&
            ImGuiSelectableFlags_AllowDoubleClick,szero)) then
            if (igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)) then
               w%sm%editsys = isys
               w%sm%editbuf = trim(w%sm%names(i))
               w%sm%editfocus = .true.
            end if
         end if
         call igSameLine(0._c_float,-1._c_float)
         call igSetCursorPosX(pos)
         call iw_text(trim(w%sm%names(i)),danger=w%sm%nameexists(i))
      end if

    end subroutine draw_file_cell

    !> Draw the help mark for the file name pattern: one tooltip line per
    !> substitution, laid out in two columns like the view-mode tooltip.
    subroutine draw_pattern_help()
      use gui_main, only: tooltip_enabled, tooltip_wrap_factor

      integer, parameter :: nsub = 6
      character(len=*), parameter :: subkey(nsub) = (/character(len=2) :: &
         "%n","%f","%i","%d","%%","*"/)
      character(len=*), parameter :: subtxt(nsub) = (/character(len=38) :: &
         "Name of the system                    ",&
         "Base name of the file                 ",&
         "Position in the tree list, zero-padded",&
         "ID of the system in the tree          ",&
         "A literal percent sign                ",&
         "Same as %i                            "/)

      integer :: i

      call iw_text("(?)",sameline=.true.)
      if (.not.tooltip_enabled) return
      if (.not.igIsItemHovered(ImGuiHoveredFlags_None)) return

      call igBeginTooltip()
      call iw_text("File Name Pattern",highlight=.true.)
      call igPushTextWrapPos(tooltip_wrap_factor * fontsize%x)
      call iw_text("The name of each file: the extension is added automatically. Double-click&
         & a name in the table to type it by hand.")
      call igPopTextWrapPos()
      call iw_text("")
      do i = 1, nsub
         call iw_text(string(trim(subkey(i)),length=len(subkey)+1),highlight=.true.)
         call iw_text(trim(subtxt(i)),sameline=.true.)
      end do
      call igEndTooltip()

    end subroutine draw_pattern_help

  end subroutine draw_save_multiple

  !xx! private procedures

  !> Rebuild the list of systems and target file w%sm%names for the current
  !> window settings and write format ifmt, and count the repeated w%sm%names
  !> and the files that exist already. All of this is recalculated only
  !> when the settings it depends on change; changed reports whether it
  !> was recalculated in this pass.
  subroutine rebuild_names(w,ifmt,changed)
    use systems, only: sysc, sys_init, sys_empty, ok_system
    use tools_io, only: string
    class(window), intent(inout) :: w
    integer, intent(in) :: ifmt
    logical, intent(out) :: changed

    character(len=:), allocatable :: sig, root, ext
    integer, allocatable :: idx(:)
    integer :: i, j, k, n, npad
    integer*8 :: ihash

    ! overrides and the edited cell are keyed by system slot, and a closed
    ! slot is handed to the next structure loaded: drop them before that
    do i = w%sm%novr, 1, -1
       if (sysc(w%sm%ovrsys(i))%status == sys_empty) call override_set(w%sm,w%sm%ovrsys(i),"")
    end do
    if (w%sm%editsys > 0) then
       if (sysc(w%sm%editsys)%status == sys_empty) w%sm%editsys = 0
    end if

    ! the list of systems: a system that is not loaded yet cannot be
    ! written. The members of a collapsed group stay unloaded until
    ! the user expands the group in the tree
    call tree_system_list(idx,(w%savemult_scope == 0))
    n = 0
    w%sm%nloading = 0
    w%sm%nhidden = 0
    do i = 1, size(idx,1)
       if (.not.ok_system(idx(i),sys_init)) then
          if (sysc(idx(i))%hidden) then
             w%sm%nhidden = w%sm%nhidden + 1
          else
             w%sm%nloading = w%sm%nloading + 1
          end if
          cycle
       end if
       n = n + 1
       idx(n) = idx(i)
    end do

    ! signature of everything the file w%sm%names depend on: the settings, and
    ! a hash of the list of systems and their w%sm%names (hashing avoids the
    ! string allocations that building the list itself would cost)
    ihash = n
    do i = 1, n
       ihash = mod(ihash * 31 + idx(i),sighashmod)
       do j = 1, len_trim(sysc(idx(i))%seed%name)
          ihash = mod(ihash * 31 + iachar(sysc(idx(i))%seed%name(j:j)),sighashmod)
       end do
    end do
    ihash = mod(ihash * 31 + w%sm%novr,sighashmod)
    do i = 1, w%sm%novr
       ihash = mod(ihash * 31 + w%sm%ovrsys(i),sighashmod)
       do j = 1, len_trim(w%sm%ovrname(i))
          ihash = mod(ihash * 31 + iachar(w%sm%ovrname(i)(j:j)),sighashmod)
       end do
    end do
    sig = string(ifmt) // "|" // trim(w%okfile) // "|" //&
       trim(w%savemult_pattern) // "|" // string(ihash)
    changed = .true.
    if (allocated(w%sm%lastsig)) then
       if (sig == w%sm%lastsig) then
          changed = .false.
          return
       end if
    end if
    w%sm%lastsig = sig

    ! build the file w%sm%names
    w%sm%nnames = n
    if (allocated(w%sm%names)) deallocate(w%sm%names)
    if (allocated(w%sm%nameidx)) deallocate(w%sm%nameidx)
    if (allocated(w%sm%nameexists)) deallocate(w%sm%nameexists)
    allocate(w%sm%names(n),w%sm%nameidx(n),w%sm%nameexists(n))
    w%sm%nameexists = .false.
    w%sm%nameidx(1:n) = idx(1:n)
    npad = len_trim(string(max(n,1)))
    ext = "." // trim(fmtext(ifmt))
    do i = 1, n
       k = override_find(w%sm,idx(i))
       if (k > 0) then
          root = trim(w%sm%ovrname(k))
       else
          root = expand_pattern(w%savemult_pattern,idx(i),i,npad)
       end if
       if (len_trim(root) == 0) root = "structure" // string(i,npad,pad0=.true.)
       ! do not repeat the extension if the name already carries it
       w%sm%names(i) = strip_ext(root,ext) // ext
    end do

    ! repeated w%sm%names: only possible if the pattern carries no index. The
    ! count is of the files involved in a collision, not of the repeats
    w%sm%ndup = 0
    if (w%sm%novr > 0 .or. .not.pattern_hasindex(w%savemult_pattern)) then
       do i = 1, n
          do j = 1, n
             if (j /= i .and. w%sm%names(i) == w%sm%names(j)) then
                w%sm%ndup = w%sm%ndup + 1
                exit
             end if
          end do
       end do
    end if

    ! files that exist already
    call update_existence(w)

  end subroutine rebuild_names

  !> Whether system isys is written as a molecule when the write format
  !> is ifmt. Some formats are inherently molecular, some always write a
  !> cell, and the rest follow the system.
  function writes_molecule(ifmt,isys)
    use systems, only: sys
    integer, intent(in) :: ifmt, isys
    logical :: writes_molecule

    if (ifmt == isformat_w_xyz .or. ifmt == isformat_w_gjf .or. ifmt == isformat_w_cml .or.&
       ifmt == isformat_w_obj .or. ifmt == isformat_w_ply .or. ifmt == isformat_w_off) then
       writes_molecule = .true.
    elseif (ifmt == isformat_w_aimsin .or. ifmt == isformat_w_pdb .or.&
       ifmt == isformat_w_pyscf .or. ifmt == isformat_w_dftbp_gen .or.&
       ifmt == isformat_w_dftbp_hsd) then
       writes_molecule = sys(isys)%c%ismolecule
    else
       writes_molecule = .false.
    end if

  end function writes_molecule

  !> The output directory of window w, without a trailing separator (so
  !> that a user-typed "dir/" does not give "dir//name").
  function savedir(w)
    use param, only: dirsep
    class(window), intent(in) :: w
    character(len=:), allocatable :: savedir

    savedir = trim(w%okfile)
    do while (len(savedir) > 1)
       if (savedir(len(savedir):) /= dirsep) exit
       savedir = savedir(1:len(savedir)-1)
    end do

  end function savedir

  !> Check which of the target file names exist in the output directory
  !> of window w, and count them.
  subroutine update_existence(w)
    use param, only: dirsep
    class(window), intent(inout) :: w

    character(len=:), allocatable :: dir
    integer :: i

    dir = savedir(w) // dirsep
    do i = 1, w%sm%nnames
       inquire(file=dir // trim(w%sm%names(i)),exist=w%sm%nameexists(i))
    end do
    w%sm%nexist = count(w%sm%nameexists(1:w%sm%nnames))

  end subroutine update_existence

  !> Index of the hand-typed file name of system isys in the override
  !> list, or 0 if it does not have one.
  function override_find(sm,isys)
    type(savemult_state), intent(in) :: sm
    integer, intent(in) :: isys
    integer :: override_find

    integer :: i

    override_find = 0
    do i = 1, sm%novr
       if (sm%ovrsys(i) == isys) then
          override_find = i
          return
       end if
    end do

  end function override_find

  !> Set the hand-typed file name root of system isys. An empty root
  !> drops the override, so that the pattern sm%names the row again.
  subroutine override_set(sm,isys,root)
    use types, only: realloc
    type(savemult_state), intent(inout) :: sm
    integer, intent(in) :: isys
    character(len=*), intent(in) :: root

    integer :: i, k

    k = override_find(sm,isys)

    ! an empty name: forget the override, closing the gap it leaves
    if (len_trim(root) == 0) then
       if (k == 0) return
       do i = k, sm%novr-1
          sm%ovrsys(i) = sm%ovrsys(i+1)
          sm%ovrname(i) = sm%ovrname(i+1)
       end do
       sm%novr = sm%novr - 1
       return
    end if

    ! replace the one that is there, or make room for a new one
    if (k == 0) then
       sm%novr = sm%novr + 1
       call realloc(sm%ovrsys,sm%novr)
       call realloc(sm%ovrname,sm%novr)
       k = sm%novr
       sm%ovrsys(k) = isys
    end if
    sm%ovrname(k) = root

  end subroutine override_set

  !> A file name typed by hand, as a root: with the characters that have
  !> no business in a file name replaced, and with the extension of format
  !> ifmt removed if it is there, so that changing the format re-extends
  !> the name instead of stacking a second extension on it.
  function typed_root(str,ifmt)
    character(len=*), intent(in) :: str
    integer, intent(in) :: ifmt
    character(len=:), allocatable :: typed_root

    typed_root = trim(adjustl(str))
    call sanitize(typed_root)
    typed_root = strip_ext(typed_root,"." // trim(fmtext(ifmt)))

  end function typed_root

  !> Whether the file name pattern carries a token that makes every
  !> expansion unique (%i, %d or *).
  function pattern_hasindex(pattern)
    character(len=*), intent(in) :: pattern
    logical :: pattern_hasindex

    integer :: k

    pattern_hasindex = .true.
    k = 1
    do while (k <= len_trim(pattern))
       if (pattern(k:k) == "*") return
       if (pattern(k:k) == "%" .and. k < len_trim(pattern)) then
          if (pattern(k+1:k+1) == "i" .or. pattern(k+1:k+1) == "d") return
          k = k + 2
       else
          k = k + 1
       end if
    end do
    pattern_hasindex = .false.

  end function pattern_hasindex

  !> Expand the file name pattern for system isys, which is number i in
  !> the list of systems being written, with the index zero-padded to
  !> npad digits.
  function expand_pattern(pattern,isys,i,npad) result(root)
    use systems, only: sysc
    use utils, only: file_name_root
    use tools_io, only: string
    use param, only: dirsep
    character(len=*), intent(in) :: pattern
    integer, intent(in) :: isys, i, npad
    character(len=:), allocatable :: root

    character(len=:), allocatable :: name, fbase, aux
    integer :: k

    ! %n: the name of the system, without the extension
    name = trim(sysc(isys)%seed%name)
    k = index(name,"|")
    if (k > 0) then
       name = file_name_root(name(1:k-1)) // "_" // name(k+1:)
    else
       name = file_name_root(name)
    end if
    call sanitize(name)
    if (len_trim(name) == 0) name = "structure"

    ! %f: the base name of the file the system was read from
    fbase = trim(sysc(isys)%seed%file)
    k = index(fbase,dirsep,back=.true.)
    if (k > 0) fbase = fbase(k+1:)
    fbase = file_name_root(fbase)
    call sanitize(fbase)
    if (len_trim(fbase) == 0) fbase = name

    ! substitute the tokens, left to right
    root = ""
    aux = trim(pattern)
    do while (len(aux) > 0)
       if (aux(1:1) == "*") then
          root = root // string(i,npad,pad0=.true.)
          aux = aux(2:)
       elseif (len(aux) >= 2 .and. aux(1:1) == "%") then
          if (aux(2:2) == "n") then
             root = root // trim(name)
          elseif (aux(2:2) == "f") then
             root = root // trim(fbase)
          elseif (aux(2:2) == "i") then
             root = root // string(i,npad,pad0=.true.)
          elseif (aux(2:2) == "d") then
             root = root // string(isys)
          elseif (aux(2:2) == "%") then
             root = root // "%"
          else
             root = root // aux(1:2)
          end if
          aux = aux(3:)
       else
          root = root // aux(1:1)
          aux = aux(2:)
       end if
    end do
    ! the literal part of the pattern is sanitized too: a separator in it
    ! would write outside the chosen directory
    call sanitize(root)

  end function expand_pattern

  !> Replace anything that has no business in a file name by an underscore.
  subroutine sanitize(str)
    character(len=:), allocatable, intent(inout) :: str

    integer :: j

    do j = 1, len(str)
       if (index(badchars,str(j:j)) > 0) str(j:j) = "_"
    end do

  end subroutine sanitize

  !> str without a trailing ext, if it carries one.
  function strip_ext(str,ext)
    character(len=*), intent(in) :: str, ext
    character(len=:), allocatable :: strip_ext

    strip_ext = str
    if (len(strip_ext) <= len(ext)) return
    if (strip_ext(len(strip_ext)-len(ext)+1:) == ext) &
       strip_ext = strip_ext(1:len(strip_ext)-len(ext))

  end function strip_ext

end submodule save_multiple
