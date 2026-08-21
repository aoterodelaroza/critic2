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

  ! maximum number of file names listed in the preview
  integer, parameter :: maxpreview = 100

  ! characters replaced by an underscore in the name of a system (%n)
  character(len=*), parameter :: badchars = " |/:*?<>" //&
     achar(92) // achar(34) // achar(39) // char(9)

  ! modulus for the hash of the system list (keeps it well inside int8)
  integer*8, parameter :: sighashmod = 2147483647_8

  ! The target file names for the current settings, recalculated only
  ! when the settings they depend on change (lastsig); the signature
  ! spares the per-frame rebuild of the names and, above all, the
  ! inquire() on every one of them. This state is module-level because
  ! only one save-multiple window can exist at a time: all three entry
  ! points parent it to the tree, and stack_create_window(orraise=-1)
  ! raises the open one instead of making a second.
  character(len=:), allocatable :: lastsig
  character(len=mlen), allocatable :: names(:)
  integer, allocatable :: nameidx(:)
  logical, allocatable :: nameexists(:)
  integer :: nnames = 0 ! number of files that will be written
  integer :: nloading = 0 ! systems in the list that are still being loaded
  integer :: nhidden = 0 ! systems in the list held back by a collapsed group
  integer :: nexist = 0 ! how many of them exist on disk already
  integer :: ndup = 0 ! how many of them are repeated in the list

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
       iw_setpos_bottomright, iw_inputtext, get_current_working_dir
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG
    use tools_io, only: string, uout
    use param, only: dirsep
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str1
    character(len=:), allocatable :: file, firsterr
    integer :: i, isys, ifmt, iaux, nfail, nwritten, nskip
    logical :: doquit, ok, okvalid, ismol, userk, usecell, changed, dirok, direxists
    logical :: anycrystal, thisrk
    logical(c_bool) :: ldum
    type(ImVec2) :: sz

    logical, save :: ttshown = .false. ! tooltip flag

    ! build the format combo options on first use; unlike the save-as
    ! window there is no auto-detect entry, as the pattern carries no
    ! extension, so skip the first entry of the shared combo string
    call build_write_format_combo()

    ! this window works on the tree, so it closes when the tree is gone
    doquit = .true.
    if (w%idparent >= 1 .and. w%idparent <= nwin) doquit = .not.win(w%idparent)%isinit

    ! initialize state; the cache is stale on reopen (the files on disk
    ! may have changed in the meantime), so drop it
    if (w%firstpass) then
       if (allocated(lastsig)) deallocate(lastsig)
       w%errmsg = ""
       w%okmsg = ""
       w%savemult_format = lastfmt
       w%savemult_rk = lastrk
       w%savemult_nosym = lastnosym
       w%savemult_cartesian = lastcartesian
       w%savemult_docell = lastdocell
       w%savemult_overwrite = .false.
       if (allocated(lastpattern)) then
          w%savemult_pattern = lastpattern
       else
          w%savemult_pattern = "%n"
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
       & tree with control-click and shift-click, or with the tree options button",ttshown)

    ! the output directory
    call iw_text("Directory",highlight=.true.)
    ldum = iw_inputtext("##savemultdir",bufsize=1023,texta=w%okfile,width=36)
    call iw_tooltip("Directory where the structure files will be written",ttshown)
    if (iw_button("Browse...",sameline=.true.)) &
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_selectdir,idparent=w%id,orraise=-1)
    call iw_tooltip("Choose the directory with a file browser",ttshown)
    dirok = (len_trim(w%okfile) > 0)
    if (dirok) then
       ! a missing directory is only a warning: inquire() on a directory is
       ! not portable, and the write itself reports the real error
       inquire(file=savedir(w),exist=direxists)
       if (.not.direxists) &
          call iw_text("This directory does not exist",danger=.true.)
    else
       call iw_text("Choose a directory to write the files to",danger=.true.)
    end if

    ! format combo
    call iw_text("Format",highlight=.true.)
    call igPushItemWidth(iw_calcwidth(45,1))
    call iw_combo_simple("##savemultformat",combostr(icombo_fmt1:),w%savemult_format,startsatone=.true.)
    call igPopItemWidth()
    call iw_tooltip("File format for the written structures",ttshown)
    ifmt = fmtperm(max(min(w%savemult_format,isformat_w_max),1))

    ! the file name pattern
    call iw_text("File Names",highlight=.true.)
    ldum = iw_inputtext("##savemultpattern",bufsize=1023,texta=w%savemult_pattern,width=36)
    call iw_tooltip("Pattern for the file names, without extension. %n is replaced by the&
       & name of the system, %i by its position in the list (zero-padded), and %s by its&
       & ID in the tree. A * is the same as %i (as in the ROOT option to WRITE BULK).&
       & The extension for the chosen format is added automatically.",ttshown)

    ! whether the write will produce molecules or crystals: some formats
    ! are inherently molecular, some follow the system, and the rest
    ! always write a cell
    ismol = (ifmt == isformat_w_xyz .or. ifmt == isformat_w_gjf .or. ifmt == isformat_w_cml .or.&
       ifmt == isformat_w_obj .or. ifmt == isformat_w_ply .or. ifmt == isformat_w_off)

    ! rebuild the list of target file names, if the settings changed
    call rebuild_names(w,ifmt,changed)
    if (changed) then
       w%okmsg = ""
       w%errmsg = ""
       ! the user asked for these systems, so get the ones that are only
       ! waiting for a worker moving
       if (nloading > 0) call launch_initialization_thread()
    end if

    ! the list can mix molecules and crystals, so the options that only
    ! apply to a crystal are offered when at least one of the systems is
    ! one; they are then used only for those systems
    anycrystal = .false.
    do i = 1, nnames
       if (.not.sys(nameidx(i))%c%ismolecule) then
          anycrystal = .true.
          exit
       end if
    end do
    if (.not.ismol) then
       if (ifmt == isformat_w_aimsin .or. ifmt == isformat_w_pdb .or. ifmt == isformat_w_pyscf .or.&
          ifmt == isformat_w_dftbp_gen .or. ifmt == isformat_w_dftbp_hsd) &
          ismol = .not.anycrystal
    end if

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
       if (ifmt == isformat_w_aimsin) &
          call iw_text("K-point grid written to a companion _control file")
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
       call iw_tooltip("Draw the unit cell edges in the written 3D model files&
          & (same as the CELL option to WRITE)",ttshown)
    end if

    ! the list of files that will be written
    call iw_text("Files (" // string(nnames) // ")",highlight=.true.)
    str1 = "##savemultpreview" // c_null_char
    sz%x = 0._c_float
    sz%y = iw_calcheight(6,0,.false.)
    if (igBeginChild_Str(c_loc(str1),sz,.true._c_bool,ImGuiWindowFlags_HorizontalScrollbar)) then
       if (nnames == 0) then
          call iw_text("No systems to write")
       else
          do i = 1, min(nnames,maxpreview)
             isys = nameidx(i)
             call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name))
             call iw_text(" -> " // trim(names(i)),sameline=.true.,danger=nameexists(i))
          end do
          if (nnames > maxpreview) &
             call iw_text("... and " // string(nnames-maxpreview) // " more")
       end if
    end if
    call igEndChild()

    ! warn about repeated names and files that would be overwritten
    if (nloading > 0) then
       call iw_text(string(nloading) // " more systems are still loading",danger=.true.)
    elseif (nhidden > 0) then
       call iw_text(string(nhidden) // " systems in collapsed groups are skipped",danger=.true.)
       call iw_tooltip("The systems inside a collapsed group are not loaded until the group&
          & is expanded in the tree, so they cannot be written",ttshown)
    elseif (ndup > 0) then
       call iw_text(string(ndup) // " file names are repeated: add %i to the pattern",danger=.true.)
    elseif (nexist > 0) then
       call iw_text(string(nexist) // " of these files exist and will be overwritten",danger=.true.)
       ldum = iw_checkbox("Overwrite existing files##savemultoverwrite",w%savemult_overwrite)
       call iw_tooltip("Allow writing over the files that already exist in this directory",ttshown)
    end if

    ! maybe the error or result message
    if (len_trim(w%errmsg) > 0) then
       call iw_text(w%errmsg,danger=.true.)
    elseif (len_trim(w%okmsg) > 0) then
       call iw_text(w%okmsg,highlight=.true.)
    end if

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(8,2)

    ! final buttons: save
    okvalid = (nnames > 0) .and. (ndup == 0) .and. dirok .and. .not.doquit .and.&
       .not.are_threads_running()
    if (okvalid) okvalid = (nexist == 0) .or. w%savemult_overwrite
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) .and. okvalid
    ok = ok .or. iw_button("Save",disabled=.not.okvalid)
    if (ok) then
       w%errmsg = ""
       w%okmsg = ""
       nfail = 0
       nwritten = 0
       nskip = 0
       firsterr = ""
       do i = 1, nnames
          isys = nameidx(i)
          ! a system may have been closed since the list was built
          if (.not.ok_system(isys,sys_init)) then
             nskip = nskip + 1
             cycle
          end if
          file = savedir(w) // dirsep // trim(names(i))
          ! FHIaims writes the k-point grid to a companion control file,
          ! which is meaningless for a molecule: decide per system
          thisrk = (ifmt == isformat_w_qein .or. ifmt == isformat_w_castepcell .or.&
             (ifmt == isformat_w_aimsin .and. .not.sys(isys)%c%ismolecule))
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
          if (nskip > 0) &
             w%okmsg = w%okmsg // " (" // string(nskip) // " systems skipped)"
          write (uout,'("Saved ",A," structure files to: ",A)') string(nwritten), savedir(w)
       else
          file = string(nfail) // " of " // string(nnames) // " files could not be written"
          if (len_trim(firsterr) > 0) file = file // ": " // trim(firsterr)
          w%errmsg = file
       end if
       ! the files on disk changed: refresh the existence check in place
       ! (invalidating the whole cache would wipe the message above). The
       ! files just written are known to the user, so the overwrite guard
       ! starts out acknowledged instead of disabling the button
       call update_existence(w)
       w%savemult_overwrite = .true.
    end if

    ! final buttons: close
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_save_multiple

  !xx! private procedures

  !> Rebuild the list of systems and target file names for the current
  !> window settings and write format ifmt, and count the repeated names
  !> and the files that exist already. All of this is recalculated only
  !> when the settings it depends on change; changed reports whether it
  !> was recalculated in this pass.
  subroutine rebuild_names(w,ifmt,changed)
    use systems, only: sysc, sys_init, ok_system
    use tools_io, only: string
    class(window), intent(in) :: w
    integer, intent(in) :: ifmt
    logical, intent(out) :: changed

    character(len=:), allocatable :: sig, root, ext
    integer, allocatable :: idx(:)
    integer :: i, j, n, npad
    integer*8 :: ihash

    ! the list of systems: a system that is not loaded yet cannot be
    ! written. The members of a collapsed group are a case of their own:
    ! the initialization threads skip hidden systems, so they stay
    ! unloaded until the user expands the group in the tree
    call tree_system_list(idx,(w%savemult_scope == 0))
    n = 0
    nloading = 0
    nhidden = 0
    do i = 1, size(idx,1)
       if (.not.ok_system(idx(i),sys_init)) then
          if (sysc(idx(i))%hidden) then
             nhidden = nhidden + 1
          else
             nloading = nloading + 1
          end if
          cycle
       end if
       n = n + 1
       idx(n) = idx(i)
    end do

    ! signature of everything the file names depend on: the settings, and
    ! a hash of the list of systems and their names (hashing avoids the
    ! string allocations that building the list itself would cost)
    ihash = n
    do i = 1, n
       ihash = mod(ihash * 31 + idx(i),sighashmod)
       do j = 1, len_trim(sysc(idx(i))%seed%name)
          ihash = mod(ihash * 31 + iachar(sysc(idx(i))%seed%name(j:j)),sighashmod)
       end do
    end do
    sig = string(ifmt) // "|" // trim(w%okfile) // "|" //&
       trim(w%savemult_pattern) // "|" // string(ihash)
    changed = .true.
    if (allocated(lastsig)) then
       if (sig == lastsig) then
          changed = .false.
          return
       end if
    end if
    lastsig = sig

    ! build the file names
    nnames = n
    if (allocated(names)) deallocate(names)
    if (allocated(nameidx)) deallocate(nameidx)
    if (allocated(nameexists)) deallocate(nameexists)
    allocate(names(n),nameidx(n),nameexists(n))
    nameexists = .false.
    nameidx(1:n) = idx(1:n)
    npad = len_trim(string(max(n,1)))
    ext = "." // trim(fmtext(ifmt))
    do i = 1, n
       root = expand_pattern(w%savemult_pattern,idx(i),i,npad)
       if (len_trim(root) == 0) root = "structure" // string(i,npad,pad0=.true.)
       ! do not repeat the extension if the pattern already carries it
       if (len(root) > len(ext)) then
          if (root(len(root)-len(ext)+1:) == ext) then
             names(i) = root
             cycle
          end if
       end if
       names(i) = root // ext
    end do

    ! repeated names: only possible if the pattern carries no index. The
    ! count is of the files involved in a collision, not of the repeats
    ndup = 0
    if (.not.pattern_hasindex(w%savemult_pattern)) then
       do i = 1, n
          do j = 1, n
             if (j /= i .and. names(i) == names(j)) then
                ndup = ndup + 1
                exit
             end if
          end do
       end do
    end if

    ! files that exist already
    call update_existence(w)

  end subroutine rebuild_names

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
    class(window), intent(in) :: w

    integer :: i

    do i = 1, nnames
       inquire(file=savedir(w) // dirsep // trim(names(i)),exist=nameexists(i))
    end do
    nexist = count(nameexists(1:nnames))

  end subroutine update_existence

  !> Whether the file name pattern carries a token that makes every
  !> expansion unique (%i, %s or *). The pattern is walked the way
  !> expand_pattern does it, so that a literal %%i does not count.
  function pattern_hasindex(pattern)
    character(len=*), intent(in) :: pattern
    logical :: pattern_hasindex

    integer :: k

    pattern_hasindex = .true.
    k = 1
    do while (k <= len_trim(pattern))
       if (pattern(k:k) == "*") return
       if (pattern(k:k) == "%" .and. k < len_trim(pattern)) then
          if (pattern(k+1:k+1) == "i" .or. pattern(k+1:k+1) == "s") return
          k = k + 2
       else
          k = k + 1
       end if
    end do
    pattern_hasindex = .false.

  end function pattern_hasindex

  !> Expand the file name pattern for system isys, which is number i in
  !> the list of systems being written, with the index zero-padded to
  !> npad digits. %n is the name of the system, %i the index in the list,
  !> %s the ID of the system in the tree, and * the same as %i.
  function expand_pattern(pattern,isys,i,npad) result(root)
    use systems, only: sysc
    use utils, only: file_name_root
    use tools_io, only: string
    character(len=*), intent(in) :: pattern
    integer, intent(in) :: isys, i, npad
    character(len=:), allocatable :: root

    character(len=:), allocatable :: name, aux
    integer :: k

    ! the name of the system, without the extension and with anything
    ! that does not belong in a file name replaced by an underscore. The
    ! whole tree name is used, because it is already the shortest string
    ! that tells the systems apart: the directory prefix is there only
    ! when it was needed, and a system read from one of several data
    ! blocks in a file is named "file.ext|block", where the block names
    ! alone are often not unique (many CIFs use data_global or data_I)
    name = trim(sysc(isys)%seed%name)
    k = index(name,"|")
    if (k > 0) then
       name = file_name_root(name(1:k-1)) // "_" // name(k+1:)
    else
       name = file_name_root(name)
    end if
    do k = 1, len(name)
       if (index(badchars,name(k:k)) > 0) name(k:k) = "_"
    end do
    if (len_trim(name) == 0) name = "structure"

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
          elseif (aux(2:2) == "i") then
             root = root // string(i,npad,pad0=.true.)
          elseif (aux(2:2) == "s") then
             root = root // string(isys)
          else
             root = root // aux(1:2)
          end if
          aux = aux(3:)
       else
          root = root // aux(1:1)
          aux = aux(2:)
       end if
    end do

  end function expand_pattern

end submodule save_multiple
