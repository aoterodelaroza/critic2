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

! Windows, new structure and new structure from library dialogs.
submodule (windows) new_struct
  use interfaces_cimgui
  implicit none

  !xx! private procedures
  ! subroutine get_library_structure_list(libfile,nst,st,ismol)
  ! subroutine draw_spg_table(ispg)

contains

  !> Draw the contents of the new structure window.
  module subroutine draw_new_struct(w)
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use gui_main, only: g, add_systems_from_seeds,&
       launch_initialization_thread, system_shorten_names
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text, iw_calcheight,&
       iw_calcwidth, buffer_to_string_array, iw_radiobutton, iw_combo_simple, iw_checkbox
    use crystalseedmod, only: crystalseed, realloc_crystalseed
    use global, only: rborder_def
    use tools_io, only: string, fopen_scratch, fclose, stripchar, deblank
    use types, only: vstring
    use param, only: newline, bohrtoa
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, str2, stropt
    logical(c_bool) :: ldum, doquit
    logical :: ok
    type(ImVec2) :: szero, sz, szavail
    integer :: i, idx, lu
    type(crystalseed), allocatable :: seed_(:)
    integer(c_int) :: iunitat

    ! window state
    logical, save :: ttshown = .false. ! tooltip flags
    integer(c_size_t), parameter :: maxatposbuf = 100000 !! atomic positions buffer
    character(kind=c_char,len=:), allocatable, target, save :: atposbuf
    integer(c_size_t), parameter :: maxsymopbuf = 10000 !! symmetry operations buffer
    character(kind=c_char,len=:), allocatable, target, save :: symopbuf
    integer(c_size_t), parameter :: maxlatvecbuf = 5000 !! lattice vectors buffer
    character(kind=c_char,len=:), allocatable, target, save :: latvecbuf
    integer(c_size_t), parameter :: maxnamebuf = 1024 !! name buffer
    character(len=:,kind=c_char), allocatable, target, save :: namebuf
    logical, save :: ismolecule = .false. ! whether the system is a molecule or a crystal
    integer, save :: iunitat_c = 2 ! crystal atpos units (0 = bohr, 1 = angstrom, 2 = fractional)
    integer, save :: iunitat_m = 2 ! mol atpos units (0 = bohr, 1 = angstrom)
    integer, save :: symopt = 1 ! symmetry option (1 = detect, 2 = spg, 3 = manual)
    integer, save :: ispg_selected = 1 ! selected space group
    integer, save :: cellopt = 1 ! lattice option (1 = parameters, 2 = lattice-vectors)
    integer, save :: iunitcel = 0 ! units for lattice (0 = bohr, 1 = angstrom)
    real(c_float), save :: scale = 1._c_float ! scale factor for lattice vectors
    real(c_float), save :: aa(3) = 10._c_float ! cell lengths
    real(c_float), save :: bb(3) = 90._c_float ! cell angles
    real(c_float), save :: rborder = rborder_def*bohrtoa ! cell border, molecule (ang)
    logical, save :: molcubic = .false. ! cell for molecule is cubic

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! crystal or molecule
    ldum = iw_radiobutton("Crystal",bool=ismolecule,boolval=.false.)
    call iw_tooltip("The new structure is a periodic crystal",ttshown)
    ldum = iw_radiobutton("Molecule",bool=ismolecule,boolval=.true.,sameline=.true.)
    call iw_tooltip("The new structure is a molecule",ttshown)

    ! name
    call iw_text("Name",highlight=.true.)
    str2 = "##name"
    ldum = igInputText(c_loc(str2),c_loc(namebuf),maxnamebuf-1,ImGuiInputTextFlags_None,c_null_funptr,c_null_ptr)

    if (.not.ismolecule) then
       ! symmetry
       call iw_text("Symmetry",highlight=.true.)

       ldum = iw_radiobutton("Detect",int=symopt,intval=1_c_int)
       call iw_tooltip("Calculate the symmetry operations from the complete list of unit cell atoms",ttshown)
       ldum = iw_radiobutton("Space group",int=symopt,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Choose a space group from the list",ttshown)
       ldum = iw_radiobutton("Manual",int=symopt,intval=3_c_int,sameline=.true.)
       call iw_tooltip("Enter symmetry operations in cif-file format (e.g. -1/2+x,z,-y)",ttshown)

       ! symmetry options
       if (symopt == 2) then
          call draw_spg_table(ispg_selected)
       elseif (symopt == 3) then
          ! manual input of symmetry operations
          call igGetContentRegionAvail(szavail)
          sz%x = szavail%x
          sz%y = 9 * igGetTextLineHeightWithSpacing()
          str = "##symopsmanual" // c_null_char
          ldum = igInputTextMultiline(c_loc(str),c_loc(symopbuf),maxsymopbuf,sz,&
             ImGuiInputTextFlags_AllowTabInput,c_null_funptr,c_null_ptr)
       end if

       ! lattice: header
       call iw_text("Lattice",highlight=.true.)

       ldum = iw_radiobutton("Cell parameters",int=cellopt,intval=1_c_int)
       call iw_tooltip("Lattice is calculated from the cell lengths and angles",ttshown)
       ldum = iw_radiobutton("Lattice vectors",int=cellopt,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Lattice vectors are given explicitly",ttshown)

       ! units and scale factor only if lattice vectors
       if (cellopt == 2) then
          call iw_combo_simple("Units##cel","Bohr" // c_null_char // "Angstrom" // c_null_char,&
             iunitcel,sameline=.true.)
          call iw_tooltip("Units for the cell parameters/lattice vectors",ttshown)

          call igSameLine(0._c_float,-1._c_float)
          str = "Scale" // c_null_char
          stropt = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(7,1))
          ldum = igInputFloat(c_loc(str),scale,0._c_float,0._c_float,&
             c_loc(stropt),ImGuiInputTextFlags_None)
          call iw_tooltip("Scale factor to multiply the lattice vectors",ttshown)
          call igPopItemWidth()
       end if

       ! lattice: body
       if (cellopt == 1) then
          ! cell lengths
          call iw_text("Cell lengths (Å): ")
          call igSameLine(0._c_float,-1._c_float)
          str = "##celllength3" // c_null_char
          stropt = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(3*8,2))
          ldum = igInputFloat3(c_loc(str),aa,c_loc(stropt),ImGuiInputTextFlags_None)

          ! cell angles
          call iw_text("Cell angles (°):  ")
          call igSameLine(0._c_float,-1._c_float)
          str = "##cellangle3" // c_null_char
          stropt = "%.3f" // c_null_char
          call igPushItemWidth(iw_calcwidth(3*8,2))
          ldum = igInputFloat3(c_loc(str),bb,c_loc(stropt),ImGuiInputTextFlags_None)
       else
          ! lattice vectors
          call igGetContentRegionAvail(szavail)
          sz%x = szavail%x
          sz%y = 4 * igGetTextLineHeightWithSpacing()
          str = "##latvecmanual" // c_null_char
          ldum = igInputTextMultiline(c_loc(str),c_loc(latvecbuf),maxlatvecbuf,sz,&
             ImGuiInputTextFlags_AllowTabInput,c_null_funptr,c_null_ptr)
       end if
    end if ! .not.ismolecule

    ! atomic positions: header
    call iw_text("Atomic positions",highlight=.true.) ! molecule
    call iw_text("(?)",sameline=.true.)
    call iw_tooltip("Give the atomic positions for this system as:" // newline //&
       "  <Sy> <x> <y> <z>" // newline //&
       "where Sy is the atomic symbol and x,y,z are the atomic coordinates.")

    ! units, does not apply in the case of crystal/cell parameters
    if (ismolecule .or. cellopt == 2) then
       call igSameLine(0._c_float,-1._c_float)
       call igSetCursorPosX(igGetCursorPosX() + 2 * g%Style%ItemSpacing%x)
       if (ismolecule) then
          call iw_combo_simple("Units","Bohr" // c_null_char // "Angstrom" // c_null_char,&
             iunitat_m,sameline=.true.)
       else
          call iw_combo_simple("Units","Bohr" // c_null_char // "Angstrom" // c_null_char // "Fractional" // c_null_char,&
             iunitat_c,sameline=.true.)
       end if
       call iw_tooltip("Units for the atomic coordinates",ttshown)
    end if

    ! atomic positions: body
    call igGetContentRegionAvail(szavail)
    sz%x = szavail%x
    if (ismolecule) then
       sz%y = max(szavail%y - iw_calcheight(2,1) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    else
       sz%y = max(szavail%y - iw_calcheight(1,0) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    end if
    str = "##atomicpositions" // c_null_char
    ldum = igInputTextMultiline(c_loc(str),c_loc(atposbuf),maxatposbuf,sz,&
       ImGuiInputTextFlags_AllowTabInput,c_null_funptr,c_null_ptr)

    ! options line
    if (ismolecule) then
       call iw_text("Structure options",highlight=.true.)

       ! cell border
       call igIndent(0._c_float)
       str = "Cell border (Å)" // c_null_char
       stropt = "%.3f" // c_null_char
       call igPushItemWidth(iw_calcwidth(7,1))
       ldum = igInputFloat(c_loc(str),rborder,0._c_float,0._c_float,&
          c_loc(stropt),ImGuiInputTextFlags_None)
       call iw_tooltip("Size of the periodic cell border around a new molecule",ttshown)
       call igPopItemWidth()
       call igSameLine(0._c_float,-1._c_float)

       ! cubic cell
       call igSetCursorPosX(igGetCursorPosX() + 2 * g%Style%ItemSpacing%x)
       ldum = iw_checkbox("Cubic cell",molcubic)
       call iw_tooltip("New molecules are created inside a cubic periodic cell",ttshown)
       call igUnindent(0._c_float)
    end if

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! final buttons: ok
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. iw_button("OK")
    if (ok) then
       ! build the input
       lu = fopen_scratch("formatted")

       ! symmetry and cell
       if (.not.ismolecule) then
          ! symmetry
          if (symopt == 2) then
             ! pass the spg string
             write (lu,'("spg ",A)') string(ispg_selected)
          elseif (symopt == 3) then
             ! pass the symm keywords
             call buffer_to_string_array(symopbuf,lu,prefix="symm ")
          end if

          ! cell
          if (cellopt == 1) then
             ! cell parameters
             write (lu,'("cell ",6(A," "),"ang")') (string(aa(i),'f',decimal=10),i=1,3), &
                (string(bb(i),'f',decimal=10),i=1,3)
          else
             ! lattice vectors
             write (lu,'("cartesian ",A)') string(scale,'f',decimal=10)
             if (iunitcel == 1) then
                write (lu,'("ang")')
             else
                write (lu,'("bohr")')
             end if
             call buffer_to_string_array(latvecbuf,lu)
             write (lu,'("endcartesian")')
          end if
       end if

       ! atomic positions
       if (ismolecule) then
          iunitat = iunitat_m
       else
          iunitat = iunitat_c
       end if
       if (iunitat == 0) then
          call buffer_to_string_array(atposbuf,lu,suffix=" bohr")
       elseif (iunitat == 1) then
          call buffer_to_string_array(atposbuf,lu,suffix=" ang")
       else
          call buffer_to_string_array(atposbuf,lu)
       end if

       ! molecular options
       if (ismolecule) then
          if (molcubic) write (lu,'("cubic")')
          write (lu,'("border ",A)') string(rborder/bohrtoa,'f',decimal=10)
       end if

       write (lu,'("end")')

       ! generate the seed
       rewind(lu)
       if (allocated(seed_)) deallocate(seed_)
       allocate(seed_(1))
       if (ismolecule) then
          call seed_(1)%parse_molecule_env(lu,ok)
       else
          call seed_(1)%parse_crystal_env(lu,ok)
       end if
       call fclose(lu)

       ! load the system and initialize
       if (ok.and.seed_(1)%isused) then
          idx = index(namebuf,c_null_char)
          seed_(1)%name = namebuf(1:idx-1)
          call add_systems_from_seeds(1,seed_)
          call launch_initialization_thread()
          doquit = .true.
       else
          deallocate(seed_)
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      ttshown = .false.
      ismolecule = .false.
      iunitat_c = 2
      iunitat_m = 1
      symopt = 1
      ispg_selected = 1
      cellopt = 1
      iunitcel = 0
      scale = 1._c_float
      aa = 10._c_float
      bb = 90._c_float
      rborder = rborder_def*bohrtoa
      molcubic = .false.

      if (allocated(atposbuf)) deallocate(atposbuf)
      allocate(character(len=maxatposbuf+1) :: atposbuf)
      atposbuf(1:26) = "H 0.00000 0.00000 0.00000" // c_null_char

      if (allocated(symopbuf)) deallocate(symopbuf)
      allocate(character(len=maxsymopbuf+1) :: symopbuf)
      symopbuf(1:70) = "## Each line is a symmetry operation of the type: " // newline //&
         "## -1/2+y,z,x" // c_null_char

      if (allocated(latvecbuf)) deallocate(latvecbuf)
      allocate(character(len=maxlatvecbuf+1) :: latvecbuf)
      latvecbuf(1:50) = "## Enter the three lattice vectors line by line" // newline // c_null_char

      if (allocated(namebuf)) deallocate(namebuf)
      allocate(character(len=maxnamebuf+1) :: namebuf)
      namebuf(1:20) = "Input structure" // c_null_char

    end subroutine init_state

    ! terminate the state for this window
    subroutine end_state()
      if (allocated(atposbuf)) deallocate(atposbuf)
      if (allocated(symopbuf)) deallocate(symopbuf)
      if (allocated(latvecbuf)) deallocate(latvecbuf)
      if (allocated(namebuf)) deallocate(namebuf)
    end subroutine end_state

  end subroutine draw_new_struct

  !xx! private procedures

  !> Draw a space group table. Entry ispg (Hall number) is selected.
  subroutine draw_spg_table(ispg)
    use utils, only: iw_text
    use spglib, only: SpglibSpaceGroupType, spg_get_spacegroup_type
    use tools_io, only: ioj_left, string, deblank, stripchar
    integer, intent(inout) :: ispg

    integer :: i
    character(kind=c_char,len=:), allocatable, target :: str
    integer(c_int) :: flags
    type(ImVec2) :: sz, szero
    type(SpglibSpaceGroupType) :: sa

    ! initialize
    szero%x = 0
    szero%y = 0

    ! space group table
    str = "##spacegrouptable" // c_null_char
    flags = ImGuiTableFlags_Borders
    flags = ior(flags,ImGuiTableFlags_RowBg)
    flags = ior(flags,ImGuiTableFlags_ScrollY)
    flags = ior(flags,ImGuiTableFlags_ScrollX)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)

    call igGetContentRegionAvail(sz)
    sz%y = 8 * igGetTextLineHeightWithSpacing()
    if (igBeginTable(c_loc(str),7,flags,sz,0._c_float)) then
       ! set up columns
       str = "Hall" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "ITA" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,1)

       str = "HM short" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,2)

       str = "HM long" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,3)

       str = "choice" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,4)

       str = "system" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,5)

       str = "Hall symbol" // c_null_char
       flags = ImGuiTableColumnFlags_WidthFixed
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,6)
       call igTableSetupScrollFreeze(0, 1) ! top row always visible
       call igTableHeadersRow()

       ! table body
       do i = 1, 530
          sa = spg_get_spacegroup_type(i)

          ! Hall
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
          if (igTableSetColumnIndex(0)) then
             str = string(i) // c_null_char
             flags = ImGuiSelectableFlags_SpanAllColumns
             flags = ior(flags,ImGuiSelectableFlags_SelectOnNav)
             if (igSelectable_Bool(c_loc(str),logical(ispg == i,c_bool),flags,szero)) &
                ispg = i
          end if

          ! ITA
          if (igTableNextColumn()) call iw_text(string(sa%number,5,ioj_left))

          ! HM-short
          if (igTableNextColumn()) then
             str = deblank(sa%international_short)
             call iw_text(trim(stripchar(str,"_")))
          end if

          ! HM-long
          if (igTableNextColumn()) then
             str = deblank(sa%international_full)
             call iw_text(trim(stripchar(str,"_")))
          end if

          ! choice
          if (igTableNextColumn()) call iw_text(string(sa%choice,6,ioj_left))

          ! system
          if (igTableNextColumn()) then
             if (sa%number >= 1 .and. sa%number <= 2) then
                str = "triclinic"
             elseif (sa%number >= 3 .and. sa%number <= 15) then
                str = "monoclinic"
             elseif (sa%number >= 16 .and. sa%number <= 74) then
                str = "orthorhombic"
             elseif (sa%number >= 75 .and. sa%number <= 142) then
                str = "tetragonal"
             elseif (sa%number >= 143 .and. sa%number <= 167) then
                str = "trigonal"
             elseif (sa%number >= 168 .and. sa%number <= 194) then
                str = "hexagonal"
             elseif (sa%number >= 195 .and. sa%number <= 230) then
                str = "cubic"
             else
                str = "??"
             end if
             call iw_text(str)
          end if

          ! hall symbol
          if (igTableNextColumn()) call iw_text(string(sa%hall_symbol))
       end do

       call igTableSetColumnWidthAutoAll(igGetCurrentTable())
       call igEndTable()
    end if

    ! current choice
    sa = spg_get_spacegroup_type(ispg)
    str = deblank(sa%international_full)
    str = "Current choice: " // trim(stripchar(str,"_"))
    if (len_trim(sa%choice) > 0) str = str // " (" // trim(sa%choice) // ")"
    call iw_text(str)

  end subroutine draw_spg_table

end submodule new_struct
