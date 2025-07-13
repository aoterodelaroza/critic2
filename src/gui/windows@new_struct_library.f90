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

! Routines for the new structure from library window.
submodule (windows) new_struct_library
  use interfaces_cimgui
  implicit none
contains

  !> Draw the contents of the new structure from library window.
  module subroutine draw_new_struct_library(w)
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use gui_main, only: g
    use systems, only: add_systems_from_seeds, launch_initialization_thread, system_shorten_names
    use utils, only: igIsItemHovered_delayed, iw_tooltip, iw_button, iw_text, iw_calcheight,&
       iw_calcwidth, iw_radiobutton, iw_checkbox
    use crystalseedmod, only: crystalseed, realloc_crystalseed
    use spglib, only: SpglibSpaceGroupType, spg_get_spacegroup_type
    use global, only: clib_file, mlib_file, rborder_def
    use tools_io, only: string, fopen_scratch, fclose, stripchar, deblank
    use types, only: realloc, vstring
    use param, only: bohrtoa
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str, stropt
    logical(c_bool) :: ldum, doquit
    logical :: ok, saveismol, doubleclicked, isset
    type(ImVec2) :: szero, sz, szavail
    integer :: i, nseed, idum
    type(crystalseed) :: seed
    type(crystalseed), allocatable :: seed_(:)

    ! window state
    logical, save :: ttshown = .false. ! tooltip flag
    logical, save :: ismolecule = .false. ! molecule/crystal choice at the top
    integer, save :: nst = 0 ! number of library structures
    type(vstring), allocatable, save :: st(:) ! library structures
    logical(c_bool), allocatable, save :: lst(:) ! selected structures (1->nst)
    integer, save :: lastselected = 0 ! last selected structure (for shift/ctrl selection)
    real(c_float), save :: rborder = rborder_def*bohrtoa
    logical, save :: molcubic = .false.

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! crystal or molecule, from library
    saveismol = ismolecule
    if (iw_radiobutton("Crystal",bool=ismolecule,boolval=.false.)) then
       if (saveismol.and..not.w%okfile_set) w%okfile = trim(clib_file)
       w%okfile_read = .true.
    end if
    call iw_tooltip("The new structure is a periodic crystal",ttshown)
    if (iw_radiobutton("Molecule",bool=ismolecule,boolval=.true.,sameline=.true.)) then
       if (.not.saveismol.and..not.w%okfile_set) w%okfile = trim(mlib_file)
       w%okfile_read = .true.
    end if
    call iw_tooltip("The new structure is a molecule",ttshown)

    ! render the rest of the window
    ! library file
    call iw_text("Source",highlight=.true.)
    if (iw_button("Library file")) then
       idum = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openlibraryfile,&
          idparent=w%id,orraise=-1)
    end if
    call iw_tooltip("Library file from where the structures are read",ttshown)
    call iw_text(w%okfile,sameline=.true.)

    ! list box
    call iw_text("Structures to load from the library file",highlight=.true.)
    str = "##listbox" // c_null_char
    call igGetContentRegionAvail(sz)
    if (ismolecule) then
       sz%y = max(sz%y - iw_calcheight(2,1) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    else
       sz%y = max(sz%y - iw_calcheight(1,0) - g%Style%ItemSpacing%y,igGetTextLineHeightWithSpacing())
    end if

    doubleclicked = .false.
    isset = .false.
    if (igBeginListBox(c_loc(str),sz)) then
       do i = 1, nst
          str = st(i)%s // c_null_char
          if (igSelectable_Bool(c_loc(str), lst(i), ImGuiSelectableFlags_AllowDoubleClick, szero)) then

             ! implement selection range with shift and control
             if (igIsKeyDown(ImGuiKey_ModShift).and.lastselected /= 0.and.lastselected /= i) then
                ! selecte a whole range
                lst = .false.
                if (lastselected > i) then
                   lst(i:lastselected) = .true.
                else
                   lst(lastselected:i) = .true.
                end if
             elseif (igIsKeyDown(ImGuiKey_ModCtrl)) then
                ! select and individual item and accumulate
                lst(i) = .true.
             else
                ! select one item, start range, remove all others
                lst = .false.
                lst(i) = .true.
                lastselected = i
                doubleclicked = igIsMouseDoubleClicked(ImGuiPopupFlags_MouseButtonLeft)
                isset = .true.
             end if
          end if
       end do
       call igEndListBox()
    end if

    if (ismolecule) then
       ! options line
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
    end if

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! final buttons: ok
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. doubleclicked
    ok = ok .or. iw_button("OK")
    if (ok) then
       nseed = count(lst(1:nst))
       if (nseed > 0) then
          ! we have systems (potentially), allocate the seed array
          allocate(seed_(nseed))
          nseed = 0
          do i = 1, nst
             if (lst(i)) then
                ! read the selected seeds from the library, check OK
                call seed%read_library(st(i)%s,ok,file=w%okfile)
                if (ok .and. ismolecule.eqv.seed%ismolecule) then
                   nseed = nseed + 1
                   seed_(nseed) = seed
                   seed_(nseed)%border = rborder/bohrtoa
                   seed_(nseed)%cubic = molcubic
                end if
             end if
          end do
          if (nseed > 0) then
             ! load the seeds as new systems
             call realloc_crystalseed(seed_,nseed)
             call add_systems_from_seeds(nseed,seed_)
             call launch_initialization_thread()
             call system_shorten_names()
          end if
       end if
       doquit = .true.
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! read the library file
    if (w%okfile_read) then
       call get_library_structure_list(w%okfile,nst,st,ismolecule)
       if (allocated(lst)) deallocate(lst)
       allocate(lst(nst))
       lst = .false.
       lastselected = 0
       w%okfile_read = .false.
    end if

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
      nst = 0
      lastselected = 0
      rborder = rborder_def*bohrtoa
      molcubic = .false.
      w%okfile = trim(clib_file)
      w%okfile_set = .false.
      w%okfile_read = .true.
    end subroutine init_state

    ! terminate the state for this window
    subroutine end_state()
      nst = 0
      if (allocated(st)) deallocate(st)
      if (allocated(lst)) deallocate(lst)
    end subroutine end_state

  end subroutine draw_new_struct_library

  !xx! private procedures

  !> Get the structure list from the library file. If ismol = .true.,
  !> read only the moelcules. If ismol = .false., read only the
  !> crystals.
  subroutine get_library_structure_list(libfile,nst,st,ismol)
    use tools_io, only: fopen_read, fclose, lgetword, getword, equal, getline
    use types, only: vstring, realloc
    character(len=:), allocatable, intent(in) :: libfile
    integer, intent(out) :: nst
    type(vstring), allocatable, intent(inout) :: st(:)
    logical, intent(in) :: ismol

    integer :: lu, lp
    character(len=:), allocatable :: word, name, line
    logical :: ok

    ! open the file
    nst = 0
    inquire(file=libfile,exist=ok)
    if (.not.ok) return
    lu = fopen_read(libfile)
    if (lu < 0) return

    ! preallocate
    if (allocated(st)) deallocate(st)
    allocate(st(10))

    ! read the structures
    main: do while (getline(lu,line))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,'structure')) then
          name = getword(line,lp)

          ok = getline(lu,line)
          if (.not.ok) exit main
          lp = 1
          word = lgetword(line,lp)
          ok = equal(word,'crystal').and..not.ismol
          ok = ok .or. (equal(word,'molecule').and.ismol)
          if (.not.ok) cycle main

          nst = nst + 1
          if (nst > size(st,1)) call realloc(st,2*nst)
          st(nst)%s = name
       endif
    end do main
    if (nst > 0) &
       call realloc(st,nst)
    ! clean up
    call fclose(lu)

  end subroutine get_library_structure_list


end submodule new_struct_library
