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

! Routines for the scfplot window.
submodule (windows) load_field
  use interfaces_cimgui
  implicit none
contains

  !> Draw the contents of the load field window.
  module subroutine draw_load_field(w)
    use fieldseedmod, only: field_detect_format
    use param, only: ifformat_wien, ifformat_elk, ifformat_pi,&
       ifformat_cube, ifformat_bincube, ifformat_abinit, ifformat_vasp,&
       ifformat_vaspnov, ifformat_qub, ifformat_xsf, ifformat_elkgrid,&
       ifformat_siestagrid, ifformat_dftb, ifformat_pwc, ifformat_wfn,&
       ifformat_wfx, ifformat_fchk, ifformat_molden, dirsep
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use systems, only: nsys, sysc, sys, sys_init, ok_system
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_radiobutton, iw_button,&
       iw_calcwidth
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: oksys, ok, doquit, disabled
    integer :: isys, i, j, ll, idx, iff, iaux
    type(ImVec2) :: szavail, szero
    real(c_float) :: combowidth
    logical(c_bool) :: is_selected, ldum
    character(len=:,kind=c_char), allocatable, target :: str1, str2, loadstr
    logical :: isgrid

    integer, parameter :: itoken_file1 = 1
    integer, parameter :: itoken_file2 = 2
    integer, parameter :: itoken_file3 = 3

    ! window state
    logical, save :: ttshown = .false. ! tooltip flag
    integer(c_int), save :: sourceopt = 0_c_int ! 0 = file, 1 = expression
    character(len=:,kind=c_char), allocatable, target, save :: file1 ! first (main) file
    character(len=:,kind=c_char), allocatable, target, save :: file1_fmtstr ! first file format string
    logical, save :: file1_set = .false.
    integer, save :: file1_format = 0
    character(len=:,kind=c_char), allocatable, target, save :: file2 ! first auxiliary file
    logical, save :: file2_set = .false.
    character(len=:,kind=c_char), allocatable, target, save :: file3 ! second auxiliaryy file
    logical, save :: file3_set = .false.
    integer, save :: iginterp = 3 ! 0 = nearest, 1 = trilinear, 2 = trispline, 3 = tricubic, 4 = smoothrho

    ! DFTB+ hsd name candidates
    character*15, parameter :: hsdnames(4) = (/&
       "wfc-3ob-3-1.hsd","wfc.3ob-3-1.hsd","wfc.pbc-0-3.hsd","wfc.mio-1-1.hsd"/)

    ! first pass when opened, reset the state
    if (w%firstpass) call init_state()

    ! initialize
    szero%x = 0
    szero%y = 0
    doquit = .false.

    !! make sure we get a system that exists
    ! check if the system still exists
    isys = w%isys
    oksys = ok_system(isys,sys_init)
    if (.not.oksys) then
       ! reset to the table selected
       isys = win(iwin_tree)%tree_selected
       oksys = ok_system(isys,sys_init)
    end if
    if (.not.oksys) then
       ! reset to the input console selected
       isys = win(iwin_tree)%inpcon_selected
       oksys = ok_system(isys,sys_init)
    end if
    if (.not.oksys) then
       ! check all systems and try to find one that is OK
       do i = 1, nsys
          oksys = ok_system(isys,sys_init)
          if (oksys) then
             isys = i
             exit
          end if
       end do
    end if
    if (.not.oksys) then
       ! this dialog does not make sense anymore, close it and exit
       call end_state()
       call w%end()
       return
    end if

    ! system combo
    call iw_text("System",highlight=.true.)
    call igSameLine(0._c_float,-1._c_float)
    call igGetContentRegionAvail(szavail)
    combowidth = max(szavail%x - g%Style%ItemSpacing%x,0._c_float)
    str1 = "##systemcombo" // c_null_char
    call igSetNextItemWidth(combowidth)
    if (oksys) &
       str2 = string(isys) // ": " // trim(sysc(isys)%seed%name) // c_null_char
    if (igBeginCombo(c_loc(str1),c_loc(str2),ImGuiComboFlags_None)) then
       do i = 1, nsys
          if (sysc(i)%status == sys_init) then
             is_selected = (isys == i)
             str2 = string(i) // ": " // trim(sysc(i)%seed%name) // c_null_char
             if (igSelectable_Bool(c_loc(str2),is_selected,ImGuiSelectableFlags_None,szero)) then
                w%isys = i
                isys = w%isys
             end if
             if (is_selected) &
                call igSetItemDefaultFocus()
          end if
       end do
       call igEndCombo()
    end if
    call iw_tooltip("Load a field for this system",ttshown)

    ! select the source
    ldum = iw_radiobutton("From File",int=sourceopt,intval=0_c_int)
    call iw_tooltip("Load the field from an external file",ttshown)
    ! ldum = iw_radiobutton("Expression",int=sourceopt,intval=1_c_int,sameline=.true.)
    ! call iw_tooltip("Load the field from an arithmetic expression involving other fields",ttshown)

    if (sourceopt == 0) then
       ! if a new file is read, detect the format if necessary
       if (w%okfile_set) then
          w%errmsg = ""
          if (w%itoken == itoken_file1) then
             file1 = w%okfile
             file1_format = w%okfile_format
             if (file1_format == 0) &
                file1_format = -field_detect_format(file=w%okfile)
             if (abs(file1_format) == ifformat_pi) file1_format = 0 ! aiPI deactivated for now
             if (file1_format == 0) then
                w%errmsg = "Unknown field file extension"
                w%okfile = ""
             end if
             file1_set = .true.
             file2 = ""
             file2_set = .false.
             file3 = ""
             file3_set = .false.
          elseif (w%itoken == itoken_file2) then
             file2 = w%okfile
             file2_set = .true.
          elseif (w%itoken == itoken_file3) then
             file3 = w%okfile
             file3_set = .true.
          end if
       end if
       ! source file and format
       call iw_text("Source",highlight=.true.)
       if (iw_button("File",danger=(file1_format==0))) then
          iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfieldfile,&
             idparent=w%id,itoken=itoken_file1)
       end if
       call iw_tooltip("File from where the field is read",ttshown)
       if (len(file1) > 0) then
          if (iw_button("Clear##clearfile1",sameline=.true.,danger=.true.)) then
             file1_format = 0
             file1 = ""
             file1_set = .false.
             file1_fmtstr = ""
             file2 = ""
             file2_set = .false.
             file3 = ""
             file3_set = .false.
          end if
       end if
       call iw_text(file1,sameline=.true.)

       isgrid = .false.
       if (file1_format /= 0) then
          call iw_text("Format: ")
          select case (abs(file1_format))
          case(ifformat_wien)
             file1_fmtstr = "WIEN"
             call iw_text("WIEN2k clmsum-style file (LAPW)",sameline=.true.)
             ! default struct file
             if (.not.file2_set) then
                ll = len_trim(sysc(isys)%seed%file)
                if (sysc(isys)%seed%file(ll-6:ll) == ".struct") &
                   file2 = sysc(isys)%seed%file
                file2_set = .true.
             end if
          case(ifformat_elk)
             file1_fmtstr = "ELK"
             call iw_text("elk STATE.OUT file (LAPW)",sameline=.true.)
             ! default GEOMETRY.OUT file
             if (.not.file2_set) then
                ll = len_trim(sysc(isys)%seed%file)
                if (sysc(isys)%seed%file(ll-3:ll) == ".OUT") &
                   file2 = sysc(isys)%seed%file
                file2_set = .true.
             end if
          case(ifformat_pi)
             file1_fmtstr = "PI"
             call iw_text("aiPI ion files (LCAO)",sameline=.true.)
          case(ifformat_cube)
             file1_fmtstr = "CUBE"
             call iw_text("Cube file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_bincube)
             file1_fmtstr = "BINCUBE"
             call iw_text("Binary cube file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_abinit)
             file1_fmtstr = "ABINIT"
             call iw_text("Abinit DEN-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_vasp)
             file1_fmtstr = "VASP"
             call iw_text("VASP CHGCAR-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_vaspnov)
             file1_fmtstr = "VASPNOV"
             call iw_text("VASP ELFCAR-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_qub)
             file1_fmtstr = "QUB"
             call iw_text("Aimpac qub file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_xsf)
             file1_fmtstr = "XSF"
             call iw_text("Xcrysden xsf-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_elkgrid)
             file1_fmtstr = "ELKGRID"
             call iw_text("elk grid file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_siestagrid)
             file1_fmtstr = "SIESTA"
             call iw_text("SIESTA RHO-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_dftb)
             file1_fmtstr = "DFTB"
             call iw_text("DFTB+ wavefunction (LCAO)",sameline=.true.)
             ! default eigenvec.bin file
             if (.not.file2_set) then
                idx = index(sysc(isys)%fullname,dirsep,.true.)
                file2 = sysc(isys)%fullname(1:idx) // "eigenvec.bin"
                inquire(file=file2,exist=ok)
                if (.not.ok) file2 = ""
                file2_set = .true.
             end if
             ! default hsd file (from the list of candidate names)
             if (.not.file3_set) then
                idx = index(sysc(isys)%fullname,dirsep,.true.)
                ok = .false.
                do j = 1, size(hsdnames,1)
                   file3 = sysc(isys)%fullname(1:idx) // trim(hsdnames(j))
                   inquire(file=file3,exist=ok)
                   if (ok) exit
                end do
                if (.not.ok) file3 = ""
                file3_set = .true.
             end if
          case(ifformat_pwc)
             file1_fmtstr = "PWC"
             call iw_text("Quantum ESPRESSO pwc file (grid+wfn)",sameline=.true.)
             isgrid = .true.
          case(ifformat_wfn)
             file1_fmtstr = "WFN"
             call iw_text("Gaussian wfn wavefunction file (molecular wfn)",sameline=.true.)
          case(ifformat_wfx)
             file1_fmtstr = "WFX"
             call iw_text("Gaussian wfx wavefunction file (molecular wfn)",sameline=.true.)
          case(ifformat_fchk)
             file1_fmtstr = "FCHK"
             call iw_text("Gaussian fchk wavefunction file (molecular wfn+basis set)",sameline=.true.)
          case(ifformat_molden)
             file1_fmtstr = "MOLDEN"
             call iw_text("molden wavefunction file (molecular wfn+basis set)",sameline=.true.)
          end select
          if (file1_format < 0) &
             call iw_text(" [auto-detected]",sameline=.true.)
       end if

       ! format-specific options and additional files
       select case (abs(file1_format))
       case(ifformat_wien)
          if (iw_button("File (.struct)",danger=.not.file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("WIEN2k structure file (.struct) used to interpret the clmsum-style file",ttshown)
          if (len(file2) > 0) then
             if (iw_button("Clear##clearfile2",sameline=.true.,danger=.true.)) &
                file2 = ""
          end if
          call iw_text(file2,sameline=.true.)
       case(ifformat_elk)
          if (iw_button("File (GEOMETRY.OUT)",danger=.not.file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("GEOMETRY.OUT structure file used to interpret the STATE.OUT",ttshown)
          if (len(file2) > 0) then
             if (iw_button("Clear##clearfile3",sameline=.true.,danger=.true.)) &
                file2 = ""
          end if
          call iw_text(file2,sameline=.true.)
       case(ifformat_dftb)
          if (iw_button("File (.bin)",danger=.not.file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("eigenvec.bin file for reading the DFTB+ wavefunction",ttshown)
          if (len(file2) > 0) then
             if (iw_button("Clear##clearfile4",sameline=.true.,danger=.true.)) &
                file2 = ""
          end if
          call iw_text(file2,sameline=.true.)

          if (iw_button("File (.hsd)",danger=.not.file3_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file3)
          end if
          call iw_tooltip("hsd file for reading the DFTB+ wavefunction",ttshown)
          if (len(file3) > 0) then
             if (iw_button("Clear##clearfile5",sameline=.true.,danger=.true.)) &
                file3 = ""
          end if
          call iw_text(file3,sameline=.true.)
       end select

       ! error message, if applicable
       if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    elseif (sourceopt == 1) then
       ! from expression
       w%errmsg = ""
    end if
    w%okfile_set = .false.

    if (sourceopt == 0 .and. isgrid) then
       call iw_text("Options",highlight=.true.)
       call iw_text("Interpolation method")
       call iw_tooltip("Choose the interpolation method for the grid",ttshown)
       ldum = iw_radiobutton("Nearest",int=iginterp,intval=0_c_int)
       call iw_tooltip("Value of the nearest grid point",ttshown)
       ldum = iw_radiobutton("Tri-linear",int=iginterp,intval=1_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional linear interpolation",ttshown)
       ldum = iw_radiobutton("Tri-spline",int=iginterp,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional spline interpolation",ttshown)
       ldum = iw_radiobutton("Tri-cubic",int=iginterp,intval=3_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional cubic polynomial interpolation",ttshown)
       ldum = iw_radiobutton("Smoothrho",int=iginterp,intval=4_c_int,sameline=.true.)
       call iw_tooltip("Polyharmonic splines + smoothing, all-electron densities only (recommended for them)",ttshown)
    end if

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! calculated whether we have enough info to continue to ok
    disabled = (len(file1) == 0)
    if (.not.disabled) then
       iff = abs(file1_format)
       if (iff == ifformat_wien .or. iff == ifformat_elk .or. iff == ifformat_dftb) &
          disabled = disabled .or. (len(file2) == 0)
       if (iff == ifformat_dftb) &
          disabled = disabled .or. (len(file3) == 0)
    end if
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. iw_button("OK",disabled=disabled)
    if (ok) then
       if (sourceopt == 0) then
          loadstr = file1_fmtstr // " " // trim(file1) // " " // trim(file2) // " " // trim(file3)
          loadstr = loadstr // " notestmt"
          call sys(isys)%load_field_string(loadstr,.false.,iff,w%errmsg)
          if (len_trim(w%errmsg) == 0) doquit = .true.
       elseif (sourceopt == 1) then
          ! todo
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) &
       doquit = .true.

    ! quit the window
    if (doquit) then
       call end_state()
       call w%end()
    end if

  contains
    ! initialize the state for this window
    subroutine init_state()
      w%errmsg = ""
      sourceopt = 0_c_int
      file1 = ""
      file1_set = .false.
      file1_format = 0
      file1_fmtstr = ""
      file2 = ""
      file2_set = .false.
      file3 = ""
      file3_set = .false.
    end subroutine init_state
    ! terminate the state for this window
    subroutine end_state()
      file1 = ""
      file1_set = .false.
      file1_format = 0
      file1_fmtstr = ""
      file2 = ""
      file2_set = .false.
      file3 = ""
      file3_set = .false.
      if (allocated(file1)) deallocate(file1)
      if (allocated(file1_fmtstr)) deallocate(file1_fmtstr)
      if (allocated(file2)) deallocate(file2)
      if (allocated(file3)) deallocate(file3)
    end subroutine end_state
  end subroutine draw_load_field

end submodule load_field
