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
    use utils, only: iw_text, iw_tooltip, iw_radiobutton, iw_button, iw_setpos_bottomright,&
       iw_helpermark
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
       isys = win(iwin_tree)%isys
       oksys = ok_system(isys,sys_init)
    end if
    if (.not.oksys) then
       ! reset to the input console selected
       isys = win(iwin_tree)%isys
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
    ldum = iw_radiobutton("From File",int=w%lf%sourceopt,intval=0_c_int)
    call iw_tooltip("Load the field from an external file",ttshown)
    ! ldum = iw_radiobutton("Expression",int=w%lf%sourceopt,intval=1_c_int,sameline=.true.)
    ! call iw_tooltip("Load the field from an arithmetic expression involving other fields",ttshown)

    if (w%lf%sourceopt == 0) then
       ! if a new file is read, detect the format if necessary
       if (w%okfile_set) then
          w%errmsg = ""
          if (w%itoken == itoken_file1) then
             w%lf%file1 = w%okfile
             w%lf%file1_format = w%okfile_format
             if (w%lf%file1_format == 0) &
                w%lf%file1_format = -field_detect_format(file=w%okfile)
             if (abs(w%lf%file1_format) == ifformat_pi) w%lf%file1_format = 0 ! aiPI deactivated for now
             if (w%lf%file1_format == 0) then
                w%errmsg = "Unknown field file extension"
                w%okfile = ""
             end if
             w%lf%file1_set = .true.
             w%lf%file2 = ""
             w%lf%file2_set = .false.
             w%lf%file3 = ""
             w%lf%file3_set = .false.
          elseif (w%itoken == itoken_file2) then
             w%lf%file2 = w%okfile
             w%lf%file2_set = .true.
          elseif (w%itoken == itoken_file3) then
             w%lf%file3 = w%okfile
             w%lf%file3_set = .true.
          end if
       end if
       ! source file and format
       call iw_text("Source",highlight=.true.)
       if (iw_button("File",danger=(w%lf%file1_format==0))) then
          iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfieldfile,&
             idparent=w%id,itoken=itoken_file1)
       end if
       call iw_tooltip("File from where the field is read",ttshown)
       if (len(w%lf%file1) > 0) then
          if (iw_button("Clear##clearfile1",sameline=.true.,danger=.true.)) then
             w%lf%file1_format = 0
             w%lf%file1 = ""
             w%lf%file1_set = .false.
             w%lf%file1_fmtstr = ""
             w%lf%file2 = ""
             w%lf%file2_set = .false.
             w%lf%file3 = ""
             w%lf%file3_set = .false.
          end if
       end if
       call iw_text(w%lf%file1,sameline=.true.)

       isgrid = .false.
       if (w%lf%file1_format /= 0) then
          call iw_text("Format: ")
          select case (abs(w%lf%file1_format))
          case(ifformat_wien)
             w%lf%file1_fmtstr = "WIEN"
             call iw_text("WIEN2k clmsum-style file (LAPW)",sameline=.true.)
             ! default struct file
             if (.not.w%lf%file2_set) then
                ll = len_trim(sysc(isys)%seed%file)
                if (ll >= 7) then
                   if (sysc(isys)%seed%file(ll-6:ll) == ".struct") &
                      w%lf%file2 = sysc(isys)%seed%file
                end if
                w%lf%file2_set = .true.
             end if
          case(ifformat_elk)
             w%lf%file1_fmtstr = "ELK"
             call iw_text("elk STATE.OUT file (LAPW)",sameline=.true.)
             ! default GEOMETRY.OUT file
             if (.not.w%lf%file2_set) then
                ll = len_trim(sysc(isys)%seed%file)
                if (ll >= 4) then
                   if (sysc(isys)%seed%file(ll-3:ll) == ".OUT") &
                      w%lf%file2 = sysc(isys)%seed%file
                end if
                w%lf%file2_set = .true.
             end if
          case(ifformat_pi)
             w%lf%file1_fmtstr = "PI"
             call iw_text("aiPI ion files (LCAO)",sameline=.true.)
          case(ifformat_cube)
             w%lf%file1_fmtstr = "CUBE"
             call iw_text("Cube file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_bincube)
             w%lf%file1_fmtstr = "BINCUBE"
             call iw_text("Binary cube file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_abinit)
             w%lf%file1_fmtstr = "ABINIT"
             call iw_text("Abinit DEN-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_vasp)
             w%lf%file1_fmtstr = "VASP"
             call iw_text("VASP CHGCAR-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_vaspnov)
             w%lf%file1_fmtstr = "VASPNOV"
             call iw_text("VASP ELFCAR-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_qub)
             w%lf%file1_fmtstr = "QUB"
             call iw_text("Aimpac qub file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_xsf)
             w%lf%file1_fmtstr = "XSF"
             call iw_text("Xcrysden xsf-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_elkgrid)
             w%lf%file1_fmtstr = "ELKGRID"
             call iw_text("elk grid file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_siestagrid)
             w%lf%file1_fmtstr = "SIESTA"
             call iw_text("SIESTA RHO-style file (grid)",sameline=.true.)
             isgrid = .true.
          case(ifformat_dftb)
             w%lf%file1_fmtstr = "DFTB"
             call iw_text("DFTB+ wavefunction (LCAO)",sameline=.true.)
             ! default eigenvec.bin file
             if (.not.w%lf%file2_set) then
                idx = index(sysc(isys)%fullname,dirsep,.true.)
                w%lf%file2 = sysc(isys)%fullname(1:idx) // "eigenvec.bin"
                inquire(file=w%lf%file2,exist=ok)
                if (.not.ok) w%lf%file2 = ""
                w%lf%file2_set = .true.
             end if
             ! default hsd file (from the list of candidate names)
             if (.not.w%lf%file3_set) then
                idx = index(sysc(isys)%fullname,dirsep,.true.)
                ok = .false.
                do j = 1, size(hsdnames,1)
                   w%lf%file3 = sysc(isys)%fullname(1:idx) // trim(hsdnames(j))
                   inquire(file=w%lf%file3,exist=ok)
                   if (ok) exit
                end do
                if (.not.ok) w%lf%file3 = ""
                w%lf%file3_set = .true.
             end if
          case(ifformat_pwc)
             w%lf%file1_fmtstr = "PWC"
             call iw_text("Quantum ESPRESSO pwc file (grid+wfn)",sameline=.true.)
             isgrid = .true.
          case(ifformat_wfn)
             w%lf%file1_fmtstr = "WFN"
             call iw_text("Gaussian wfn wavefunction file (molecular wfn)",sameline=.true.)
          case(ifformat_wfx)
             w%lf%file1_fmtstr = "WFX"
             call iw_text("Gaussian wfx wavefunction file (molecular wfn)",sameline=.true.)
          case(ifformat_fchk)
             w%lf%file1_fmtstr = "FCHK"
             call iw_text("Gaussian fchk wavefunction file (molecular wfn+basis set)",sameline=.true.)
          case(ifformat_molden)
             w%lf%file1_fmtstr = "MOLDEN"
             call iw_text("molden wavefunction file (molecular wfn+basis set)",sameline=.true.)
          end select
          if (w%lf%file1_format < 0) &
             call iw_text(" [auto-detected]",sameline=.true.)
       end if

       ! format-specific options and additional files
       select case (abs(w%lf%file1_format))
       case(ifformat_wien)
          if (iw_button("File (.struct)",danger=.not.w%lf%file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("WIEN2k structure file (.struct) used to interpret the clmsum-style file",ttshown)
          if (len(w%lf%file2) > 0) then
             if (iw_button("Clear##clearfile2",sameline=.true.,danger=.true.)) &
                w%lf%file2 = ""
          end if
          call iw_text(w%lf%file2,sameline=.true.)
       case(ifformat_elk)
          if (iw_button("File (GEOMETRY.OUT)",danger=.not.w%lf%file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("GEOMETRY.OUT structure file used to interpret the STATE.OUT",ttshown)
          if (len(w%lf%file2) > 0) then
             if (iw_button("Clear##clearfile3",sameline=.true.,danger=.true.)) &
                w%lf%file2 = ""
          end if
          call iw_text(w%lf%file2,sameline=.true.)
       case(ifformat_dftb)
          if (iw_button("File (.bin)",danger=.not.w%lf%file2_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file2)
          end if
          call iw_tooltip("eigenvec.bin file for reading the DFTB+ wavefunction",ttshown)
          if (len(w%lf%file2) > 0) then
             if (iw_button("Clear##clearfile4",sameline=.true.,danger=.true.)) &
                w%lf%file2 = ""
          end if
          call iw_text(w%lf%file2,sameline=.true.)

          if (iw_button("File (.hsd)",danger=.not.w%lf%file3_set,disabled=disabled)) then
             iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
                idparent=w%id,itoken=itoken_file3)
          end if
          call iw_tooltip("hsd file for reading the DFTB+ wavefunction",ttshown)
          if (len(w%lf%file3) > 0) then
             if (iw_button("Clear##clearfile5",sameline=.true.,danger=.true.)) &
                w%lf%file3 = ""
          end if
          call iw_text(w%lf%file3,sameline=.true.)
       end select

       ! error message, if applicable
       if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    elseif (w%lf%sourceopt == 1) then
       ! from expression
       w%errmsg = ""
    end if
    w%okfile_set = .false.

    if (w%lf%sourceopt == 0 .and. isgrid) then
       call iw_text("Options",highlight=.true.)
       ! a help marker, not a tooltip on the text: a plain text has no item id,
       ! so a delayed (ttshown) tooltip attached to it never fires
       call iw_text("Interpolation method")
       call iw_helpermark("Choose the interpolation method for the grid")
       ldum = iw_radiobutton("Nearest",int=w%lf%iginterp,intval=0_c_int)
       call iw_tooltip("Value of the nearest grid point",ttshown)
       ldum = iw_radiobutton("Tri-linear",int=w%lf%iginterp,intval=1_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional linear interpolation",ttshown)
       ldum = iw_radiobutton("Tri-spline",int=w%lf%iginterp,intval=2_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional spline interpolation",ttshown)
       ldum = iw_radiobutton("Tri-cubic",int=w%lf%iginterp,intval=3_c_int,sameline=.true.)
       call iw_tooltip("Three-dimensional cubic polynomial interpolation",ttshown)
       ldum = iw_radiobutton("Smoothrho",int=w%lf%iginterp,intval=4_c_int,sameline=.true.)
       call iw_tooltip("Polyharmonic splines + smoothing, all-electron densities only (recommended for them)",ttshown)
    end if

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(8,2)

    ! calculated whether we have enough info to continue to ok
    disabled = (len(w%lf%file1) == 0)
    if (.not.disabled) then
       iff = abs(w%lf%file1_format)
       if (iff == ifformat_wien .or. iff == ifformat_elk .or. iff == ifformat_dftb) &
          disabled = disabled .or. (len(w%lf%file2) == 0)
       if (iff == ifformat_dftb) &
          disabled = disabled .or. (len(w%lf%file3) == 0)
    end if
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. iw_button("OK",disabled=disabled)
    if (ok) then
       if (w%lf%sourceopt == 0) then
          loadstr = w%lf%file1_fmtstr // " " // trim(w%lf%file1) // " " // trim(w%lf%file2) // " " // trim(w%lf%file3)
          loadstr = loadstr // " notestmt"
          call sys(isys)%load_field_string(loadstr,.false.,iff,w%errmsg)
          if (len_trim(w%errmsg) == 0) doquit = .true.
       elseif (w%lf%sourceopt == 1) then
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
      w%lf%sourceopt = 0_c_int
      w%lf%file1 = ""
      w%lf%file1_set = .false.
      w%lf%file1_format = 0
      w%lf%file1_fmtstr = ""
      w%lf%file2 = ""
      w%lf%file2_set = .false.
      w%lf%file3 = ""
      w%lf%file3_set = .false.
    end subroutine init_state
    ! terminate the state for this window
    subroutine end_state()
      w%lf%file1 = ""
      w%lf%file1_set = .false.
      w%lf%file1_format = 0
      w%lf%file1_fmtstr = ""
      w%lf%file2 = ""
      w%lf%file2_set = .false.
      w%lf%file3 = ""
      w%lf%file3_set = .false.
      if (allocated(w%lf%file1)) deallocate(w%lf%file1)
      if (allocated(w%lf%file1_fmtstr)) deallocate(w%lf%file1_fmtstr)
      if (allocated(w%lf%file2)) deallocate(w%lf%file2)
      if (allocated(w%lf%file3)) deallocate(w%lf%file3)
    end subroutine end_state
  end subroutine draw_load_field

end submodule load_field
