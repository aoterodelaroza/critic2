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

! Routines for the load field window.
submodule (windows) load_field
  use interfaces_cimgui
  implicit none
contains

  !> Draw the contents of the load field window.
  module subroutine draw_load_field(w)
    use fieldseedmod, only: field_detect_format, field_format_string
    use fieldmod, only: type_grid
    use param, only: ifformat_wien, ifformat_elk, ifformat_pi, ifformat_vasp,&
       ifformat_vaspnov, ifformat_dftb, ifformat_pwc, ifformat_fchk, ifformat_molden,&
       dirsep, newline
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG,&
       BIND_OK_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS
    use systems, only: nsys, sysc, sys, sys_init, ok_system
    use gui_main, only: g
    use utils, only: iw_text, iw_tooltip, iw_radiobutton, iw_button, iw_setpos_bottomright,&
       iw_helpermark, iw_begintabitem, iw_inputtext, iw_inputint, iw_inputint3,&
       iw_inputfloat, iw_checkbox, iw_combo_simple, iw_field_combo, iw_arith_help,&
       iw_arith_help_button, iw_calcwidth
    use tools_io, only: string, nameguess, quoteword
    class(window), intent(inout), target :: w

    logical :: oksys, ok, doquit, disabled, resetsys, gotofile
    integer :: isys, i, j, ll, idx, iff, iaux, sourceopt
    character(len=:,kind=c_char), allocatable :: reason
    integer(c_int) :: tabflags
    type(ImVec2) :: szavail, szero
    real(c_float) :: combowidth
    logical(c_bool) :: is_selected
    logical :: ldum
    character(len=:,kind=c_char), allocatable, target :: str1, str2, loadstr

    integer, parameter :: itoken_file1 = 1
    integer, parameter :: itoken_file2 = 2
    integer, parameter :: itoken_file3 = 3

    ! FFT operations, in the order of the lf%ifftop combo
    character*4, parameter :: fftopname(0:11) = (/"gx  ","gy  ","gz  ","hxx ","hxy ",&
       "hxz ","hyy ","hyz ","hzz ","gmod","lap ","pot "/)
    ! TYPNUC values, in the order of the lf%itypnuc combo
    integer, parameter :: typnucval(0:3) = (/-3,-1,1,3/)

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
    sourceopt = 0

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
       ! check all systems and try to find one that is OK
       do i = 1, nsys
          oksys = ok_system(i,sys_init)
          if (oksys) then
             isys = i
             exit
          end if
       end do
    end if
    if (.not.oksys) then
       ! this dialog does not make sense anymore, close it and exit
       call w%end()
       return
    end if

    ! system combo
    call iw_text("System",highlight=.true.,alignframe=.true.)
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

    ! reset the per-system widget state when the system changes (or its
    ! species count does, e.g. from the geometry editor)
    resetsys = (w%lf%isys_last /= isys) .or. .not.allocated(w%lf%zpsp)
    if (.not.resetsys) resetsys = (size(w%lf%zpsp) /= sys(isys)%c%nspc)
    if (resetsys) then
       w%lf%file1 = ""
       w%lf%file1_format = 0
       w%lf%file2 = ""
       w%lf%file2_set = .false.
       w%lf%file3 = ""
       w%lf%file3_set = .false.
       w%lf%sourcefid = -1
       w%lf%sizeoffid = -1
       if (allocated(w%lf%zpsp)) deallocate(w%lf%zpsp)
       allocate(w%lf%zpsp(sys(isys)%c%nspc))
       w%lf%zpsp = -1
       w%lf%isys_last = isys
       call validate_expr()
    end if

    ! if a new file arrived from a file dialog, process it
    if (w%okfile_set) then
       w%errmsg = ""
       if (w%itoken == itoken_file1) then
          w%lf%file1 = w%okfile
          call intake_file1(w%okfile_format)
       elseif (w%itoken == itoken_file2) then
          w%lf%file2 = w%okfile
          w%lf%file2_set = .true.
       elseif (w%itoken == itoken_file3) then
          w%lf%file3 = w%okfile
          w%lf%file3_set = .true.
       end if
    end if
    w%okfile_set = .false.

    ! an external request (Load Field From File...) selects the From File
    ! tab and pops the file dialog, as if the File button had been pressed
    gotofile = w%lf%gotofile
    w%lf%gotofile = .false.
    if (gotofile) then
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfieldfile,&
          idparent=w%id,itoken=itoken_file1)
    end if

    ! the source tabs
    str1 = "##loadfieldtabbar" // c_null_char
    tabflags = ImGuiTabItemFlags_None
    if (gotofile .or. w%firstpass) tabflags = ImGuiTabItemFlags_SetSelected
    if (igBeginTabBar(c_loc(str1),ImGuiTabBarFlags_None)) then
       if (iw_begintabitem("From File##loadfieldtab_file",flags=tabflags)) then
          sourceopt = 0
          call draw_source_file()
          call igEndTabItem()
       end if
       if (iw_begintabitem("Expression##loadfieldtab_expr")) then
          sourceopt = 1
          call draw_source_expression()
          call igEndTabItem()
       end if
       if (iw_begintabitem("Promolecular##loadfieldtab_promol")) then
          sourceopt = 2
          call draw_source_promol()
          call igEndTabItem()
       end if
       if (iw_begintabitem("Grid Transform##loadfieldtab_transform")) then
          sourceopt = 3
          call draw_source_transform()
          call igEndTabItem()
       end if
       call igEndTabBar()
    end if

    ! common and advanced options
    call draw_common_options()
    if (sourceopt /= 3) &
       call draw_advanced_options()

    ! error message, if applicable
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    ! calculate whether we have enough info to continue to ok, and say
    ! why not
    reason = ok_reason()
    disabled = (len_trim(reason) > 0)
    if (disabled) call iw_text(reason,disabled=.true.,wrap=.true.)

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(8,2)
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG))
    ok = ok .or. iw_button("OK",disabled=disabled)
    if (ok .and..not.disabled) then
       w%errmsg = ""
       loadstr = build_load_string()
       if (len_trim(w%errmsg) == 0) then
          call sys(isys)%load_field_string(loadstr,.false.,iff,w%errmsg)
          if (len_trim(w%errmsg) == 0) then
             ! expand the field list so the new field is visible
             sysc(isys)%showfields = .true.
             doquit = .true.
          end if
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) &
       doquit = .true.

    ! quit the window
    if (doquit) &
       call w%end()

  contains
    !> Draw the contents of the from-a-file tab.
    subroutine draw_source_file()
      character(len=:), allocatable :: fmtkey, fmtdesc
      logical :: fmtgrid

      ! source file and format
      call iw_text("Source",highlight=.true.)
      if (iw_button("File",danger=(w%lf%file1_format==0))) then
         iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openfieldfile,&
            idparent=w%id,itoken=itoken_file1)
      end if
      call iw_tooltip("File from where the field is read",ttshown)
      if (len(w%lf%file1) > 0) then
         if (iw_button("Clear##clearfile1",sameline=.true.,danger=.true.)) then
            w%lf%file1 = ""
            call intake_file1()
         end if
         call iw_tooltip("Clear the file",ttshown)
      end if
      ! editable path: a typed or pasted path is accepted on enter
      if (iw_inputtext("##loadfieldpath",bufsize=4095,texta=w%lf%file1,width=45,&
         sameline=.true.,notlive=.true.)) &
         call intake_file1()
      if (len_trim(w%lf%file1) > 0) then
         call iw_tooltip(trim(w%lf%file1),ttshown)
      else
         call iw_tooltip("Path to the field file (type it or use the File button)",ttshown)
      end if

      if (w%lf%file1_format /= 0) then
         call field_format_string(abs(w%lf%file1_format),fmtkey,fmtdesc,fmtgrid)
         call iw_text("Format: ")
         call iw_text(fmtdesc,sameline=.true.)
         if (w%lf%file1_format < 0) &
            call iw_text(" [auto-detected]",sameline=.true.)

         ! propose default auxiliary files
         select case (abs(w%lf%file1_format))
         case(ifformat_wien)
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
            ! default GEOMETRY.OUT file
            if (.not.w%lf%file2_set) then
               ll = len_trim(sysc(isys)%seed%file)
               if (ll >= 4) then
                  if (sysc(isys)%seed%file(ll-3:ll) == ".OUT") &
                     w%lf%file2 = sysc(isys)%seed%file
               end if
               w%lf%file2_set = .true.
            end if
         case(ifformat_dftb)
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
         end select
      end if

      ! format-specific options and additional files
      select case (abs(w%lf%file1_format))
      case(ifformat_wien)
         call aux_file_row("File (.struct)","WIEN2k structure file (.struct) used to&
            & interpret the clmsum-style file",itoken_file2,w%lf%file2,w%lf%file2_set,&
            danger=.not.w%lf%file2_set,filter="WIEN2k struct (struct){.struct}")
      case(ifformat_elk)
         call aux_file_row("File (GEOMETRY.OUT)","GEOMETRY.OUT structure file used to&
            & interpret the STATE.OUT",itoken_file2,w%lf%file2,w%lf%file2_set,&
            danger=.not.w%lf%file2_set,filter="elk (OUT){.OUT}")
         call aux_file_row("File (extra .OUT)","Additional elk .OUT file with more&
            & scalar field data (optional)",itoken_file3,w%lf%file3,w%lf%file3_set,&
            filter="elk (OUT){.OUT}")
      case(ifformat_dftb)
         call aux_file_row("File (.bin)","eigenvec.bin file for reading the DFTB+&
            & wavefunction",itoken_file2,w%lf%file2,w%lf%file2_set,&
            danger=.not.w%lf%file2_set,filter="eigenvector file (bin){.bin}")
         call aux_file_row("File (.hsd)","hsd file for reading the DFTB+ wavefunction",&
            itoken_file3,w%lf%file3,w%lf%file3_set,danger=.not.w%lf%file3_set,&
            filter="hsd file (hsd){.hsd}")
      case(ifformat_pwc)
         call aux_file_row("File (.chk, up)","wannier90 checkpoint file for the spin-up&
            & (or only) channel (optional)",itoken_file2,w%lf%file2,w%lf%file2_set,&
            filter="wannier90 checkpoint (chk){.chk}")
         if (len(w%lf%file2) > 0) then
            call aux_file_row("File (.chk, down)","wannier90 checkpoint file for the&
               & spin-down channel (optional)",itoken_file3,w%lf%file3,w%lf%file3_set,&
               filter="wannier90 checkpoint (chk){.chk}")
         elseif (len(w%lf%file3) > 0) then
            ! the down chk makes no sense without the up chk
            w%lf%file3 = ""
            w%lf%file3_set = .false.
         end if

         ! pwc options
         call iw_combo_simple("Spin channel##loadfieldpwcspin","Both"//c_null_char//"Up"//c_null_char//&
            "Down"//c_null_char,w%lf%pwcspin)
         call iw_tooltip("Spin channel used to build the density",ttshown)
         ldum = iw_inputtext("k-points##loadfieldpwckpt",bufsize=255,texta=w%lf%pwckpt,width=25)
         call iw_tooltip("Space-separated list of k-point indices used to build the density (empty = all)",ttshown)
         ldum = iw_inputtext("Bands##loadfieldpwcband",bufsize=255,texta=w%lf%pwcband,width=25)
         call iw_tooltip("Space-separated list of band indices used to build the density (empty = all)",ttshown)
         ldum = iw_checkbox("Energy range (eV)##loadfieldpwcuserange",w%lf%pwcuserange)
         call iw_tooltip("Use only the states in this energy range to build the density",ttshown)
         if (w%lf%pwcuserange) then
            ldum = iw_inputfloat("Min.##loadfieldpwcemin",w%lf%pwcemin,width=10,sameline=.true.)
            ldum = iw_inputfloat("Max.##loadfieldpwcemax",w%lf%pwcemax,width=10,sameline=.true.)
         end if
      case(ifformat_vasp,ifformat_vaspnov)
         call iw_combo_simple("Block##loadfieldvaspblk","1: density"//c_null_char//&
            "2: spin / x-magnetization"//c_null_char//"3: y-magnetization"//c_null_char//&
            "4: z-magnetization"//c_null_char,w%lf%vaspblk)
         call iw_tooltip("Data block read from the VASP file",ttshown)
      case(ifformat_molden,ifformat_fchk)
         if (abs(w%lf%file1_format) == ifformat_molden) then
            call iw_combo_simple("Dialect##loadfieldmoldendialect","Auto-detect"//c_null_char//&
               "PSI4"//c_null_char//"ORCA"//c_null_char,w%lf%moldendialect)
            call iw_tooltip("Program that generated the molden file",ttshown)
         end if
         ldum = iw_checkbox("Read virtual orbitals##loadfieldreadvirtual",w%lf%readvirtual)
         call iw_tooltip("Read the unoccupied (virtual) molecular orbitals as well",ttshown)
      end select

    end subroutine draw_source_file

    !> Draw the contents of the expression tab.
    subroutine draw_source_expression()

      call iw_text("Expression",highlight=.true.,alignframe=.true.)
      call iw_helpermark("Arithmetic expression used to compute the new field, typically&
         & involving the other fields of this system (e.g. $1+$2)."//newline//&
         iw_arith_help//newline//"Click the Help button for more info.")
      call iw_arith_help_button("##loadfieldexprhelp",ttshown)
      if (iw_inputtext("##loadfieldexprtext",bufsize=2047,texta=w%lf%expr)) &
         call validate_expr()
      if (len_trim(w%lf%expr) > 0) then
         if (len_trim(w%lf%exprerr) > 0) then
            call iw_text(w%lf%exprerr,danger=.true.,wrap=.true.)
         else
            call iw_text("valid expression",disabled=.true.)
         end if
      end if

      call iw_text("Domain",highlight=.true.)
      ldum = iw_radiobutton("Automatic",int=w%lf%exprdom,intval=0_c_int)
      call iw_tooltip("New grid if all fields in the expression are grids of the same size;&
         & otherwise a ghost field evaluated on the fly",ttshown)
      ldum = iw_radiobutton("Ghost field",int=w%lf%exprdom,intval=1_c_int,sameline=.true.)
      call iw_tooltip("Evaluate the expression on the fly every time the field value is needed",ttshown)
      ldum = iw_radiobutton("Grid",int=w%lf%exprdom,intval=2_c_int,sameline=.true.)
      call iw_tooltip("Evaluate the expression once on a grid",ttshown)
      if (w%lf%exprdom == 2) &
         call draw_grid_size()

    end subroutine draw_source_expression

    !> Draw the contents of the promolecular tab.
    subroutine draw_source_promol()

      call iw_text("Density",highlight=.true.)
      ldum = iw_radiobutton("Promolecular",int=w%lf%promopt,intval=0_c_int)
      call iw_tooltip("Sum of in-vacuo atomic densities, evaluated analytically",ttshown)
      ldum = iw_radiobutton("Promolecular on a grid",int=w%lf%promopt,intval=1_c_int,sameline=.true.)
      call iw_tooltip("Sum of in-vacuo atomic densities, evaluated once on a grid",ttshown)
      ldum = iw_radiobutton("Core on a grid",int=w%lf%promopt,intval=2_c_int,sameline=.true.)
      call iw_tooltip("Sum of atomic core densities on a grid; requires activating the&
         & core augmentation (Advanced options)",ttshown)

      if (w%lf%promopt >= 1) &
         call draw_grid_size()
      if (w%lf%promopt == 2 .and. .not.goodzpsp()) &
         call iw_text("Core densities require activating the core augmentation&
            & (Advanced options)",danger=.true.,wrap=.true.)

    end subroutine draw_source_promol

    !> Draw the contents of the grid transform tab. The resulting field
    !> is derived from existing fields of this system. Note the AS
    !> FFT/RESAMPLE/CLM code paths ignore the trailing LOAD options
    !> (except the name), so the common/advanced options are not shown
    !> for this tab.
    subroutine draw_source_transform()

      call iw_combo_simple("Transform##loadfieldgtransf","Fourier-transform derivative"//c_null_char//&
         "Resample grid"//c_null_char,w%lf%gtransf)
      call iw_tooltip("Operation that generates the new field",ttshown)

      select case (w%lf%gtransf)
      case (0)
         ! fft
         call iw_combo_simple("Operation##loadfieldfftop","Gradient x-component"//c_null_char//&
            "Gradient y-component"//c_null_char//"Gradient z-component"//c_null_char//&
            "Hessian xx"//c_null_char//"Hessian xy"//c_null_char//"Hessian xz"//c_null_char//&
            "Hessian yy"//c_null_char//"Hessian yz"//c_null_char//"Hessian zz"//c_null_char//&
            "Gradient norm"//c_null_char//"Laplacian"//c_null_char//&
            "Potential (Poisson)"//c_null_char,w%lf%ifftop)
         call iw_tooltip("Derivative of the source field calculated via fast Fourier transform",ttshown)
         call field_row("Source field","##loadfieldfftsource",w%lf%sourcefid,&
            "Grid field the operation applies to",onlygrid=.true.)
         if (w%lf%ifftop == 11) then
            ldum = iw_checkbox("Potential in Rydberg units##loadfieldfftry",w%lf%fftry)
            call iw_tooltip("Return the potential in Rydberg instead of Hartree units",ttshown)
         end if
      case (1)
         ! resample
         call field_row("Source field","##loadfieldresamplesource",w%lf%sourcefid,&
            "Grid field to be resampled",onlygrid=.true.)
         call draw_grid_size()
      end select

    end subroutine draw_source_transform

    !> Draw a labeled row with a field selector combo and a tooltip.
    subroutine field_row(label,strid,fid,tip,onlygrid)
      character(len=*,kind=c_char), intent(in) :: label, strid, tip
      integer, intent(inout) :: fid
      logical, intent(in) :: onlygrid

      call iw_text(label,alignframe=.true.)
      call igSameLine(0._c_float,-1._c_float)
      ldum = iw_field_combo(strid,isys,fid,width=iw_calcwidth(30,1),onlygrid=onlygrid)
      call iw_tooltip(tip,ttshown)

    end subroutine field_row

    !> Draw a file button + Clear button + file name row for an
    !> auxiliary file. Opens a modal file dialog with token itok whose
    !> result lands in file. Clearing resets fileset (if present) so a
    !> default is re-proposed.
    subroutine aux_file_row(label,tip,itok,file,fileset,danger,filter)
      character(len=*,kind=c_char), intent(in) :: label, tip
      integer, intent(in) :: itok
      character(len=:,kind=c_char), allocatable, intent(inout) :: file
      logical, intent(inout), optional :: fileset
      logical, intent(in), optional :: danger
      character(len=*,kind=c_char), intent(in), optional :: filter

      logical :: danger_

      danger_ = .false.
      if (present(danger)) danger_ = danger
      if (iw_button(label,danger=danger_)) then
         iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_openonefilemodal,&
            idparent=w%id,itoken=itok,dialog_filter=filter)
      end if
      call iw_tooltip(tip,ttshown)
      if (len(file) > 0) then
         if (iw_button("Clear##clearfile" // string(itok),sameline=.true.,danger=.true.)) then
            file = ""
            if (present(fileset)) fileset = .false.
         end if
         call iw_tooltip("Clear the file",ttshown)
      end if
      ! long paths show their tail; the full path goes in a tooltip
      call iw_text(pathtail(file),sameline=.true.)
      if (len_trim(file) > 0) &
         call iw_tooltip(trim(file),ttshown)

    end subroutine aux_file_row

    !> Tail of a long path, prefixed with an ellipsis, for display.
    function pathtail(file) result(str)
      character(len=*,kind=c_char), intent(in) :: file
      character(len=:,kind=c_char), allocatable :: str

      integer, parameter :: maxlen = 40

      if (len_trim(file) <= maxlen) then
         str = trim(file)
      else
         str = "..." // file(len_trim(file)-maxlen+1:len_trim(file))
      end if

    end function pathtail

    !> Process a new main file (from the file dialog, the path input, or
    !> the Clear button): set the format (iforce, if given and nonzero,
    !> else by extension) and reset the auxiliary files.
    subroutine intake_file1(iforce)
      integer, intent(in), optional :: iforce

      w%errmsg = ""
      w%lf%file1_format = 0
      if (present(iforce)) w%lf%file1_format = iforce
      if (w%lf%file1_format == 0 .and. len_trim(w%lf%file1) > 0) then
         w%lf%file1_format = -field_detect_format(file=w%lf%file1)
         if (abs(w%lf%file1_format) == ifformat_pi) w%lf%file1_format = 0 ! aiPI deactivated for now
         if (w%lf%file1_format == 0) w%errmsg = "Unknown field file extension"
      end if
      w%lf%file2 = ""
      w%lf%file2_set = .false.
      w%lf%file3 = ""
      w%lf%file3_set = .false.

    end subroutine intake_file1

    !> Probe-evaluate the expression against the current system and
    !> store the error message (empty = valid) in the state.
    subroutine validate_expr()
      character(len=:), allocatable :: errm

      w%lf%exprerr = ""
      if (len_trim(w%lf%expr) == 0) return
      call sys(isys)%check_expression(w%lf%expr,errm)
      if (len_trim(errm) > 0) w%lf%exprerr = trim(errm)

    end subroutine validate_expr

    !> Draw the shared grid-size widget (explicit size or same size as
    !> an existing grid field).
    subroutine draw_grid_size()

      ldum = iw_radiobutton("Grid size",int=w%lf%ngridopt,intval=0_c_int)
      call iw_tooltip("Number of grid points along each lattice vector",ttshown)
      if (w%lf%ngridopt == 0) &
         ldum = iw_inputint3("##loadfieldngrid",w%lf%ngrid,width=4*3,sameline=.true.)
      ldum = iw_radiobutton("Same size as field",int=w%lf%ngridopt,intval=1_c_int)
      call iw_tooltip("Use the same grid size as an existing grid field",ttshown)
      if (w%lf%ngridopt == 1) then
         call igSameLine(0._c_float,-1._c_float)
         ldum = iw_field_combo("##loadfieldsizeof",isys,w%lf%sizeoffid,width=iw_calcwidth(25,1),&
            onlygrid=.true.)
      end if

    end subroutine draw_grid_size

    !> Draw the common options (name, grid interpolation, normalization).
    subroutine draw_common_options()

      call iw_text("Options",highlight=.true.)

      ! field name
      ldum = iw_inputtext("Name##loadfieldname",bufsize=255,texta=w%lf%name,width=25,&
         flags=ImGuiInputTextFlags_CharsNoBlank)
      call iw_tooltip("Identifier for the new field, usable in expressions and keywords (optional)",ttshown)

      ! grid options, only if the load produces a grid and honors them
      if (producesgrid()) then
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

         ! normalization
         ldum = iw_checkbox("Normalize##loadfielddonorm",w%lf%donorm)
         call iw_tooltip("Rescale the grid values so its cell integral equals this value",ttshown)
         if (w%lf%donorm) then
            ldum = iw_inputfloat("##loadfieldnormval",w%lf%normval,width=10,sameline=.true.)
            call iw_tooltip("Target value of the cell integral of the grid",ttshown)
         end if
      end if

    end subroutine draw_common_options

    !> Draw the advanced options inside a collapsed tree node.
    subroutine draw_advanced_options()
      character(len=:,kind=c_char), allocatable, target :: strad

      strad = "Advanced options##loadfieldadvanced" // c_null_char
      if (igTreeNodeEx_Str(c_loc(strad),ImGuiTreeNodeFlags_None)) then
         ! typnuc
         call iw_combo_simple("Nuclear critical point signature##loadfieldtypnuc",&
            "-3 (maximum)"//c_null_char//"-1"//c_null_char//"+1"//c_null_char//&
            "+3 (minimum)"//c_null_char,w%lf%itypnuc)
         call iw_tooltip("Signature of the critical points at the nuclear positions&
            & (-3 for fields with maxima at the nuclei, like the density)",ttshown)

         ! numerical/analytical derivatives
         call iw_combo_simple("Derivatives##loadfieldnumana","Default"//c_null_char//&
            "Numerical"//c_null_char//"Analytical"//c_null_char,w%lf%inumana)
         call iw_tooltip("Method for calculating the field derivatives",ttshown)

         ! dftb: exact/approximate
         if (fileformat() == ifformat_dftb) then
            call iw_combo_simple("Sum evaluation##loadfieldiexact","Default"//c_null_char//&
               "Exact"//c_null_char//"Approximate"//c_null_char,w%lf%iexact)
            call iw_tooltip("Evaluate the lattice sums exactly or approximately",ttshown)
         end if

         ! wien2k normalization
         if (fileformat() == ifformat_wien) then
            call iw_combo_simple("Normalization##loadfieldwiennorm","Default"//c_null_char//&
               "Density (RHONORM)"//c_null_char//"Potential (VNORM)"//c_null_char,w%lf%iwiennorm)
            call iw_tooltip("Normalization convention of the WIEN2k file&
               & (density-style or potential-style)",ttshown)
         end if

         ! muffin-tin discontinuity test
         if (fileformat() == ifformat_wien .or. fileformat() == ifformat_elk) then
            ldum = iw_checkbox("Test the muffin-tin discontinuity##loadfieldtestrmt",w%lf%testrmt)
            call iw_tooltip("Check the size of the discontinuity of the field&
               & at the muffin-tin surfaces (takes time)",ttshown)
         end if

         ! smoothrho options
         if (producesgrid() .and. w%lf%iginterp == 4) then
            ldum = iw_inputint("Environment points (NENV)##loadfieldnenv",w%lf%nenv,width=6)
            call iw_tooltip("Number of environment grid points for the smoothrho fit (0 = default)",ttshown)
            ldum = iw_inputfloat("Distance factor (FDMAX)##loadfieldfdmax",w%lf%fdmax,width=8)
            call iw_tooltip("Maximum distance factor for the smoothrho environment (0 = default)",ttshown)
         end if

         ! core augmentation and pseudopotential charges
         ldum = iw_checkbox("Activate core augmentation##loadfielddocore",w%lf%docore)
         call iw_tooltip("Augment the field with core density contributions, built from&
            & the pseudopotential charges (ZPSP) below: the number of valence electrons&
            & of each species (leave at -1 to skip a species)",ttshown)
         if (w%lf%docore) then
            do i = 1, sys(isys)%c%nspc
               ldum = iw_inputint(trim(sys(isys)%c%spc(i)%name) // "##loadfieldzpsp" // string(i),&
                  w%lf%zpsp(i),width=5)
            end do
         end if

         call igTreePop()
      end if

    end subroutine draw_advanced_options

    !> Format of the main file if the file tab is active, 0 otherwise.
    function fileformat() result(jff)
      integer :: jff

      jff = 0
      if (sourceopt == 0) jff = abs(w%lf%file1_format)

    end function fileformat

    !> Whether the current load produces a grid AND its load path
    !> honors the grid options (the grid-transform paths ignore them).
    function producesgrid() result(isg)
      logical :: isg

      character(len=:), allocatable :: fmtkey, fmtdesc

      isg = .false.
      select case (sourceopt)
      case (0)
         call field_format_string(fileformat(),fmtkey,fmtdesc,isg)
      case (1)
         isg = (w%lf%exprdom == 2)
      case (2)
         isg = (w%lf%promopt >= 1)
      end select

    end function producesgrid

    !> Why the OK button must be disabled (not enough information to
    !> build the LOAD command); empty = OK is enabled.
    function ok_reason() result(reason)
      character(len=:,kind=c_char), allocatable :: reason

      integer :: jff

      character(len=*,kind=c_char), parameter :: zpspmsg =&
         "Set at least one ZPSP charge (Advanced options) to enable OK"

      reason = ""
      select case (sourceopt)
      case (0)
         ! from a file
         if (len(w%lf%file1) == 0) then
            reason = "Choose a source file to enable OK"
         elseif (w%lf%file1_format == 0) then
            reason = "The format of the source file is unknown"
         else
            jff = fileformat()
            if ((jff == ifformat_wien .or. jff == ifformat_elk .or. jff == ifformat_dftb) .and.&
               len(w%lf%file2) == 0 .or. jff == ifformat_dftb .and. len(w%lf%file3) == 0) &
               reason = "Choose the auxiliary file(s) to enable OK"
         end if
      case (1)
         ! expression; a failed probe blocks OK only for ghost fields,
         ! where the load would be refused anyway (a grid can still be
         ! valid if the expression is only singular at the probe point)
         if (len_trim(w%lf%expr) == 0) then
            reason = "Enter an expression to enable OK"
         elseif (w%lf%exprdom == 1 .and. len_trim(w%lf%exprerr) > 0) then
            reason = "The expression is invalid"
         elseif (w%lf%exprdom == 2) then
            reason = sizereason()
         end if
      case (2)
         ! promolecular / core
         if (w%lf%promopt >= 1) reason = sizereason()
         if (len_trim(reason) == 0 .and. w%lf%promopt == 2 .and. .not.goodzpsp()) &
            reason = zpspmsg
      case (3)
         ! grid transform
         if (.not.sys(isys)%goodfield(w%lf%sourcefid,type=type_grid)) then
            reason = "Select a source field to enable OK"
         elseif (w%lf%gtransf == 1) then
            reason = sizereason()
         end if
      end select

      ! an activated core augmentation needs at least one ZPSP charge (the
      ! grid-transform tab does not use the advanced options)
      if (len_trim(reason) == 0 .and. sourceopt /= 3 .and. w%lf%docore .and. .not.goodzpsp()) &
         reason = zpspmsg

    end function ok_reason

    !> Why the shared grid-size widget holds an invalid size; empty =
    !> the size is valid.
    function sizereason() result(reason)
      character(len=:,kind=c_char), allocatable :: reason

      reason = ""
      if (w%lf%ngridopt == 0) then
         if (any(w%lf%ngrid <= 0)) reason = "The grid size must be positive"
      elseif (.not.sys(isys)%goodfield(w%lf%sizeoffid,type=type_grid)) then
         reason = "Select a grid field for the size to enable OK"
      end if

    end function sizereason

    !> Whether the core augmentation is active with usable ZPSP charges.
    function goodzpsp() result(good)
      logical :: good

      good = w%lf%docore .and. allocated(w%lf%zpsp)
      if (good) good = any(w%lf%zpsp >= 0)

    end function goodzpsp

    !> Build the LOAD command from the window state (without the LOAD
    !> keyword). Sets w%errmsg and returns an empty string on error.
    function build_load_string() result(lstr)
      character(len=:,kind=c_char), allocatable :: lstr
      character(len=:,kind=c_char), allocatable :: qc
      character(len=:), allocatable :: fmtkey, fmtdesc
      logical :: fmtgrid
      integer :: k, jff

      lstr = ""
      select case (sourceopt)
      case (0)
         ! from a file; paths with blanks are quoted (getword strips a
         ! leading-quote-delimited token)
         jff = fileformat()
         call field_format_string(jff,fmtkey,fmtdesc,fmtgrid)
         lstr = fmtkey // " " // qfile(w%lf%file1) // " " // qfile(w%lf%file2) //&
            " " // qfile(w%lf%file3)

         ! per-format options
         if (jff == ifformat_vasp .or. jff == ifformat_vaspnov) then
            if (w%lf%vaspblk > 0) lstr = lstr // " " // string(w%lf%vaspblk+1)
         elseif (jff == ifformat_pwc) then
            if (w%lf%pwcspin > 0) lstr = lstr // " spin " // string(w%lf%pwcspin)
            if (len_trim(w%lf%pwckpt) > 0) lstr = lstr // " kpt " // trim(w%lf%pwckpt)
            if (len_trim(w%lf%pwcband) > 0) lstr = lstr // " band " // trim(w%lf%pwcband)
            if (w%lf%pwcuserange) &
               lstr = lstr // " erange " // string(w%lf%pwcemin,'f',decimal=6) // " " //&
               string(w%lf%pwcemax,'f',decimal=6)
         elseif (jff == ifformat_molden) then
            if (w%lf%moldendialect == 1) then
               lstr = lstr // " psi4"
            elseif (w%lf%moldendialect == 2) then
               lstr = lstr // " orca"
            end if
         elseif (jff == ifformat_dftb) then
            if (w%lf%iexact == 1) then
               lstr = lstr // " exact"
            elseif (w%lf%iexact == 2) then
               lstr = lstr // " approximate"
            end if
         elseif (jff == ifformat_wien) then
            if (w%lf%iwiennorm == 1) then
               lstr = lstr // " rhonorm"
            elseif (w%lf%iwiennorm == 2) then
               lstr = lstr // " vnorm"
            end if
         end if
         if ((jff == ifformat_molden .or. jff == ifformat_fchk) .and. w%lf%readvirtual) &
            lstr = lstr // " readvirtual"
         if (jff == ifformat_wien .or. jff == ifformat_elk) then
            if (.not.w%lf%testrmt) lstr = lstr // " notestmt"
         end if
      case (1)
         ! expression: quote it with whichever quote character it does not contain
         if (index(w%lf%expr,'"') == 0) then
            qc = '"'
         elseif (index(w%lf%expr,"'") == 0) then
            qc = "'"
         else
            w%errmsg = "The expression cannot contain both single and double quotes"
            return
         end if
         lstr = "as " // qc // trim(w%lf%expr) // qc
         if (w%lf%exprdom == 1) then
            lstr = lstr // " ghost"
         elseif (w%lf%exprdom == 2) then
            lstr = lstr // " " // sizestring()
         end if
      case (2)
         ! promolecular / core
         if (w%lf%promopt == 0) then
            lstr = "promolecular"
         elseif (w%lf%promopt == 1) then
            lstr = "as promolecular " // sizestring()
         else
            lstr = "as core " // sizestring()
         end if
      case (3)
         ! grid transform
         select case (w%lf%gtransf)
         case (0)
            lstr = "as fft " // trim(fftopname(w%lf%ifftop)) // " " // string(w%lf%sourcefid)
            if (w%lf%ifftop == 11 .and. w%lf%fftry) lstr = lstr // " ry"
         case (1)
            lstr = "as resample " // string(w%lf%sourcefid) // " " // sizestring()
         end select
      end select

      ! field name (works in all load paths)
      if (len_trim(w%lf%name) > 0) lstr = lstr // " name " // trim(w%lf%name)

      ! the grid-transform paths ignore the rest of the options
      if (sourceopt == 3) return

      ! grid options
      if (producesgrid()) then
         select case (w%lf%iginterp)
         case (0)
            lstr = lstr // " nearest"
         case (1)
            lstr = lstr // " trilinear"
         case (2)
            lstr = lstr // " trispline"
         case (3)
            lstr = lstr // " tricubic"
         case (4)
            lstr = lstr // " smoothrho"
            if (w%lf%nenv > 0) lstr = lstr // " nenv " // string(w%lf%nenv)
            if (w%lf%fdmax > 0._c_float) &
               lstr = lstr // " fdmax " // string(w%lf%fdmax,'f',decimal=4)
         end select
         if (w%lf%donorm) lstr = lstr // " normalize " // string(w%lf%normval,'f',decimal=6)
      end if

      ! advanced options
      if (w%lf%itypnuc /= 0) lstr = lstr // " typnuc " // string(typnucval(w%lf%itypnuc))
      if (w%lf%inumana == 1) then
         lstr = lstr // " numerical"
      elseif (w%lf%inumana == 2) then
         lstr = lstr // " analytical"
      end if
      if (goodzpsp()) then
         ! the LOAD parser only accepts short element symbols in ZPSP, not
         ! species labels; emit one pair per element (first species wins)
         lstr = lstr // " zpsp"
         do k = 1, sys(isys)%c%nspc
            if (w%lf%zpsp(k) < 0) cycle
            if (any(sys(isys)%c%spc(1:k-1)%z == sys(isys)%c%spc(k)%z .and. w%lf%zpsp(1:k-1) >= 0)) cycle
            lstr = lstr // " " // trim(nameguess(sys(isys)%c%spc(k)%z,.true.)) // " " //&
               string(w%lf%zpsp(k))
         end do
      else
         ! core augmentation is deactivated by default
         lstr = lstr // " nocore"
      end if

    end function build_load_string

    !> Quote a file path for the LOAD command if it contains blanks.
    !> Sets w%errmsg if the path cannot be represented.
    function qfile(file) result(qf)
      character(len=*,kind=c_char), intent(in) :: file
      character(len=:,kind=c_char), allocatable :: qf

      qf = quoteword(file)
      if (index(qf," ") > 0 .and. qf(1:1) /= '"' .and. qf(1:1) /= "'") &
         w%errmsg = "A file path cannot contain blanks and both quote characters: " // qf

    end function qfile

    !> String for the shared grid-size widget (explicit numbers or a
    !> SIZEOF clause).
    function sizestring() result(sstr)
      character(len=:,kind=c_char), allocatable :: sstr

      if (w%lf%ngridopt == 0) then
         sstr = string(w%lf%ngrid(1)) // " " // string(w%lf%ngrid(2)) // " " //&
            string(w%lf%ngrid(3))
      else
         sstr = "sizeof " // string(w%lf%sizeoffid)
      end if

    end function sizestring

    ! initialize the state for this window; window_init already reset
    ! w%lf, but allocatable-string components cannot have default
    ! initializers, so allocate them here
    subroutine init_state()
      w%lf%file1 = ""
      w%lf%file2 = ""
      w%lf%file3 = ""
      w%lf%expr = ""
      w%lf%exprerr = ""
      w%lf%name = ""
      w%lf%pwckpt = ""
      w%lf%pwcband = ""
    end subroutine init_state
  end subroutine draw_load_field

end submodule load_field
