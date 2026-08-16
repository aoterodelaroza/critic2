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

! Routines for the save as... window.
submodule (windows) saveas
  use interfaces_cimgui
  use param, only: isformat_w_xyz, isformat_w_gjf, isformat_w_cml, isformat_w_obj,&
     isformat_w_ply, isformat_w_off, isformat_w_gaussian_periodic, isformat_w_qein,&
     isformat_w_aimsin, isformat_w_vasp, isformat_w_abinit, isformat_w_elk,&
     isformat_w_tessel, isformat_w_critic, isformat_w_cif, isformat_w_crystal,&
     isformat_w_shelx, isformat_w_octave, isformat_w_dcpdb, isformat_w_gulp,&
     isformat_w_lammps, isformat_w_siesta_fdf, isformat_w_siesta_struct,&
     isformat_w_dftbp_hsd, isformat_w_dftbp_gen, isformat_w_pyscf,&
     isformat_w_tinkerfrac, isformat_w_pdb, isformat_w_castepcell, isformat_w_alamode,&
     isformat_w_unknown, isformat_w_max
  implicit none

  ! format names, indexed by the isformat_w constants
  character(len=32), parameter :: fmtnames(isformat_w_max) = (/&
     "xyz file (.xyz)                 ",& ! isformat_w_xyz
     "Gaussian input (.gjf)           ",& ! isformat_w_gjf
     "CML file (.cml)                 ",& ! isformat_w_cml
     "Wavefront obj (.obj)            ",& ! isformat_w_obj
     "PLY file (.ply)                 ",& ! isformat_w_ply
     "OFF file (.off)                 ",& ! isformat_w_off
     "Gaussian input, periodic (.gau) ",& ! isformat_w_gaussian_periodic
     "Quantum ESPRESSO input (.pwi)   ",& ! isformat_w_qein
     "FHIaims input (.in)             ",& ! isformat_w_aimsin
     "VASP POSCAR (.POSCAR)           ",& ! isformat_w_vasp
     "abinit input (.abin)            ",& ! isformat_w_abinit
     "elk GEOMETRY.OUT (.elk)         ",& ! isformat_w_elk
     "tessel input (.tess)            ",& ! isformat_w_tessel
     "critic2 input (.cri)            ",& ! isformat_w_critic
     "CIF file (.cif)                 ",& ! isformat_w_cif
     "CRYSTAL input (.d12)            ",& ! isformat_w_crystal
     "SHELX res (.res)                ",& ! isformat_w_shelx
     "escher/octave (.m)              ",& ! isformat_w_octave
     "dcp database (.db)              ",& ! isformat_w_dcpdb
     "GULP input (.gin)               ",& ! isformat_w_gulp
     "LAMMPS data (.lammps)           ",& ! isformat_w_lammps
     "SIESTA fdf (.fdf)               ",& ! isformat_w_siesta_fdf
     "SIESTA STRUCT_IN (.struct_in)   ",& ! isformat_w_siesta_struct
     "DFTB+ hsd (.hsd)                ",& ! isformat_w_dftbp_hsd
     "DFTB+ gen (.gen)                ",& ! isformat_w_dftbp_gen
     "pyscf script (.pyscf)           ",& ! isformat_w_pyscf
     "TINKER frac (.frac)             ",& ! isformat_w_tinkerfrac
     "pdb file (.pdb)                 ",& ! isformat_w_pdb
     "CASTEP cell (.cell)             ",& ! isformat_w_castepcell
     "alamode input (.alm.in)         "/) ! isformat_w_alamode

  ! canonical extension for each format, indexed by the isformat_w constants
  character(len=9), parameter :: fmtext(isformat_w_max) = (/&
     "xyz      ","gjf      ","cml      ","obj      ","ply      ",&
     "off      ","gau      ","pwi      ","in       ","POSCAR   ",&
     "abin     ","elk      ","tess     ","cri      ","cif      ",&
     "d12      ","res      ","m        ","db       ","gin      ",&
     "lammps   ","fdf      ","struct_in","hsd      ","gen      ",&
     "pyscf    ","frac     ","pdb      ","cell     ","alm.in   "/)

  ! order of the formats in the combo
  integer, parameter :: fmtperm(isformat_w_max) = (/isformat_w_aimsin,isformat_w_cif,&
     isformat_w_qein,isformat_w_vasp,isformat_w_castepcell,isformat_w_shelx,&
     isformat_w_crystal,isformat_w_abinit,isformat_w_elk,isformat_w_gjf,&
     isformat_w_gaussian_periodic,isformat_w_xyz,isformat_w_cml,isformat_w_alamode,&
     isformat_w_critic,isformat_w_octave,isformat_w_dcpdb,isformat_w_gulp,&
     isformat_w_lammps,isformat_w_siesta_fdf,isformat_w_siesta_struct,&
     isformat_w_dftbp_hsd,isformat_w_dftbp_gen,isformat_w_pyscf,isformat_w_tinkerfrac,&
     isformat_w_tessel,isformat_w_pdb,isformat_w_obj,isformat_w_ply,isformat_w_off/)

  ! options for the format combo (built once on first use)
  character(kind=c_char,len=:), allocatable, target :: combostr

  ! cached existence check for the FHIaims companion control file
  character(len=:), allocatable :: lastcontrol
  logical :: lastcontrol_exists = .false.

  ! settings remembered from the last successful save (this session).
  ! The 50 rklength default is a deliberate GUI choice (the CLI writers
  ! default to 40 when no rklength is given).
  real(c_float) :: lastrk = 50._c_float
  logical :: lastnosym = .false.
  logical :: lastcartesian = .false.

contains

  !> Draw the save as window
  module subroutine draw_saveas(w)
    use windows, only: wintype_dialog, wpurp_dialog_savefile
    use systems, only: sys, sysc, sys_init, ok_system
    use crystalmod, only: struct_detect_write_format
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip, iw_checkbox,&
       iw_combo_simple, iw_dragfloat_realc, iw_close_event, iw_setpos_bottomright,&
       iw_inputtext, file_name_root
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG
    use tools_io, only: string, uout
    class(window), intent(inout), target :: w

    logical :: doquit, ok, okvalid, syschanged, changed, ismol, userk
    integer :: i, isys, iview, iaux, ifmt, idetect
    logical(c_bool) :: ldum

    logical, save :: ttshown = .false. ! tooltip flag

    ! build the format combo options on first use
    if (.not.allocated(combostr)) then
       combostr = "Auto-detect" // c_null_char
       do i = 1, isformat_w_max
          combostr = combostr // trim(fmtnames(fmtperm(i))) // c_null_char
       end do
    end if

    ! this window writes the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    if (.not.doquit) doquit = .not.ok_system(isys,sys_init)

    ! initialize state
    if (w%firstpass) then
       w%okfile = okfile_default(isys,"structure","in")
       w%saveas_format = 0
       w%errmsg = ""
       w%saveas_rk = lastrk
       w%saveas_nosym = lastnosym
       w%saveas_cartesian = lastcartesian
    end if

    ! the system being written
    if (.not.doquit) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)
    end if

    ! file name: editable field plus browse button
    call iw_text("File Name",highlight=.true.)
    ldum = iw_inputtext("##saveasfile",bufsize=1023,texta=w%okfile,width=36)
    call iw_tooltip("File the structure will be written to (with auto-detect,&
       & the extension selects the format)",ttshown)
    if (iw_button("Browse...",sameline=.true.)) &
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savefile,idparent=w%id,orraise=-1)
    call iw_tooltip("Choose the file with a file browser",ttshown)
    call w%okfile_warn_overwrite()

    ! format combo
    call igPushItemWidth(iw_calcwidth(45,1))
    call iw_combo_simple("Format##saveasformat",combostr,w%saveas_format,changed=changed)
    call igPopItemWidth()
    call iw_tooltip("File format for the written structure. If auto-detect,&
       & the format is chosen based on the file extension",ttshown)

    ! if the user set a format, change the file extension to match it
    if (changed .and. w%saveas_format > 0) &
       w%okfile = file_name_root(w%okfile) // "." // trim(fmtext(fmtperm(w%saveas_format)))

    ! the write format resulting from the combo selection; the
    ! auto-detected format is FHIaims input if the extension is unknown
    call struct_detect_write_format(w%okfile,idetect)
    if (w%saveas_format == 0) then
       ifmt = idetect
       if (ifmt == isformat_w_unknown) ifmt = isformat_w_aimsin
       ! show the detected format below the combo
       call iw_text("Detected:",highlight=.true.)
       call iw_text(trim(fmtnames(ifmt)),sameline=.true.)
    else
       ifmt = fmtperm(w%saveas_format)
       ! warn if the file extension does not match the selected format
       if (idetect /= ifmt) &
          call iw_text("Extension does not match the selected format",danger=.true.)
    end if

    ! whether the write will produce a molecule or a crystal: some
    ! formats are inherently molecular, some follow the system, the
    ! rest always write a cell
    ismol = .false.
    if (.not.doquit) then
       if (ifmt == isformat_w_xyz .or. ifmt == isformat_w_gjf .or. ifmt == isformat_w_cml .or.&
          ifmt == isformat_w_obj .or. ifmt == isformat_w_ply .or. ifmt == isformat_w_off) then
          ismol = .true.
       elseif (ifmt == isformat_w_aimsin .or. ifmt == isformat_w_pdb .or. ifmt == isformat_w_pyscf .or.&
          ifmt == isformat_w_dftbp_gen .or. ifmt == isformat_w_dftbp_hsd) then
          ismol = sys(isys)%c%ismolecule
       end if
       call iw_text("Format:",highlight=.true.)
       call iw_text(trim(merge("molecule","crystal ",ismol)),sameline=.true.)
    end if

    ! the k-point grid (rklength) applies to espresso, CASTEP, and to
    ! FHIaims crystals (for a molecule, FHIaims would write a
    ! meaningless k-point grid to the companion control file)
    userk = (ifmt == isformat_w_qein .or. ifmt == isformat_w_castepcell .or.&
       (ifmt == isformat_w_aimsin .and. .not.ismol))

    ! options, only for the formats that support them
    if (userk .or. ifmt == isformat_w_cif .or. ifmt == isformat_w_crystal .or.&
       ifmt == isformat_w_shelx) &
       call iw_text("Options",highlight=.true.)

    ! k-point grid (rklength)
    if (userk) then
       ldum = iw_dragfloat_realc("RKlength##saveasrk",x1=w%saveas_rk,speed=1._c_float,&
          min=1._c_float,max=200._c_float,decimal=1,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Length parameter for calculating the k-point grid",ttshown)
       if (ifmt == isformat_w_aimsin) then
          call iw_text("K-point grid written to a companion _control file")
          ! warn if the companion file exists (check cached on okfile change)
          if (.not.allocated(lastcontrol)) lastcontrol = ""
          if (trim(w%okfile) // "_control" /= lastcontrol) then
             lastcontrol = trim(w%okfile) // "_control"
             inquire(file=lastcontrol,exist=lastcontrol_exists)
          end if
          if (lastcontrol_exists) &
             call iw_text("The _control file exists and will be overwritten",danger=.true.)
       end if
    end if

    ! do not use symmetry
    if (ifmt == isformat_w_cif .or. ifmt == isformat_w_crystal .or.&
       ifmt == isformat_w_shelx) then
       ldum = iw_checkbox("Do not use symmetry##saveasnosym",w%saveas_nosym)
       call iw_tooltip("Write the structure without using the space group symmetry (P1)",ttshown)
    end if

    ! Cartesian coordinates (crystals only; ismol is .false. when doquit)
    if (ifmt == isformat_w_aimsin .and. .not.ismol .and. .not.doquit) then
       ldum = iw_checkbox("Cartesian coordinates##saveascartesian",w%saveas_cartesian)
       call iw_tooltip("Write Cartesian atomic coordinates (atom) instead of&
          & fractional coordinates (atom_frac)",ttshown)
    end if

    ! maybe the error message
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(8,2)

    ! final buttons: OK
    okvalid = (len_trim(w%okfile) > 0) .and. .not.doquit
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) .and. okvalid
    ok = ok .or. iw_button("OK",disabled=.not.okvalid)
    if (ok) then
       ! write the structure; pass the rklength only if it applies
       if (userk) then
          call sys(isys)%c%write_any_file(w%okfile,w%errmsg,iwformat=ifmt,&
             rklength=real(w%saveas_rk,8),nosym=w%saveas_nosym,cartesian=w%saveas_cartesian)
       else
          call sys(isys)%c%write_any_file(w%okfile,w%errmsg,iwformat=ifmt,&
             nosym=w%saveas_nosym,cartesian=w%saveas_cartesian)
       end if

       ! on success, remember the settings, report, and quit
       if (len_trim(w%errmsg) == 0) then
          lastrk = w%saveas_rk
          lastnosym = w%saveas_nosym
          lastcartesian = w%saveas_cartesian
          call okfile_save_dir(w%okfile)
          write (uout,'("Saved structure file: ",A)') trim(w%okfile)
          doquit = .true.
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_saveas

end submodule saveas
