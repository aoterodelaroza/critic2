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

! Routines for the extract cluster as molecule window.
submodule (windows) extract
  use interfaces_cimgui
  implicit none

  ! region shapes
  integer(c_int), parameter :: er_sphere = 0
  integer(c_int), parameter :: er_cube = 1
  integer(c_int), parameter :: er_cells = 2

  ! settings remembered from the last successful extraction (this session)
  integer(c_int) :: lastregion = er_sphere
  real(c_float) :: lastrsph = 5._c_float
  real(c_float) :: lastrcub = 5._c_float
  integer(c_int) :: lastnx(3) = (/1_c_int,1_c_int,1_c_int/)
  logical :: lastborder = .false.
  logical :: lastmolmotif = .false.
  logical :: lastenviron = .false.
  logical :: lastcloseafter = .true.

contains

  !> Draw the extract cluster as molecule window
  module subroutine draw_extract(w)
    use systems, only: sys, sysc, sys_init, ok_system, add_systems_from_seeds,&
       launch_initialization_thread
    use crystalseedmod, only: crystalseed
    use fragmentmod, only: fragment
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip, iw_checkbox,&
       iw_radiobutton, iw_dragfloat_realc, iw_inputint3,&
       iw_close_event, iw_setpos_bottomright
    use tools_io, only: string, uout
    use param, only: bohrtoa, pi
    class(window), intent(inout), target :: w

    logical :: doquit, okvalid, syschanged, ismol, environok, useenv
    integer :: i, isys, iview, nmol
    real*8 :: rad, x0(3), xlo(3), xhi(3)
    logical(c_bool) :: ldum
    type(fragment) :: fr
    type(fragment), allocatable :: fr0(:)
    logical, allocatable :: isdiscrete(:)
    type(crystalseed), allocatable :: seed(:)

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    if (.not.doquit) doquit = .not.ok_system(isys,sys_init)
    ismol = .false.
    if (.not.doquit) ismol = sys(isys)%c%ismolecule

    ! initialize state
    if (w%firstpass) then
       w%extract_region = lastregion
       w%extract_rsph = lastrsph
       w%extract_rcub = lastrcub
       w%extract_nx = lastnx
       w%extract_border = lastborder
       w%extract_molmotif = lastmolmotif
       w%extract_environ = lastenviron
       w%extract_closeafter = lastcloseafter
    end if

    ! default center (cell center for crystals, origin for molecules);
    ! re-applied when the anchor view changes system. Also compute the
    ! rough atom number density (bohr^-3) for the atom-count estimate:
    ! atoms per cell volume for crystals, per bounding box for molecules
    if (w%firstpass .or. syschanged) then
       w%extract_center = merge(0._c_float,0.5_c_float,ismol)
       w%errmsg = ""
       w%extract_atdens = 0._c_float
       if (.not.doquit) then
          if (ismol) then
             xlo = huge(1d0)
             xhi = -huge(1d0)
             do i = 1, sys(isys)%c%ncel
                xlo = min(xlo,sys(isys)%c%atcel(i)%r)
                xhi = max(xhi,sys(isys)%c%atcel(i)%r)
             end do
             ! pad thin dimensions (e.g. planar molecules)
             xhi = max(xhi - xlo,4d0)
             w%extract_atdens = real(sys(isys)%c%ncel / (xhi(1)*xhi(2)*xhi(3)),c_float)
          elseif (sys(isys)%c%omega > 0d0) then
             w%extract_atdens = real(sys(isys)%c%ncel / sys(isys)%c%omega,c_float)
          end if
       end if
    end if

    ! the unit cells region is only available for crystals
    if (ismol .and. w%extract_region == er_cells) &
       w%extract_region = er_sphere

    ! whether the whole-molecules (environment) option applies: molecular
    ! crystals, or molecular systems with more than one molecule
    environok = .false.
    if (.not.doquit) then
       if (ismol) then
          environok = (sys(isys)%c%nmol > 1)
       else
          environok = sys(isys)%c%ismol3d
       end if
    end if
    useenv = (w%extract_region == er_sphere .and. w%extract_environ .and. environok)

    ! the system being cut
    if (.not.doquit) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)
    end if

    ! close after extracting
    ldum = iw_checkbox("Close after extracting##extractcloseafter",w%extract_closeafter)
    call iw_tooltip("If checked, extracting closes this window and selects the new&
       & system in the view. Otherwise, the window stays open and the view does not&
       & change",ttshown)

    ! region radio buttons
    call iw_text("Region",highlight=.true.)
    ldum = iw_radiobutton("Sphere##extractregion",int=w%extract_region,intval=er_sphere,&
       sameline=.true.)
    call iw_tooltip("Cut all atoms in a sphere around the center",ttshown)
    ldum = iw_radiobutton("Cube##extractregion",int=w%extract_region,intval=er_cube,&
       sameline=.true.)
    call iw_tooltip("Cut all atoms in a cube around the center",ttshown)
    if (.not.ismol) then
       ldum = iw_radiobutton("Unit cell##extractregion",int=w%extract_region,intval=er_cells,&
          sameline=.true.)
       call iw_tooltip("Cut whole unit cells",ttshown)
    end if

    ! center (all regions except cells)
    if (w%extract_region /= er_cells) then
       if (ismol) then
          call iw_text("Center (Å)")
       else
          call iw_text("Center (fractional)")
       end if
       call igPushItemWidth(iw_calcwidth(27,1))
       ldum = iw_dragfloat_realc("##extractcenter",x3=w%extract_center,speed=0.01_c_float,&
          decimal=4,sameline=.true.)
       call igPopItemWidth()
       call iw_tooltip("Center of the region, in fractional coordinates for crystals&
          & and Å for molecules",ttshown)
    end if

    ! region parameters
    if (w%extract_region == er_sphere) then
       ldum = iw_dragfloat_realc("Radius (Å)##extractrsph",x1=w%extract_rsph,speed=0.1_c_float,&
          min=0.1_c_float,max=500._c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Radius of the sphere",ttshown)
       ! rough number of atoms in the sphere from the atom density
       if (w%extract_atdens > 0._c_float) then
          rad = real(w%extract_rsph,8) / bohrtoa
          call iw_text("(~" // string(nint(4d0*pi/3d0 * rad**3 * real(w%extract_atdens,8))) //&
             " atoms)",sameline=.true.)
       end if
       if (environok) then
          ldum = iw_checkbox("Whole molecules##extractenviron",w%extract_environ)
          call iw_tooltip("Cut all whole molecules whose center of mass is within the&
             & radius (ENVIRON option to WRITE) instead of the atoms in the sphere",ttshown)
       end if
    elseif (w%extract_region == er_cube) then
       ldum = iw_dragfloat_realc("Half-edge (Å)##extractrcub",x1=w%extract_rcub,speed=0.1_c_float,&
          min=0.1_c_float,max=500._c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
       call iw_tooltip("Half the edge length of the cube",ttshown)
    elseif (w%extract_region == er_cells) then
       call iw_text("Cells along a, b, c")
       ldum = iw_inputint3("##extractnx",w%extract_nx,width=3*5,sameline=.true.)
       w%extract_nx = max(min(w%extract_nx,50_c_int),1_c_int)
       call iw_tooltip("Number of unit cells along each crystallographic axis",ttshown)
       ldum = iw_checkbox("Include border atoms##extractborder",w%extract_border)
       call iw_tooltip("Include atoms near the boundary of the chosen cells (BORDER&
          & option to WRITE)",ttshown)
    end if

    ! complete boundary molecules (not applicable to whole-molecule cuts)
    if (.not.useenv) then
       ldum = iw_checkbox("Complete boundary molecules##extractmolmotif",w%extract_molmotif)
       call iw_tooltip("Add the atoms required to complete the molecules cut by the&
          & region boundary (MOLMOTIF option to WRITE)",ttshown)
    end if

    ! maybe the error message
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(10,2)

    ! final buttons: Extract
    okvalid = .not.doquit
    if (iw_button("Extract",disabled=.not.okvalid)) then
       w%errmsg = ""

       ! center of the region: fractional directly for crystals; for
       ! molecules, Å in the molecular frame -> internal Cartesian ->
       ! fractional (same conversion as the WRITE keyword)
       x0 = real(w%extract_center,8)
       if (ismol) then
          x0 = x0 / bohrtoa - sys(isys)%c%molx0
          x0 = sys(isys)%c%c2x(x0)
       end if

       ! build the fragment for the selected region (radii in bohr)
       if (w%extract_region == er_sphere) then
          rad = real(w%extract_rsph,8) / bohrtoa
          if (useenv) then
             fr = sys(isys)%c%listatoms_environ(rad,x0)
          else
             fr = sys(isys)%c%listatoms_sphcub(rsph=rad,xsph=x0)
          end if
       elseif (w%extract_region == er_cube) then
          rad = real(w%extract_rcub,8) / bohrtoa
          fr = sys(isys)%c%listatoms_sphcub(rcub=rad,xcub=x0)
       else
          fr = sys(isys)%c%listatoms_cells(int(w%extract_nx),w%extract_border)
       end if

       if (fr%nat == 0) then
          w%errmsg = "No atoms in the region"
       else
          ! complete the boundary molecules; refuse non-finite (periodic)
          ! components, which from_fragment cannot handle
          if (.not.useenv .and. w%extract_molmotif) then
             call sys(isys)%c%listmolecules(fr,nmol,fr0,isdiscrete)
             if (.not.all(isdiscrete(1:nmol))) then
                w%errmsg = "The boundary molecules form a periodic network"
             else
                call fr%merge_array(fr0,.false.)
             end if
          end if
       end if
       if (len_trim(w%errmsg) == 0) then
          ! for molecules, shift back to the user (absolute) Cartesian frame
          if (ismol) then
             do i = 1, fr%nat
                fr%at(i)%r = fr%at(i)%r + sys(isys)%c%molx0
             end do
          end if

          ! build the seed and add the new system; select it in the
          ! view only if the window closes after extracting
          allocate(seed(1))
          call seed(1)%from_fragment(fr)
          seed(1)%name = trim(sysc(isys)%seed%name) // " (cluster)"
          call add_systems_from_seeds(1,seed,noselect=.not.w%extract_closeafter)
          call launch_initialization_thread()
          write (uout,'("Extracted cluster (",A," atoms) from system ",A,": ",A)') &
             string(fr%nat), string(isys), trim(sysc(isys)%seed%name)

          ! remember the settings
          lastregion = w%extract_region
          lastrsph = w%extract_rsph
          lastrcub = w%extract_rcub
          lastnx = w%extract_nx
          lastborder = w%extract_border
          lastmolmotif = w%extract_molmotif
          lastenviron = w%extract_environ
          lastcloseafter = w%extract_closeafter

          ! maybe close the window
          if (w%extract_closeafter) doquit = .true.
       end if
    end if
    call iw_tooltip("Cut the region from the system and add it as a new molecular system",ttshown)

    ! final buttons: close
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_extract

end submodule extract
