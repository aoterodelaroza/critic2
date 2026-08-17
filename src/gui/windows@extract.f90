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

  ! what is taken from the region (mutually exclusive, as in the WRITE
  ! keyword, which refuses MOLMOTIF and ENVIRON together)
  integer(c_int), parameter :: ei_atoms = 0 ! the atoms inside the region
  integer(c_int), parameter :: ei_molmotif = 1 ! whole molecules with any atom inside (MOLMOTIF)
  integer(c_int), parameter :: ei_environ = 2 ! whole molecules with COM inside (ENVIRON)

  ! settings remembered from the last successful extraction (this session)
  integer(c_int) :: lastregion = er_sphere
  real(c_float) :: lastrsph = 5._c_float
  real(c_float) :: lastrcub = 5._c_float
  integer(c_int) :: lastnx(3) = (/1_c_int,1_c_int,1_c_int/)
  logical :: lastborder = .false.
  integer(c_int) :: lastinclude = ei_atoms
  logical :: lastcloseafter = .true.

  ! transient center marker: color, opacity, and radius limits (bohr)
  real(c_float), parameter :: center_rgb(3) = (/1.0_c_float,0.6_c_float,0.1_c_float/)
  real(c_float), parameter :: center_alpha = 0.6_c_float
  real*8, parameter :: center_radmin = 1.7d0 ! larger than most drawn atoms (0.7*rcov), so it survives burial
  real*8, parameter :: center_radmax = 3.0d0

contains

  !> Draw the extract cluster as molecule window
  module subroutine draw_extract(w)
    use systems, only: sys, sysc, sys_init, ok_system, add_systems_from_seeds,&
       launch_initialization_thread
    use crystalseedmod, only: crystalseed
    use fragmentmod, only: fragment
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip, iw_checkbox,&
       iw_radiobutton, iw_dragfloat_realc, iw_intstepper,&
       iw_close_event, iw_setpos_bottomright
    use tools_io, only: string, uout
    use param, only: bohrtoa, pi
    class(window), intent(inout), target :: w

    character(len=*,kind=c_char), parameter :: ttnx = &
       "Number of unit cells along each crystallographic axis"

    logical :: doquit, syschanged, ismol, molmotifok, environok, okpick, hovered
    integer :: i, isys, iview, nmol, ipad
    integer(c_int) :: inc
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
       w%extract_include = lastinclude
       w%extract_closeafter = lastcloseafter
       w%extract_picking = .false.
    end if

    ! handle a pending center pick commanded to the parent view; the
    ! anchor-gone and pick-taken-over cases only clear the pending flag
    ! (the view's owner watchdog releases the mode if this window dies)
    if (w%extract_picking) then
       if (doquit) then
          w%extract_picking = .false.
       elseif (win(iview)%vmdata%owner /= w%id) then
          w%extract_picking = .false.
       elseif (syschanged .or. w%extract_pick%is_stale(sysc(isys)%timelastchange_geometry)) then
          ! the view moved to another system or the geometry changed (stale
          ! cell-atom ids): cancel the pick
          call win(iview)%viewmode_release_forced(w%id)
          w%extract_picking = .false.
       elseif (win(iview)%viewmode >= 0) then
          ! the pick finished (flag = 0 means cancelled)
          if (win(iview)%vmdata%flag == 1) then
             if (win(iview)%vmdata%idx(1) > 0) then
                ! an atom was clicked: use its position
                i = win(iview)%vmdata%idx(1)
                call center_from_frac(sys(isys)%c%atcel(i)%x + win(iview)%vmdata%idx(2:4))
             else
                ! empty space: unproject the click onto the plane of the scene center
                call view_click_frame(iview,win(iview)%vmdata%xpos,x0,okpick)
                if (okpick) call center_from_frac(sys(isys)%c%c2x(x0))
             end if
          end if
          w%extract_picking = .false.
          win(iview)%vmdata%idx = 0
          win(iview)%vmdata%flag = 0
       end if
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

    ! which whole-molecule rules this system admits: completing the
    ! molecules does something as soon as one fragment of the motif is
    ! discrete (only those are completed), while selecting molecules by
    ! their center needs every fragment to be discrete. Both can change
    ! without the system changing (recalculating the bonds redoes the
    ! molecular motif), so they only hide options; the chosen rule is
    ! not overwritten, and is read through ei_effective() instead
    molmotifok = .false.
    environok = .false.
    if (.not.doquit) then
       molmotifok = (sys(isys)%c%nmoldiscrete > 0)
       environok = molmotifok .and. (sys(isys)%c%nmoldiscrete == sys(isys)%c%nmol)
    end if

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
    call iw_text("Region",highlight=.true.,alignframe=.true.)
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
       call iw_text("Center",highlight=.true.,alignframe=.true.)
       if (iw_button("Pick##extractpick",sameline=.true.,disabled=(w%extract_picking .or. doquit))) then
          w%extract_picking = .true.
          call w%extract_pick%arm()
          call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the region center",w%id,&
             acceptempty=.true.)
       end if
       hovered = logical(igIsItemHovered(ImGuiHoveredFlags_None))
       call iw_tooltip("Pick the center in the view: click an atom to use its position,&
          & or empty space to use the clicked point",ttshown)
       ! call igPushItemWidth(iw_calcwidth(24,1))
       ldum = iw_dragfloat_realc("##extractcenter",x3=w%extract_center,speed=0.01_c_float,&
          decimal=4,sameline=.true.)
       ! call igPopItemWidth()
       hovered = hovered .or. logical(igIsItemHovered(ImGuiHoveredFlags_None)) .or. logical(igIsItemActive())
       call iw_tooltip("Center of the region",ttshown)
       if (ismol) then
          call iw_text("(Å)",sameline=.true.)
       else
          call iw_text("(frac)",sameline=.true.)
       end if

       ! transient marker at the center while hovering the Pick button,
       ! hovering/dragging the coordinates, or picking; show_transient_sphere
       ! takes the absolute Cartesian frame (with molx0 for molecules); radius
       ! scaled to the scene so it is visible
       if ((hovered .or. w%extract_picking) .and. .not.doquit) then
          if (sysc(isys)%sc%isinit /= 0) then
             x0 = sys(isys)%c%x2c(center_to_frac())
             if (ismol) x0 = x0 + sys(isys)%c%molx0
             rad = min(max(0.05d0 * real(sysc(isys)%sc%scenerad,8),center_radmin),center_radmax)
             call sysc(isys)%sc%show_transient_sphere(w%id,1,x0,rad,center_rgb,center_alpha)
          end if
       end if
    end if

    ! region parameters
    if (w%extract_region == er_sphere) then
       call iw_text("Radius",highlight=.true.,alignframe=.true.)
       ldum = iw_dragfloat_realc("(Å)##extractrsph",x1=w%extract_rsph,speed=0.1_c_float,&
          min=0.1_c_float,max=500._c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
       call iw_tooltip("Radius of the sphere",ttshown)
       rad = real(w%extract_rsph,8) / bohrtoa
       call show_atom_estimate(4d0*pi/3d0 * rad**3)
    elseif (w%extract_region == er_cube) then
       call iw_text("Half-edge",highlight=.true.,alignframe=.true.)
       ldum = iw_dragfloat_realc("(Å)##extractrcub",x1=w%extract_rcub,speed=0.1_c_float,&
          min=0.1_c_float,max=500._c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
       call iw_tooltip("Half the edge length of the cube",ttshown)
       rad = real(w%extract_rcub,8) / bohrtoa
       call show_atom_estimate(8d0 * rad**3)
    elseif (w%extract_region == er_cells) then
       ! same form as the periodicity widget in the representation editor
       ! (all three fields share a digit width, and clamp themselves)
       call iw_text("Cells",highlight=.true.,alignframe=.true.)
       ipad = ceiling(log10(max(maxval(w%extract_nx),1) + 0.1))
       ldum = iw_intstepper("aaxis##extractnx",w%extract_nx(1),label="a:",minval=1_c_int,&
          maxval=50_c_int,ndigit=ipad,notlive=.true.,sameline=.true.,tooltip=ttnx)
       ldum = iw_intstepper("baxis##extractnx",w%extract_nx(2),label="b:",minval=1_c_int,&
          maxval=50_c_int,ndigit=ipad,notlive=.true.,sameline=.true.,tooltip=ttnx)
       ldum = iw_intstepper("caxis##extractnx",w%extract_nx(3),label="c:",minval=1_c_int,&
          maxval=50_c_int,ndigit=ipad,notlive=.true.,sameline=.true.,tooltip=ttnx)
       ! the border only adds atoms to an atom cut
       if (ei_effective() /= ei_environ) then
          ldum = iw_checkbox("Include border atoms##extractborder",w%extract_border)
          call iw_tooltip("Include atoms near the boundary of the chosen cells (BORDER&
             & option to WRITE)",ttshown)
       end if
    end if

    ! what is taken from the region? Only the rules this system admits
    ! are offered, and if that leaves the atoms rule alone the one-entry
    ! group would read as broken, so it goes away. The buttons run on the
    ! rule in force, and only a click writes the choice back
    if (molmotifok) then
       inc = ei_effective()
       call iw_text("Include",highlight=.true.,alignframe=.true.)
       if (iw_radiobutton("Atoms in the region##extractinclude",int=inc,&
          intval=ei_atoms)) w%extract_include = inc
       call iw_tooltip("Extract only atoms inside the region. Molecules crossing the boundary&
          & are left incomplete",ttshown)
       if (iw_radiobutton("Whole molecules with any atom in the region##extractinclude",&
          int=inc,intval=ei_molmotif)) w%extract_include = inc
       call iw_tooltip("Extract every molecule that has at least one atom inside the&
          & region, completing it with atoms from outside",ttshown)
       if (environok) then
          if (iw_radiobutton("Whole molecules centered in the region##extractinclude",&
             int=inc,intval=ei_environ)) w%extract_include = inc
          call iw_tooltip("Extract every molecule whose center of mass&
             & is in the region, even if some of their atoms stick out.",ttshown)
       end if
    end if

    ! maybe the error message
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(12,2)

    ! final buttons: Extract
    if (iw_button("Extract",disabled=doquit)) then
       w%errmsg = ""
       inc = ei_effective()

       ! center of the region in crystallographic coordinates (same
       ! conversion as the WRITE keyword)
       x0 = center_to_frac()

       ! build the fragment for the selected region (radii in bohr); the
       ! center-of-mass rule picks whole molecules by itself, replacing
       ! the atom cut
       if (w%extract_region == er_sphere) then
          rad = real(w%extract_rsph,8) / bohrtoa
          if (inc == ei_environ) then
             fr = sys(isys)%c%listatoms_molcenter(rsph=rad,xsph=x0)
          else
             fr = sys(isys)%c%listatoms_sphcub(rsph=rad,xsph=x0)
          end if
       elseif (w%extract_region == er_cube) then
          rad = real(w%extract_rcub,8) / bohrtoa
          if (inc == ei_environ) then
             fr = sys(isys)%c%listatoms_molcenter(rcub=rad,xcub=x0)
          else
             fr = sys(isys)%c%listatoms_sphcub(rcub=rad,xcub=x0)
          end if
       else
          if (inc == ei_environ) then
             fr = sys(isys)%c%listatoms_molcenter(nx=int(w%extract_nx))
          else
             fr = sys(isys)%c%listatoms_cells(int(w%extract_nx),w%extract_border)
          end if
       end if

       if (fr%nat == 0) then
          w%errmsg = "No atoms in the region"
       else
          ! complete the molecules cut by the region boundary. A component
          ! that extends indefinitely (a framework) cannot be completed
          ! with a finite number of atoms, so it is left as it was cut:
          ! merge only the discrete ones, on top of the atoms already in
          ! the fragment (merge_array discards repeats)
          if (inc == ei_molmotif) then
             call sys(isys)%c%listmolecules(fr,nmol,fr0,isdiscrete)
             if (any(isdiscrete(1:nmol))) &
                call fr%merge_array(pack(fr0(1:nmol),isdiscrete(1:nmol)),.true.)
          end if

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
          lastinclude = w%extract_include
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

  contains
    ! The region center is stored in extract_center as displayed:
    ! fractional for crystals, Å in the absolute molecular frame for
    ! molecules. These convert to/from crystallographic coordinates.
    function center_to_frac() result(x)
      real*8 :: x(3)
      x = real(w%extract_center,8)
      if (ismol) x = sys(isys)%c%c2x(x / bohrtoa - sys(isys)%c%molx0)
    end function center_to_frac

    subroutine center_from_frac(x)
      real*8, intent(in) :: x(3)
      if (ismol) then
         w%extract_center = real((sys(isys)%c%x2c(x) + sys(isys)%c%molx0) * bohrtoa,c_float)
      else
         w%extract_center = real(x,c_float)
      end if
    end subroutine center_from_frac

    ! The include rule actually in force. A rule this system does not
    ! admit falls back to the next one down, but the stored choice is
    ! kept (the molecular motif may come back on the next rebond)
    ! Write, next to the region size, a rough estimate of the number of
    ! atoms in a region of volume vol (bohr^3), from the atom number
    ! density of the system. A molecule has no periodic images, so the
    ! cut cannot contain more atoms than the system itself
    subroutine show_atom_estimate(vol)
      real*8, intent(in) :: vol

      integer :: nest

      if (w%extract_atdens <= 0._c_float) return
      nest = nint(vol * real(w%extract_atdens,8))
      if (ismol) nest = min(nest,sys(isys)%c%ncel)
      call iw_text("(~" // string(nest) // " atoms)",sameline=.true.)

    end subroutine show_atom_estimate

    function ei_effective() result(ei)
      integer(c_int) :: ei
      ei = w%extract_include
      if (ei == ei_environ .and. .not.environok) ei = ei_molmotif
      if (ei == ei_molmotif .and. .not.molmotifok) ei = ei_atoms
    end function ei_effective
  end subroutine draw_extract

end submodule extract
