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
  use param, only: bohrtoa
  implicit none

  ! region shapes
  integer(c_int), parameter :: er_sphere = 0
  integer(c_int), parameter :: er_cube = 1
  integer(c_int), parameter :: er_cells = 2

  ! how many molecules come out of the region
  integer(c_int), parameter :: em_one = 0 ! a single molecular system
  integer(c_int), parameter :: em_several = 1 ! the n-mers it contains (NMER)

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
  integer(c_int) :: lastmode = em_one
  integer(c_int) :: lastinclude = ei_atoms
  integer(c_int) :: lastnmer = 2
  logical :: lastnmer_do(nmer_max) = .true.
  integer(c_int) :: lastnmer_any(nmer_max) = 0
  real(c_float) :: lastnmer_dist(nmer_max) = 1e6_c_float
  logical :: lastcloseafter = .true.

  ! the region preview style (region_rgb/region_alpha/region_edgerad) is
  ! shared with the isosurface editor and lives in the windows module

  ! transient center marker: color, opacity, and radius limits (bohr)
  real(c_float), parameter :: center_rgb(3) = (/1.0_c_float,0.6_c_float,0.1_c_float/)
  real(c_float), parameter :: center_alpha = 0.6_c_float
  real*8, parameter :: center_radmin = 1.7d0 ! larger than most drawn atoms (0.7*rcov), so it survives burial
  real*8, parameter :: center_radmax = 3.0d0

contains

  !> Draw the extract cluster as molecule window
  module subroutine draw_extract(w)
    use systems, only: sys, sysc, sys_init, ok_system, add_systems_from_seeds,&
       add_group_master, launch_initialization_thread
    use crystalmod, only: nmer_name
    use crystalseedmod, only: crystalseed
    use fragmentmod, only: fragment
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip, iw_checkbox,&
       iw_radiobutton, iw_dragfloat_realc, iw_intstepper, iw_combo_simple,&
       iw_table_column, iw_close_event, iw_setpos_bottomright
    use tools_io, only: string, uout
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG
    use param, only: pi, eye
    class(window), intent(inout), target :: w

    character(len=*,kind=c_char), parameter :: ttnx = &
       "Number of unit cells along each crystallographic axis"

    logical :: doquit, syschanged, ismol, molmotifok, environok, severalok, okpick, hovered
    integer :: i, isys, iview, nmol, ipad
    integer(c_int) :: inc, emode, tflags
    real*8 :: rad, x0(3), xlo(3), xhi(3), vbox(3,3)
    real(c_float) :: dmax
    type(ImVec2) :: sz0
    character(len=:,kind=c_char), allocatable, target :: str1
    logical(c_bool) :: ldum
    logical :: ok
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
       w%extract_mode = lastmode
       w%extract_include = lastinclude
       w%extract_nmer = lastnmer
       w%extract_nmer_do = lastnmer_do
       w%extract_nmer_any = lastnmer_any
       w%extract_nmer_dist = lastnmer_dist
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
       w%okmsg = ""
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
    severalok = .false.
    if (.not.doquit) then
       molmotifok = (sys(isys)%c%nmoldiscrete > 0)
       environok = molmotifok .and. (sys(isys)%c%nmoldiscrete == sys(isys)%c%nmol)
       ! with a single molecule in the system there is nothing to break into n-mers
       severalok = (sys(isys)%c%nmoldiscrete > 1)
    end if

    ! the system being cut
    if (.not.doquit) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)
    end if

    ! how many molecules come out of the region: the region as a single
    ! molecule, or the n-mers it contains
    if (severalok) then
       emode = em_effective()
       call iw_text("Extract",highlight=.true.)
       if (iw_radiobutton("Single molecule##extractmode",int=emode,intval=em_one)) &
          w%extract_mode = emode
       call iw_tooltip("Extract the region as a single molecular system",ttshown)
       if (iw_radiobutton("Monomers, dimers, etc.##extractmode",int=emode,intval=em_several,&
          sameline=.true.)) w%extract_mode = emode
       call iw_tooltip("Extract the molecules in the region as separate monomers, dimers, etc.",ttshown)
    end if

    ! region radio buttons
    call iw_text("Region",highlight=.true.)
    ldum = iw_radiobutton("Sphere##extractregion",int=w%extract_region,intval=er_sphere)
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
       call iw_text("Center",alignframe=.true.)
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
             rad = min(max(0.05d0 * real(sysc(isys)%sc%scenerad,8),center_radmin),center_radmax)
             call sysc(isys)%sc%show_transient_sphere(w%id,1,center_to_cart(),rad,center_rgb,center_alpha)
          end if
       end if
    end if

    ! show the chosen region in the view for as long as the window is open:
    ! a translucent sphere, or the wireframe of the cube or cell block
    if (.not.doquit) then
       if (sysc(isys)%sc%isinit /= 0) then
          if (w%extract_region == er_sphere) then
             call sysc(isys)%sc%show_transient_sphere(w%id,2,center_to_cart(),&
                real(w%extract_rsph,8)/bohrtoa,region_rgb,region_alpha)
          else
             if (w%extract_region == er_cube) then
                ! axis-aligned box of side 2*half-edge, centered on the region center
                rad = real(w%extract_rcub,8) / bohrtoa
                vbox = 2d0 * rad * eye
                x0 = center_to_cart() - rad
             else
                ! the cell block, anchored at the origin of the main cell; the
                ! columns of m_x2c are the lattice vectors in cartesian
                do i = 1, 3
                   vbox(:,i) = sys(isys)%c%m_x2c(:,i) * int(w%extract_nx(i))
                end do
                x0 = 0d0
                if (ismol) x0 = sys(isys)%c%molx0
             end if
             call sysc(isys)%sc%show_transient_box(w%id,2,x0,vbox,region_edgerad,region_rgb,region_alpha)
          end if
       end if
    end if

    ! region parameters
    if (w%extract_region == er_sphere) then
       call iw_text("Radius",alignframe=.true.)
       ldum = iw_dragfloat_realc("(Å)##extractrsph",x1=w%extract_rsph,speed=0.1_c_float,&
          min=0.1_c_float,max=500._c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
       call iw_tooltip("Radius of the sphere",ttshown)
       call show_atom_estimate(region_volume())
    elseif (w%extract_region == er_cube) then
       call iw_text("Half-edge",alignframe=.true.)
       ldum = iw_dragfloat_realc("(Å)##extractrcub",x1=w%extract_rcub,speed=0.1_c_float,&
          min=0.1_c_float,max=500._c_float,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
       call iw_tooltip("Half the edge length of the cube",ttshown)
       call show_atom_estimate(region_volume())
    elseif (w%extract_region == er_cells) then
       ! same form as the periodicity widget in the representation editor
       ! (all three fields share a digit width, and clamp themselves)
       call iw_text("Cells",alignframe=.true.)
       ipad = ceiling(log10(max(maxval(w%extract_nx),1) + 0.1))
       ldum = iw_intstepper("aaxis##extractnx",w%extract_nx(1),label="a:",minval=1_c_int,&
          maxval=50_c_int,ndigit=ipad,notlive=.true.,sameline=.true.,tooltip=ttnx)
       ldum = iw_intstepper("baxis##extractnx",w%extract_nx(2),label="b:",minval=1_c_int,&
          maxval=50_c_int,ndigit=ipad,notlive=.true.,sameline=.true.,tooltip=ttnx)
       ldum = iw_intstepper("caxis##extractnx",w%extract_nx(3),label="c:",minval=1_c_int,&
          maxval=50_c_int,ndigit=ipad,notlive=.true.,sameline=.true.,tooltip=ttnx)
       ! a whole number of cells has an exact atom count, but only for a plain
       ! atom cut with no border shell: the border and both molecular rules
       ! bring in atoms from outside the chosen cells
       if (.not.doquit) then
          if (ei_effective() == ei_atoms .and. .not.w%extract_border) then
             call show_atom_estimate(nexact=sys(isys)%c%ncel * product(int(w%extract_nx)))
          else
             call show_atom_estimate(region_volume())
          end if
       end if
       ! the border only adds atoms to an atom cut
       if (ei_effective() /= ei_environ) then
          ldum = iw_checkbox("Include border atoms##extractborder",w%extract_border)
          call iw_tooltip("Include atoms near the boundary of the chosen cells",ttshown)
       end if
    end if

    ! what is taken from the region? Only the rules this system admits
    ! are offered, and if that leaves the atoms rule alone the one-entry
    ! group would read as broken, so it goes away. The buttons run on the
    ! rule in force, and only a click writes the choice back
    if (molmotifok .and. em_effective() == em_one) then
       inc = ei_effective()
       call iw_text("Include",highlight=.true.)
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

    ! n-mer options: the highest order, then one row per order with
    ! whether it is generated, how the center-of-mass distance filter is
    ! applied, and the cutoff itself (capped by the size of the region)
    if (em_effective() == em_several) then
       call iw_text("Include",highlight=.true.)
       call iw_text("Up to n =",alignframe=.true.)
       ldum = iw_intstepper("nmer##extractnmer",w%extract_nmer,minval=1_c_int,&
          maxval=int(nmer_max,c_int),notlive=.true.,sameline=.true.,&
          tooltip="Highest n-mer order: 2 gives monomers and dimers, 3 adds trimers, and so on")
       dmax = real(region_diameter() * bohrtoa,c_float)

       tflags = ImGuiTableFlags_None
       tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
       tflags = ior(tflags,ImGuiTableFlags_RowBg)
       tflags = ior(tflags,ImGuiTableFlags_Borders)
       tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
       str1 = "##extractnmertable" // c_null_char
       sz0%x = 0._c_float
       sz0%y = 0._c_float
       if (igBeginTable(c_loc(str1),5,tflags,sz0,0._c_float)) then
          call iw_table_column("Nmer",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Write",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Distances",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Cutoff (Å)",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Number ('any')",id=4_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
          call igTableHeadersRow()
          do i = 1, int(w%extract_nmer)
             call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
             if (igTableSetColumnIndex(0_c_int)) call iw_text(string(i))
             if (igTableSetColumnIndex(1_c_int)) then
                ldum = iw_checkbox("##extractnmerdo" // string(i),w%extract_nmer_do(i))
                call iw_tooltip("Generate the " // trim(nmer_name(i)) // "s",ttshown)
             end if
             ! a monomer has no pair of molecules, so it has no distances to filter
             if (igTableSetColumnIndex(2_c_int)) then
                if (i == 1) then
                   call iw_text("--",disabled=.true.)
                else
                   call igPushItemWidth(iw_calcwidth(8,1))
                   call iw_combo_simple("##extractnmerany" // string(i),"all" // c_null_char //&
                      "any" // c_null_char,w%extract_nmer_any(i))
                   call igPopItemWidth()
                   call iw_tooltip("all: keep the " // trim(nmer_name(i)) // " only if every pair of&
                      & molecules is within the cutoff. any: keep it if at least one pair is",ttshown)
                end if
             end if
             if (igTableSetColumnIndex(3_c_int)) then
                if (i == 1) then
                   call iw_text("--",disabled=.true.)
                else
                   if (dmax > 0._c_float) &
                      w%extract_nmer_dist(i) = min(w%extract_nmer_dist(i),dmax)
                   ldum = iw_dragfloat_realc("##extractnmerdist" // string(i),&
                      x1=w%extract_nmer_dist(i),speed=0.1_c_float,min=0._c_float,max=dmax,&
                      decimal=2,flags=ImGuiSliderFlags_AlwaysClamp)
                   call iw_tooltip("Maximum distance between the centers of mass of two molecules&
                      & in the " // trim(nmer_name(i)) // "; at most the size of the region",ttshown)
                end if
             end if
             if (igTableSetColumnIndex(4_c_int)) then
                if (w%extract_nmer_do(i)) then
                   call iw_text("~" // string(nmer_estimate(i)))
                else
                   call iw_text("--",disabled=.true.)
                end if
             end if
          end do
          call igEndTable()
       end if
    end if

    ! the outcome of the last extraction
    if (len_trim(w%errmsg) > 0) then
       call iw_text(w%errmsg,danger=.true.)
    elseif (len_trim(w%okmsg) > 0) then
       call iw_text(w%okmsg,highlight=.true.)
    end if

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(34,2,ncheck=1)

    ! close after extracting, next to the button it applies to
    ldum = iw_checkbox("Close after extracting##extractcloseafter",w%extract_closeafter)
    call iw_tooltip("If checked, extracting closes this window and selects the new&
       & system in the view. Otherwise, the window stays open and the view does not&
       & change",ttshown)

    ! final buttons: Extract
    ok = iw_button("Extract",disabled=doquit,sameline=.true.)
    ok = ok .or. (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG) .and. .not.doquit)
    if (ok) then
       w%errmsg = ""
       w%okmsg = ""
       if (em_effective() == em_several) then
          call extract_nmers()
       else
          call extract_one()
       end if

       if (len_trim(w%errmsg) == 0) then
          ! remember the settings
          lastregion = w%extract_region
          lastrsph = w%extract_rsph
          lastrcub = w%extract_rcub
          lastnx = w%extract_nx
          lastborder = w%extract_border
          lastmode = w%extract_mode
          lastinclude = w%extract_include
          lastnmer = w%extract_nmer
          lastnmer_do = w%extract_nmer_do
          lastnmer_any = w%extract_nmer_any
          lastnmer_dist = w%extract_nmer_dist
          lastcloseafter = w%extract_closeafter

          ! maybe close the window
          if (w%extract_closeafter) doquit = .true.
       end if
    end if
    call iw_tooltip("Cut the region from the system and add its monomers, dimers,...&
       & to the tree as new systems",ttshown)

    ! final buttons: close
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding (the OK binding is
    ! handled above by the Extract button, and closes only if requested)
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

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

    ! The region center in the absolute Cartesian frame (bohr), which is
    ! the frame the scene's transient shapes are drawn in
    function center_to_cart() result(x)
      real*8 :: x(3)
      x = sys(isys)%c%x2c(center_to_frac())
      if (ismol) x = x + sys(isys)%c%molx0
    end function center_to_cart

    subroutine center_from_frac(x)
      real*8, intent(in) :: x(3)
      if (ismol) then
         w%extract_center = real((sys(isys)%c%x2c(x) + sys(isys)%c%molx0) * bohrtoa,c_float)
      else
         w%extract_center = real(x,c_float)
      end if
    end subroutine center_from_frac

    ! Write, next to the region size, the number of atoms the region
    ! contains: nexact if it is known exactly, otherwise a rough estimate
    ! from the atom number density of the system and the region volume vol
    ! (bohr^3). A molecule has no periodic images, so the cut cannot
    ! contain more atoms than the system itself
    subroutine show_atom_estimate(vol,nexact)
      real*8, intent(in), optional :: vol
      integer, intent(in), optional :: nexact

      integer :: nest

      if (present(nexact)) then
         call iw_text("(" // string(nexact) // " atoms)",sameline=.true.)
         return
      end if
      if (w%extract_atdens <= 0._c_float) return
      nest = nint(vol * real(w%extract_atdens,8))
      if (ismol) nest = min(nest,sys(isys)%c%ncel)
      call iw_text("(~" // string(nest) // " atoms)",sameline=.true.)

    end subroutine show_atom_estimate

    ! Extract the region as a single molecular system.
    subroutine extract_one()

      inc = ei_effective()

      ! center of the region in crystallographic coordinates
      x0 = center_to_frac()

      ! build the fragment for the selected region
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
         ! complete the molecules cut by the region boundaryj
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
         w%okmsg = "Extracted one system with " // string(fr%nat) // " atoms"
      end if

    end subroutine extract_one

    ! Extract the k-mers (k up to extract_nmer) that can be made from
    ! the molecules whose center of mass is inside the region, subject
    ! to the per-order center-of-mass distance filter. Every k-mer
    ! contains at least one main-cell molecule and is built with it
    ! first, so k-mers related to each other by a lattice translation
    ! are generated once.
    subroutine extract_nmers()
      use crystalseedmod, only: realloc_crystalseed
      use tools_math, only: nchoosek, comb

      integer :: k, l, j, npool, nseed, imain, ncomb, icnt(nmer_max), ntotal, imas, ipass
      integer :: idmol, iat, lvec(3)
      integer, allocatable :: icomb(:)
      logical, allocatable :: ismain(:)
      real*8, allocatable :: com(:,:)
      type(fragment), allocatable :: frpool(:)
      type(fragment) :: frk, frdum
      real*8 :: ntot

      ! the pool: the molecules with their center of mass in the region
      x0 = center_to_frac()
      if (w%extract_region == er_sphere) then
         rad = real(w%extract_rsph,8) / bohrtoa
         frdum = sys(isys)%c%listatoms_molcenter(rsph=rad,xsph=x0,fr0=frpool)
      elseif (w%extract_region == er_cube) then
         rad = real(w%extract_rcub,8) / bohrtoa
         frdum = sys(isys)%c%listatoms_molcenter(rcub=rad,xcub=x0,fr0=frpool)
      else
         frdum = sys(isys)%c%listatoms_molcenter(nx=int(w%extract_nx),fr0=frpool)
      end if
      npool = 0
      if (allocated(frpool)) npool = size(frpool,1)
      if (npool == 0) then
         w%errmsg = "No molecule has its center of mass in the region"
         return
      end if

      ! which of them come from the main cell, and where their centers are
      allocate(ismain(npool),com(3,npool))
      do j = 1, npool
         idmol = sys(isys)%c%idatcelmol(1,frpool(j)%at(1)%cidx)
         iat = sys(isys)%c%idatcelmol(2,frpool(j)%at(1)%cidx)
         lvec = nint(frpool(j)%at(1)%x - sys(isys)%c%mol(idmol)%at(iat)%x)
         ismain(j) = all(lvec == 0)
         com(:,j) = frpool(j)%cmass()
      end do
      if (.not.any(ismain)) then
         w%errmsg = "No main-cell molecule has its center of mass in the region"
         return
      end if

      ! give up before enumerating if there are too many combinations to walk
      ntot = 0d0
      do k = 1, int(w%extract_nmer)
         if (.not.w%extract_nmer_do(k) .or. k > npool) cycle
         ntot = ntot + binomial(npool,k)
      end do

      ! walk the combinations, keeping the ones anchored on the main cell
      ! that pass the distance filter. The first pass only counts, so that a
      ! batch too large for the tree is refused before any of its groups has
      ! been created; the second pass builds them
      ntotal = 0
      icnt = 0
      do ipass = 1, 2
         do k = 1, int(w%extract_nmer)
            if (.not.w%extract_nmer_do(k) .or. k > npool) cycle
            if (ipass == 2) then
               if (icnt(k) == 0) cycle
               if (allocated(seed)) deallocate(seed)
               allocate(seed(icnt(k)))
            end if
            nseed = 0
            ncomb = nchoosek(npool,k)
            allocate(icomb(k))
            do l = 1, ncomb
               call comb(npool,k,l,icomb)

               ! anchor the k-mer on the first main-cell molecule it contains
               imain = 0
               do j = 1, k
                  if (ismain(icomb(j))) then
                     imain = j
                     exit
                  end if
               end do
               if (imain == 0) cycle
               if (k > 1) then
                  if (.not.distok(k,icomb,com)) cycle
               end if
               nseed = nseed + 1
               if (ipass == 1) then
                  ntotal = ntotal + 1
                  cycle
               end if

               ! build the k-mer with the main-cell molecule first
               frk = frpool(icomb(imain))
               do j = 1, k
                  if (j == imain) cycle
                  call frk%append(frpool(icomb(j)))
               end do
               if (ismol) then
                  do j = 1, frk%nat
                     frk%at(j)%r = frk%at(j)%r + sys(isys)%c%molx0
                  end do
               end if
               if (nseed > size(seed,1)) call realloc_crystalseed(seed,2*nseed)
               call seed(nseed)%from_fragment(frk)
               seed(nseed)%name = trim(sysc(isys)%seed%name) // " (" // nmer_name(k) //&
                  " " // string(nseed) // ")"
            end do
            deallocate(icomb)

            ! this order goes into the tree as its own collapsed group, headed
            ! by an entry that names it
            if (ipass == 1) then
               icnt(k) = nseed
            else
               call add_group_master(nmer_name(k) // "s",trim(sysc(isys)%seed%name),imas)
               call add_systems_from_seeds(nseed,seed,imaster=imas,noselect=.true.)
               deallocate(seed)
            end if
         end do

         if (ntotal == 0) then
            w%errmsg = "No n-mer passes the distance filter"
            return
         end if
      end do

      call launch_initialization_thread()
      write (uout,'("Extracted ",A," n-mers from system ",A,": ",A)') &
         string(ntotal), string(isys), trim(sysc(isys)%seed%name)
      w%okmsg = "Extracted " // string(ntotal) // " systems:"
      do k = 1, int(w%extract_nmer)
         if (icnt(k) == 0) cycle
         w%okmsg = w%okmsg // " " // string(icnt(k)) // " " // nmer_name(k) // "s"
         write (uout,'("  ",A,": ",A)') nmer_name(k) // "s", string(icnt(k))
      end do
      write (uout,*)

    end subroutine extract_nmers

    ! Whether the k-mer made of the pool molecules ic(1:k), with centers
    ! of mass com, passes the distance filter of order k: either every
    ! pair of molecules is within the cutoff, or at least one pair is
    function distok(k,ic,com) result(ok)
      integer, intent(in) :: k, ic(k)
      real*8, intent(in) :: com(:,:)
      logical :: ok

      integer :: a, b
      real*8 :: dcut, d
      logical :: needall

      dcut = real(w%extract_nmer_dist(k),8) / bohrtoa
      needall = (w%extract_nmer_any(k) == 0)
      ok = needall
      do a = 1, k-1
         do b = a+1, k
            d = norm2(com(:,ic(a)) - com(:,ic(b)))
            if (needall) then
               if (d > dcut) then
                  ok = .false.
                  return
               end if
            elseif (d <= dcut) then
               ok = .true.
               return
            end if
         end do
      end do

    end function distok

    ! The binomial coefficient in floating point, to size up the
    ! enumeration without overflowing an integer
    function binomial(n,k) result(c)
      integer, intent(in) :: n, k
      real*8 :: c

      integer :: i

      c = 1d0
      do i = 1, k
         c = c * real(n-k+i,8) / real(i,8)
      end do

    end function binomial

    ! Rough number of n-mers of order k the current settings would make (any cutoff).
    function nmer_estimate(k) result(nest)
      integer, intent(in) :: k
      integer :: nest

      real*8 :: moldens, npool, nmain, nnear, c
      integer :: i

      nest = 0
      if (doquit .or. w%extract_atdens <= 0._c_float) return
      if (sys(isys)%c%ncel == 0 .or. sys(isys)%c%nmol == 0) return

      ! molecules per bohr^3, and how many of them the region holds
      moldens = real(w%extract_atdens,8) * sys(isys)%c%nmol / sys(isys)%c%ncel
      npool = region_volume() * moldens
      nmain = min(real(sys(isys)%c%nmoldiscrete,8),npool)
      if (nmain < 1d0) return
      if (k == 1) then
         nest = nint(nmain)
         return
      end if

      ! partners of an anchor: those within the cutoff, at most the pool
      nnear = 4d0*pi/3d0 * (real(w%extract_nmer_dist(k),8) / bohrtoa)**3 * moldens
      nnear = min(nnear,npool) - 1d0
      if (nnear < real(k-1,8)) return

      ! nmain anchors times the ways of choosing k-1 partners
      c = 1d0
      do i = 1, k-1
         c = c * (nnear - real(k-1-i,8)) / real(i,8)
         if (c > 1d9) exit
      end do
      nest = nint(min(nmain * c,1d9))

    end function nmer_estimate

    ! Volume of the chosen region (bohr^3)
    function region_volume() result(v)
      real*8 :: v

      real*8 :: rad0

      v = 0d0
      if (doquit) return
      if (w%extract_region == er_sphere) then
         rad0 = real(w%extract_rsph,8) / bohrtoa
         v = 4d0*pi/3d0 * rad0**3
      elseif (w%extract_region == er_cube) then
         rad0 = real(w%extract_rcub,8) / bohrtoa
         v = 8d0 * rad0**3
      else
         v = sys(isys)%c%omega * product(real(int(w%extract_nx),8))
      end if

    end function region_volume

    ! The largest distance between two points of the chosen region
    ! (bohr): the largest center-of-mass distance a k-mer drawn from it
    ! can possibly have
    function region_diameter() result(d)
      real*8 :: d

      real*8 :: v(3)
      integer :: i2, i3

      d = 0d0
      if (doquit) return
      if (w%extract_region == er_sphere) then
         d = 2d0 * real(w%extract_rsph,8) / bohrtoa
      elseif (w%extract_region == er_cube) then
         d = 2d0 * sqrt(3d0) * real(w%extract_rcub,8) / bohrtoa
      else
         ! the longest body diagonal of the cell block
         ! the eight body diagonals come in +/- pairs of equal length
         do i2 = -1, 1, 2
            do i3 = -1, 1, 2
               v = sys(isys)%c%x2c(real((/1,i2,i3/) * int(w%extract_nx),8))
               d = max(d,norm2(v))
            end do
         end do
      end if

    end function region_diameter

    function em_effective() result(em)
      integer(c_int) :: em
      em = w%extract_mode
      if (.not.severalok) em = em_one
    end function em_effective

    function ei_effective() result(ei)
      integer(c_int) :: ei
      ei = w%extract_include
      if (ei == ei_environ .and. .not.environok) ei = ei_molmotif
      if (ei == ei_molmotif .and. .not.molmotifok) ei = ei_atoms
    end function ei_effective
  end subroutine draw_extract

end submodule extract
