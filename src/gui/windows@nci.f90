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

! The NCI plot window (the NCIPLOT keyword): the reduced density
! gradient and the signed density (the density times the sign of the
! second eigenvalue of its Hessian) of a field on a grid over a
! region, between fragments, loaded as fields of the system.
submodule (windows) nci
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! color and opacity of the region between the fragments (Show region)
  real(c_float), parameter :: nci_region_rgb(3) = (/0.55_c_float,0.75_c_float,1.0_c_float/)
  real(c_float), parameter :: nci_region_alpha = 0.25_c_float

  ! default buffer of the bounding-box region (molecules), angstrom;
  ! tighter than iso_bbox_buffer_def because the NCI regions are
  ! short-range and the grid points far from the atoms are wasted
  real*8, parameter :: nci_bbox_buffer = 2d0

  ! color bins of the points of the s vs. sign(lambda2)*rho plot
  integer, parameter :: nscbin = 64

contains

  !> Draw the NCI plot window: the density (a field of the system
  !> shown by the anchor view), the region and the grid over it, the
  !> fragments between which the plot is made, and the cutoffs of the
  !> plot. Calculate runs the NCIPLOT kernel as a blocking job
  !> (run_nci) and loads the grids as two new fields of the system.
  module subroutine draw_nci(w)
    use systems, only: sys, sysc, sys_init, ok_system, field_cost, atlisttype_ncel_bohr
    use representations, only: iso_region_cell, iso_region_bbox, iso_region_to_box, iso_region_seed, iso_rgb_hl,&
       iso_region_point_from_cart, iso_isgridfield, iso_region_parallel
    use nci, only: nci_options, nci_options_set_field, nci_frag_region
    use types, only: field_evaluation_avail
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_setpos_bottomright,&
       iw_separatortext, iw_combo_simple, iw_dragfloat_real8, iw_checkbox, iw_field_combo,&
       iw_calcwidth, iw_calcheight, iw_table_column, iw_table_headers_row, iw_highlight_selectable,&
       iw_cell_right, iw_close_button, duration_string, file_name_root, file_name_base,&
       iw_colormap_lut
    use gui_main, only: ColorHighlightScene, g
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer, parameter :: itable_domains = 1 ! the domain table, for table_write_*
    ! colors of the lines of the scatter plot and of the labels of their values
    real(c_float), parameter :: rgb_range(3) = (/0.3_c_float,0.6_c_float,1._c_float/)
    real(c_float), parameter :: rgb_isoval(3) = (/0.3_c_float,0.8_c_float,0.4_c_float/)
    real(c_float), parameter :: rgb_cscale(3) = (/0.95_c_float,0.45_c_float,0.85_c_float/)
    ! limits of the isovalue and of the sign(lambda2)*rho range
    real*8, parameter :: dimplot_min = 0.001d0, dimplot_max = 10d0
    real*8, parameter :: srhorange_max = 10d0

    logical :: doquit, goodsys, syschanged, ismol, isgrid, okbox, okfield, okcalc, ch, navail, capped, ok
    logical :: maskchanged, objchanged, lcommit
    integer :: isys, iview, istat, irow, ifield, nn(3), ihover, nmolsel
    real*8 :: xc(3), box(3,0:3), cscale(2)
    real(c_float) :: xbutton

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)

    ! initialize the state on the first pass
    if (w%firstpass) then
       w%nc = nci_state()
       w%errmsg = ""
    end if

    ! a new system: the reference field, on the cell (crystals) or
    ! around the atoms (molecules), no fragments; the same system with
    ! a new geometry: the fragments and the results refer to the old one
    if (goodsys) then
       ismol = sys(isys)%c%ismolecule
       if (.not.w%nc%sys%is(isys) .or. syschanged) then
          call w%nc%sys%set(isys)
          call reset_system()
          call drop_results()
       elseif (w%nc%timelast /= sysc(isys)%timelastchange_geometry) then
          call clear_fragments()
          call drop_results()
       end if

       ! the reference field, if the field is gone
       if (.not.w%nc%fref%ok(sys(isys))) call set_field(sys(isys)%iref)

       ! the molecule selection follows the molecule list; by default,
       ! every molecule of a molecular system is a fragment if there is
       ! more than one (not in a crystal: the fragment densities have
       ! no lattice sum, so the criterion would fail near the cell faces)
       ! (with one molecule, there are no molecule rows: size zero)
       nmolsel = merge(sys(isys)%c%nmol,0,sys(isys)%c%nmol > 1)
       if (allocated(w%nc%molsel)) then
          if (size(w%nc%molsel,1) /= nmolsel) deallocate(w%nc%molsel)
       end if
       if (.not.allocated(w%nc%molsel)) then
          allocate(w%nc%molsel(nmolsel))
          w%nc%molsel = ismol
       end if
    end if

    ! a pending pick of a region point in the view and the region mode
    ! it was armed on; dropped, releasing the view, if the system is
    ! not ready, the window moved to another view, or the region mode
    ! changed
    if (w%nc%picking) then
       if (goodsys .and. iview == w%nc%pickview .and. w%nc%iregion == w%nc%pickmode) then
          call view_pick_result(iview,w%id,isys,w%nc%pick,istat,xc)
          if (istat == ipick_point) &
             call iso_region_point_from_cart(isys,w%nc%iregion,xc,w%nc%rgn_x(:,w%nc%pickrow))
          if (istat /= ipick_pending) w%nc%picking = .false.
       else
          call cancel_pick()
       end if
    end if

    okbox = .false.
    okfield = .false.
    nn = 0
    if (goodsys) then
       ! the system
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! the density
       call draw_density()
       ifield = w%nc%fref%id
       okfield = sys(isys)%goodfield(ifield)
       isgrid = iso_isgridfield(isys,ifield)
       if (.not.isgrid) w%nc%opts%fourierint = .false.

       ! the region
       if (region_form(isys,w%nc%iregion,w%nc%rgn_x,w%nc%sec(nci_sec_region),"nci","the NCI plot",&
          ttshown,irow,pickoff=(w%nc%picking),bboxbuf=nci_bbox_buffer)) w%errmsg = ""
       if (irow >= 0) then
          w%nc%picking = .true.
          w%nc%pickrow = irow
          w%nc%pickmode = w%nc%iregion
          w%nc%pickview = iview
          call w%nc%pick%arm()
          call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the region point: an atom,&
             & a critical point, or a point in space",w%id,acceptempty=.true.)
       end if
       call iso_region_to_box(isys,w%nc%iregion,w%nc%rgn_x,box,okbox)
       if (.not.okbox) &
          call iw_text("The region is degenerate (zero volume)",danger=.true.,wrap=.true.)

       ! the grid: the native grid of a grid field that spans the whole
       ! cell (a partial grid does not cover the cell it is placed on)
       navail = isgrid .and. (w%nc%iregion == iso_region_cell)
       if (navail) navail = .not.sys(isys)%f(ifield)%grid%partial
       ! (leaving the native grid: the NCI default, not the isosurface one)
       if (.not.navail .and. w%nc%ilevel == 0) w%nc%ilevel = nci_level_def
       ch = grid_level_form(isys,ifield,w%nc%iregion,box,okbox,navail,w%nc%ilevel,w%nc%icustom,&
          w%nc%nptscustom,w%nc%ptsangcustom,w%nc%sec(nci_sec_grid),"nci","the NCI plot",ttshown,nn,capped)
       if (w%nc%ilevel == 0) nn = sys(isys)%f(ifield)%grid%n
       if (w%nc%sec(nci_sec_grid)%open) call draw_cost()

       ! the fragments and the cutoffs; a committed change of a mask
       ! cutoff masks the result again and updates the field and the
       ! objects, one of the plot settings only the objects
       maskchanged = .false.
       objchanged = .false.
       call draw_fragments()
       call draw_plot()

       ! the domains (the row under the mouse is highlighted on the
       ! isosurface), then redo what the changes in the Plot section
       ! require (the cutoffs and the range lines of the scatter plot)
       ihover = 0
       if (associated(w%nc%res)) then
          call draw_domains()
          if (maskchanged) then
             call nci_refresh(w,isys,iview)
          else
             if (objchanged) then
                if (w%nc%fsref%ok(sys(isys))) &
                   call nci_find_domains(w,sys(isys)%f(w%nc%fsref%id)%grid%f)
                call nci_objects(w,isys,iview)
             end if
          end if
          call nci_highlight(w,isys,iview,ihover)
       end if
       ! (for the scatter plot of the next frame, drawn above the table)
       w%nc%hoverdom = ihover

       ! the region between the fragments
       call draw_region()

       ! the output files
       call draw_output()
       if (w%okfile_set) then
          if (associated(w%nc%res)) then
             call request_write()
          else
             w%errmsg = "There is no result to write"
          end if
       end if

       ! the region in the anchor view, for as long as the window is open
       ! (not the cell of a crystal, which the cell shows)
       if (okbox .and. .not.(w%nc%iregion == iso_region_cell .and..not.ismol)) &
          call region_show_box(iview,w%uid,isys,box(:,0),box(:,1:3))

       if (len_trim(w%errmsg) > 0) then
          call iw_text(w%errmsg,danger=.true.,wrap=.true.)
       elseif (len_trim(w%okmsg) > 0) then
          call iw_text(w%okmsg,wrap=.true.)
       end if
    end if
    w%okfile_set = .false.

    ! bottom row: reset, calculate, close; on the left, whether the
    ! result is still that of the form
    okcalc = goodsys
    if (okcalc) okcalc = okbox .and. okfield .and. all(nn > 0) .and. .not.w%nc%picking
    call iw_setpos_bottomright(len("Reset") + len("Calculate") + len("Close"),3)
    if (okcalc) then
       if (result_stale()) then
          xbutton = igGetCursorPosX()
          call igSetCursorPosX(g%Style%WindowPadding%x)
          call iw_text("settings changed",danger=.true.,alignframe=.true.)
          call iw_tooltip("The density, the derivatives, the region, the grid, the fragments,&
             & or the voids changed since the last calculation: press Calculate to update the plot",&
             ttshown)
          call igSameLine(xbutton,-1._c_float)
       end if
    end if
    if (iw_button("Reset##ncireset",danger=.true.,disabled=.not.goodsys)) then
       call reset_form()
       if (associated(w%nc%res)) then
          call nci_refresh(w,isys,iview)
          ok = nci_cscale(w,cscale,reset=.true.)
       end if
    end if
    call iw_tooltip("Reset all the settings of the window to their defaults (the result of the&
       & last calculation is kept)",ttshown,whendisabled=.true.)
    if (iw_button("Calculate##ncicalc",danger=.true.,disabled=.not.okcalc,sameline=.true.)) &
       call request_job()
    call iw_tooltip("Calculate the reduced density gradient and the signed density on the grid,&
       & and load them as fields of the system",ttshown,whendisabled=.true.)
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) then
       call cancel_pick()
       call w%end()
    end if

  contains

    !> Use field id of the system as the density: the defaults of
    !> NCIPLOT that depend on it (nci_options_set_field), the signed
    !> density range following the new default density cutoff, and the
    !> native grid of a grid field on the whole cell (the CLI default).
    !> The window has no density cutoff (the range does its job): the
    !> one in the options is the default of the field, which gives the
    !> default range and the initial color scale of the isosurface.
    subroutine set_field(id)
      integer, intent(in) :: id

      integer :: id0
      real*8 :: rhoprev

      id0 = id
      if (.not.sys(isys)%goodfield(id0)) id0 = 0
      call w%nc%fref%set(sys(isys),id0)
      rhoprev = w%nc%opts%rhoplot
      call nci_options_set_field(w%nc%opts,sys(isys),id0)
      ! a bound of the range at or beyond the default of the previous
      ! field (the default range, or none at all) moves to that of the
      ! new one, unless that would invert a range the user set beyond it
      associate(o => w%nc%opts)
        if (o%srhorange(1) <= -rhoprev .and. -o%rhoplot <= o%srhorange(2)) o%srhorange(1) = -o%rhoplot
        if (o%srhorange(2) >= rhoprev .and. o%rhoplot >= o%srhorange(1)) o%srhorange(2) = o%rhoplot
      end associate
      if (iso_isgridfield(isys,id0) .and. w%nc%iregion == iso_region_cell) then
         if (.not.sys(isys)%f(id0)%grid%partial) w%nc%ilevel = 0
      end if

    end subroutine set_field

    !> Set every setting of the form to its default for system isys:
    !> the default cutoffs, grid and output options, and those of
    !> reset_system. The results of the last calculation and the state
    !> of the sections are kept, and so are the grid and the fragments
    !> of the options, which describe the result (filled by request_job).
    subroutine reset_form()

      type(nci_state) :: d
      type(nci_options) :: o0

      o0 = w%nc%opts
      w%nc%opts = nci_options()
      if (associated(w%nc%res)) then
         w%nc%opts%periodic = o0%periodic
         w%nc%opts%x0 = o0%x0
         w%nc%opts%xmat = o0%xmat
         w%nc%opts%nstep = o0%nstep
         w%nc%opts%nfrag = o0%nfrag
         w%nc%opts%frcell = o0%frcell
         call move_alloc(o0%fr,w%nc%opts%fr)
      end if
      w%nc%ilevel = d%ilevel
      w%nc%icustom = d%icustom
      w%nc%nptscustom = d%nptscustom
      w%nc%ptsangcustom = d%ptsangcustom
      w%nc%molmotif = d%molmotif
      call reset_system()

    end subroutine reset_form

    !> Set the settings of the form that depend on system isys to their
    !> defaults: the reference field, the cell (crystals) or the box
    !> around the atoms (molecules), no fragments, and the default
    !> output file. Cancels a pending pick and clears the messages.
    subroutine reset_system()

      call cancel_pick()
      w%errmsg = ""
      w%okmsg = ""
      if (sys(isys)%c%ismolecule) then
         w%nc%iregion = iso_region_bbox
      else
         w%nc%iregion = iso_region_cell
      end if
      call iso_region_seed(isys,w%nc%iregion,w%nc%rgn_x,bboxbuf=nci_bbox_buffer)
      call set_field(sys(isys)%iref)
      call clear_fragments()
      w%okfile = okfile_default(isys,"nci","dat")

    end subroutine reset_system

    !> Forget the fragments: the molecule selection and the user-made
    !> fragments.
    subroutine clear_fragments()

      if (allocated(w%nc%molsel)) deallocate(w%nc%molsel)
      if (allocated(w%nc%userfr)) deallocate(w%nc%userfr)

    end subroutine clear_fragments

    !> Forget the results of the last calculation (the fields it
    !> loaded stay in the system).
    subroutine drop_results()

      call w%end_nci()
      w%nc%dragline = 0
      w%nc%fsref = field_ref()
      w%nc%frref = field_ref()
      w%nc%isoref = rep_ref()
      w%nc%cbref = rep_ref()
      w%nc%timelast = sysc(isys)%timelastchange_geometry
      w%errmsg = ""

    end subroutine drop_results

    !> The number of fragments the form defines: the selected
    !> molecules and the selected user-made fragments.
    function nfrag_total() result(nfrag)
      integer :: nfrag

      nfrag = 0
      if (allocated(w%nc%molsel)) nfrag = count(w%nc%molsel)
      if (allocated(w%nc%userfr)) nfrag = nfrag + count(w%nc%userfr%sel)

    end function nfrag_total

    !> The Density section: the field, and (More) how a grid field is
    !> differentiated.
    subroutine draw_density()
      integer :: id, ic
      logical :: isgridf

      isgridf = iso_isgridfield(isys,w%nc%fref%id)
      call iw_separatortext("Density",w%nc%sec(nci_sec_density),hasmore=isgridf)
      if (.not.w%nc%sec(nci_sec_density)%open) return

      call iw_text("Field",alignframe=.true.)
      id = w%nc%fref%id
      if (iw_field_combo("##ncifield",isys,id,width=iw_calcwidth(40,0),sameline=.true.)) then
         call set_field(id)
         w%errmsg = ""
      end if
      call iw_tooltip("Field whose density the NCI plot is calculated from. Field 0 is the&
         & promolecular density (the sum of the in-vacuo atomic densities), as in NCIPLOT&
         & without a wavefunction",ttshown)

      if (w%nc%sec(nci_sec_density)%more .and. isgridf) then
         call iw_text("Derivatives",alignframe=.true.)
         ic = merge(1,0,w%nc%opts%fourierint)
         call iw_combo_simple("##nciderivs","Interpolant" // c_null_char // "Fourier transform" //&
            c_null_char,ic,sameline=.true.,changed=ch)
         if (ch) w%nc%opts%fourierint = (ic == 1)
         call iw_tooltip("How the gradient and the Hessian of a grid field are calculated: by&
            & differentiating the interpolant (INTERPOLATION in NCIPLOT) or by Fourier transform&
            & of the grid (the default for grids without core densities)",ttshown)
      end if

      ! derivatives of a grid field from its interpolant: unreliable
      ! unless smoothrho (shown also with the Derivatives combo hidden)
      if (.not.w%nc%opts%fourierint) &
         call warn_grid_interp(isys,w%nc%fref%id,"derivatives (gradient and Hessian)")

    end subroutine draw_density

    !> The Estimate Cost row of the Grid section: the cost of the
    !> density, its gradient and its Hessian (kept with the field,
    !> shared by all windows); the Fourier transform of a grid has no
    !> per-point cost.
    subroutine draw_cost()
      logical :: measure
      real*8 :: secs
      type(field_evaluation_avail) :: request

      measure = iw_button("Estimate Cost##nciestimate",disabled=.not.(okbox .and. okfield) .or.&
         w%nc%opts%fourierint)
      call iw_tooltip("Measure the cost of the density and its derivatives, to estimate the time the&
         & grid takes (the measurement is kept with the field and shared by all windows; not for&
         & derivatives by Fourier transform)",ttshown,whendisabled=.true.)
      if (.not.okfield .or. w%nc%opts%fourierint) return
      if (measure) w%errmsg = ""
      call request%field_nder2()
      secs = field_cost(isys,ifield,measure,request)
      if (measure .and. secs < 0d0) w%errmsg = "Could not measure the cost of the field"
      if (secs >= 0d0 .and. all(nn > 0)) &
         call iw_text("~" // duration_string(secs * real(product(int(nn,8)),8)),sameline=.true.,&
            disabled=.true.)

    end subroutine draw_cost

    !> The Fragments section: a table with the molecules of the system
    !> (if more than one) and the fragments the user made from the
    !> atoms selected in the view, each selected or not as a fragment
    !> of the plot (the user-made ones can be removed), and, with any
    !> fragment selected, the fragment criterion parameters.
    subroutine draw_fragments()
      integer :: i, icol, nmolrow, nrow, irm, nfrag
      integer, allocatable :: idx(:)
      integer(c_int) :: flags
      logical :: havesel
      type(fragment), pointer :: fr
      logical, pointer :: sel
      character(kind=c_char,len=:), allocatable, target :: str1
      type(ImVec2) :: sz0
      type(c_ptr), target :: clipper
      type(ImGuiListClipper), pointer :: clipper_f

      call iw_separatortext("Fragments",w%nc%sec(nci_sec_fragments))
      if (.not.w%nc%sec(nci_sec_fragments)%open) return

      ! what is plotted: the fragment criterion needs two fragments
      nfrag = nfrag_total()
      if (nfrag > 1) then
         call iw_text("Only regions between fragments are plotted",disabled=.true.,wrap=.true.)
      elseif (nfrag == 1) then
         call iw_text("All NCI domains are plotted (select at least two fragments)",&
            disabled=.true.,wrap=.true.)
      else
         call iw_text("All NCI domains are plotted",disabled=.true.,wrap=.true.)
      end if

      ! the table: the molecules (if more than one), then the
      ! user-made fragments; a removal is applied after the table so
      ! the rows keep their IDs for the frame
      nmolrow = 0
      if (allocated(w%nc%molsel)) nmolrow = size(w%nc%molsel,1)
      nrow = nmolrow
      if (allocated(w%nc%userfr)) nrow = nrow + size(w%nc%userfr,1)
      irm = 0
      if (nrow > 0) then
         flags = ImGuiTableFlags_None
         flags = ior(flags,ImGuiTableFlags_Resizable)
         flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
         flags = ior(flags,ImGuiTableFlags_Borders)
         flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
         flags = ior(flags,ImGuiTableFlags_ScrollY)
         str1 = "##ncifragments" // c_null_char
         sz0%x = 0
         sz0%y = iw_calcheight(min(5,nrow)+1,0,.false.)
         if (igBeginTable(c_loc(str1),5,flags,sz0,0._c_float)) then
            icol = -1
            call iw_table_column("(remove)",icol=icol,flags=ImGuiTableColumnFlags_NoHeaderLabel)
            call iw_table_column("Id",icol=icol)
            call iw_table_column("Nat",icol=icol)
            call iw_table_column("Select",icol=icol)
            call iw_table_column(cmass_header(isys),icol=icol,flags=ImGuiTableColumnFlags_WidthStretch)
            call iw_table_headers_row(freezetop=.true.,autofit=.true.)

            clipper = ImGuiListClipper_ImGuiListClipper()
            call ImGuiListClipper_Begin(clipper,nrow,-1._c_float)
            do while(ImGuiListClipper_Step(clipper))
               call c_f_pointer(clipper,clipper_f)
               do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                  call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)

                  ! remove (user-made fragments only); the close icon is one
                  ! line of text high: centered on the row frames
                  if (i > nmolrow .and. igTableSetColumnIndex(0)) then
                     call igSetCursorPosY(igGetCursorPosY() + g%Style%FramePadding%y)
                     if (iw_close_button("##ncifragrm" // string(i))) irm = i - nmolrow
                     call iw_tooltip("Remove this fragment",ttshown)
                  end if

                  ! the fragment of the row and its selection
                  if (i <= nmolrow) then
                     fr => sys(isys)%c%mol(i)
                     sel => w%nc%molsel(i)
                  else
                     fr => w%nc%userfr(i-nmolrow)%fr
                     sel => w%nc%userfr(i-nmolrow)%sel
                  end if

                  ! id, with the row highlight: the atoms of the fragment
                  if (igTableSetColumnIndex(1)) then
                     call iw_text(string(i),alignframe=.true.)
                     if (iw_highlight_selectable("##ncifragrow_" // string(i))) then
                        ! (copied: the strided component array would be a temporary)
                        idx = fr%at(1:fr%nat)%cidx
                        call sysc(isys)%highlight_atoms(.true.,idx,atlisttype_ncel_bohr,&
                           spread(ColorHighlightScene,2,fr%nat))
                     end if
                  end if

                  ! number of atoms, whether it is a fragment of the plot,
                  ! center of mass
                  if (igTableSetColumnIndex(2)) call iw_text(string(fr%nat),alignframe=.true.)
                  if (igTableSetColumnIndex(3)) then
                     ch = iw_checkbox("##ncifragsel" // string(i),sel)
                     call iw_tooltip("Use this as a fragment of the plot",ttshown)
                  end if
                  if (igTableSetColumnIndex(4)) &
                     call iw_text(cmass_string(isys,fr%cmass(.false.)),alignframe=.true.)
               end do
            end do
            call ImGuiListClipper_End(clipper)
            call ImGuiListClipper_destroy(clipper)
            call igEndTable()
         end if
      end if
      if (irm > 0) w%nc%userfr = (/w%nc%userfr(1:irm-1),w%nc%userfr(irm+1:)/)

      ! select all, clear, add from the selection in the view
      if (nrow > 0) then
         if (iw_button("All##ncifragall")) call select_all_fragments(1)
         call iw_tooltip("Use every row of the table as a fragment",ttshown)
         if (iw_button("None##ncifragnone",sameline=.true.)) call select_all_fragments(0)
         call iw_tooltip("Use no row of the table as a fragment",ttshown)
         if (iw_button("Toggle##ncifragtoggle",sameline=.true.)) call select_all_fragments(-1)
         call iw_tooltip("Invert the selection of the rows of the table",ttshown)
      end if
      havesel = .false.
      if (allocated(sysc(isys)%highlight_rgba)) havesel = any(sysc(isys)%highlight_rgba >= 0._c_float)
      if (iw_button("Add from Selection##ncifragadd",sameline=(nrow > 0),disabled=.not.havesel)) &
         call add_selfrag()
      call iw_tooltip("Add the atoms selected in the view, together, as a new fragment",ttshown,&
         whendisabled=.true.)

      ! the fragment criterion, with two or more fragments selected
      ! (after the table, so a checkbox changed this frame counts)
      if (nfrag_total() > 1) then
         call iw_text("Dominance",alignframe=.true.)
         ch = iw_dragfloat_real8("##ncirhoparam",x1=w%nc%opts%rhoparam,speed=0.005d0,decimal=3,&
            min=0d0,max=1d0,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,notlive=.true.,&
            committed=lcommit)
         maskchanged = maskchanged .or. lcommit
         call iw_tooltip("A point is between fragments if no single fragment contributes more&
            & than this fraction of the promolecular density of all the fragments (RHOPARAM)",ttshown)
         call iw_text("Share",alignframe=.true.,gap=.true.)
         ch = iw_dragfloat_real8("##ncirhoparam2",x1=w%nc%opts%rhoparam2,speed=0.005d0,decimal=3,&
            min=0d0,max=1d0,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,notlive=.true.,&
            committed=lcommit)
         maskchanged = maskchanged .or. lcommit
         call iw_tooltip("and the fragments together contribute at least this fraction of the&
            & promolecular density of the whole system (RHOPARAM2)",ttshown)
         ch = iw_checkbox("Show region##ncishowregion",w%nc%showregion,labelleft=.true.)
         call iw_tooltip("Show the region between the fragments, where the NCI plot is drawn, as a&
            & translucent surface in the view (calculates the promolecular densities of the&
            & fragments on the grid, which Calculate then reuses)",ttshown)
      end if

    end subroutine draw_fragments

    !> Select (mode = 1), deselect (0) or invert the selection (-1) of
    !> every row of the fragment table.
    subroutine select_all_fragments(mode)
      integer, intent(in) :: mode

      if (allocated(w%nc%molsel)) w%nc%molsel = merge(.not.w%nc%molsel,mode == 1,mode < 0)
      if (allocated(w%nc%userfr)) w%nc%userfr%sel = merge(.not.w%nc%userfr%sel,mode == 1,mode < 0)

    end subroutine select_all_fragments

    !> The Plot section: the s vs. signed density plot of the result,
    !> the cutoffs of the isosurface plot (RDG isovalue, sign(lambda2)*rho
    !> range) and, More, the void criterion. The density cutoff
    !> (CUTPLOT) and ONLYNEG of the CLI are ranges too, so the range
    !> does their job; the color scale belongs to the isosurface
    !> object. The dat-file cutoffs (CUTOFFS) stay at their defaults:
    !> the plot zooms instead.
    subroutine draw_plot()

      call iw_separatortext("Plot",w%nc%sec(nci_sec_plot),hasmore=.true.)
      if (.not.w%nc%sec(nci_sec_plot)%open) return
      call draw_scatter()

      associate(o => w%nc%opts)
        call iw_text("Isovalue (s)",alignframe=.true.,rgb=rgb_isoval)
        ch = iw_dragfloat_real8("##ncidimplot",x1=o%dimplot,speed=0.005d0,decimal=3,min=dimplot_min,&
           max=dimplot_max,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,notlive=.true.,committed=lcommit)
        objchanged = objchanged .or. lcommit
        call iw_tooltip("Value of the reduced density gradient (s) of the isosurface (DIMPLOT)",ttshown)
        call iw_text("sign(λ2)ρ range (a.u.)",alignframe=.true.,gap=.true.,rgb=rgb_range)
        ch = iw_dragfloat_real8("##ncisrhorange",x2=o%srhorange,speed=0.001d0,decimal=4,&
           min=-srhorange_max,max=srhorange_max,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,&
           notlive=.true.,committed=lcommit)
        if (o%srhorange(2) < o%srhorange(1)) o%srhorange(2) = o%srhorange(1)
        maskchanged = maskchanged .or. lcommit
        call iw_tooltip("Points with signed density (the density times the sign of the second&
           & eigenvalue of its Hessian) outside this range are left out of the isosurface&
           & (SRHORANGE). An upper bound of zero keeps only the attractive interactions.",ttshown)

        ! the color scale of the isosurface, with the plot (its values
        ! are those of the isosurface object and its color bar)
        if (associated(w%nc%res) .and. associated(w%nc%drv)) then
           if (nci_cscale(w,cscale)) then
              call iw_text("Color scale (a.u.)",alignframe=.true.,rgb=rgb_cscale)
              ch = iw_dragfloat_real8("##ncicscale",x2=cscale,speed=0.001d0,decimal=4,&
                 min=-srhorange_max,max=srhorange_max,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
              if (ch) then
                 if (cscale(2) < cscale(1)) cscale(2) = cscale(1)
                 ok = nci_cscale(w,cscale,set=.true.)
              end if
              call iw_tooltip("Values of sign(λ2)ρ at the two ends of the color scale of the&
                 & isosurface (blue to red) and its color bar",ttshown)
           end if
        end if

        if (w%nc%sec(nci_sec_plot)%more) then
           ch = iw_checkbox("Voids only##ncidovoid",o%dovoid,labelleft=.true.)
           maskchanged = maskchanged .or. ch
           call iw_tooltip("Plot only the regions where the promolecular density is below a&
              & value (VOID)",ttshown)
           if (o%dovoid) then
              call iw_text("Void density (a.u.)",alignframe=.true.,gap=.true.)
              ch = iw_dragfloat_real8("##ncirhovoid",x1=o%rho_void,speed=0.0005d0,decimal=4,&
                 min=0.0001d0,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,notlive=.true.,&
                 committed=lcommit)
              maskchanged = maskchanged .or. lcommit
              call iw_tooltip("Promolecular density below which a point belongs to a void",ttshown)
           end if
        end if
      end associate

    end subroutine draw_plot

    !> The s vs. signed density plot of the result (top of the Plot
    !> section): the points under the dat-file cutoffs, with draggable
    !> lines for the sign(lambda2)*rho range (two), the isovalue (one,
    !> horizontal) and the ends of the color scale of the isosurface
    !> (two). Letting go of a line masks the result again (range) or
    !> redoes the objects (isovalue); the color scale changes as its
    !> lines move.
    subroutine draw_scatter()
      integer :: i, held, icmap, nb
      real(c_double) :: xdrag
      real(c_float) :: lut(3,nscbin)
      real(c_double), allocatable, target :: hx(:), hy(:)
      character(kind=c_char,len=:), allocatable, target :: str1, str2
      type(ImVec2) :: szplot
      type(ImVec4) :: col

      ! a line let go while the plot is not drawn lands nowhere
      if (.not.associated(w%nc%res) .or. .not.associated(w%nc%drv)) then
         w%nc%dragline = 0
         return
      end if

      associate(o => w%nc%opts, d => w%nc%drv)
        str1 = "##nciscatter" // c_null_char
        szplot%x = -1._c_float
        szplot%y = iw_calcheight(1,10,.false.)
        held = 0
        if (ipBeginPlot(c_loc(str1),szplot,ior(ImPlotFlags_NoTitle,ImPlotFlags_NoLegend))) then
           str1 = "sign(λ2)ρ (a.u.)" // c_null_char
           str2 = "s" // c_null_char
           call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_None,ImPlotAxisFlags_None)
           ! the axes span the points (the cutoffs of the dat file of
           ! the result) when a result arrives; the user can zoom and
           ! pan them afterwards
           if (w%nc%fitplot) then
              xdrag = real(w%nc%datcut(1),c_double)
              call ipSetupAxisLimits(ImAxis_X1,-xdrag,xdrag,ImPlotCond_Always)
              xdrag = real(w%nc%datcut(2),c_double)
              call ipSetupAxisLimits(ImAxis_Y1,0._c_double,xdrag,ImPlotCond_Always)
              w%nc%fitplot = .false.
           end if

           ! the points, colored by their signed density as the
           ! isosurface (its colormap and color scale), gray outside the
           ! sign(lambda2)*rho range (drawn first, under the others):
           ! one scatter per color bin, as implot has no per-point colors
           if (allocated(d%scx)) then
              if (size(d%scx,1) > 0) then
                 ok = nci_cscale(w,cscale,icmap=icmap)
                 call nci_scatter_bins(d,cscale,o%srhorange)
                 call iw_colormap_lut(icmap,lut)
                 do i = 0, nscbin
                    nb = d%scbin(i) - d%scbin(i-1)
                    if (nb == 0) cycle
                    str1 = "##ncipoints" // string(i) // c_null_char
                    if (i == 0 .or. .not.ok) then
                       ! (all gray if the isosurface is not colored by
                       ! the signed density of the window)
                       col = ImVec4(0.55_c_float,0.55_c_float,0.55_c_float,0.6_c_float)
                    else
                       col = ImVec4(lut(1,i),lut(2,i),lut(3,i),0.7_c_float)
                    end if
                    call ipSetNextMarkerStyle(ImPlotMarker_Square,1.5_c_float,col,0._c_float,col)
                    call ipPlotScatter(c_loc(str1),c_loc(d%scx(d%scbin(i-1)+1)),c_loc(d%scy(d%scbin(i-1)+1)),&
                       int(nb,c_int),0_c_int,0_c_int)
                 end do

                 ! the points of the domain under the mouse in the
                 ! domain table, on top, in the highlight color of the
                 ! isosurface
                 if (w%nc%hoverdom > 0 .and. allocated(d%scl)) then
                    hx = pack(d%scx,d%scl == w%nc%hoverdom)
                    hy = pack(d%scy,d%scl == w%nc%hoverdom)
                    if (size(hx) > 0) then
                       str1 = "##ncipointshl" // c_null_char
                       col = ImVec4(iso_rgb_hl(1),iso_rgb_hl(2),iso_rgb_hl(3),1._c_float)
                       call ipSetNextMarkerStyle(ImPlotMarker_Square,2.5_c_float,col,0._c_float,col)
                       call ipPlotScatter(c_loc(str1),c_loc(hx),c_loc(hy),int(size(hx),c_int),0_c_int,0_c_int)
                    end if
                 end if
              end if
           end if

           ! the sign(lambda2)*rho range, the isovalue and the color
           ! scale, inside the span of the points; not counted in a fit
           ! (double click), which is to the points. The range lines
           ! are submitted first, so they take the click where they
           ! sit on the color scale lines (the default), and drawn
           ! wider, so both show there.
           col = ImVec4(rgb_range(1),rgb_range(2),rgb_range(3),1._c_float)
           do i = 1, 2
              xdrag = real(o%srhorange(i),c_double)
              if (logical(ipDragLineX(int(i-1,c_int),xdrag,col,5._c_float,ImPlotDragToolFlags_NoFit))) then
                 held = i
                 o%srhorange(i) = max(min(real(xdrag,8),w%nc%datcut(1)),-w%nc%datcut(1))
              end if
           end do

           ! the color scale
           if (nci_cscale(w,cscale)) then
              col = ImVec4(rgb_cscale(1),rgb_cscale(2),rgb_cscale(3),1._c_float)
              do i = 1, 2
                 xdrag = real(cscale(i),c_double)
                 if (logical(ipDragLineX(int(i+2,c_int),xdrag,col,2._c_float,ImPlotDragToolFlags_NoFit))) then
                    held = i + 3
                    cscale(i) = max(min(real(xdrag,8),w%nc%datcut(1)),-w%nc%datcut(1))
                    if (i == 1) cscale(2) = max(cscale(2),cscale(1))
                    if (i == 2) cscale(1) = min(cscale(1),cscale(2))
                    ok = nci_cscale(w,cscale,set=.true.)
                 end if
              end do
           end if

           ! the isovalue
           col = ImVec4(rgb_isoval(1),rgb_isoval(2),rgb_isoval(3),1._c_float)
           xdrag = real(o%dimplot,c_double)
           if (logical(ipDragLineY(2_c_int,xdrag,col,2._c_float,ImPlotDragToolFlags_NoFit))) then
              held = 3
              o%dimplot = max(min(real(xdrag,8),w%nc%datcut(2),dimplot_max),dimplot_min)
           end if
           if (o%srhorange(2) < o%srhorange(1)) o%srhorange(2) = o%srhorange(1)
           call ipEndPlot()
        end if
        call iw_tooltip("Reduced density gradient against the signed density at the grid&
           & points under the default cutoffs of the dat file of NCIPLOT, colored as the&
           & isosurface (gray outside the sign(λ2)ρ range). Drag the lines to&
           & change the isovalue (green), the sign(λ2)ρ range (blue) and the color scale&
           & (pink) of the isosurface (where a pink line is under a blue one, move the blue&
           & one or use the field); the mouse wheel zooms, a drag pans, a right-drag zooms&
           & into a box and a double click fits the points",ttshown)
        ! the values land when their line is let go (the color scale
        ! has been set while it moved)
        if (w%nc%dragline == 3 .and. held == 0) then
           objchanged = .true.
        elseif ((w%nc%dragline == 1 .or. w%nc%dragline == 2) .and. held == 0) then
           maskchanged = .true.
        end if
        w%nc%dragline = held
        if (allocated(d%scx)) call iw_text(string(size(d%scx,1)) // " points shown",disabled=.true.)
      end associate

    end subroutine draw_scatter

    !> The Domains section: the connected regions of the grid under
    !> the RDG isovalue (the pieces of the isosurface plot), by
    !> decreasing volume, with their signed density range; the row
    !> under the mouse is highlighted on the isosurface.
    subroutine draw_domains()
      integer :: i, icol, nrow, iwrite
      integer(c_int) :: flags
      character(kind=c_char,len=:), allocatable, target :: str1
      type(ImVec2) :: sz0
      type(c_ptr), target :: clipper
      type(ImGuiListClipper), pointer :: clipper_f

      if (.not.associated(w%nc%drv)) return
      call iw_separatortext("Domains (" // string(w%nc%drv%ndom) // ")",w%nc%sec(nci_sec_domains))
      if (.not.w%nc%sec(nci_sec_domains)%open) return
      associate(d => w%nc%drv)
        if (d%ndom == 0) then
           call iw_text("No points of the grid are under the RDG isovalue",disabled=.true.,wrap=.true.)
           return
        end if
        call iw_text("Total volume (Å³): " // string(d%voltot * bohrtoa**3,'f',decimal=3),wrap=.true.)

        flags = ImGuiTableFlags_None
        flags = ior(flags,ImGuiTableFlags_Resizable)
        flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
        flags = ior(flags,ImGuiTableFlags_Borders)
        flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
        flags = ior(flags,ImGuiTableFlags_ScrollY)
        str1 = "##ncidomains" // c_null_char
        nrow = min(8,d%ndom)
        sz0%x = 0
        sz0%y = iw_calcheight(nrow+1,0,.false.)
        if (igBeginTable(c_loc(str1),6,flags,sz0,0._c_float)) then
           icol = -1
           call iw_table_column("Id",icol=icol)
           call iw_table_column("Volume (Å³)",icol=icol)
           call iw_table_column("Electrons",icol=icol)
           call iw_table_column("Points",icol=icol)
           call iw_table_column("Min sign(λ2)ρ",icol=icol)
           call iw_table_column("Max sign(λ2)ρ",icol=icol,flags=ImGuiTableColumnFlags_WidthStretch)
           call iw_table_headers_row(freezetop=.true.,autofit=.true.,writemenu=iwrite,ttshown=ttshown)

           clipper = ImGuiListClipper_ImGuiListClipper()
           call ImGuiListClipper_Begin(clipper,d%ndom,-1._c_float)
           call w%table_write_begin(itable_domains,iwrite,clipper)
           do while(ImGuiListClipper_Step(clipper))
              call c_f_pointer(clipper,clipper_f)
              do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                 call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
                 if (igTableSetColumnIndex(0)) then
                    call iw_text(string(i))
                    if (iw_highlight_selectable("##ncidomain_" // string(i))) ihover = i
                 end if
                 if (igTableSetColumnIndex(1)) call iw_cell_right(string(d%dom_vol(i) * bohrtoa**3,'f',decimal=3))
                 if (igTableSetColumnIndex(2)) call iw_cell_right(string(d%dom_q(i),'f',decimal=5))
                 if (igTableSetColumnIndex(3)) call iw_cell_right(string(d%dom_npts(i)))
                 if (igTableSetColumnIndex(4)) call iw_cell_right(string(d%dom_rmin(i),'f',decimal=5))
                 if (igTableSetColumnIndex(5)) call iw_cell_right(string(d%dom_rmax(i),'f',decimal=5))
              end do
           end do
           call ImGuiListClipper_End(clipper)
           call ImGuiListClipper_destroy(clipper)
           call w%table_write_end(itable_domains,table_title("NCI domains",isys))
           call igEndTable()
        end if
        call iw_tooltip("Connected regions of the grid where the masked reduced density gradient is&
           & below the isovalue, by decreasing volume, with the electron density integrated over&
           & each (electrons) and the range of the signed density in each. The region under the&
           & mouse is highlighted in the view and in the plot",ttshown)
      end associate

    end subroutine draw_domains

    !> The Output section: the files NCIPLOT writes, from a root name
    !> chosen in a dialog, and whether the molecules of the xyz file
    !> are completed.
    subroutine draw_output()
      integer :: idum
      logical :: okwrite

      call iw_separatortext("Output",w%nc%sec(nci_sec_output))
      if (.not.w%nc%sec(nci_sec_output)%open) return

      ! a masked result to write
      okwrite = associated(w%nc%res) .and. associated(w%nc%drv)
      if (okwrite) okwrite = allocated(w%nc%drv%inter)

      ch = iw_checkbox("Complete molecules##ncimolmotif",w%nc%molmotif,labelleft=.true.)
      call iw_tooltip("Complete the molecules that cross the boundary of the cell or the region&
         & in the xyz file written next to the cube files (MOLMOTIF)",ttshown)
      if (iw_button("Write Files...##nciwrite",sameline=.true.,disabled=.not.okwrite)) &
         idum = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savencifiles,idparent=w%id,orraise=-1)
      call iw_tooltip("Write the files of the NCIPLOT keyword for the last calculation: the cube files&
         & of the signed density and the masked RDG (and the promolecular density, with fragments or voids),&
         & the s vs. density data, the VMD script with its xyz file, and the checkpoint, named&
         & from the file chosen without its extension (existing files with those names are&
         & replaced)",ttshown,whendisabled=.true.)

    end subroutine draw_output

    !> Request the blocking job that writes the output files
    !> (run_nci), named from the file chosen in the dialog.
    subroutine request_write()

      w%errmsg = ""
      w%okmsg = ""
      if (.not.w%request_block()) then
         w%errmsg = "Another calculation is starting: press the button again"
         return
      end if
      w%nc%jwrite = .true.
      w%nc%jpromol = .false.
      w%nc%jsys = isys

    end subroutine request_write

    !> Add the atoms selected in the view (the cell atoms with a
    !> highlight, at their cell positions) as a new fragment, selected.
    subroutine add_selfrag()
      integer :: i, n
      integer, allocatable :: iat(:)
      real*8, allocatable :: x(:,:)
      type(nci_userfrag) :: uf

      call sysc(isys)%highlighted_atom_list(n,iat)
      if (n == 0) return
      allocate(x(3,n))
      do i = 1, n
         x(:,i) = sys(isys)%c%atcel(iat(i))%r
      end do
      uf%fr = sys(isys)%c%identify_fragment(n,x)
      if (allocated(w%nc%userfr)) then
         w%nc%userfr = (/w%nc%userfr,uf/)
      else
         w%nc%userfr = (/uf/)
      end if

    end subroutine add_selfrag

    !> Request the blocking job that calculates the NCI plot (run_nci)
    !> on the grid and fragments of the form (form_grid_frags).
    subroutine request_job()

      w%errmsg = ""
      w%okmsg = ""
      if (.not.w%request_block()) then
         w%errmsg = "Another calculation is starting: press the button again"
         return
      end if
      call form_grid_frags(w%nc%opts)
      w%nc%jwrite = .false.
      w%nc%jpromol = .false.
      w%nc%jsys = isys

    end subroutine request_job

    !> Set the grid and the fragments of options o from the form: the
    !> grid is the nodes of the density field (native level) or nn
    !> points over the region (origin and steps, cell-frame Cartesian
    !> bohr; periodic over the whole cell of a crystal), the fragments
    !> are the selected molecules and the selected user-made fragments
    !> (none if fewer than two: the criterion needs two).
    subroutine form_grid_frags(o)
      type(nci_options), intent(inout) :: o

      integer :: j, k, nfrag

      associate(c => sys(isys)%c)
        o%periodic = (.not.ismol .and. w%nc%iregion == iso_region_cell)
        o%nstep = nn
        if (w%nc%ilevel == 0) then
           o%x0 = 0d0
           do j = 1, 3
              o%xmat(:,j) = c%m_x2c(:,j) / real(nn(j),8)
           end do
        else
           o%x0 = box(:,0)
           do j = 1, 3
              o%xmat(:,j) = box(:,j) / real(nn(j),8)
           end do
        end if

        ! (crystal fragments are cell molecules or sets of cell atoms:
        ! all their lattice translations are in the fragment)
        o%frcell = .not.ismol
        nfrag = nfrag_total()
        if (allocated(o%fr)) deallocate(o%fr)
        allocate(o%fr(max(nfrag,1)))
        o%nfrag = 0
        if (nfrag > 1) then
           if (allocated(w%nc%molsel)) then
              do j = 1, size(w%nc%molsel,1)
                 if (.not.w%nc%molsel(j)) cycle
                 o%nfrag = o%nfrag + 1
                 o%fr(o%nfrag) = c%mol(j)
              end do
           end if
           if (allocated(w%nc%userfr)) then
              do k = 1, size(w%nc%userfr,1)
                 if (.not.w%nc%userfr(k)%sel) cycle
                 o%nfrag = o%nfrag + 1
                 o%fr(o%nfrag) = w%nc%userfr(k)%fr
              end do
           end if
        end if
      end associate

    end subroutine form_grid_frags


    !> The region between the fragments (Show region): with two or
    !> more fragments and a valid grid, find promolecular grids for the
    !> grid and fragments of the form (the result, or those calculated
    !> for the region), calculating them in a blocking job if there are
    !> none, and show the fragment criterion as a translucent transient
    !> isosurface in the anchor view. Runs every frame, outside the
    !> sections.
    subroutine draw_region()
      type(nci_options) :: o
      type(nci_result), pointer :: src
      integer :: isrc, itrep, i, m(3)
      logical :: found, changed, parchanged
      real*8 :: xr(3,0:3), xsh(3), x0(3)

      ! (no fragment criterion with fewer than two fragments: the
      ! region is turned off, not just hidden)
      if (nfrag_total() < 2) w%nc%showregion = .false.
      if (.not.w%nc%showregion) then
         ! (off: the grids are not needed any more)
         if (associated(w%nc%pre)) deallocate(w%nc%pre)
         if (allocated(w%nc%region_g)) deallocate(w%nc%region_g)
         return
      end if
      if (.not.okbox .or. any(nn < 2) .or. w%nc%picking) return
      call form_grid_frags(o)

      ! the grids: the result if it is for this grid and these
      ! fragments, else those calculated for the region, else
      ! calculate them (unless a job is already on its way, or a
      ! widget is being dragged: the region form updates live)
      src => null()
      isrc = 0
      if (associated(w%nc%res)) then
         if (allocated(w%nc%res%rhofrag) .and. same_grid_frags(o,w%nc%opts)) then
            src => w%nc%res
            isrc = 1
         end if
      end if
      if (isrc == 0 .and. associated(w%nc%pre)) then
         if (same_grid_frags(o,w%nc%prekey)) then
            src => w%nc%pre
            isrc = 2
         end if
      end if
      if (isrc == 0) then
         if (w%nc%jsys == 0 .and. .not.igIsAnyItemActive()) then
            if (w%request_block()) then
               w%nc%prekey = o
               w%nc%jwrite = .false.
               w%nc%jpromol = .true.
               w%nc%jsys = isys
            end if
         end if
         return
      end if

      ! the criterion field, if its source, the grids, or the
      ! parameters changed (Dominance and Share: once the drag is
      ! released, not live)
      parchanged = any(w%nc%region_par /= (/w%nc%opts%rhoparam,w%nc%opts%rhoparam2/))
      if (parchanged) parchanged = .not.igIsAnyItemActive()
      m = region_roll(o,src%n)
      if (.not.allocated(w%nc%region_g) .or. w%nc%region_src /= isrc .or. w%nc%region_gen /= w%nc%gen .or.&
         parchanged .or. any(w%nc%region_shift /= m)) then
         call nci_frag_region(w%nc%opts,src,w%nc%region_g)
         do i = 1, 3
            if (m(i) /= 0) w%nc%region_g = cshift(w%nc%region_g,m(i),i)
         end do
         w%nc%region_src = isrc
         w%nc%region_gen = w%nc%gen
         w%nc%region_par = (/w%nc%opts%rhoparam,w%nc%opts%rhoparam2/)
         w%nc%region_shift = m
         changed = .true.
      else
         changed = .false.
      end if

      ! the transient isosurface: the parallelepiped with corners on
      ! the first and last nodes of the (rolled) grid, which the
      ! isosurface samples endpoint-inclusively (user-frame angstrom)
      if (iview < 1 .or. iview > nwin) return
      if (.not.associated(win(iview)%sc)) return
      call win(iview)%sc%show_transient_iso(w%uid,2,itrep,found)
      if (itrep == 0) return
      xsh = 0d0
      if (ismol) xsh = sys(isys)%c%molx0
      x0 = src%x0 + matmul(src%xmat,real(m,8))
      xr(:,0) = (x0 + xsh) * bohrtoa
      do i = 1, 3
         xr(:,i) = (x0 + src%xmat(:,i) * real(src%n(i)-1,8) + xsh) * bohrtoa
      end do
      associate(r => win(iview)%sc%reptrans(itrep))
        if (.not.found .or. r%iso%fref%id /= 0 .or. .not.r%iso%isvalue()) then
           call r%iso%set_field(isys,0,value=.true.)
           r%name = "NCI fragment region"
           r%iso%slot(1)%isoval = 0d0
           r%iso%slot(1)%rgb = nci_region_rgb
           r%iso%slot(1)%alpha = nci_region_alpha
           changed = .true.
        end if
        if (changed .or. .not.r%iso%grid_isapplied(src%n,iso_region_parallel,xr)) then
           call r%iso%apply_grid(src%n,iso_region_parallel,xr)
           call r%iso%set_samples(isys,w%nc%region_g,.false.)
           win(iview)%sc%forcebuildlists = .true.
        end if
      end associate

    end subroutine draw_region

    !> The nodes the region of a periodic grid of options o (n nodes)
    !> is rolled by along each axis, so it covers the cell in which the
    !> anchor view draws the atoms (display_cell_window: the display
    !> origin, or the middle of the vacuum of a slab), like the
    !> periodic isosurfaces.
    function region_roll(o,n) result(m)
      use representations, only: display_cell_window
      use display, only: rep_display
      type(nci_options), intent(in) :: o
      integer, intent(in) :: n(3)
      integer :: m(3)

      m = 0
      if (.not.o%periodic) return
      if (iview < 1 .or. iview > nwin) return
      if (.not.associated(win(iview)%sc)) return
      m = nint(display_cell_window(isys,win(iview)%sc%disp,rep_display()) * real(n,8))

    end function region_roll

    !> Whether there is a result and the form asks for another
    !> calculation: a different density field or derivatives, grid
    !> (region and level) or fragments (form_grid_frags), or voids
    !> on a result without the promolecular density. The cutoffs, the
    !> plot settings and the fragment criterion only mask the result
    !> again.
    function result_stale() result(stale)
      logical :: stale

      type(nci_options) :: o

      stale = .false.
      if (.not.associated(w%nc%res)) return
      associate(r => w%nc%res, opts => w%nc%opts)
        stale = (opts%ifield /= r%ifield) .or. (opts%fourierint .neqv. r%fourierint)
        if (.not.stale) stale = opts%dovoid .and. .not.allocated(r%rhoat)
        if (.not.stale) then
           ! (opts holds the grid and the fragments of the result)
           call form_grid_frags(o)
           stale = .not.same_grid_frags(o,opts)
        end if
      end associate

    end function result_stale

    !> Cancel a pending pick, releasing the view.
    subroutine cancel_pick()
      if (.not.w%nc%picking) return
      if (w%nc%pickview >= 1 .and. w%nc%pickview <= nwin) &
         call win(w%nc%pickview)%viewmode_release_forced(w%id)
      w%nc%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_nci

  !> Draw the overlay of the pending job of the NCI plot window.
  module subroutine block_nci(w)
    use systems, only: sys, sysc, ok_system, sys_init
    use utils, only: iw_wait_overlay, file_name_root
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    character(len=:), allocatable :: info
    integer :: isys, ifield
    type(nci_options), pointer :: o

    info = ""
    isys = w%nc%jsys
    ! (writing: the field of the result, which the form may no longer show)
    ifield = w%nc%opts%ifield
    if (w%nc%jwrite .and. associated(w%nc%res)) ifield = w%nc%res%ifield
    if (ok_system(isys,sys_init)) then
       info = "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (sys(isys)%goodfield(ifield)) &
          info = info // newline // "Field:  " // string(ifield) // ": " //&
          trim(sys(isys)%f(ifield)%name)
       o => w%nc%opts
       if (w%nc%jpromol) o => w%nc%prekey
       info = info // newline // "Grid:   " // string(o%nstep(1)) // " x " //&
          string(o%nstep(2)) // " x " // string(o%nstep(3))
       if (w%nc%jpromol) info = info // newline // "Fragments: " // string(o%nfrag)
    end if
    if (w%nc%jpromol) then
       call iw_wait_overlay("Calculating the region between the fragments...",info,.true.)
    elseif (w%nc%jwrite) then
       info = info // newline // "Files:  " // file_name_root(w%okfile) // "*"
       call iw_wait_overlay("Writing the NCI files...",info,.false.)
    else
       call iw_wait_overlay("Calculating the NCI plot...",info,.true.)
    end if

  end subroutine block_nci

  !> Run the pending job of the NCI plot window, called by the main
  !> loop after the frame with the overlay: the reduced density
  !> gradient and the signed density on the grid of the options
  !> (request_job), then the masks, the fields and the objects
  !> (nci_refresh). The calculation can be cancelled with Esc, and
  !> then the previous results are gone and no field is loaded.
  module subroutine run_nci(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_cplist
    use nci, only: nci_compute, nci_compute_promol
    use gui_main, only: begin_cancellable, end_cancellable
    use interfaces_glfw, only: glfwGetTime
    use utils, only: duration_string
    use tools_io, only: uout, string
    class(window), intent(inout), target :: w

    integer :: isys
    real*8 :: t
    logical :: cancelled
    character(len=:), allocatable :: errmsg
    type(nci_result), pointer :: pre

    isys = w%nc%jsys
    w%nc%jsys = 0
    if (.not.ok_system(isys,sys_init)) return
    if (w%nc%jwrite) then
       call nci_write_files(w,isys)
       return
    end if

    ! the promolecular grids of the region between the fragments; a
    ! cancel or an error turns the region off, so it is not requested
    ! again on the next frame
    if (w%nc%jpromol) then
       w%nc%jpromol = .false.
       if (.not.associated(w%nc%pre)) allocate(w%nc%pre)
       call begin_cancellable()
       call nci_compute_promol(w%nc%prekey,sys(isys),w%nc%pre,errmsg)
       cancelled = end_cancellable()
       if (cancelled) errmsg = "The calculation of the region was cancelled"
       if (len_trim(errmsg) > 0) then
          w%errmsg = errmsg
          deallocate(w%nc%pre)
          w%nc%showregion = .false.
          return
       end if
       w%nc%gen = w%nc%gen + 1
       return
    end if

    if (.not.sys(isys)%goodfield(w%nc%opts%ifield)) then
       w%errmsg = "The density field is gone"
       return
    end if
    w%errmsg = ""

    ! the grids; a cancel or an error leaves no result
    if (.not.associated(w%nc%res)) allocate(w%nc%res)
    t = glfwGetTime()
    call begin_cancellable()
    ! (the promolecular grids of the region, if they are for this grid
    ! and these fragments, are taken over instead of calculated again;
    ! a disassociated pointer is an absent argument)
    pre => null()
    if (associated(w%nc%pre)) then
       if (same_grid_frags(w%nc%opts,w%nc%prekey)) pre => w%nc%pre
    end if
    call nci_compute(w%nc%opts,sys(isys),w%nc%res,errmsg,pre=pre)
    cancelled = end_cancellable()
    ! (the grids of the region are now in the result, or for another
    ! grid than the form's)
    if (associated(w%nc%pre)) deallocate(w%nc%pre)
    t = glfwGetTime() - t
    if (cancelled) errmsg = "The calculation was cancelled: no field was loaded"
    if (len_trim(errmsg) > 0) then
       w%errmsg = errmsg
       deallocate(w%nc%res)
       return
    end if
    w%nc%timelast = sysc(isys)%timelastchange_geometry
    w%nc%gen = w%nc%gen + 1
    ! the dat-file cutoffs of this result (the defaults of its field),
    ! for its scatter plot and its files, and fit the plot to them
    w%nc%datcut = (/w%nc%opts%rhocut,w%nc%opts%dimcut/)
    w%nc%fitplot = .true.

    ! the masks, the fields and the objects
    call nci_refresh(w,isys,w%anchor_view())
    call sysc(isys)%post_event(lastchange_cplist)
    sysc(isys)%showfields = .true.
    if (w%nc%fsref%ok(sys(isys)) .and. w%nc%frref%ok(sys(isys))) &
       write (uout,'("* NCI plot grids (",A," x ",A," x ",A," points from field ",A,", calculated in ",A,&
       &") loaded as fields ",A," (RDG) and ",A," (signed density)"/)') &
       string(w%nc%res%n(1)), string(w%nc%res%n(2)), string(w%nc%res%n(3)), string(w%nc%res%ifield),&
       duration_string(t), string(w%nc%fsref%id), string(w%nc%frref%id)

  end subroutine run_nci

  !> Drop the heap results of the NCI plot window w, and the domain
  !> labels its isosurface object carries (the object outlives the
  !> results).
  module subroutine end_nci(w)
    class(window), intent(inout), target :: w

    integer :: iview

    if (associated(w%nc%res)) deallocate(w%nc%res)
    if (associated(w%nc%drv)) deallocate(w%nc%drv)
    if (associated(w%nc%pre)) deallocate(w%nc%pre)
    if (allocated(w%nc%region_g)) deallocate(w%nc%region_g)
    iview = w%anchor_view()
    if (iview >= 1 .and. iview <= nwin) then
       if (win(iview)%isopen .and. associated(win(iview)%sc)) then
          if (w%nc%isoref%ok(win(iview)%sc)) then
             call win(iview)%sc%rep(w%nc%isoref%id)%iso%clear_groups()
             win(iview)%sc%forcebuildlists = .true.
          end if
       end if
    end if

  end subroutine end_nci

  !xx! private procedures

  !> Whether options a and b have the same grid and the same
  !> fragments (atoms, species and positions).
  function same_grid_frags(a,b) result(same)
    type(nci_options), intent(in) :: a, b
    logical :: same

    integer :: i, j

    same = all(a%nstep == b%nstep) .and. all(a%x0 == b%x0) .and. all(a%xmat == b%xmat) .and.&
       a%nfrag == b%nfrag .and. (a%frcell .eqv. b%frcell)
    if (.not.same .or. a%nfrag == 0) return
    if (.not.allocated(a%fr) .or. .not.allocated(b%fr)) then
       same = .false.
       return
    end if
    do i = 1, a%nfrag
       if (a%fr(i)%nat /= b%fr(i)%nat) then
          same = .false.
          return
       end if
       do j = 1, a%fr(i)%nat
          if (any(a%fr(i)%at(j)%r /= b%fr(i)%at(j)%r) .or. a%fr(i)%spc(a%fr(i)%at(j)%is)%z /=&
             b%fr(i)%spc(b%fr(i)%at(j)%is)%z) then
             same = .false.
             return
          end if
       end do
    end do

  end function same_grid_frags

  !> Write the files of the NCIPLOT keyword for the result of the NCI
  !> plot window w of system isys, named from the file chosen in its
  !> dialog (w%okfile) without the extension: the checkpoint, then
  !> the dat, cube, and VMD files (nci_write_files; the VMD script
  !> names the files next to it by their base name).
  subroutine nci_write_files(w,isys)
    use systems, only: sys
    use nci, only: nci_chk_write, nci_write_files_ => nci_write_files, nci_options
    use utils, only: file_name_root, file_name_base
    use tools_io, only: uout
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys

    type(nci_options) :: o
    real*8 :: cs(2)
    logical :: ok
    integer :: id
    character(len=:), allocatable :: root, errmsg

    w%okmsg = ""
    root = file_name_root(w%okfile)
    if (.not.associated(w%nc%res) .or. .not.associated(w%nc%drv)) then
       w%errmsg = "There is no result to write"
       return
    end if
    if (.not.allocated(w%nc%drv%inter)) then
       w%errmsg = "The result is not masked: change a cutoff or recalculate"
       return
    end if
    id = w%nc%fsref%id
    if (.not.sys(isys)%goodfield(id,uid=w%nc%fsref%uid,n=w%nc%res%n)) then
       w%errmsg = "The masked RDG field is gone: recalculate"
       return
    end if

    ! the density cutoff (CUTPLOT) only sets the color scale of the
    ! VMD script here: the masked grid is that of the window; the
    ! dat file has the cutoffs the result was calculated with
    o = w%nc%opts
    ok = nci_cscale(w,cs)
    o%rhoplot = maxval(abs(cs))
    o%rhocut = w%nc%datcut(1)
    o%dimcut = w%nc%datcut(2)
    call nci_chk_write(root // ".chk_nci",o,sys(isys),w%nc%res,errmsg)
    if (len_trim(errmsg) == 0) &
       call nci_write_files_(root,o,sys(isys),w%nc%res,w%nc%drv%inter,sys(isys)%f(id)%grid%f,&
       w%nc%molmotif,errmsg,refroot=file_name_base(root))
    w%errmsg = errmsg
    if (len_trim(errmsg) > 0) return
    w%okmsg = "NCI plot files written with root " // root
    write (uout,'("* NCI plot files written with root: ",A/)') root
    call okfile_save_dir(w%okfile)

  end subroutine nci_write_files

  !> Mask the result of the NCI plot window w for system isys with the
  !> cutoffs of its options (nci_masks) and put the grids in their
  !> fields: the masked reduced density gradient replaces the data of
  !> its field (loaded if it is gone or on another grid), the signed
  !> density is loaded if its field is gone; then the objects of the
  !> scene of view iview (nci_objects). An error from the masks (a
  !> criterion that needs a grid the result does not have) is
  !> reported in the window and leaves the fields and the objects as
  !> they were.
  subroutine nci_refresh(w,isys,iview)
    use systems, only: sys, replace_grid_field
    use nci, only: nci_masks, nci_options
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    integer, intent(in) :: iview

    type(nci_options) :: o
    logical :: ok
    real*8, allocatable :: splot(:,:,:)
    character(len=:), allocatable :: errmsg

    if (.not.associated(w%nc%res)) return
    if (.not.associated(w%nc%drv)) allocate(w%nc%drv)

    ! the masks (the scatter and the domains of a failed pass are not
    ! shown)
    ! (no density cutoff in the window: the range is the cutoff)
    o = w%nc%opts
    o%rhoplot = huge(1d0)
    call nci_masks(o,w%nc%res,inter=w%nc%drv%inter,splot=splot,errmsg=errmsg)
    if (len_trim(errmsg) > 0) then
       w%errmsg = errmsg
       w%nc%drv = nci_derived()
       return
    end if
    w%errmsg = ""

    ! the scatter and the domains
    call nci_scatter(w)
    call nci_find_domains(w,splot)

    ! the fields
    associate(r => w%nc%res, o => w%nc%opts)
      if (.not.sys(isys)%goodfield(w%nc%frref%id,uid=w%nc%frref%uid,n=r%n)) &
         call nci_load_field(w,isys,w%nc%frref,"NCI sign(l2)*rho of field " // string(r%ifield),r%srho)
      ok = .false.
      if (sys(isys)%goodfield(w%nc%fsref%id,uid=w%nc%fsref%uid,n=r%n)) &
         ok = replace_grid_field(isys,w%nc%fsref%id,splot)
      if (.not.ok) call nci_load_field(w,isys,w%nc%fsref,"NCI RDG of field " // string(r%ifield),splot)
    end associate

    call nci_objects(w,isys,iview)

  end subroutine nci_refresh

  !> Make or update the objects of the NCI plot window w for system
  !> isys in the scene of view iview (nothing if it is not an open
  !> view of real space with a scene): an isosurface of the masked RDG
  !> field at the RDG isovalue, colored by the signed density field
  !> on the blue-green-red scale (blue = attractive, red = repulsive:
  !> the NCI convention), and its color bar. The color scale is set
  !> when the isosurface starts being colored by that field, to the
  !> default density cutoff of the field (or the widest bound of the
  !> range); after that it belongs to the object, and changes made to
  !> it in the object editor stay.
  subroutine nci_objects(w,isys,iview)
    use systems, only: sys
    use scenes, only: space_real
    use representations, only: reptype_isosurface, repflavor_isosurface, reptype_colorbar,&
       repflavor_colorbar, iso_map_field
    use utils, only: iw_cmap_jet
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    integer, intent(in) :: iview

    integer :: irep, icb
    real*8 :: cs(2)
    logical :: ok

    if (.not.w%nc%fsref%ok(sys(isys)) .or. .not.w%nc%frref%ok(sys(isys))) return
    if (iview < 1 .or. iview > nwin) return
    if (.not.win(iview)%isopen .or. win(iview)%type /= wintype_view) return
    if (.not.associated(win(iview)%sc)) return
    if (win(iview)%sc%isinit == 0 .or. win(iview)%sc%space /= space_real) return

    associate(sc => win(iview)%sc)
      ! the isosurface, on the masked RDG field at its own grid
      if (.not.w%nc%isoref%ok(sc)) then
         call sc%add_representation(reptype_isosurface,repflavor_isosurface,id=irep)
         if (irep <= 0) return
         sc%rep(irep)%name = "NCI plot"
         call w%nc%isoref%set(sc,irep)
         call sc%undo_note("Add " // trim(sc%rep(irep)%name))
      end if
      irep = w%nc%isoref%id
      associate(r => sc%rep(irep))
        if (.not.r%iso%fref%same(w%nc%fsref)) call r%iso%set_field(isys,w%nc%fsref%id,value=.true.)
        if (r%iso%niso < 1) call r%iso%add_iso(isoval=w%nc%opts%dimplot)
        r%iso%slot(1)%isoval = w%nc%opts%dimplot
        r%iso%slot(1)%shown = .true.
        if (r%iso%slot(1)%imap_mode /= iso_map_field .or. .not.r%iso%slot(1)%mapref%same(w%nc%frref)) then
           ok = nci_cscale(w,cs)
           call r%iso%slot(1)%set_map_field(isys,w%nc%frref%id,icmap=iw_cmap_jet,maprange=cs)
        end if
        r%shown = .true.
        ! the domains of the grid points, for the highlight
        if (associated(w%nc%drv)) then
           if (allocated(w%nc%drv%lbl)) call r%iso%set_groups(w%nc%drv%lbl)
        end if
      end associate

      ! the color bar of that isosurface
      if (.not.w%nc%cbref%ok(sc)) then
         call sc%add_representation(reptype_colorbar,repflavor_colorbar,id=icb)
         if (icb <= 0) return
         sc%rep(icb)%name = "NCI color bar"
         call w%nc%cbref%set(sc,icb)
         call sc%undo_note("Add " // trim(sc%rep(icb)%name))
      end if
      icb = w%nc%cbref%id
      sc%rep(icb)%colorbar%isouid = sc%rep(irep)%uid
      sc%rep(icb)%colorbar%slotuid = sc%rep(irep)%iso%slot(1)%uid
      sc%rep(icb)%colorbar%customtitle = .true.
      sc%rep(icb)%colorbar%title = "sign(λ2)ρ (a.u.)"
      sc%rep(icb)%shown = .true.
      sc%forcebuildlists = .true.
    end associate

  end subroutine nci_objects

  !> Sample the points of the result of the NCI plot window w under
  !> the dat-file cutoffs of the result (nci_in_dat, datcut) into the
  !> scatter of its derived data: at most nscmax of them, evenly
  !> spaced over the grid.
  subroutine nci_scatter(w)
    use nci, only: nci_in_dat, nci_options
    type(window), intent(inout), target :: w

    integer, parameter :: nscmax = 20000

    type(nci_options) :: o
    integer :: i, j, k, m, nsc, stride, n(3)
    logical, allocatable :: msk(:,:,:)

    if (.not.associated(w%nc%res) .or. .not.associated(w%nc%drv)) return
    o = w%nc%opts
    o%rhocut = w%nc%datcut(1)
    o%dimcut = w%nc%datcut(2)
    associate(r => w%nc%res, d => w%nc%drv)
      if (allocated(d%scx)) deallocate(d%scx)
      if (allocated(d%scy)) deallocate(d%scy)
      if (allocated(d%scidx)) deallocate(d%scidx)
      if (allocated(d%scl)) deallocate(d%scl)
      if (allocated(d%scbin)) deallocate(d%scbin)
      if (.not.allocated(d%inter)) return
      n = r%n
      allocate(msk(n(1),n(2),n(3)))
      do k = 1, n(3)
         do j = 1, n(2)
            do i = 1, n(1)
               msk(i,j,k) = nci_in_dat(o,r,d%inter,i,j,k)
            end do
         end do
      end do
      ! every stride-th point under the cutoffs, with its index in the
      ! flattened grid
      stride = max((count(msk) + nscmax - 1) / nscmax,1)
      nsc = (count(msk) + stride - 1) / stride
      allocate(d%scx(nsc),d%scy(nsc),d%scidx(nsc))
      m = 0
      nsc = 0
      do k = 1, n(3)
         do j = 1, n(2)
            do i = 1, n(1)
               if (.not.msk(i,j,k)) cycle
               m = m + 1
               if (mod(m-1,stride) /= 0) cycle
               nsc = nsc + 1
               d%scx(nsc) = real(r%srho(i,j,k),c_double)
               d%scy(nsc) = real(r%rdg(i,j,k),c_double)
               d%scidx(nsc) = i + n(1) * (j - 1 + n(2) * (k - 1))
            end do
         end do
      end do
    end associate

  end subroutine nci_scatter

  !> Sort the scatter points of the derived data d in place into the
  !> color bins: bin 0 for the points outside the sign(lambda2)*rho
  !> range rng, bins 1..nscbin over the color scale cs for the others
  !> (clamped at its ends, as the isosurface); unless they are sorted
  !> for those already.
  subroutine nci_scatter_bins(d,cs,rng)
    type(nci_derived), intent(inout), target :: d
    real*8, intent(in) :: cs(2), rng(2)

    integer :: i, n
    integer, allocatable :: ib(:), pos(:)
    real(c_double), allocatable :: tx(:), ty(:)
    integer, allocatable :: tidx(:), tl(:)
    real*8 :: lo, hi

    if (.not.allocated(d%scx)) return
    if (allocated(d%scbin)) then
       if (all(d%scb_cs == cs) .and. all(d%scb_rng == rng)) return
    end if
    n = size(d%scx,1)
    lo = min(cs(1),cs(2))
    hi = max(cs(1),cs(2))
    allocate(ib(n))
    do i = 1, n
       if (d%scx(i) < rng(1) .or. d%scx(i) > rng(2)) then
          ib(i) = 0
       else
          ib(i) = 1 + nint(min(max((d%scx(i) - lo) / max(hi - lo,1d-30),0d0),1d0) * (nscbin-1))
       end if
    end do

    ! counting sort into the bins
    if (allocated(d%scbin)) deallocate(d%scbin)
    allocate(d%scbin(-1:nscbin),pos(0:nscbin))
    d%scbin = 0
    do i = 1, n
       d%scbin(ib(i)) = d%scbin(ib(i)) + 1
    end do
    do i = 0, nscbin
       d%scbin(i) = d%scbin(i) + d%scbin(i-1)
    end do
    pos = d%scbin(-1:nscbin-1)
    ! (the grid indices and domains of the points go with them)
    tx = d%scx
    ty = d%scy
    tidx = d%scidx
    if (allocated(d%scl)) tl = d%scl
    do i = 1, n
       pos(ib(i)) = pos(ib(i)) + 1
       d%scx(pos(ib(i))) = tx(i)
       d%scy(pos(ib(i))) = ty(i)
       d%scidx(pos(ib(i))) = tidx(i)
       if (allocated(d%scl)) d%scl(pos(ib(i))) = tl(i)
    end do
    d%scb_cs = cs
    d%scb_rng = rng

  end subroutine nci_scatter_bins

  !> Find the domains of the isosurface plot of the NCI plot window w
  !> (nci_domains on the masked RDG f) into its derived data; the
  !> labels wait there for the isosurface object (nci_objects).
  subroutine nci_find_domains(w,f)
    use nci, only: nci_domains
    type(window), intent(inout), target :: w
    real*8, intent(in) :: f(:,:,:)

    character(len=:), allocatable :: errmsg
    integer :: i, l, n(3)

    if (.not.associated(w%nc%res) .or. .not.associated(w%nc%drv)) return
    associate(d => w%nc%drv)
      call nci_domains(w%nc%opts,w%nc%res,f,d%ndom,d%lbl,d%dom_npts,d%dom_vol,d%dom_rmin,d%dom_rmax,&
         errmsg,q=d%dom_q)
      if (len_trim(errmsg) > 0) then
         w%errmsg = errmsg
         d%ndom = 0
      end if
      d%voltot = 0d0
      if (allocated(d%dom_vol)) d%voltot = sum(d%dom_vol)

      ! the domain of each scatter point, for its highlight in the plot
      if (allocated(d%scl)) deallocate(d%scl)
      if (allocated(d%scidx) .and. allocated(d%lbl) .and. d%ndom > 0) then
         allocate(d%scl(size(d%scidx)))
         n = shape(d%lbl)
         do i = 1, size(d%scidx)
            l = d%scidx(i) - 1
            d%scl(i) = d%lbl(mod(l,n(1))+1,mod(l/n(1),n(2))+1,l/(n(1)*n(2))+1)
         end do
      end if
    end associate

  end subroutine nci_find_domains

  !> Highlight domain ihover (0 = none) of the NCI plot window w
  !> (system isys, view iview) on its isosurface object, if the object
  !> is still the one the window made: on the masked RDG field at its
  !> native grid, which the domain labels are on.
  subroutine nci_highlight(w,isys,iview,ihover)
    use systems, only: sys
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    integer, intent(in) :: iview
    integer, intent(in) :: ihover

    if (iview < 1 .or. iview > nwin) return
    if (.not.win(iview)%isopen .or. .not.associated(win(iview)%sc)) return
    if (.not.w%nc%isoref%ok(win(iview)%sc)) return
    associate(r => win(iview)%sc%rep(w%nc%isoref%id))
      if (.not.allocated(r%iso%lbl)) return
      if (.not.r%iso%fref%same(w%nc%fsref) .or. any(r%iso%nptsxyz /= 0) .or.&
         .not.sys(isys)%goodfield(w%nc%fsref%id,uid=w%nc%fsref%uid,n=shape(r%iso%lbl))) then
         ! the object was changed in its editor: the labels mean nothing to it
         call r%iso%clear_groups()
         win(iview)%sc%forcebuildlists = .true.
         return
      end if
      if (r%iso%ihighlight /= ihover) then
         r%iso%ihighlight = ihover
         win(iview)%sc%forcebuildlists = .true.
      end if
    end associate

  end subroutine nci_highlight

  !> Load grid g (on the grid of the result of the NCI plot window w)
  !> as a field of system isys named name, in place of the field ref
  !> points to if it is still there and not the reference or the
  !> density of the form, and point ref to it.
  subroutine nci_load_field(w,isys,ref,name,g)
    use systems, only: sys, load_grid_field
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    type(field_ref), intent(inout) :: ref
    character(len=*), intent(in) :: name
    real*8, intent(in) :: g(:,:,:)

    integer :: j
    real*8 :: x2cl(3,3)

    if (ref%ok(sys(isys))) then
       if (ref%id /= sys(isys)%iref .and. .not.ref%same(w%nc%fref)) call sys(isys)%unload_field(ref%id)
    end if
    do j = 1, 3
       x2cl(:,j) = w%nc%res%xmat(:,j) * real(w%nc%res%n(j),8)
    end do
    call ref%set(sys(isys),load_grid_field(isys,name,g,sys(isys)%c%c2x(w%nc%res%x0),x2cl))

  end subroutine nci_load_field

  !> The color scale of the signed density of the isosurface object
  !> of the NCI plot window w, if it is colored by the signed density
  !> field of the window (returns .true.): the values cs at its two
  !> ends are read or, if set, written (as a fixed range, which the
  !> color bar of the isosurface follows); with reset, the default is
  !> written. If there is no such object, cs is the default: minus to
  !> plus the default density cutoff of the field, or the widest bound
  !> of the sign(lambda2)*rho range. icmap is the colormap of the
  !> isosurface (jet, the default, if there is no such object).
  function nci_cscale(w,cs,set,reset,icmap) result(ok)
    use representations, only: iso_map_field
    use utils, only: iw_cmap_jet
    type(window), intent(inout), target :: w
    real*8, intent(inout) :: cs(2)
    logical, intent(in), optional :: set, reset
    integer, intent(out), optional :: icmap
    logical :: ok

    integer :: iview, irep
    logical :: set_, reset_

    ok = .false.
    set_ = .false.
    if (present(set)) set_ = set
    reset_ = .false.
    if (present(reset)) reset_ = reset
    if (.not.set_) then
       cs(2) = max(w%nc%opts%rhoplot,maxval(abs(w%nc%opts%srhorange)))
       cs(1) = -cs(2)
    end if
    if (present(icmap)) icmap = iw_cmap_jet
    iview = w%anchor_view()
    if (iview < 1 .or. iview > nwin) return
    if (win(iview)%type /= wintype_view .or. .not.associated(win(iview)%sc)) return
    associate(sc => win(iview)%sc)
      irep = w%nc%isoref%get(sc)
      if (irep <= 0) return
      associate(r => sc%rep(irep))
        if (r%iso%niso < 1) return
        if (r%iso%slot(1)%imap_mode /= iso_map_field .or. .not.r%iso%slot(1)%mapref%same(w%nc%frref)) return
        if (set_ .or. reset_) then
           r%iso%slot(1)%maprange = cs
           r%iso%slot(1)%maprange_auto = .false.
           sc%forcebuildlists = .true.
           call sc%undo_note("Edit " // trim(r%name))
        else
           cs = r%iso%slot(1)%maprange
        end if
        if (present(icmap)) icmap = r%iso%slot(1)%icmap
        ok = .true.
      end associate
    end associate

  end function nci_cscale

end submodule nci
