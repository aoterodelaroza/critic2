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
  implicit none

  ! pending pick of the NCI plot window (nci_state%pickkind)
  integer, parameter :: ncipick_none = 0
  integer, parameter :: ncipick_region = 1 ! a region point (pickrow, in region mode pickmode)
  integer, parameter :: ncipick_molecule = 2 ! a fragment: the whole molecule of the clicked atom
  integer, parameter :: ncipick_atom = 3 ! a fragment: the clicked atom only

contains

  !> Draw the NCI plot window: the density (a field of the system
  !> shown by the anchor view), the region and the grid over it, the
  !> fragments between which the plot is made, and the cutoffs of the
  !> plot. Calculate runs the NCIPLOT kernel as a blocking job
  !> (run_nci) and loads the grids as two new fields of the system.
  module subroutine draw_nci(w)
    use systems, only: sys, sysc, sys_init, ok_system, field_cost
    use representations, only: iso_region_cell, iso_region_bbox, iso_region_to_box, iso_region_seed,&
       iso_region_point_from_cart, iso_isgridfield, iso_defaultlevel
    use nci, only: nci_options_set_field
    use types, only: field_evaluation_avail
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_setpos_bottomright,&
       iw_separatortext, iw_combo_simple, iw_dragfloat_real8, iw_checkbox, iw_field_combo,&
       iw_calcwidth, duration_string
    use gui_main, only: ColorHighlightScene
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, ismol, isgrid, okbox, okfield, okcalc, ch, navail, capped
    logical :: maskchanged, objchanged, lcommit
    integer :: isys, iview, istat, irow, ifield, nn(3), ihl, hltype, atidx(4)
    real*8 :: xc(3), box(3,0:3)

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)

    ! initialize the state on the first pass
    if (w%firstpass) then
       w%nc = nci_state()
       w%nc%ilevel = iso_defaultlevel
       w%errmsg = ""
    end if

    ! a new system: the reference field, on the cell (crystals) or
    ! around the atoms (molecules), no fragments; the same system with
    ! a new geometry: the fragments and the results refer to the old one
    if (goodsys) then
       ismol = sys(isys)%c%ismolecule
       if (.not.w%nc%sys%is(isys) .or. syschanged) then
          call w%nc%sys%set(isys)
          call cancel_pick()
          w%errmsg = ""
          if (ismol) then
             w%nc%iregion = iso_region_bbox
          else
             w%nc%iregion = iso_region_cell
          end if
          call iso_region_seed(isys,w%nc%iregion,w%nc%rgn_x)
          call set_field(sys(isys)%iref)
          call clear_fragments()
          call drop_results()
       elseif (w%nc%timelast /= sysc(isys)%timelastchange_geometry) then
          call clear_fragments()
          call drop_results()
       end if

       ! the reference field, if the field is gone
       if (.not.w%nc%fref%ok(sys(isys))) call set_field(sys(isys)%iref)

       ! the molecule selection follows the molecule list
       if (allocated(w%nc%molsel)) then
          if (size(w%nc%molsel,1) /= sys(isys)%c%nmol) deallocate(w%nc%molsel)
       end if
       if (.not.allocated(w%nc%molsel)) then
          allocate(w%nc%molsel(sys(isys)%c%nmol))
          w%nc%molsel = .false.
       end if
    end if

    ! a pending pick in the view it was armed on: a fragment, or a
    ! region point in the mode it was armed in; dropped, releasing the
    ! view, if the system is not ready, the window moved to another
    ! view, or the region mode changed
    if (w%nc%pickkind /= ncipick_none) then
       if (goodsys .and. iview == w%nc%pickview .and.&
          (w%nc%pickkind /= ncipick_region .or. w%nc%iregion == w%nc%pickmode)) then
          call view_pick_result(iview,w%id,isys,w%nc%pick,istat,xc,atidx=atidx)
          if (istat == ipick_point) then
             if (w%nc%pickkind == ncipick_region) then
                call iso_region_point_from_cart(isys,w%nc%iregion,xc,w%nc%rgn_x(:,w%nc%pickrow))
             else
                call add_pickfrag(atidx)
             end if
          end if
          if (istat /= ipick_pending) w%nc%pickkind = ncipick_none
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
          ttshown,irow,pickoff=(w%nc%pickkind /= ncipick_none))) w%errmsg = ""
       if (irow >= 0) then
          w%nc%pickkind = ncipick_region
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
       if (associated(w%nc%res)) then
          if (maskchanged) then
             call nci_refresh(w,isys,iview)
          elseif (objchanged) then
             call nci_objects(w,isys,iview)
          end if
       end if

       ! the results
       if (associated(w%nc%res)) call draw_results()

       ! the region in the anchor view, for as long as the window is open
       ! (not the cell of a crystal, which the cell shows)
       if (okbox .and. .not.(w%nc%iregion == iso_region_cell .and..not.ismol)) &
          call region_show_box(iview,w%uid,isys,box(:,0),box(:,1:3))

       if (allocated(w%errmsg)) then
          if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
       end if
    end if

    ! bottom row: calculate, close
    okcalc = goodsys
    if (okcalc) okcalc = okbox .and. okfield .and. all(nn > 0) .and. w%nc%pickkind == ncipick_none
    call iw_setpos_bottomright(len("Calculate") + len("Close"),2)
    if (iw_button("Calculate##ncicalc",danger=.true.,disabled=.not.okcalc)) call request_job()
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
    !> density range following the new density cutoff, and the native
    !> grid of a grid field on the whole cell (the CLI default).
    subroutine set_field(id)
      integer, intent(in) :: id

      integer :: id0
      real*8 :: rhoprev

      id0 = id
      if (.not.sys(isys)%goodfield(id0)) id0 = 0
      call w%nc%fref%set(sys(isys),id0)
      rhoprev = w%nc%opts%rhoplot
      call nci_options_set_field(w%nc%opts,sys(isys),id0)
      call srhorange_follow(rhoprev)
      if (iso_isgridfield(isys,id0) .and. w%nc%iregion == iso_region_cell) then
         if (.not.sys(isys)%f(id0)%grid%partial) w%nc%ilevel = 0
      end if

    end subroutine set_field

    !> Keep the signed density range inside the density cutoff after
    !> the cutoff changed from rhoprev: a bound that was at the old
    !> cutoff (no restriction) follows it to the new one.
    subroutine srhorange_follow(rhoprev)
      real*8, intent(in) :: rhoprev

      associate(o => w%nc%opts)
        if (o%srhorange(1) <= -rhoprev) o%srhorange(1) = -o%rhoplot
        if (o%srhorange(2) >= rhoprev) o%srhorange(2) = o%rhoplot
        o%srhorange = min(max(o%srhorange,-o%rhoplot),o%rhoplot)
      end associate

    end subroutine srhorange_follow

    !> Forget the fragments: the molecule selection and the picks.
    subroutine clear_fragments()

      if (allocated(w%nc%molsel)) deallocate(w%nc%molsel)
      if (allocated(w%nc%pickfr)) deallocate(w%nc%pickfr)

    end subroutine clear_fragments

    !> Forget the results of the last calculation (the fields it
    !> loaded stay in the system).
    subroutine drop_results()

      if (associated(w%nc%res)) deallocate(w%nc%res)
      w%nc%fsref = field_ref()
      w%nc%frref = field_ref()
      w%nc%isoref = rep_ref()
      w%nc%cbref = rep_ref()
      w%nc%voidok = .false.
      w%nc%timelast = sysc(isys)%timelastchange_geometry
      w%errmsg = ""

    end subroutine drop_results

    !> The number of fragments the form defines: the selected
    !> molecules and the picked fragments.
    function nfrag_total() result(nfrag)
      integer :: nfrag

      nfrag = 0
      if (allocated(w%nc%molsel)) nfrag = count(w%nc%molsel)
      if (allocated(w%nc%pickfr)) nfrag = nfrag + size(w%nc%pickfr,1)

    end function nfrag_total

    !> The Density section: the field, and (More) how a grid field is
    !> differentiated and whether the promolecular density is also
    !> calculated.
    subroutine draw_density()
      integer :: id, ic

      call iw_separatortext("Density",w%nc%sec(nci_sec_density),hasmore=.true.)
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

      if (w%nc%sec(nci_sec_density)%more) then
         if (iso_isgridfield(isys,w%nc%fref%id)) then
            call iw_text("Derivatives",alignframe=.true.)
            ic = merge(1,0,w%nc%opts%fourierint)
            call iw_combo_simple("##nciderivs","Interpolant" // c_null_char // "Fourier transform" //&
               c_null_char,ic,sameline=.true.,changed=ch)
            if (ch) w%nc%opts%fourierint = (ic == 1)
            call iw_tooltip("How the gradient and the Hessian of a grid field are calculated: by&
               & differentiating the interpolant (INTERPOLATION in NCIPLOT) or by Fourier transform&
               & of the grid (the default for grids without core densities)",ttshown)
         end if
         ch = iw_checkbox("Promolecular density##ncipromol",w%nc%opts%dopromol,labelleft=.true.)
         call iw_tooltip("Also calculate the promolecular density on the grid (the -pdens cube&
            & file of NCIPLOT)",ttshown)
      end if

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

    !> The Fragments section: the molecules of the system (selected in
    !> a table) and the fragments picked in the view; More, the
    !> fragment criterion parameters.
    subroutine draw_fragments()
      integer :: j, irm, nfrag

      call iw_separatortext("Fragments",w%nc%sec(nci_sec_fragments),hasmore=.true.)
      if (.not.w%nc%sec(nci_sec_fragments)%open) return

      ! the molecules of the system
      if (sys(isys)%c%nmol > 1) then
         ihl = 0
         hltype = 0
         ch = mol_table_widget(isys,ihl,hltype,shown=w%nc%molsel,select=.true.)
         if (ihl > 0) &
            call sysc(isys)%highlight_atoms(.true.,(/ihl/),hltype,reshape(ColorHighlightScene,(/4,1/)))
         if (iw_button("Select All##ncimolall")) w%nc%molsel = .true.
         call iw_tooltip("Make every molecule of the system a fragment",ttshown)
         if (iw_button("Clear##ncimolnone",sameline=.true.)) w%nc%molsel = .false.
         call iw_tooltip("Deselect all the molecules",ttshown)
      end if

      ! fragments picked in the view
      call iw_text("Pick",alignframe=.true.)
      if (iw_button("Molecule##ncipickmol",sameline=.true.,disabled=(w%nc%pickkind /= ncipick_none))) &
         call arm_fragpick(ncipick_molecule)
      call iw_tooltip("Click an atom in the view to add its whole molecule, at the lattice&
         & translation shown, as a fragment",ttshown,whendisabled=.true.)
      if (iw_button("Atom##ncipickatom",sameline=.true.,disabled=(w%nc%pickkind /= ncipick_none))) &
         call arm_fragpick(ncipick_atom)
      call iw_tooltip("Click an atom in the view to add it, on its own, as a fragment",ttshown,&
         whendisabled=.true.)
      if (allocated(w%nc%pickfr)) then
         if (iw_button("Remove All##ncipickclear",sameline=.true.)) deallocate(w%nc%pickfr)
         call iw_tooltip("Remove all the picked fragments",ttshown)
      end if
      ! the list, one row per fragment; a removal is applied after the
      ! loop so the rows keep their IDs for the frame
      irm = 0
      if (allocated(w%nc%pickfr)) then
         do j = 1, size(w%nc%pickfr,1)
            if (iw_button("Remove##ncipickrm" // string(j))) irm = j
            call iw_tooltip("Remove this fragment",ttshown)
            call iw_text(w%nc%pickfr(j)%label,sameline=.true.)
         end do
         if (irm > 0) w%nc%pickfr = (/w%nc%pickfr(1:irm-1),w%nc%pickfr(irm+1:)/)
      end if

      ! how many
      nfrag = nfrag_total()
      if (nfrag > 0) then
         call iw_text(string(nfrag) // " fragments: only the regions between them are plotted",&
            disabled=.true.,wrap=.true.)
      else
         call iw_text("No fragments: all the regions are plotted",disabled=.true.,wrap=.true.)
      end if

      if (w%nc%sec(nci_sec_fragments)%more) then
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
      end if

    end subroutine draw_fragments

    !> The Plot section: the cutoffs of the isosurface plot (RDG
    !> isovalue, density cutoff, signed density range, only negative)
    !> and, More, those of the s vs. rho plot and the void criterion.
    subroutine draw_plot()
      real*8 :: rhoprev

      call iw_separatortext("Plot",w%nc%sec(nci_sec_plot),hasmore=.true.)
      if (.not.w%nc%sec(nci_sec_plot)%open) return

      associate(o => w%nc%opts)
        call iw_text("Isovalue (s)",alignframe=.true.)
        ch = iw_dragfloat_real8("##ncidimplot",x1=o%dimplot,speed=0.005d0,decimal=3,min=0.001d0,&
           max=10d0,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,notlive=.true.,committed=lcommit)
        objchanged = objchanged .or. lcommit
        call iw_tooltip("Value of the reduced density gradient (s) of the isosurface (DIMPLOT)",ttshown)
        call iw_text("Density cutoff (a.u.)",alignframe=.true.,gap=.true.)
        rhoprev = o%rhoplot
        ch = iw_dragfloat_real8("##ncirhoplot",x1=o%rhoplot,speed=0.001d0,decimal=4,min=0.0001d0,&
           max=10d0,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,notlive=.true.,committed=lcommit)
        call iw_tooltip("Points with density above this value are left out of the isosurface,&
           & and the color scale of the signed density spans minus to plus this value (CUTPLOT)",&
           ttshown)
        if (ch) call srhorange_follow(rhoprev)
        maskchanged = maskchanged .or. lcommit

        call iw_text("Signed density range (a.u.)",alignframe=.true.)
        ch = iw_dragfloat_real8("##ncisrhorange",x2=o%srhorange,speed=0.001d0,decimal=4,&
           min=-o%rhoplot,max=o%rhoplot,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,&
           notlive=.true.,committed=lcommit)
        if (o%srhorange(2) < o%srhorange(1)) o%srhorange(2) = o%srhorange(1)
        maskchanged = maskchanged .or. lcommit
        call iw_tooltip("Points with signed density (the density times the sign of the second&
           & eigenvalue of its Hessian) outside this range are left out of the isosurface&
           & (SRHORANGE)",ttshown)
        ch = iw_checkbox("Only negative##ncionlyneg",o%onlyneg,gap=.true.,labelleft=.true.)
        maskchanged = maskchanged .or. ch
        call iw_tooltip("Keep only the points where the second eigenvalue of the Hessian is&
           & negative (ONLYNEG)",ttshown)

        if (w%nc%sec(nci_sec_plot)%more) then
           call iw_text("Plot cutoffs: density",alignframe=.true.)
           ! (the dat-file cutoffs do not touch the isosurface)
           ch = iw_dragfloat_real8("##ncirhocut",x1=o%rhocut,speed=0.005d0,decimal=3,&
              min=0d0,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
           call iw_tooltip("Points with density above this value are left out of the s vs.&
              & density plot (the dat file; CUTOFFS)",ttshown)
           call iw_text("s",alignframe=.true.,gap=.true.)
           ch = iw_dragfloat_real8("##ncidimcut",x1=o%dimcut,speed=0.005d0,decimal=3,&
              min=0d0,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
           call iw_tooltip("Points with reduced density gradient above this value are left out&
              & of the s vs. density plot (the dat file; CUTOFFS)",ttshown)
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

    !> The Results section: the grid of the last calculation, the
    !> time it took, and the fields it was loaded as.
    subroutine draw_results()
      character(len=:), allocatable :: str

      call iw_separatortext("Results",w%nc%sec(nci_sec_results))
      if (.not.w%nc%sec(nci_sec_results)%open) return

      associate(r => w%nc%res)
        str = "Grid: " // string(r%n(1)) // " x " // string(r%n(2)) // " x " // string(r%n(3)) //&
           " points, calculated in " // duration_string(w%nc%secs)
        if (sys(isys)%goodfield(r%ifield)) str = str // " from field " // string(r%ifield)
        call iw_text(str,wrap=.true.)
      end associate
      if (w%nc%fsref%ok(sys(isys))) &
         call iw_text("Reduced density gradient (masked by the cutoffs): field " // string(w%nc%fsref%id))
      if (w%nc%frref%ok(sys(isys))) &
         call iw_text("Signed density: field " // string(w%nc%frref%id))
      if (associated(win(iview)%sc)) then
         if (w%nc%isoref%ok(win(iview)%sc)) &
            call iw_text("Isosurface object: " // trim(win(iview)%sc%rep(w%nc%isoref%id)%name))
      end if
      if (w%nc%opts%dovoid .and. w%nc%voidok) then
         call iw_text("Void charge (a.u.): " // string(w%nc%voidq,'f',decimal=6) //&
            ", promolecular charge (a.u.): " // string(w%nc%voidqp,'f',decimal=6) //&
            ", volume (bohr^3): " // string(w%nc%voidv,'f',decimal=6),wrap=.true.)
      end if

    end subroutine draw_results

    !> Arm a fragment pick of kind (ncipick_molecule or ncipick_atom)
    !> in the anchor view.
    subroutine arm_fragpick(kind)
      integer, intent(in) :: kind

      w%nc%pickkind = kind
      w%nc%pickview = iview
      call w%nc%pick%arm()
      if (kind == ncipick_molecule) then
         call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick an atom: its molecule is added as a&
            & fragment",w%id)
      else
         call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick an atom: it is added as a fragment",&
            w%id)
      end if

    end subroutine arm_fragpick

    !> Add the fragment of the picked atom idx (cell index and lattice
    !> vector): the atom, or its whole molecule at that lattice
    !> translation (the pending pick kind).
    subroutine add_pickfrag(idx)
      integer, intent(in) :: idx(4)

      integer :: idc, im, k, dl(3)
      real*8 :: x(3,1)
      type(nci_pickfrag) :: pf

      idc = idx(1)
      if (idc == 0) then
         w%errmsg = "The picked point is not an atom"
         return
      end if
      associate(c => sys(isys)%c)
         if (w%nc%pickkind == ncipick_atom .or. .not.allocated(c%idatcelmol)) then
            x(:,1) = c%atcel(idc)%r + c%x2c(real(idx(2:4),8))
            pf%fr = c%identify_fragment(1,x)
            pf%label = "Atom " // string(idc) // " (" // trim(c%spc(c%atcel(idc)%is)%name) // ")"
         else
            ! the molecule, displaced so its copy of the atom is the one picked
            im = c%idatcelmol(1,idc)
            k = c%idatcelmol(2,idc)
            pf%fr = c%mol(im)
            dl = idx(2:4) - c%mol(im)%at(k)%lvec
            call pf%fr%translate(dl)
            pf%label = "Molecule " // string(im) // " (" // string(pf%fr%nat) // " atoms)"
            if (any(dl /= 0)) pf%label = pf%label // " at " // string(dl(1)) // "," //&
               string(dl(2)) // "," // string(dl(3))
         end if
      end associate
      if (allocated(w%nc%pickfr)) then
         w%nc%pickfr = (/w%nc%pickfr,pf/)
      else
         w%nc%pickfr = (/pf/)
      end if
      w%errmsg = ""

    end subroutine add_pickfrag

    !> Request the blocking job that calculates the NCI plot (run_nci)
    !> on the grid and fragments of the form: the grid is the nodes of
    !> the density field (native level) or nn points over the region
    !> (origin and steps, cell-frame Cartesian bohr; periodic over the
    !> whole cell of a crystal), the fragments are the selected
    !> molecules and the picked fragments. The far corner x1 of the
    !> options (the atoms the output files list) is set when the files
    !> are written.
    subroutine request_job()
      integer :: j, k, nfrag

      w%errmsg = ""
      if (.not.w%request_block()) then
         w%errmsg = "Another calculation is starting: press the button again"
         return
      end if
      associate(o => w%nc%opts, c => sys(isys)%c)
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

        ! the fragments
        nfrag = nfrag_total()
        if (allocated(o%fr)) deallocate(o%fr)
        allocate(o%fr(max(nfrag,1)))
        o%nfrag = 0
        if (allocated(w%nc%molsel)) then
           do j = 1, size(w%nc%molsel,1)
              if (.not.w%nc%molsel(j)) cycle
              o%nfrag = o%nfrag + 1
              o%fr(o%nfrag) = c%mol(j)
           end do
        end if
        if (allocated(w%nc%pickfr)) then
           do k = 1, size(w%nc%pickfr,1)
              o%nfrag = o%nfrag + 1
              o%fr(o%nfrag) = w%nc%pickfr(k)%fr
           end do
        end if
      end associate
      w%nc%jsys = isys

    end subroutine request_job

    !> Cancel a pending pick, releasing the view.
    subroutine cancel_pick()
      if (w%nc%pickkind == ncipick_none) return
      if (w%nc%pickview >= 1 .and. w%nc%pickview <= nwin) &
         call win(w%nc%pickview)%viewmode_release_forced(w%id)
      w%nc%pickkind = ncipick_none
    end subroutine cancel_pick

  end subroutine draw_nci

  !> Draw the overlay of the pending job of the NCI plot window.
  module subroutine block_nci(w)
    use systems, only: sys, sysc, ok_system, sys_init
    use utils, only: iw_wait_overlay
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    character(len=:), allocatable :: info
    integer :: isys

    info = ""
    isys = w%nc%jsys
    if (ok_system(isys,sys_init)) then
       info = "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (sys(isys)%goodfield(w%nc%opts%ifield)) &
          info = info // newline // "Field:  " // string(w%nc%opts%ifield) // ": " //&
          trim(sys(isys)%f(w%nc%opts%ifield)%name)
       info = info // newline // "Grid:   " // string(w%nc%opts%nstep(1)) // " x " //&
          string(w%nc%opts%nstep(2)) // " x " // string(w%nc%opts%nstep(3))
    end if
    call iw_wait_overlay("Calculating the NCI plot...",info,.true.)

  end subroutine block_nci

  !> Run the pending job of the NCI plot window, called by the main
  !> loop after the frame with the overlay: the reduced density
  !> gradient and the signed density on the grid of the options
  !> (request_job), then the masks, the fields and the objects
  !> (nci_refresh). The calculation can be cancelled with Esc, and
  !> then the previous results are gone and no field is loaded.
  module subroutine run_nci(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_cplist
    use nci, only: nci_compute
    use gui_main, only: begin_cancellable, end_cancellable
    use interfaces_glfw, only: glfwGetTime
    use tools_io, only: uout, string
    class(window), intent(inout), target :: w

    integer :: isys
    real*8 :: t
    logical :: cancelled
    character(len=:), allocatable :: errmsg

    isys = w%nc%jsys
    w%nc%jsys = 0
    if (.not.ok_system(isys,sys_init)) return
    if (.not.sys(isys)%goodfield(w%nc%opts%ifield)) then
       w%errmsg = "The density field is gone"
       return
    end if
    w%errmsg = ""

    ! the grids; a cancel or an error leaves no result
    if (.not.associated(w%nc%res)) allocate(w%nc%res)
    t = glfwGetTime()
    call begin_cancellable()
    call nci_compute(w%nc%opts,sys(isys),w%nc%res,errmsg)
    cancelled = end_cancellable()
    t = glfwGetTime() - t
    if (cancelled) errmsg = "The calculation was cancelled: no field was loaded"
    if (len_trim(errmsg) > 0) then
       w%errmsg = errmsg
       deallocate(w%nc%res)
       return
    end if
    w%nc%secs = t
    w%nc%timelast = sysc(isys)%timelastchange_geometry

    ! the masks, the fields and the objects
    call nci_refresh(w,isys,w%anchor_view())
    call sysc(isys)%post_event(lastchange_cplist)
    sysc(isys)%showfields = .true.
    if (w%nc%fsref%ok(sys(isys)) .and. w%nc%frref%ok(sys(isys))) &
       write (uout,'("* NCI plot grids loaded as fields ",A," (RDG) and ",A," (signed density)"/)') &
       string(w%nc%fsref%id), string(w%nc%frref%id)

  end subroutine run_nci

  !xx! private procedures

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
    use nci, only: nci_masks
    use tools_io, only: string
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    integer, intent(in) :: iview

    logical :: ok
    real*8, allocatable :: splot(:,:,:)
    character(len=:), allocatable :: errmsg

    if (.not.associated(w%nc%res)) return

    ! the masks (the void integrals of a failed pass are not shown)
    w%nc%voidok = .false.
    call nci_masks(w%nc%opts,w%nc%res,splot=splot,voidq=w%nc%voidq,voidqp=w%nc%voidqp,&
       voidv=w%nc%voidv,errmsg=errmsg)
    if (len_trim(errmsg) > 0) then
       w%errmsg = errmsg
       return
    end if
    w%errmsg = ""
    w%nc%voidok = .true.

    ! the fields
    associate(r => w%nc%res, o => w%nc%opts)
      if (.not.sys(isys)%goodfield(w%nc%frref%id,uid=w%nc%frref%uid,n=r%n)) &
         call nci_load_field(w,isys,w%nc%frref,"NCI sign(l2)*rho of field " // string(o%ifield),r%srho)
      ok = .false.
      if (sys(isys)%goodfield(w%nc%fsref%id,uid=w%nc%fsref%uid,n=r%n)) &
         ok = replace_grid_field(isys,w%nc%fsref%id,splot)
      if (.not.ok) call nci_load_field(w,isys,w%nc%fsref,"NCI RDG of field " // string(o%ifield),splot)
    end associate

    call nci_objects(w,isys,iview)

  end subroutine nci_refresh

  !> Make or update the objects of the NCI plot window w for system
  !> isys in the scene of view iview (nothing if it is not an open
  !> view of real space with a scene): an isosurface of the masked RDG
  !> field at the RDG isovalue, colored by the signed density field
  !> over minus to plus the density cutoff on the blue-green-red scale
  !> (blue = attractive, red = repulsive: the NCI convention), and
  !> its color bar.
  subroutine nci_objects(w,isys,iview)
    use systems, only: sys
    use scenes, only: space_real
    use representations, only: reptype_isosurface, repflavor_isosurface, reptype_colorbar,&
       repflavor_colorbar
    use utils, only: iw_cmap_jet
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    integer, intent(in) :: iview

    integer :: irep, icb

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
        call r%iso%slot(1)%set_map_field(isys,w%nc%frref%id,icmap=iw_cmap_jet,&
           maprange=(/-w%nc%opts%rhoplot,w%nc%opts%rhoplot/))
        r%shown = .true.
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

end submodule nci
