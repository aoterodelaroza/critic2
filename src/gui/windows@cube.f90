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

! The cube window (the CUBE keyword): the values of a field, one of
! its derivatives, an arithmetic expression, or a wavefunction, on a
! grid over a region, written to a cube, xsf, or VASP file, or
! loaded as a new field of the system.
submodule (windows) cube
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! more points than this are warned about (a large file, a long wait)
  integer*8, parameter :: npts_warn = 50000000_8

contains

  !> Draw the cube window: the values of a field of the system shown
  !> by the anchor view, one of its derivatives, an arithmetic
  !> expression, or a wavefunction of a QE grid field, on a grid over
  !> a region (the isosurface regions, or the nodes of a grid field),
  !> as CUBE calculates them, written to the file chosen in a dialog
  !> or loaded as a new field of the system.
  module subroutine draw_cube(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use fieldmod, only: type_grid
    use representations, only: iso_region_cell, iso_region_bbox, iso_region_to_box, iso_region_seed,&
       iso_region_point_from_cart, rep_shape, shapekind_box
    use rhoplot, only: cubepsi_none, cubepsi_psink, cubepsi_unk, cubepsi_wannier, cubepsi_mlwf
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_setpos_bottomright,&
       iw_separatortext, iw_combo_simple, iw_inputint, iw_inputint3, iw_dragfloat_real8, iw_checkbox,&
       duration_string
    use global, only: precisecube
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, syschanged, ismol, isgrid, isqe, okbox, okfield, ch, okcalc, okload, reqwrite
    integer :: isys, iview, istat, irow, ifield, i, nn(3)
    real*8 :: xc(3), box(3,0:3), prev0(3)
    type(field_ref) :: req

    logical, save :: ttshown = .false. ! tooltip flag

    ! the wavefunctions in the order of the combo (MLWF only if the
    ! field has a wannier90 checkpoint)
    integer, parameter :: psi_of_combo(0:4) = (/cubepsi_none,cubepsi_psink,cubepsi_unk,cubepsi_wannier,&
       cubepsi_mlwf/)

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)

    ! initialize the state on the first pass
    ! (but for a field asked for from the tree)
    if (w%firstpass) then
       req = w%cu%req
       reqwrite = w%cu%reqwrite
       w%cu = cube_state()
       w%cu%req = req
       w%cu%reqwrite = reqwrite
       w%errmsg = ""
    end if

    ! a new system: the reference field, on its own grid if it has one,
    ! or on the cell (crystals) or around the atoms (molecules)
    if (goodsys) then
       ismol = sys(isys)%c%ismolecule
       if (.not.w%cu%sys%is(isys) .or. syschanged) then
          call w%cu%sys%set(isys)
          call cancel_pick()
          w%errmsg = ""
          if (ismol) then
             w%cu%iregion = iso_region_bbox
          else
             w%cu%iregion = iso_region_cell
          end if
          call set_field(sys(isys)%iref)
          call iso_region_seed(isys,merge(iso_region_cell,w%cu%iregion,w%cu%iregion == region_native),&
             w%cu%rgn_x)
          w%okfile = okfile_default(isys,"grid","cube")
       end if

       ! a field asked for from the tree, on its own grid if it has one
       if (w%cu%req%uid > 0) then
          if (w%cu%req%ok(sys(isys))) then
             call set_field(w%cu%req%id)
             w%errmsg = ""
             if (w%cu%reqwrite) &
                i = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savegridfile,idparent=w%id,orraise=-1)
          end if
          w%cu%req = field_ref()
          w%cu%reqwrite = .false.
       end if

       ! the reference field, if the field is gone
       if (.not.w%cu%fref%ok(sys(isys))) then
          call w%cu%fref%set(sys(isys),sys(isys)%iref)
          w%cu%costest = -1d0
       end if
    end if

    ! a pending pick in the view it was armed on; dropped, releasing the
    ! view, if the system is not ready, the window moved to another
    ! view, or the region mode changed
    if (w%cu%picking) then
       if (goodsys .and. iview == w%cu%pickview .and. w%cu%iregion == w%cu%pickmode) then
          call view_pick_result(iview,w%id,isys,w%cu%pick,istat,xc)
          if (istat == ipick_point) &
             call iso_region_point_from_cart(isys,w%cu%iregion,xc,w%cu%rgn_x(:,w%cu%pickrow))
          if (istat /= ipick_pending) w%cu%picking = .false.
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

       ! the quantity
       call iw_separatortext("Property",w%cu%sec(cube_sec_property),hasmore=.not.w%cu%isexpr)
       if (w%cu%sec(cube_sec_property)%open) then
          ch = quantity_form(isys,w%cu%isexpr,w%cu%fref,w%cu%nti,w%cu%expr,w%cu%exprerr,&
             w%cu%sec(cube_sec_property)%more,"cube","the grid","the grid",ttshown)
          ! another quantity costs differently
          if (ch) w%cu%costest = -1d0
       end if
       ifield = -1
       if (.not.w%cu%isexpr) ifield = w%cu%fref%id
       isgrid = nativefield(ifield)
       isqe = .false.
       if (isgrid) isqe = sys(isys)%f(ifield)%grid%isqe
       if (.not.isqe) w%cu%ipsi = cubepsi_none
       if (w%cu%isexpr) then
          okfield = (len_trim(w%cu%expr) > 0 .and. len_trim(w%cu%exprerr) == 0)
       else
          okfield = sys(isys)%goodfield(ifield)
       end if

       ! the wavefunction, for a QE grid field; it is written as it is
       ! (no derivatives, as in CUBE)
       if (isqe) call draw_wavefunction()
       if (w%cu%ipsi /= cubepsi_none) w%cu%nti = 0

       ! the region: a wavefunction is on the nodes of the grid
       if (w%cu%ipsi /= cubepsi_none) w%cu%iregion = region_native
       if (region_form(isys,w%cu%iregion,w%cu%rgn_x,w%cu%sec(cube_sec_region),"cube","the grid",ttshown,&
          irow,pickoff=w%cu%picking,nativeok=(isgrid .and. w%cu%nti == 0))) then
          w%errmsg = ""
          w%cu%costest = -1d0
       end if
       if (w%cu%ipsi /= cubepsi_none .and. w%cu%sec(cube_sec_region)%open) &
          call iw_text("A wavefunction is written on the nodes of the grid of its field",disabled=.true.,&
             wrap=.true.)
       if (irow >= 0) then
          w%cu%picking = .true.
          w%cu%pickrow = irow
          w%cu%pickmode = w%cu%iregion
          w%cu%pickview = iview
          call w%cu%pick%arm()
          call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the region point: an atom,&
             & a critical point, or a point in space",w%id,acceptempty=.true.)
       end if
       if (w%cu%iregion == region_native) then
          okbox = .true.
       else
          call iso_region_to_box(isys,w%cu%iregion,w%cu%rgn_x,box,okbox)
          if (.not.okbox) &
             call iw_text("The region is degenerate (zero volume)",danger=.true.,wrap=.true.)
       end if

       ! the grid
       if (okbox) nn = grid_size()
       call draw_grid()

       ! the region in the anchor view, for as long as the window is open
       ! (not the cell of a crystal or the grid, which the cell shows)
       if (okbox .and. w%cu%iregion /= region_native .and.&
          .not.(w%cu%iregion == iso_region_cell .and..not.ismol)) then
          if (associated(win(iview)%sc)) then
             if (win(iview)%sc%isinit /= 0) then
                prev0 = box(:,0)
                if (ismol) prev0 = prev0 + sys(isys)%c%molx0
                call win(iview)%sc%show_transient_shapes(w%uid,1,(/rep_shape(kind=shapekind_box,&
                   x1=prev0,v=box(:,1:3),rad=region_edgerad,rgb=region_rgb,alpha=region_alpha)/))
             end if
          end if
       end if

       if (allocated(w%errmsg)) then
          if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
       end if

       ! the file chosen in the dialog
       if (w%okfile_set .and. okbox .and. okfield) call request_job(.false.)
    end if
    w%okfile_set = .false.

    ! bottom row: write the file, load the field, close; a field
    ! cannot be on a supercell of the system or have only the header,
    ! and in a crystal it spans the cell (a grid field on part of it is
    ! only evaluated inside the main cell)
    okcalc = goodsys
    if (okcalc) okcalc = okbox .and. okfield
    okload = okcalc .and..not.w%cu%header
    if (okload .and. w%cu%iregion == region_native) then
       okload = .not.(w%cu%dosuper .and. w%cu%ipsi == cubepsi_none .and. any(w%cu%nsuper > 1))
    elseif (okload) then
       okload = sys(isys)%c%ismolecule .or. w%cu%iregion == iso_region_cell
    end if
    call iw_setpos_bottomright(len("Write to File") + len("Load as Field") + len("Close"),3)
    if (iw_button("Write to File##cubewrite",danger=.true.,disabled=.not.okcalc)) &
       i = stack_create_window(wintype_dialog,.true.,wpurp_dialog_savegridfile,idparent=w%id,orraise=-1)
    call iw_tooltip("Calculate the grid and write it to a file: a cube file for any region,&
       & or an XSF or VASP file for the whole cell",ttshown,whendisabled=.true.)
    if (iw_button("Load as Field##cubeload",danger=.true.,sameline=.true.,disabled=.not.okload)) &
       call request_job(.true.)
    call iw_tooltip("Calculate the grid and load it as a new field of the system (not with Header&
       & only or a supercell; in a crystal, only on the unit cell or the native grid)",ttshown,&
       whendisabled=.true.)
    if (iw_button("Close",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) then
       call cancel_pick()
       call w%end()
    end if

  contains

    !> Whether field id of the system is a periodic grid field, whose
    !> nodes can be written as they are (CUBE GRID).
    function nativefield(id) result(ok)
      integer, intent(in) :: id
      logical :: ok

      ok = sys(isys)%goodfield(id)
      if (ok) ok = (sys(isys)%f(id)%type == type_grid .and. allocated(sys(isys)%f(id)%grid))
      if (ok) ok = .not.sys(isys)%f(id)%grid%partial

    end function nativefield

    !> Write field id of the system as it is (no derivative, no
    !> wavefunction), on its own grid if it has one.
    subroutine set_field(id)
      integer, intent(in) :: id

      call w%cu%fref%set(sys(isys),id)
      w%cu%isexpr = .false.
      w%cu%nti = 0
      w%cu%ipsi = cubepsi_none
      w%cu%costest = -1d0
      if (nativefield(id)) w%cu%iregion = region_native

    end subroutine set_field

    !> The Wavefunction section: the field itself, or one of its Bloch
    !> or Wannier functions (CUBE MLWF, WANNIER, UNK, PSINK), its band,
    !> k-point, spin, lattice translation, and part.
    subroutine draw_wavefunction()
      integer :: ic, nopt, jpsi
      character(len=:), allocatable :: str

      associate(qe => sys(isys)%f(ifield)%grid%qe)
        ! the band, k-point, and spin within this field's, also if the
        ! section is closed (the field may have changed)
        w%cu%ibnd = max(min(w%cu%ibnd,qe%nbnd),1)
        w%cu%ik = max(min(w%cu%ik,qe%nks),1)
        if (qe%nspin /= 2) w%cu%ispin = 1

        call iw_separatortext("Wavefunction",w%cu%sec(cube_sec_wfn))
        if (.not.w%cu%sec(cube_sec_wfn)%open) return
        str = "None (the field)" // c_null_char // "Bloch function" // c_null_char //&
           "Periodic part of the Bloch function" // c_null_char // "Wannier function" // c_null_char
        nopt = 4
        if (sys(isys)%f(ifield)%grid%iswan) then
           str = str // "Maximally localized Wannier function" // c_null_char
           nopt = 5
        end if
        ic = 0
        do jpsi = 0, nopt-1
           if (psi_of_combo(jpsi) == w%cu%ipsi) ic = jpsi
        end do
        call iw_text("Kind",alignframe=.true.)
        call iw_combo_simple("##cubepsi",str,ic,sameline=.true.)
        w%cu%ipsi = psi_of_combo(max(min(ic,nopt-1),0))
        call iw_tooltip("Write the field itself, or one of the wavefunctions in its QE file: a Bloch&
           & function (psi_nk), its periodic part (u_nk), a Wannier function, or a maximally&
           & localized Wannier function (with a wannier90 checkpoint)",ttshown)
        if (w%cu%ipsi == cubepsi_none) return

        call iw_text("Band",alignframe=.true.)
        ch = iw_inputint("##cubebnd",w%cu%ibnd,width=5,sameline=.true.)
        w%cu%ibnd = max(min(w%cu%ibnd,qe%nbnd),1)
        call iw_tooltip("Band of the wavefunction (1 to " // string(qe%nbnd) // ")",ttshown)
        if (w%cu%ipsi == cubepsi_psink .or. w%cu%ipsi == cubepsi_unk) then
           call iw_text("k-point",alignframe=.true.,gap=.true.)
           ch = iw_inputint("##cubek",w%cu%ik,width=5,sameline=.true.)
           w%cu%ik = max(min(w%cu%ik,qe%nks),1)
           call iw_tooltip("k-point of the Bloch function (1 to " // string(qe%nks) // ")",ttshown)
        end if
        if (qe%nspin == 2) then
           call iw_text("Spin",alignframe=.true.,gap=.true.)
           ch = iw_inputint("##cubespin",w%cu%ispin,width=3,sameline=.true.)
           w%cu%ispin = max(min(w%cu%ispin,2),1)
           call iw_tooltip("Spin of the wavefunction (1 = up, 2 = down)",ttshown)
        end if
        if (w%cu%ipsi /= cubepsi_unk) then
           call iw_text("Lattice translation",alignframe=.true.)
           ch = iw_inputint3("##cubeinr",w%cu%inr,width=3*3,sameline=.true.)
           call iw_tooltip("Lattice vector (in crystallographic coordinates) the wavefunction is&
              & displaced by",ttshown)
        end if
        call iw_text("Part",alignframe=.true.)
        call iw_combo_simple("##cubepart","Absolute value" // c_null_char // "Real part" // c_null_char //&
           "Imaginary part" // c_null_char,w%cu%ipart,sameline=.true.)
        call iw_tooltip("The part of the complex wavefunction that is written",ttshown)
      end associate

    end subroutine draw_wavefunction

    !> The number of points of the grid along each axis: the nodes of
    !> the field (times the supercell), the number given, or the length
    !> of each side over the spacing (as CUBE does).
    function grid_size() result(n)
      integer :: n(3)

      integer :: j

      if (w%cu%iregion == region_native) then
         n = sys(isys)%f(ifield)%grid%n
         if (w%cu%dosuper .and. w%cu%ipsi == cubepsi_none) n = n * w%cu%nsuper
      elseif (w%cu%igrid == igrid_npts) then
         n = max(w%cu%nn,2)
      else
         do j = 1, 3
            n(j) = nint(norm2(box(:,j)) / (w%cu%spacing / bohrtoa)) + 1
         end do
      end if

    end function grid_size

    !> The Grid section: the number of points or their spacing (the
    !> nodes of the field and, More, their shift and a supercell, for
    !> the native grid), how many points that is and the size of the
    !> file. More: only the header.
    subroutine draw_grid()
      integer*8 :: ntot
      real*8 :: mb

      call iw_separatortext("Grid",w%cu%sec(cube_sec_grid),hasmore=.true.)
      if (.not.w%cu%sec(cube_sec_grid)%open) return

      if (w%cu%iregion == region_native) then
         if (w%cu%sec(cube_sec_grid)%more) then
            call iw_text("Shift",alignframe=.true.)
            ch = iw_inputint3("##cubeshift",w%cu%ishift,width=4*3,sameline=.true.)
            call iw_tooltip("Displace the nodes by these numbers of points along each lattice vector&
               & (cube files only)",ttshown)
            if (w%cu%ipsi /= cubepsi_none) call igBeginDisabled(.true._c_bool)
            ch = iw_checkbox("Supercell##cubesuper",w%cu%dosuper,gap=.true.,labelleft=.true.)
            call iw_tooltip("Repeat the grid along the lattice vectors, in a supercell (not for&
               & wavefunctions)",ttshown,whendisabled=.true.)
            if (w%cu%dosuper) then
               ch = iw_inputint3("##cubensuper",w%cu%nsuper,width=3*3,sameline=.true.)
               w%cu%nsuper = max(w%cu%nsuper,1)
               call iw_tooltip("Number of copies of the cell along each lattice vector",ttshown,&
                  whendisabled=.true.)
            end if
            if (w%cu%ipsi /= cubepsi_none) call igEndDisabled()
         end if
      else
         call iw_text("Given as",alignframe=.true.)
         call iw_combo_simple("##cubegrid","Number of points" // c_null_char // "Spacing" // c_null_char,&
            w%cu%igrid,sameline=.true.)
         call iw_tooltip("Give the number of points along each side of the region, or the distance&
            & between neighboring points",ttshown)
         if (w%cu%igrid == igrid_npts) then
            call iw_text("Points",alignframe=.true.,gap=.true.)
            ch = iw_inputint3("##cubenn",w%cu%nn,width=5*3,sameline=.true.)
            w%cu%nn = max(w%cu%nn,2)
            call iw_tooltip("Number of points along each side of the region",ttshown)
         else
            call iw_text("Spacing (Å)",alignframe=.true.,gap=.true.)
            ch = iw_dragfloat_real8("##cubespacing",x1=w%cu%spacing,speed=0.002d0,decimal=3,min=0.005d0,&
               sameline=.true.,flags=ImGuiSliderFlags_AlwaysClamp)
            call iw_tooltip("Distance between neighboring points, approximately: each side gets the&
               & whole number of points closest to it",ttshown)
         end if
      end if

      ! the size of the grid and of the file
      if (all(nn > 0)) then
         ntot = product(int(nn,8))
         mb = real(ntot,8) * merge(23.2d0,13.2d0,precisecube) / 1d6
         call iw_text(string(nn(1)) // " x " // string(nn(2)) // " x " // string(nn(3)) // " = " //&
            string(ntot) // " points, about " // string(mb,'f',decimal=1) // " MB as a cube file",&
            disabled=.true.,wrap=.true.)
         if (ntot > npts_warn) &
            call iw_text("This is a very large grid: it may take long to calculate and write",danger=.true.,&
               wrap=.true.)
      end if

      ! the time the grid takes (not for the nodes of a field or for the
      ! header only, which are not calculated)
      if (w%cu%iregion /= region_native .and..not.w%cu%header) then
         if (iw_button("Estimate Cost##cubeestimate",disabled=.not.(okbox .and. okfield))) &
            call estimate_cost()
         call iw_tooltip("Estimate the time needed to calculate the grid with the selected options",ttshown,&
            whendisabled=.true.)
         if (w%cu%costest >= 0d0 .and. all(nn > 0)) &
            call iw_text("~" // duration_string(w%cu%costest * real(product(int(nn,8)),8)),sameline=.true.,&
               disabled=.true.)
      end if

      if (w%cu%sec(cube_sec_grid)%more) then
         ch = iw_checkbox("Header only##cubeheader",w%cu%header,labelleft=.true.)
         call iw_tooltip("Write only the header of the file (the cell and the atoms), not the values&
            & (not for binary cube files)",ttshown)
      end if

    end subroutine draw_grid

    !> Measure the cost (seconds per point) of calculating the quantity
    !> on the region, into w%cu%costest: the grid is calculated on a
    !> coarse grid over the same box, as the job does, in batches of
    !> increasing size until the timing is meaningful. An expression
    !> on the unit cell is first tried on the grids of its fields, as
    !> the job does; that is timed on the whole grid.
    subroutine estimate_cost()
      use interfaces_glfw, only: glfwGetTime
      use rhoplot, only: rhoplot_cube_values
      use grid3mod, only: grid3
      use iso_c_binding, only: c_loc
      !$ use omp_lib, only: omp_get_max_threads
      integer :: j, nb, np(3), nthr
      real*8 :: t, xdp(3,3)
      real*8, allocatable :: lf(:,:,:)
      character(len=:), allocatable :: errmsg
      type(grid3) :: faux

      real*8, parameter :: timetarget = 0.1d0 ! measure until a batch takes this long
      integer, parameter :: nbmax = 40 ! largest batch: nbmax**3 points

      w%cu%costest = -1d0
      w%errmsg = ""
      if (w%cu%isexpr .and. w%cu%iregion == iso_region_cell) then
         t = glfwGetTime()
         call faux%new_eval(c_loc(sys(isys)),c_loc(sys(isys)%c),nn,trim(w%cu%expr),sys(isys)%c%m_x2c)
         t = glfwGetTime() - t
         if (faux%isinit) then
            w%cu%costest = t / real(product(int(nn,8)),8)
            return
         end if
      end if

      ! the kernel runs its threads over the third axis: give it as
      ! many planes as the job has threads
      nthr = 1
      !$ nthr = omp_get_max_threads()
      nb = 4
      do while (.true.)
         np = min(nb,nn)
         np(3) = min(max(nb,nthr),nn(3))
         do j = 1, 3
            xdp(:,j) = box(:,j) / real(np(j),8)
         end do
         t = glfwGetTime()
         call rhoplot_cube_values(sys(isys),box(:,0),xdp,np,ifield,trim(w%cu%expr),w%cu%nti,.false.,lf,errmsg)
         t = glfwGetTime() - t
         if (len_trim(errmsg) > 0) then
            w%errmsg = errmsg
            return
         end if
         if (t >= timetarget .or. nb >= nbmax .or. all(np == nn)) exit
         nb = min(2 * nb,nbmax)
      end do
      w%cu%costest = t / real(product(np),8)

    end subroutine estimate_cost

    !> Request the blocking job that calculates the grid and loads it
    !> as a field (load) or writes it to the file chosen in the dialog
    !> (w%okfile, checked against the region) (run_cube), with the
    !> system, the field (-1 = the expression) and the points of the
    !> region (the native grid is set up by the job).
    subroutine request_job(load)
      use rhoplot, only: rhoplot_cube_outform, cubefmt_cube, cubefmt_bincube
      logical, intent(in) :: load
      integer :: outform, j

      w%errmsg = ""
      if (.not.load) then
         outform = rhoplot_cube_outform(w%okfile)
         if (w%cu%iregion /= region_native .and. w%cu%iregion /= iso_region_cell .and.&
            outform /= cubefmt_cube .and. outform /= cubefmt_bincube) then
            w%errmsg = "XSF and VASP files hold grids on the whole cell: choose the unit cell or the native&
               & grid as the region, or write a cube file"
            return
         end if
      end if
      if (.not.w%request_block()) return
      w%cu%jload = load
      w%cu%jsys = isys
      w%cu%jfield = ifield
      w%cu%jnn = nn
      if (w%cu%iregion /= region_native) then
         w%cu%jx0 = box(:,0)
         do j = 1, 3
            w%cu%jxd(:,j) = box(:,j) / real(nn(j),8)
         end do
      end if

    end subroutine request_job

    !> Cancel a pending pick, releasing the view.
    subroutine cancel_pick()
      if (.not.w%cu%picking) return
      if (w%cu%pickview >= 1 .and. w%cu%pickview <= nwin) &
         call win(w%cu%pickview)%viewmode_release_forced(w%id)
      w%cu%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_cube

  !> Draw the overlay of the pending job of the cube window.
  module subroutine block_cube(w)
    use systems, only: sys, sysc, ok_system, sys_init
    use utils, only: iw_wait_overlay
    use tools_io, only: string
    use param, only: newline
    class(window), intent(inout), target :: w

    character(len=:), allocatable :: info
    integer :: isys

    info = ""
    isys = w%cu%jsys
    if (ok_system(isys,sys_init)) then
       info = "System: " // string(isys) // ": " // trim(sysc(isys)%seed%name)
       if (sys(isys)%goodfield(w%cu%jfield)) &
          info = info // newline // "Field:  " // string(w%cu%jfield) // ": " //&
          trim(sys(isys)%f(w%cu%jfield)%name)
    end if
    call iw_wait_overlay("Calculating the grid...",info,.not.w%cu%header)

  end subroutine block_cube

  !> Run the pending job of the cube window, called by the main loop
  !> after the frame with the overlay: calculate the grid (the nodes
  !> of the field, maybe in a supercell, or the points of the region
  !> given by request_job) and load it as a new field of the system
  !> (jload) or write it to w%okfile, in the format of its
  !> extension. The calculation can be cancelled with Esc, and then
  !> no field is loaded and no file is written.
  module subroutine run_cube(w)
    use systems, only: sys, sysc, ok_system, sys_init, lastchange_cplist
    use iso_c_binding, only: c_loc
    use crystalmod, only: crystal
    use rhoplot, only: rhoplot_cube_values, rhoplot_cube_gridcheck, rhoplot_cube_gridvalues,&
       rhoplot_cube_write, rhoplot_cube_outform, cubepsi_none, cubepsi_unk
    use representations, only: iso_region_cell
    use gui_main, only: begin_cancellable, end_cancellable
    use tools_io, only: uout, fopen_write, fopen_append, fclose, string
    class(window), intent(inout), target :: w

    integer :: isys, ifield, ia(3), na(3), lu, j, ishift(3), inr(3), nn(3), id
    real*8 :: x0(3), xd(3,3), x00(3,3), x2cl(3,3)
    real*8, allocatable :: lf(:,:,:)
    logical :: native, existed, super
    character(len=:), allocatable :: errmsg, name
    type(crystal), pointer :: cr
    type(crystal), target :: nc

    isys = w%cu%jsys
    ifield = w%cu%jfield
    if (.not.ok_system(isys,sys_init)) return
    if (ifield >= 0) then
       if (.not.sys(isys)%goodfield(ifield)) return
    end if
    native = (w%cu%iregion == region_native)
    errmsg = ""

    ! the file can be written (the writers stop the program if not),
    ! checked without changing it; if it is new, it is removed again
    ! when nothing is written
    existed = .true. ! loading: no file to remove on failure
    if (.not.w%cu%jload) then
       inquire(file=trim(w%okfile),exist=existed)
       if (existed) then
          lu = fopen_append(trim(w%okfile),errstop=.false.)
       else
          lu = fopen_write(trim(w%okfile),errstop=.false.)
       end if
       if (lu < 0) then
          w%errmsg = "Could not open the file for writing: " // trim(w%okfile)
          return
       end if
       call fclose(lu)
    end if

    ! the grid: the nodes of the field (maybe in a supercell), or the
    ! points of the region
    cr => sys(isys)%c
    nn = w%cu%jnn
    x0 = w%cu%jx0
    xd = w%cu%jxd
    ia = (/1,2,3/)
    na = -1
    if (native) then
       super = w%cu%dosuper .and. w%cu%ipsi == cubepsi_none .and. any(w%cu%nsuper > 1)
       errmsg = rhoplot_cube_gridcheck(sys(isys),ifield,w%cu%ipsi,w%cu%ibnd,super)
       if (len_trim(errmsg) > 0) then
          call fail(errmsg)
          return
       end if
       if (super) then
          na = w%cu%nsuper
          x00 = 0d0
          do j = 1, 3
             x00(j,j) = na(j)
          end do
          nc = sys(isys)%c
          call nc%newcell(x00,noenv=.true.,errmsg=errmsg)
          if (len_trim(errmsg) > 0) then
             call fail(errmsg)
             return
          end if
          cr => nc
       end if
       x0 = 0d0
       do j = 1, 3
          xd(:,j) = cr%m_x2c(:,j) / real(nn(j),8)
       end do
    end if

    ! the values, unless only the header is written; the periodic part
    ! of the Bloch function takes no lattice translation (as in CUBE)
    if (.not.w%cu%header) then
       inr = w%cu%inr
       if (w%cu%ipsi == cubepsi_unk) inr = 0
       call begin_cancellable()
       if (native) then
          call rhoplot_cube_gridvalues(sys(isys),ifield,ia,na,w%cu%ipsi,w%cu%ibnd,w%cu%ik,inr,&
             w%cu%ispin,w%cu%ipart,lf)
       elseif (ifield < 0) then
          call rhoplot_cube_values(sys(isys),x0,xd,nn,-1,trim(w%cu%expr),0,&
             w%cu%iregion == iso_region_cell,lf,errmsg)
       else
          call rhoplot_cube_values(sys(isys),x0,xd,nn,ifield,"",w%cu%nti,.false.,lf,errmsg)
       end if
       if (end_cancellable()) then
          if (w%cu%jload) then
             call fail("The calculation was cancelled: no field was loaded")
          else
             call fail("The calculation was cancelled: no file was written")
          end if
          return
       end if
       if (len_trim(errmsg) > 0) then
          call fail(errmsg)
          return
       end if
    end if

    ! the field: on the cell (the nodes of the field) or the region
    ! (cryst. origin, Cartesian span), named as CLI LOAD AS names it;
    ! made in its slot, so the grid is not copied again
    if (w%cu%jload) then
       do j = 1, 3
          x2cl(:,j) = xd(:,j) * real(nn(j),8)
       end do
       if (ifield < 0) then
          name = trim(w%cu%expr)
       else
          name = "$" // string(ifield) // trim(quant_short(w%cu%nti))
          if (w%cu%ipsi /= cubepsi_none) name = name // " (band " // string(w%cu%ibnd) // ", k-point " //&
             string(w%cu%ik) // ", spin " // string(w%cu%ispin) // ")"
       end if
       name = "<generated>, grid: " // name
       id = sys(isys)%getfieldnum()
       call sys(isys)%f(id)%load_grid_from_array3(sys(isys)%c,id,name,lf,matmul(sys(isys)%c%m_c2x,x0),x2cl)
       deallocate(lf)
       sys(isys)%f(id)%file = "<generated>"
       sys(isys)%f(id)%sptr = c_loc(sys(isys))
       call sys(isys)%set_reference(id,.true.)
       w%errmsg = ""
       call sysc(isys)%post_event(lastchange_cplist)
       sysc(isys)%showfields = .true.
       write (uout,'("* Grid loaded as field ",A,": ",A/)') string(id), name
       return
    end if

    ! the file
    ishift = 0
    if (native) ishift = w%cu%ishift
    call rhoplot_cube_write(cr,rhoplot_cube_outform(w%okfile),lf,trim(w%okfile),w%cu%header,xd,x0,&
       ishift,errmsg)
    if (len_trim(errmsg) > 0) then
       call fail(errmsg)
       return
    end if
    w%errmsg = ""
    write (uout,'("* CUBE written to file: ",A/)') string(w%okfile)
    call okfile_save_dir(w%okfile)

  contains

    !> Report msg as the reason nothing was done, and remove the file
    !> if it is the empty one run_cube made.
    subroutine fail(msg)
      character(len=*), intent(in) :: msg

      integer :: lu, ios

      w%errmsg = msg
      if (existed) return
      open(newunit=lu,file=trim(w%okfile),status='old',iostat=ios)
      if (ios == 0) close(lu,status='delete')

    end subroutine fail

  end subroutine run_cube

end submodule cube
