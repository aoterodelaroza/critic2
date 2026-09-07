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

! Routines for the crystal voids window.
submodule (windows) voids
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! the atomic radii used in the packing tab
  integer, parameter :: vdrad_vdw = 0 ! van der Waals radii
  integer, parameter :: vdrad_cov = 1 ! covalent radii
  integer, parameter :: vdrad_nnm = 2 ! half the distance to the nearest neighbor

  ! bohr^3 to Å^3, for the volumes reported by every tab
  real*8, parameter :: fac3 = bohrtoa**3

  ! grids larger than this are refused
  integer, parameter :: maxgridpts = 2000000

  ! a run expected to last longer than this says so next to the button, and
  ! a grid of more than bigwarngridpts points does too when the cost of a
  ! point could not be measured
  real*8, parameter :: bigwarn_secs = 3d0
  integer, parameter :: bigwarngridpts = 500000

contains

  !> Draw the crystal voids window: measure the empty space in the system
  !> shown by the anchor view with three different definitions of a void,
  !> one per tab.
  module subroutine draw_voids(w)
    use systems, only: sysc, sys, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_close_event, iw_begintabitem,&
       iw_setpos_bottomright
    use tools_io, only: string
    class(window), intent(inout), target :: w

    logical :: doquit, goodsys, iscrystal, syschanged
    integer :: isys, iview
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: str1

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)
    iscrystal = goodsys
    if (iscrystal) iscrystal = .not.sys(isys)%c%ismolecule

    ! initialize the form on the first pass
    if (w%firstpass) then
       w%vd = voids_state()
       w%errmsg = ""
    end if

    ! The results describe one geometry of one system: drop them when either
    ! changes underneath them.
    if (goodsys) then
       if (w%vd%isys /= isys) then
          ! a different system: everything measured for the old one goes
          w%vd%pol_ic = 0
          w%vd%iso_secs = -1d0
          w%vd%iso_secs_ncel = -1
          w%vd%isys = isys
          call drop_results()
       elseif (w%vd%timelast /= sysc(isys)%timelastchange_geometry) then
          ! the same system moved: the results are stale but the cost of a
          ! grid point is not
          call drop_results()
       end if
    end if

    ! the system the voids are measured in
    if (goodsys) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)
    end if

    ! all three definitions need a periodic cell to measure the empty space in
    if (goodsys .and. .not.iscrystal) &
       call iw_text("The void analysis can only be done in crystals",danger=.true.,wrap=.true.)

    ! the three ways of measuring the empty space, one per tab
    if (iscrystal) then
       str1 = "##drawvoids_tabbar" // c_null_char
       flags = ImGuiTabBarFlags_None
       if (igBeginTabBar(c_loc(str1),flags)) then
          if (iw_begintabitem("Isosurface##drawvoids_isotab")) then
             call draw_isosurface_tab(w,isys,ttshown)
             call igEndTabItem()
          end if
          call iw_tooltip("The voids are the regions where the promolecular density is low",ttshown)

          if (iw_begintabitem("Polyhedra##drawvoids_poltab")) then
             call draw_polyhedra_tab(w,isys,ttshown)
             call igEndTabItem()
          end if
          call iw_tooltip("The empty space is whatever the coordination polyhedra do not cover",ttshown)

          if (iw_begintabitem("Packing##drawvoids_pcktab")) then
             call draw_packing_tab(w,isys,ttshown)
             call igEndTabItem()
          end if
          call iw_tooltip("The empty space is whatever falls outside the atomic spheres",ttshown)
          call igEndTabBar()
       end if
    end if

    ! the error message from the last calculation, whichever tab ran it
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    ! close button
    call iw_setpos_bottomright(5,1)
    if (iw_button("Close")) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  contains
    ! Forget what every tab calculated, and the message that came with it.
    subroutine drop_results()

      w%vd%iso_done = .false.
      w%vd%pol_done = .false.
      w%vd%pck_done = .false.
      w%vd%timelast = sysc(isys)%timelastchange_geometry
      w%errmsg = ""

    end subroutine drop_results

  end subroutine draw_voids

  !xx! private procedures

  !> Draw the isosurface tab: the voids are the connected regions where the
  !> promolecular density is lower than the isovalue.
  subroutine draw_isosurface_tab(w,isys,ttshown)
    use systems, only: sys
    use representations, only: iso_estimate_cost, iso_region_cell
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_calcheight,&
       iw_table_column, duration_string
    use tools_io, only: string, ioj_right
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    logical, intent(inout) :: ttshown

    logical :: toobig
    integer :: i, n(3)
    integer*8 :: npts
    real*8 :: tcost, rdum, xdum(3,0:3)
    integer(c_int) :: flags
    character(kind=c_char,len=:), allocatable, target :: str1, s
    type(ImVec2) :: sz0
    type(c_ptr), target :: clipper
    type(ImGuiListClipper), pointer :: clipper_f

    integer, parameter :: ic_iso_id = 0
    integer, parameter :: ic_iso_vol = 1
    integer, parameter :: ic_iso_pct = 2
    integer, parameter :: ic_iso_x = 3
    integer, parameter :: ic_iso_rho = 4

    ! the isovalue that separates a void from the rest of the cell. Changing
    ! either of the two settings makes the results on screen stale, so they
    ! go away until the user asks for the calculation again
    call iw_text("Isovalue",highlight=.true.,alignframe=.true.)
    if (iw_dragfloat_real8("(a.u.)##voidsisoval",x1=w%vd%iso_isoval,speed=1d-4,&
       min=1d-6,max=1d0,decimal=5,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)) &
       w%vd%iso_done = .false.
    call iw_tooltip("A point belongs to a void if the promolecular density (the sum of&
       & the in-vacuo atomic densities) at that point is lower than this value",ttshown)

    ! the grid the density is sampled on, chosen by its spacing
    call iw_text("Grid spacing",highlight=.true.,alignframe=.true.)
    if (iw_dragfloat_real8("(Å)##voidsspacing",x1=w%vd%iso_spacing,speed=0.005d0,&
       min=0.02d0,max=1d0,decimal=3,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)) &
       w%vd%iso_done = .false.
    call iw_tooltip("Approximate distance between two consecutive points of the grid on&
       & which the promolecular density is calculated. A finer grid resolves the shape of&
       & the voids better but takes longer",ttshown)

    ! the point count is counted in integer*8: three grid dimensions that the
    ! spacing slider can reach overflow a default integer, and a wrapped
    ! negative count would sail past the size check below
    call grid_from_spacing(sys(isys)%c%aa,w%vd%iso_spacing,n)
    npts = int(n(1),8) * int(n(2),8) * int(n(3),8)
    toobig = (npts > int(maxgridpts,8))
    call iw_text("(" // string(n(1)) // " × " // string(n(2)) // " × " // string(n(3)) //&
       " = " // string(npts) // " points)",alignframe=.true.,danger=toobig)

    ! grid point cost estimation
    if (.not.toobig) then
       ! re-measure when there is no measurement yet, or when the atom count
       ! has moved enough to matter
       if (w%vd%iso_secs <= 0d0 .or. &
          abs(sys(isys)%c%ncel - w%vd%iso_secs_ncel) > max(1,w%vd%iso_secs_ncel/4)) then
          xdum = 0d0
          rdum = iso_estimate_cost(isys,0,iso_region_cell,xdum,n)
          ! a failed measurement is not recorded, so it is retried; its
          ! failure paths return before the benchmark loop and cost nothing
          if (rdum > 0d0) then
             w%vd%iso_secs = rdum
             w%vd%iso_secs_ncel = sys(isys)%c%ncel
          end if
       end if
    end if

    ! run the calculation
    if (iw_button("Calculate##voidsisocalc",danger=.true.,disabled=toobig)) call run_isosurface()
    call iw_tooltip("Calculate the promolecular density on the grid and group the points&
       & below the isovalue into voids. The window does not respond while it runs",ttshown)
    if (toobig) then
       call iw_text("Too many grid points; increase the grid spacing",danger=.true.,&
          sameline=.true.)
    elseif (w%vd%iso_secs > 0d0) then
       tcost = real(npts,8) * w%vd%iso_secs
       s = "(~" // duration_string(tcost)
       if (tcost > bigwarn_secs) s = s // ", with the window frozen until it is done"
       call iw_text(s // ")",sameline=.true.,danger=(tcost > bigwarn_secs))
    elseif (npts > int(bigwarngridpts,8)) then
       ! the cost could not be measured: fall back on the point count, so a
       ! grid that will freeze the window for a long time still says so
       call iw_text("(a grid this size takes a while, with the window frozen&
          & until it is done)",sameline=.true.,danger=.true.)
    end if

    ! the results of the last run
    if (w%vd%iso_done) then
       call iw_text("Total void volume",highlight=.true.)
       call iw_text(string(w%vd%iso_vtot*fac3,'f',decimal=4) // " Å³ (" //&
          string(w%vd%iso_vtot/sys(isys)%c%omega*100d0,'f',decimal=2) // "% of the cell)",&
          sameline=.true.)
       call iw_text("Cell volume",highlight=.true.)
       call iw_text(string(sys(isys)%c%omega*fac3,'f',decimal=4) // " Å³",sameline=.true.)
       call iw_text("Number of voids",highlight=.true.)
       call iw_text(string(w%vd%iso_nvoid),sameline=.true.)

       ! one row per void, largest first
       flags = ImGuiTableFlags_None
       flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
       flags = ior(flags,ImGuiTableFlags_RowBg)
       flags = ior(flags,ImGuiTableFlags_Borders)
       flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
       flags = ior(flags,ImGuiTableFlags_ScrollY)
       flags = ior(flags,ImGuiTableFlags_ScrollX)
       str1 = "##tablevoidsiso" // c_null_char
       sz0%x = 0
       sz0%y = iw_calcheight(min(w%vd%iso_nvoid,10)+1,0,.false.)
       if (igBeginTable(c_loc(str1),5,flags,sz0,0._c_float)) then
          call iw_table_column("Id",id=ic_iso_id,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Volume (Å³)",id=ic_iso_vol,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("% cell",id=ic_iso_pct,flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("Deepest point (fractional)",id=ic_iso_x,&
             flags=ImGuiTableColumnFlags_WidthFixed)
          call iw_table_column("ρ (a.u.)",id=ic_iso_rho,flags=ImGuiTableColumnFlags_WidthFixed)
          call igTableSetupScrollFreeze(0,1)
          call igTableHeadersRow()
          call igTableSetColumnWidthAutoAll(igGetCurrentTable())

          ! the rows go through a clipper: there is one void per connected
          ! region and only the visible rows need to be emitted
          clipper = ImGuiListClipper_ImGuiListClipper()
          call ImGuiListClipper_Begin(clipper,w%vd%iso_nvoid,-1._c_float)
          do while (ImGuiListClipper_Step(clipper))
             call c_f_pointer(clipper,clipper_f)
             do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                if (igTableSetColumnIndex(ic_iso_id)) &
                   call iw_text(string(i))
                if (igTableSetColumnIndex(ic_iso_vol)) &
                   call iw_text(string(w%vd%iso_vol(i)*fac3,'f',length=12,decimal=5,justify=ioj_right))
                if (igTableSetColumnIndex(ic_iso_pct)) &
                   call iw_text(string(w%vd%iso_vol(i)/sys(isys)%c%omega*100d0,'f',length=8,&
                   decimal=3,justify=ioj_right))
                if (igTableSetColumnIndex(ic_iso_x)) then
                   s = string(w%vd%iso_x(1,i),'f',length=9,decimal=5,justify=ioj_right) //&
                      string(w%vd%iso_x(2,i),'f',length=9,decimal=5,justify=ioj_right) //&
                      string(w%vd%iso_x(3,i),'f',length=9,decimal=5,justify=ioj_right)
                   call iw_text(s)
                end if
                if (igTableSetColumnIndex(ic_iso_rho)) &
                   call iw_text(string(w%vd%iso_rho(i),'e',decimal=4))
             end do
          end do
          call ImGuiListClipper_End(clipper)
          call ImGuiListClipper_destroy(clipper)
          call igEndTable()
       end if
    end if

  contains
    ! Sample the promolecular density and collect the voids in it.
    subroutine run_isosurface()
      character(len=:), allocatable :: errmsg

      w%vd%iso_done = .false.
      call sys(isys)%c%promolecular_voids(n,w%vd%iso_isoval,w%vd%iso_vtot,w%vd%iso_nvoid,&
         w%vd%iso_vol,w%vd%iso_x,w%vd%iso_rho,errmsg)
      w%errmsg = errmsg
      w%vd%iso_done = (len_trim(errmsg) == 0)

    end subroutine run_isosurface

  end subroutine draw_isosurface_tab

  !> Draw the polyhedra tab: the volume of the coordination polyhedra, and
  !> the part of the cell they do not cover.
  subroutine draw_polyhedra_tab(w,isys,ttshown)
    use systems, only: sys
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_combo_simple,&
       iw_calcheight, iw_table_column
    use global, only: bondfactor
    use tools_io, only: string, ioj_right, ioj_center
    use param, only: atmcov
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    logical, intent(inout) :: ttshown

    logical :: changed
    integer :: i, j, iz1, iz2, nat, nf, ier, nspc
    integer(c_int) :: flags
    real*8 :: dmin, dmax, vol
    character(kind=c_char,len=:), allocatable, target :: str1, str2, s
    type(ImVec2) :: sz0

    integer, parameter :: ic_pol_id = 0
    integer, parameter :: ic_pol_at = 1
    integer, parameter :: ic_pol_mult = 2
    integer, parameter :: ic_pol_x = 3
    integer, parameter :: ic_pol_nv = 4
    integer, parameter :: ic_pol_d = 5
    integer, parameter :: ic_pol_nf = 6
    integer, parameter :: ic_pol_vol = 7

    nspc = sys(isys)%c%nspc

    ! the species selection, defaulted the first time this system is seen
    if (w%vd%pol_ic < 1 .or. w%vd%pol_ic > nspc .or. w%vd%pol_iv < 1 .or. w%vd%pol_iv > nspc) then
       w%vd%pol_ic = 1
       w%vd%pol_iv = min(2,nspc)
       w%vd%pol_rmax = -1d0
    end if
    iz1 = sys(isys)%c%spc(w%vd%pol_ic)%z
    iz2 = sys(isys)%c%spc(w%vd%pol_iv)%z

    ! the distance range follows the species until the user changes it
    if (w%vd%pol_rmax < 0d0) call reset_distance_range()

    ! the species at the center of the polyhedra and at their vertices
    str2 = ""
    do i = 1, nspc
       str2 = str2 // trim(sys(isys)%c%spc(i)%name) // c_null_char
    end do
    ! changing any of these makes the results on screen stale, so they go
    ! away until the user asks for the calculation again
    call iw_text("Center",highlight=.true.,alignframe=.true.)
    call iw_combo_simple("##voidspolcenter",str2,w%vd%pol_ic,sameline=.true.,&
       changed=changed,startsatone=.true.)
    if (changed) then
       w%vd%pol_rmax = -1d0
       w%vd%pol_done = .false.
    end if
    call iw_tooltip("Species of the atom at the center of the coordination polyhedra",ttshown)

    call iw_text("Vertices",highlight=.true.,alignframe=.true.,sameline=.true.)
    call iw_combo_simple("##voidspolvertex",str2,w%vd%pol_iv,sameline=.true.,&
       changed=changed,startsatone=.true.)
    if (changed) then
       w%vd%pol_rmax = -1d0
       w%vd%pol_done = .false.
    end if
    call iw_tooltip("Species of the atoms at the vertices of the coordination polyhedra",ttshown)

    ! the distance range that decides which atoms are vertices
    call iw_text("Distance range",highlight=.true.,alignframe=.true.)
    if (iw_dragfloat_real8("##voidspolrmin",x1=w%vd%pol_rmin,speed=0.01d0,min=0d0,&
       max=1d2,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)) &
       w%vd%pol_done = .false.
    call iw_tooltip("Shortest center-vertex distance (Å)",ttshown)
    if (iw_dragfloat_real8("(Å)##voidspolrmax",x1=w%vd%pol_rmax,speed=0.01d0,min=0d0,&
       max=1d2,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)) &
       w%vd%pol_done = .false.
    call iw_tooltip("Longest center-vertex distance (Å)",ttshown)
    if (iw_button("Reset##voidspolreset",sameline=.true.)) then
       call reset_distance_range()
       w%vd%pol_done = .false.
    end if
    call iw_tooltip("Restore the default distance range for these two species (from zero to&
       & the sum of covalent radii times the bond factor)",ttshown)

    ! run the calculation
    if (iw_button("Calculate##voidspolcalc",danger=.true.)) call run_polyhedra()
    call iw_tooltip("Build the coordination polyhedron of every non-equivalent atom of the&
       & center species and calculate its volume",ttshown)

    ! the results of the last run
    if (w%vd%pol_done) then
       call iw_text("Volume inside the polyhedra",highlight=.true.)
       call iw_text(string(w%vd%pol_vtot*fac3,'f',decimal=4) // " Å³ (" //&
          string(w%vd%pol_vtot/sys(isys)%c%omega*100d0,'f',decimal=2) // "% of the cell)",&
          sameline=.true.)
       call iw_text("Volume outside the polyhedra",highlight=.true.)
       call iw_text(string((sys(isys)%c%omega-w%vd%pol_vtot)*fac3,'f',decimal=4) // " Å³ (" //&
          string((1d0-w%vd%pol_vtot/sys(isys)%c%omega)*100d0,'f',decimal=2) // "% of the cell)",&
          sameline=.true.)
       call iw_text("Cell volume",highlight=.true.)
       call iw_text(string(sys(isys)%c%omega*fac3,'f',decimal=4) // " Å³",sameline=.true.)
       if (w%vd%pol_n == 0) then
          call iw_text("No coordination polyhedra were found in this distance range",&
             danger=.true.,wrap=.true.)
       elseif (w%vd%pol_vtot > sys(isys)%c%omega) then
          ! the volumes are added up one polyhedron at a time, which only
          ! measures the empty space if the polyhedra do not overlap. They
          ! add up to more than the cell when they do, and the difference
          ! below is then not a volume of anything
          call iw_text("The polyhedra add up to more than the cell volume: they overlap&
             & each other, so the volume outside them is not meaningful. Use a shorter&
             & maximum distance",danger=.true.,wrap=.true.)
       end if

       ! one row per polyhedron
       if (w%vd%pol_n > 0) then
          flags = ImGuiTableFlags_None
          flags = ior(flags,ImGuiTableFlags_NoSavedSettings)
          flags = ior(flags,ImGuiTableFlags_RowBg)
          flags = ior(flags,ImGuiTableFlags_Borders)
          flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
          flags = ior(flags,ImGuiTableFlags_ScrollY)
          flags = ior(flags,ImGuiTableFlags_ScrollX)
          str1 = "##tablevoidspol" // c_null_char
          sz0%x = 0
          sz0%y = iw_calcheight(min(w%vd%pol_n,10)+1,0,.false.)
          if (igBeginTable(c_loc(str1),8,flags,sz0,0._c_float)) then
             call iw_table_column("Id",id=ic_pol_id,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Atom",id=ic_pol_at,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Mult",id=ic_pol_mult,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Coordinates (fractional)",id=ic_pol_x,&
                flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("nv",id=ic_pol_nv,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Distances (Å)",id=ic_pol_d,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("nf",id=ic_pol_nf,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Volume (Å³)",id=ic_pol_vol,flags=ImGuiTableColumnFlags_WidthFixed)
             call igTableSetupScrollFreeze(0,1)
             call igTableHeadersRow()
             call igTableSetColumnWidthAutoAll(igGetCurrentTable())

             do i = 1, w%vd%pol_n
                j = w%vd%pol_id(i)
                call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                if (igTableSetColumnIndex(ic_pol_id)) &
                   call iw_text(string(j))
                if (igTableSetColumnIndex(ic_pol_at)) &
                   call iw_text(string(sys(isys)%c%at(j)%name,4,ioj_center))
                if (igTableSetColumnIndex(ic_pol_mult)) &
                   call iw_text(string(sys(isys)%c%at(j)%mult))
                if (igTableSetColumnIndex(ic_pol_x)) then
                   s = string(sys(isys)%c%at(j)%x(1),'f',length=9,decimal=5,justify=ioj_right) //&
                      string(sys(isys)%c%at(j)%x(2),'f',length=9,decimal=5,justify=ioj_right) //&
                      string(sys(isys)%c%at(j)%x(3),'f',length=9,decimal=5,justify=ioj_right)
                   call iw_text(s)
                end if
                if (igTableSetColumnIndex(ic_pol_nv)) &
                   call iw_text(string(w%vd%pol_nv(i)))
                if (igTableSetColumnIndex(ic_pol_d)) then
                   s = string(w%vd%pol_dmin(i)*bohrtoa,'f',length=8,decimal=4,justify=ioj_right) //&
                      " to " //&
                      string(w%vd%pol_dmax(i)*bohrtoa,'f',length=8,decimal=4,justify=ioj_right)
                   call iw_text(s)
                end if
                if (igTableSetColumnIndex(ic_pol_nf)) &
                   call iw_text(string(w%vd%pol_nf(i)))
                if (igTableSetColumnIndex(ic_pol_vol)) &
                   call iw_text(string(w%vd%pol_vol(i)*fac3,'f',length=12,decimal=5,justify=ioj_right))
             end do
             call igEndTable()
          end if
       end if
    end if

  contains
    ! The default distance range for the two selected species: from zero to
    ! the sum of covalent radii times the bond factor, as in POLYHEDRA.
    subroutine reset_distance_range()

      w%vd%pol_rmin = 0d0
      w%vd%pol_rmax = (atmcov(iz1) + atmcov(iz2)) * bondfactor * bohrtoa

    end subroutine reset_distance_range

    ! Build the coordination polyhedron of every non-equivalent atom of the
    ! center species.
    subroutine run_polyhedra()
      integer :: k

      w%errmsg = ""
      w%vd%pol_done = .false.
      w%vd%pol_n = 0
      w%vd%pol_vtot = 0d0
      if (allocated(w%vd%pol_id)) deallocate(w%vd%pol_id)
      if (allocated(w%vd%pol_nv)) deallocate(w%vd%pol_nv)
      if (allocated(w%vd%pol_nf)) deallocate(w%vd%pol_nf)
      if (allocated(w%vd%pol_dmin)) deallocate(w%vd%pol_dmin)
      if (allocated(w%vd%pol_dmax)) deallocate(w%vd%pol_dmax)
      if (allocated(w%vd%pol_vol)) deallocate(w%vd%pol_vol)
      allocate(w%vd%pol_id(sys(isys)%c%nneq),w%vd%pol_nv(sys(isys)%c%nneq))
      allocate(w%vd%pol_nf(sys(isys)%c%nneq),w%vd%pol_dmin(sys(isys)%c%nneq))
      allocate(w%vd%pol_dmax(sys(isys)%c%nneq),w%vd%pol_vol(sys(isys)%c%nneq))

      do k = 1, sys(isys)%c%nneq
         if (sys(isys)%c%at(k)%is /= w%vd%pol_ic) cycle
         call sys(isys)%c%coord_polyhedron(sys(isys)%c%at(k)%x,w%vd%pol_iv,0,&
            w%vd%pol_rmin/bohrtoa,w%vd%pol_rmax/bohrtoa,nat,dmin,dmax,nf,vol,ier)
         if (nat <= 2) cycle ! fewer than three vertices in range: no polyhedron
         if (ier /= 0) then
            w%errmsg = "Failed to triangulate the coordination polyhedron of atom " // string(k)
            cycle
         end if
         w%vd%pol_n = w%vd%pol_n + 1
         w%vd%pol_id(w%vd%pol_n) = k
         w%vd%pol_nv(w%vd%pol_n) = nat
         w%vd%pol_nf(w%vd%pol_n) = nf
         w%vd%pol_dmin(w%vd%pol_n) = dmin
         w%vd%pol_dmax(w%vd%pol_n) = dmax
         w%vd%pol_vol(w%vd%pol_n) = vol
         w%vd%pol_vtot = w%vd%pol_vtot + sys(isys)%c%at(k)%mult * vol
      end do
      w%vd%pol_done = .true.

    end subroutine run_polyhedra

  end subroutine draw_polyhedra_tab

  !> Draw the packing tab: the atoms are spheres and the void is whatever
  !> falls outside all of them.
  subroutine draw_packing_tab(w,isys,ttshown)
    use systems, only: sys
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_combo_simple
    use tools_io, only: string
    use param, only: atmvdw, atmcov
    type(window), intent(inout), target :: w
    integer, intent(in) :: isys
    logical, intent(inout) :: ttshown

    logical :: changed, ismc
    real*8 :: vvoid

    ! which spheres the atoms are. The nearest-neighbor spheres never overlap,
    ! so their volume is a sum and there is nothing to sample
    call iw_text("Atomic radii",highlight=.true.,alignframe=.true.)
    call iw_combo_simple("##voidspckradii","van der Waals" // c_null_char // "covalent" //&
       c_null_char // "half the nearest-neighbor distance" // c_null_char,w%vd%pck_radii,&
       sameline=.true.,changed=changed)
    if (changed) w%vd%pck_done = .false.
    call iw_tooltip("Radius assigned to each atom. The van der Waals and covalent radii&
       & come from critic2's internal tables and the spheres may overlap; half the distance&
       & to the nearest neighbor gives spheres that never do",ttshown)
    ismc = (w%vd%pck_radii /= vdrad_nnm)

    ! the Monte Carlo sampling stops when it reaches this relative error
    call igBeginDisabled(logical(.not.ismc,c_bool))
    call iw_text("Precision",highlight=.true.,alignframe=.true.)
    if (iw_dragfloat_real8("##voidspckprec",x1=w%vd%pck_prec,speed=1d-3,min=1d-4,&
       max=1d-1,decimal=4,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)) &
       w%vd%pck_done = .false.
    call iw_tooltip("The volume covered by the spheres is calculated by Monte Carlo sampling,&
       & which stops when the standard deviation of the volume divided by the volume itself&
       & reaches this value. A smaller value takes longer",ttshown)
    call igEndDisabled()

    ! run the calculation
    if (iw_button("Calculate##voidspckcalc",danger=.true.)) call run_packing()
    call iw_tooltip("Calculate the volume covered by the atomic spheres and the empty space&
       & left outside them",ttshown)

    ! the results of the last run
    if (w%vd%pck_done) then
       vvoid = sys(isys)%c%omega - w%vd%pck_vfill

       call iw_text("Volume inside the spheres",highlight=.true.)
       call iw_text(string(w%vd%pck_vfill*fac3,'f',decimal=4) // pmstring(w%vd%pck_err*fac3) //&
          " Å³",sameline=.true.)
       call iw_text("Void (interstitial) volume",highlight=.true.)
       call iw_text(string(vvoid*fac3,'f',decimal=4) // pmstring(w%vd%pck_err*fac3) // " Å³",&
          sameline=.true.)
       call iw_text("Cell volume",highlight=.true.)
       call iw_text(string(sys(isys)%c%omega*fac3,'f',decimal=4) // " Å³",sameline=.true.)
       call iw_text("Packing ratio",highlight=.true.)
       call iw_text(string(w%vd%pck_vfill/sys(isys)%c%omega*100d0,'f',decimal=4) //&
          pmstring(w%vd%pck_err/sys(isys)%c%omega*100d0) // "%",sameline=.true.)
    end if

  contains
    ! The "+- error" part of a result, dropped when the volume is exact.
    function pmstring(err) result(s)
      real*8, intent(in) :: err
      character(len=:), allocatable :: s

      s = ""
      if (err > 0d0) s = " ± " // string(err,'f',decimal=4)

    end function pmstring

    ! Calculate the volume covered by the atomic spheres.
    subroutine run_packing()

      w%errmsg = ""
      select case (w%vd%pck_radii)
      case (vdrad_nnm)
         ! the spheres do not overlap: the volume is the sum of their volumes
         w%vd%pck_vfill = sys(isys)%c%get_pack_ratio() / 100d0 * sys(isys)%c%omega
         w%vd%pck_err = 0d0
      case (vdrad_cov)
         w%vd%pck_vfill = sys(isys)%c%vdw_volume(w%vd%pck_prec,atmcov)
         w%vd%pck_err = w%vd%pck_prec * w%vd%pck_vfill
      case default ! vdrad_vdw
         w%vd%pck_vfill = sys(isys)%c%vdw_volume(w%vd%pck_prec,atmvdw)
         w%vd%pck_err = w%vd%pck_prec * w%vd%pck_vfill
      end select
      w%vd%pck_done = .true.

    end subroutine run_packing

  end subroutine draw_packing_tab

  !> Number of grid points along each lattice vector that comes closest to
  !> the target spacing (in Å) for a cell with lengths aa (in bohr). At
  !> least two points are used along each direction.
  subroutine grid_from_spacing(aa,spacing,n)
    real*8, intent(in) :: aa(3)
    real*8, intent(in) :: spacing
    integer, intent(out) :: n(3)

    integer :: i

    do i = 1, 3
       n(i) = max(nint(aa(i) * bohrtoa / spacing),2)
    end do

  end subroutine grid_from_spacing

end submodule voids
