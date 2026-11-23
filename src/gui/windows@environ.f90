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

! Routines for the atomic environments window (the ENVIRON, ECON and
! COORD keywords, for one atom).
submodule (windows) environ
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

  ! transient markers in the view: halos around the centre, the
  ! coordinated neighbours, and the neighbours of the hovered row
  ! (color, opacity), and the radius of a halo relative to the
  ! covalent radius of the atom (the atoms are drawn at 0.7)
  real(c_float), parameter :: center_rgb(3) = (/1.0_c_float,0.6_c_float,0.1_c_float/)
  real(c_float), parameter :: coord_rgb(3) = (/0.2_c_float,0.7_c_float,1.0_c_float/)
  real(c_float), parameter :: hover_rgb(3) = (/1.0_c_float,1.0_c_float,0.2_c_float/)
  real(c_float), parameter :: halo_alpha = 0.35_c_float
  real(c_float), parameter :: hover_alpha = 0.7_c_float
  real*8, parameter :: halo_scale = 0.85d0
  real*8, parameter :: halo_radmin = 0.5d0

  ! the columns of the neighbour table
  integer(c_int), parameter :: ic_shell = 0
  integer(c_int), parameter :: ic_name = 1
  integer(c_int), parameter :: ic_id = 2
  integer(c_int), parameter :: ic_cidx = 3
  integer(c_int), parameter :: ic_lvec = 4
  integer(c_int), parameter :: ic_count = 5
  integer(c_int), parameter :: ic_dist = 6
  integer(c_int), parameter :: ic_coord = 7
  integer(c_int), parameter :: ic_wecon = 8

contains

  !> Draw the atomic environments window: the neighbours of an atom of
  !> the system shown by the anchor view (chosen in a list or picked in
  !> the view), as ENVIRON lists them, with the coordination (the bonds
  !> of the atom) and the ECoN weight of each.
  module subroutine draw_environ(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_close_event,&
       iw_setpos_bottomright, iw_table_column, iw_table_headers_row, iw_calcheight,&
       iw_separatortext, iw_highlight_selectable, iw_cell_right, iw_combo_simple,&
       iw_intstepper, iw_checkbox
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer, parameter :: itable_neigh = 1 ! the table, for table_write_*

    logical :: doquit, goodsys, syschanged, ismol, ldum, changed
    integer :: i, isys, iview, istat, iwrite, ineq, nrow, ncol, ihover
    integer :: atpick(4)
    integer(c_int) :: icol(0:8) ! position of each column in the table (-1 = not shown)
    integer(c_int) :: tflags
    real*8 :: xc(3)
    type(ImVec2) :: sz0
    type(c_ptr) :: clipper
    type(ImGuiListClipper), pointer :: clipper_f
    character(len=:), allocatable :: s
    character(len=:,kind=c_char), allocatable, target :: str1

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window operates on the system shown by its anchor view
    doquit = .not.w%anchor(iview,isys,syschanged)
    goodsys = .not.doquit
    if (goodsys) goodsys = ok_system(isys,sys_init)
    ismol = .false.
    if (goodsys) ismol = sys(isys)%c%ismolecule

    ! initialize the form on the first pass
    if (w%firstpass) then
       w%env = environ_state()
       w%errmsg = ""
    end if

    ! a new system: the centre goes to the first atom
    if (goodsys) then
       if (.not.w%env%sys%is(isys) .or. syschanged) then
          call w%env%sys%set(isys)
          call cancel_pick()
          w%errmsg = ""
          w%env%idx = (/1,0,0,0/)
          w%env%ispc = 0
       end if
    end if

    ! a pending pick in the view it was armed on; dropped, releasing the
    ! view, if the system is not ready or the window moved to another view
    if (w%env%picking) then
       if (goodsys .and. iview == w%env%pickview) then
          call view_pick_result(iview,w%id,isys,w%env%pick,istat,xc,atidx=atpick)
          if (istat == ipick_point .and. atpick(1) > 0) w%env%idx = atpick
          if (istat /= ipick_pending) w%env%picking = .false.
       else
          call cancel_pick()
       end if
    end if

    ! a centre that is gone after a change of geometry, or a species
    ! that is gone
    if (goodsys) then
       if (w%env%idx(1) < 1 .or. w%env%idx(1) > sys(isys)%c%ncel) w%env%idx = (/1,0,0,0/)
       if (ismol) w%env%idx(2:4) = 0
       if (w%env%ispc < 0 .or. w%env%ispc > sys(isys)%c%nspc) w%env%ispc = 0
    end if

    ihover = 0
    if (goodsys) then
       ! the neighbour list, calculated again if anything it depends on changed
       call update_list()

       ! the system
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)

       ! the centre
       call iw_separatortext("Center",w%env%seccenter)
       if (w%env%seccenter%open) then
          associate(c => sys(isys)%c, ls => w%env%ls)
            ! the entries of the combo: the non-equivalent atoms
            if (ls%atomkeysys /= w%env%sys%uid .or. ls%atomkeytime /= sysc(isys)%timelastchange_geometry) then
               ls%atomkeysys = w%env%sys%uid
               ls%atomkeytime = sysc(isys)%timelastchange_geometry
               ls%atomstr = ""
               do i = 1, c%nneq
                  ls%atomstr = ls%atomstr // trim(c%at(i)%name) // " (" // string(i) // ")" // c_null_char
               end do
            end if
            call iw_text("Atom",alignframe=.true.)
            ineq = c%atcel(w%env%idx(1))%idx - 1
            call iw_combo_simple("##environatom",ls%atomstr,ineq,sameline=.true.,changed=changed)
            if (changed) then
               ! the first cell atom of that non-equivalent atom, in the main cell
               do i = 1, c%ncel
                  if (c%atcel(i)%idx == ineq + 1) then
                     w%env%idx = (/i,0,0,0/)
                     exit
                  end if
               end do
            end if
            if (ismol) then
               call iw_tooltip("Atom whose environment is shown",ttshown)
            else
               call iw_tooltip("Atom whose environment is shown (a non-equivalent atom; picking&
                  & it in the view chooses a particular image)",ttshown)
            end if
            if (iw_button("Pick##environpick",sameline=.true.,disabled=w%env%picking)) then
               w%env%picking = .true.
               w%env%pickview = iview
               call w%env%pick%arm()
               call win(iview)%viewmode_set_forced(vm_pick_atom,"Pick the central atom",w%id)
            end if
            call iw_tooltip("Pick the atom in the view",ttshown)

            if (.not.ismol) then
               call iw_text("Cell atom")
               call iw_text(string(w%env%idx(1)),sameline=.true.)
               call iw_tooltip("Index of the atom in the complete list of the unit cell",ttshown)
               call iw_text("Lattice vector",gap=.true.)
               call iw_text(lvec_string(w%env%idx(2:4)),sameline=.true.)
               call iw_tooltip("Lattice vector of the image of the atom marked in the view",ttshown)
            end if
          end associate
       end if

       ! the neighbours
       call iw_separatortext("Neighbors (" // string(w%env%ls%nat) // ")",w%env%secneigh)
       if (w%env%secneigh%open) then
          call iw_text("Up to",alignframe=.true.)
          call iw_combo_simple("##environupmode","Distance (Å)" // c_null_char // "Shells" //&
             c_null_char,w%env%upmode,sameline=.true.)
          call iw_tooltip("List the neighbors up to a distance from the atom, or up to a number&
             & of shells (atoms at the same distance)",ttshown)
          if (w%env%upmode == 0) then
             ldum = iw_dragfloat_real8("##environupdist",x1=w%env%updist,speed=0.02d0,min=0.1d0,&
                max=50d0,decimal=2,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
             call iw_tooltip("Largest distance from the atom to the neighbors in the list",ttshown)
          else
             ldum = iw_intstepper("##environupshell",w%env%upshell,minval=1_c_int,maxval=100_c_int,&
                sameline=.true.)
             call iw_tooltip("Number of shells of neighbors in the list",ttshown)
          end if

          call iw_text("Species",alignframe=.true.,gap=.true.)
          s = "All" // c_null_char
          do i = 1, sys(isys)%c%nspc
             s = s // trim(sys(isys)%c%spc(i)%name) // c_null_char
          end do
          call iw_combo_simple("##environspc",s,w%env%ispc,sameline=.true.)
          call iw_tooltip("List only the neighbors of this species",ttshown)

          ldum = iw_checkbox("Group by shells",w%env%groupshells,labelleft=.true.)
          call iw_tooltip("One row per shell and atom (a shell may hold several kinds of atom),&
             & with the number of atoms",ttshown)
          ldum = iw_checkbox("Show in view",w%env%showview,labelleft=.true.,gap=.true.)
          call iw_tooltip("Mark the atom (orange) and its coordinated neighbors (blue) in the&
             & view; the neighbors of the row under the mouse are marked in yellow",ttshown)

          ! the table
          call draw_table()
       end if

       ! the centre and its neighbours in the anchor view
       call draw_markers()

       ! the error message from writing the table
       if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)
    end if

    ! bottom row: close
    call iw_setpos_bottomright(len("Close"),1)
    if (iw_button("Close")) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) then
       call cancel_pick()
       call w%end()
    end if

  contains

    !> Calculate the neighbour list of the centre if anything it
    !> depends on changed: the system, the atom, the cutoff, the
    !> species, the geometry, or the bonds.
    subroutine update_list()
      use param, only: icrd_crys
      real*8, allocatable :: econ(:), econ1(:), nd(:), nd1(:)
      real*8, allocatable :: up2d
      integer, allocatable :: up2sh, ispc0
      integer, allocatable :: iperm(:)
      integer :: i, j, k, iat
      real*8 :: time

      iat = w%env%idx(1)
      associate(c => sys(isys)%c, ls => w%env%ls)
        ! the mean weighted distance of the ECoN of the atom (all
        ! species), for the weights; it does not depend on the cutoff
        time = sysc(isys)%timelastchange_geometry
        if (ls%ndkeysys /= w%env%sys%uid .or. ls%ndkeyiat /= iat .or. ls%ndkeytime /= time) then
           ls%ndkeysys = w%env%sys%uid
           ls%ndkeyiat = iat
           ls%ndkeytime = time
           allocate(econ(0:c%nspc),econ1(0:c%nspc),nd(0:c%nspc),nd1(0:c%nspc))
           call c%econ_atom(c%atcel(iat)%x,econ,econ1,nd,nd1)
           ls%nd0 = nd(0)
           ls%keysys = 0
        end if

        time = max(sysc(isys)%timelastchange_geometry,sysc(isys)%timelastchange_rebond)
        if (ls%keysys == w%env%sys%uid .and. ls%keyiat == iat .and. ls%keymode == w%env%upmode .and.&
           ls%keyspc == w%env%ispc .and. ls%keytime == time) then
           if (w%env%upmode == 0 .and. ls%keydist == w%env%updist) return
           if (w%env%upmode == 1 .and. ls%keyshell == w%env%upshell) return
        end if
        ls%keysys = w%env%sys%uid
        ls%keyiat = iat
        ls%keymode = w%env%upmode
        ls%keyspc = w%env%ispc
        ls%keytime = time
        ls%keydist = w%env%updist
        ls%keyshell = w%env%upshell

        ! the list, around the atom in the main cell (an unallocated
        ! criterion is an absent argument)
        if (w%env%upmode == 0) then
           up2d = w%env%updist / bohrtoa
        else
           up2sh = w%env%upshell
        end if
        if (w%env%ispc > 0) ispc0 = w%env%ispc
        call c%list_near_atoms(c%atcel(iat)%x,icrd_crys,.true.,ls%nat,ls%eid,ls%dist,ls%lvec,ls%ishell,&
           up2d=up2d,up2sh=up2sh,ispc0=ispc0,nozero=.true.)

        ! within a shell (the same distance), by non-equivalent atom
        allocate(iperm(ls%nat))
        do i = 1, ls%nat
           iperm(i) = i
        end do
        do i = 2, ls%nat
           k = iperm(i)
           j = i - 1
           do while (j >= 1)
              if (ls%ishell(iperm(j)) /= ls%ishell(k) .or.&
                 c%atcel(ls%eid(iperm(j)))%idx <= c%atcel(ls%eid(k))%idx) exit
              iperm(j+1) = iperm(j)
              j = j - 1
           end do
           iperm(j+1) = k
        end do
        ls%eid = ls%eid(iperm)
        ls%dist = ls%dist(iperm)
        ls%lvec = ls%lvec(:,iperm)
        ls%ishell = ls%ishell(iperm)

        ! the groups: runs in the same shell and non-equivalent atom
        if (allocated(ls%igroup)) deallocate(ls%igroup)
        allocate(ls%igroup(ls%nat+1))
        ls%ngroup = 0
        do i = 1, ls%nat
           if (i > 1) then
              if (ls%ishell(i) == ls%ishell(i-1) .and.&
                 c%atcel(ls%eid(i))%idx == c%atcel(ls%eid(i-1))%idx) cycle
           end if
           ls%ngroup = ls%ngroup + 1
           ls%igroup(ls%ngroup) = i
        end do
        ls%igroup(ls%ngroup+1) = ls%nat + 1

        ! the ECoN weight of each neighbour
        if (allocated(ls%wecon)) deallocate(ls%wecon)
        allocate(ls%wecon(ls%nat))
        ls%wecon = 0d0
        if (ls%nd0 > 0d0) ls%wecon = exp(1d0 - (ls%dist/ls%nd0)**6)

        ! coordinated: bonded to the atom
        if (allocated(ls%iscoord)) deallocate(ls%iscoord)
        allocate(ls%iscoord(ls%nat))
        ls%iscoord = .false.
        if (allocated(c%nstar)) then
           do i = 1, ls%nat
              do j = 1, c%nstar(iat)%ncon
                 if (c%nstar(iat)%idcon(j) == ls%eid(i) .and.&
                    all(c%nstar(iat)%lcon(:,j) == ls%lvec(:,i))) then
                    ls%iscoord(i) = .true.
                    exit
                 end if
              end do
           end do
        end if
      end associate

    end subroutine update_list

    !> The table of neighbours (the section is open): one row per
    !> neighbour, or one per group (shell and atom) if grouped. Sets
    !> ihover to the hovered row.
    subroutine draw_table()
      integer :: i, i1, i2
      logical :: grouped

      grouped = w%env%groupshells
      tflags = ImGuiTableFlags_None
      tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
      tflags = ior(tflags,ImGuiTableFlags_RowBg)
      tflags = ior(tflags,ImGuiTableFlags_Borders)
      tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
      tflags = ior(tflags,ImGuiTableFlags_ScrollY)
      ! a different table in each mode, so each keeps its own column widths
      if (grouped) then
         str1 = "##environtableshells" // c_null_char
         nrow = w%env%ls%ngroup
      else
         str1 = "##environtable" // c_null_char
         nrow = w%env%ls%nat
      end if

      ! the columns of this mode
      ncol = 0
      icol = -1
      call addcol(ic_shell)
      call addcol(ic_name)
      call addcol(ic_id)
      if (grouped) then
         call addcol(ic_count)
      else
         ! in a molecule, the cell id is the id
         if (.not.ismol) then
            call addcol(ic_cidx)
            call addcol(ic_lvec)
         end if
      end if
      call addcol(ic_dist)
      call addcol(ic_coord)
      call addcol(ic_wecon)

      ! the rest of the window, but at least a few rows (the stretch
      ! beyond them is elastic)
      call igGetContentRegionAvail(sz0)
      sz0%x = 0._c_float
      sz0%y = sz0%y - iw_calcheight(1,0,.true.)
      sz0%y = max(sz0%y,iw_calcheight(6,0,.false.))
      w%heightslack = sz0%y - iw_calcheight(6,0,.false.)
      if (.not.igBeginTable(c_loc(str1),ncol,tflags,sz0,0._c_float)) return

      ! the headers, in the order of the columns
      call column(ic_shell,"Shell")
      call column(ic_name,"Atom")
      call column(ic_id,"Id")
      call column(ic_cidx,"Cell id")
      call column(ic_lvec,"Lattice vector")
      call column(ic_count,"Count")
      call column(ic_dist,"Distance (Å)")
      call column(ic_coord,"Coordinated")
      call column(ic_wecon,"ECoN weight")
      call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)

      clipper = ImGuiListClipper_ImGuiListClipper()
      call ImGuiListClipper_Begin(clipper,nrow,igGetTextLineHeightWithSpacing())
      call w%table_write_begin(itable_neigh,iwrite,clipper)
      do while (ImGuiListClipper_Step(clipper))
         call c_f_pointer(clipper,clipper_f)
         do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
            ! the neighbours in this row: one, or the atoms of a group
            call row_range(i,i1,i2)

            call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
            associate(c => sys(isys)%c, ls => w%env%ls, j => w%env%ls%eid(i1))
              if (setcol(ic_shell)) &
                 call iw_cell_right(string(ls%ishell(i1)))
              ! selectable spans the row; hovering it marks its atoms in the view
              if (iw_highlight_selectable("##environrow" // string(i))) ihover = i
              if (setcol(ic_name)) &
                 call iw_text(trim(c%at(c%atcel(j)%idx)%name))
              if (setcol(ic_id)) &
                 call iw_cell_right(string(c%atcel(j)%idx))
              if (setcol(ic_count)) &
                 call iw_cell_right(string(i2-i1+1))
              if (setcol(ic_cidx)) &
                 call iw_cell_right(string(j))
              if (setcol(ic_lvec)) &
                 call iw_text(lvec_string(ls%lvec(:,i1)+w%env%idx(2:4)))
              if (setcol(ic_dist)) &
                 call iw_cell_right(string(ls%dist(i1)*bohrtoa,'f',decimal=5))
              if (setcol(ic_coord)) then
                 if (grouped) then
                    call iw_cell_right(string(count(ls%iscoord(i1:i2))))
                 elseif (ls%iscoord(i1)) then
                    call iw_text("yes")
                 end if
              end if
              if (setcol(ic_wecon)) &
                 call iw_cell_right(string(ls%wecon(i1),'f',decimal=4))
            end associate
         end do
      end do
      call ImGuiListClipper_End(clipper)
      call ImGuiListClipper_destroy(clipper)
      call w%table_write_end(itable_neigh,table_title("Environment of atom " //&
         string(w%env%idx(1)),isys))
      call igEndTable()

    end subroutine draw_table

    !> Add column id to the columns of the table (at the end).
    subroutine addcol(id)
      integer(c_int), intent(in) :: id
      icol(id) = int(ncol,c_int)
      ncol = ncol + 1
    end subroutine addcol

    !> Move to column id of the current row; false if the table does
    !> not have it, or if it is clipped.
    function setcol(id) result(ok)
      integer(c_int), intent(in) :: id
      logical :: ok
      ok = .false.
      if (icol(id) < 0) return
      ok = igTableSetColumnIndex(icol(id))
    end function setcol

    !> The header of column id, if the table has it.
    subroutine column(id,title)
      integer(c_int), intent(in) :: id
      character(len=*), intent(in) :: title
      if (icol(id) < 0) return
      call iw_table_column(title,id=id,flags=ImGuiTableColumnFlags_WidthFixed)
    end subroutine column

    !> The neighbours i1 to i2 in row irow of the table: one, or the
    !> atoms of a group if grouped.
    subroutine row_range(irow,i1,i2)
      integer, intent(in) :: irow
      integer, intent(out) :: i1, i2

      if (w%env%groupshells) then
         i1 = w%env%ls%igroup(irow)
         i2 = w%env%ls%igroup(irow+1) - 1
      else
         i1 = irow
         i2 = irow
      end if

    end subroutine row_range

    !> Mark the centre, its coordinated neighbours, and the neighbours
    !> of the hovered row in the anchor view, with translucent halos.
    subroutine draw_markers()
      use representations, only: rep_shape, shapekind_sphere, shapekind_cylinder
      type(rep_shape), allocatable :: shp(:)
      integer :: i, n, i1, i2
      logical :: inhover
      real*8 :: x0(3)

      if (.not.associated(win(iview)%sc)) return
      if (win(iview)%sc%isinit == 0) return
      ! not armed this frame = not shown
      if (.not.w%env%showview) return

      ! the atoms of the hovered row
      i1 = 1
      i2 = 0
      if (ihover > 0) call row_range(ihover,i1,i2)

      allocate(shp(1 + 2 * w%env%ls%nat))
      n = 1
      x0 = atpos(w%env%idx(1),(/0,0,0/))
      shp(1) = rep_shape(kind=shapekind_sphere,x1=x0,rad=halo_rad(w%env%idx(1)),&
         rgb=center_rgb,alpha=halo_alpha,rim=.true.)
      do i = 1, w%env%ls%nat
         inhover = (i >= i1 .and. i <= i2)
         if (.not.inhover .and. .not.w%env%ls%iscoord(i)) cycle
         n = n + 1
         if (inhover) then
            shp(n) = rep_shape(kind=shapekind_sphere,x1=atpos(w%env%ls%eid(i),w%env%ls%lvec(:,i)),&
               rad=halo_rad(w%env%ls%eid(i)),rgb=hover_rgb,alpha=hover_alpha,rim=.true.)
            ! and a line from the centre (it shows past the bonds)
            n = n + 1
            shp(n) = rep_shape(kind=shapekind_cylinder,x1=x0,&
               v=reshape((/atpos(w%env%ls%eid(i),w%env%ls%lvec(:,i))-x0,0d0,0d0,0d0,0d0,0d0,0d0/),(/3,3/)),&
               rad=0.05d0,rgb=hover_rgb)
         else
            shp(n) = rep_shape(kind=shapekind_sphere,x1=atpos(w%env%ls%eid(i),w%env%ls%lvec(:,i)),&
               rad=halo_rad(w%env%ls%eid(i)),rgb=coord_rgb,alpha=halo_alpha,rim=.true.)
         end if
      end do
      call win(iview)%sc%show_transient_shapes(w%uid,1,shp(1:n))

    end subroutine draw_markers

    !> Position in the scene (Cartesian, bohr; molecules: in the
    !> frame of the shapes) of cell atom iat with lattice vector lv
    !> from the centre's cell, in the image of the centre.
    function atpos(iat,lv) result(x)
      integer, intent(in) :: iat, lv(3)
      real*8 :: x(3)
      associate(c => sys(isys)%c)
        x = c%x2c(c%atcel(iat)%x + lv + w%env%idx(2:4))
        if (ismol) x = x + c%molx0
      end associate
    end function atpos

    !> Radius of the halo around cell atom iat.
    function halo_rad(iat) result(r)
      use param, only: atmcov0
      integer, intent(in) :: iat
      real*8 :: r
      associate(c => sys(isys)%c)
        r = max(halo_scale * atmcov0(c%spc(c%atcel(iat)%is)%z),halo_radmin)
      end associate
    end function halo_rad

    !> A lattice vector as text.
    function lvec_string(lv) result(str)
      integer, intent(in) :: lv(3)
      character(len=:), allocatable :: str
      str = string(lv(1)) // " " // string(lv(2)) // " " // string(lv(3))
    end function lvec_string

    !> Forget a pending pick and release the view it was armed on.
    subroutine cancel_pick()
      if (.not.w%env%picking) return
      if (w%env%pickview >= 1 .and. w%env%pickview <= nwin) &
         call win(w%env%pickview)%viewmode_release_forced(w%id)
      w%env%picking = .false.
    end subroutine cancel_pick

  end subroutine draw_environ

end submodule environ
