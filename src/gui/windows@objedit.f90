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

! Object editing in the view (vm_objedit): the handlers of the object
! types edited with the mouse in the view (2D drawings, 3D shapes), and
! their helpers. The 2D handler works in the NDC of the render buffer.
! The 3D handlers work on the screen projection of the objects:
! a point under the mouse is the atom under it (unless the no-snap bind
! is held) or the point at a given depth, and the objects and their
! handles are hit in texture pixels.
submodule (windows) objedit
  use interfaces_cimgui
  implicit none

  ! half the side of the drawn handles (pixels)
  real(c_float), parameter :: handle_px = 4.5_c_float

  ! the edges of a box (corners as in box_corners)
  integer, parameter :: box_edges(2,12) = reshape((/1,2, 1,3, 1,5, 2,4, 2,6, 3,4, 3,7, 4,8,&
     5,6, 5,7, 6,8, 7,8/),shape(box_edges))

contains

  !> The object editing view mode (vm_objedit) for the 2D drawing r,
  !> with the input of this frame inp: draw new shapes with the tool,
  !> select and edit (move, drag the handles of) or remove the existing
  !> ones, finish a polygon or polyline in progress (the exit bind, which
  !> is then marked as used, the finish key, a double click). The shape
  !> being drawn or dragged is shown over the view (ImGui draw list) and
  !> goes to the drawing when the operation ends.
  module subroutine planar_events(w,r,inp)
    use representations, only: representation,&
       planarkind_NUM, planarkind_ellipse, planarkind_rect,&
       planarkind_arrow, planarkind_curve, planarkind_freehand, planar_template, planar_append,&
       planar_curve_default_bend,&
       planar_delete, planar_handles, planar_drag_handle, planar_move,&
       planar_simplify
    use keybindings, only: is_bind_event, BIND_OBJEDIT_FINISH, BIND_OBJEDIT_DELPOINT,&
       BIND_OBJEDIT_DELETE
    class(window), intent(inout), target :: w
    type(representation), intent(inout) :: r
    type(objedit_input), intent(inout) :: inp

    integer :: ikind, nh, i, k
    real*8 :: d(2)
    real*8, allocatable :: xh(:,:)
    logical :: ok, changed

    real*8, parameter :: free_px = 2d0 ! spacing of the recorded freehand points (pixels)
    real*8, parameter :: simplify_px = 0.75d0 ! freehand simplification tolerance (pixels)

    ikind = inp%itool - objtool_kind0
    changed = .false.
    if (r%planar%isel > r%planar%nshape) r%planar%isel = 0
    if (w%oe%op == objop_handle .or. w%oe%op == objop_move) then
       if (w%oe%iitem < 1 .or. w%oe%iitem > r%planar%nshape) w%oe%op = objop_none
    end if

    ! press
    if (inp%press) then
       if (w%oe%op == objop_clicks) then
          ! polygon/polyline in progress: a double click finishes it if it
          ! has enough points (its first click placed the last one); a
          ! click adds a point
          if (inp%dbl .and. clicks_complete(w)) then
             call planar_finish_clicks(w,r)
          elseif (norm2(inp%xm - w%oe%planar%sh%x(:,w%oe%planar%sh%npt)) > objedit_drag_px * inp%pxs) then
             call planar_addpoint(w%oe%planar%sh,inp%xm)
          end if
       elseif (inp%itool == objtool_select) then
          w%oe%op = objop_none
          w%oe%x0 = inp%xm
          w%oe%moved = .false.

          ! a handle of the selected (and shown) shape
          if (r%planar%isel > 0) then
             if (r%planar%shape(r%planar%isel)%shown) then
                call planar_handles(r%planar%shape(r%planar%isel),nh,xh)
                do i = nh, 1, -1
                   if (norm2(inp%xm - xh(:,i)) <= objedit_hit_px * inp%pxs) then
                      w%oe%op = objop_handle
                      w%oe%ih = i
                      w%oe%iitem = r%planar%isel
                      exit
                   end if
                end do
             end if
          end if

          ! a shape: the topmost under the mouse
          if (w%oe%op == objop_none) then
             k = planar_shape_at(r,inp%xm,objedit_hit_px*inp%pxs)
             r%planar%isel = k
             if (k > 0) then
                w%oe%op = objop_move
                w%oe%iitem = k
             end if
          end if
          if (w%oe%op /= objop_none) then
             w%oe%planar%sh0 = r%planar%shape(w%oe%iitem)
             w%oe%planar%sh = w%oe%planar%sh0
          end if
       elseif (inp%itool == objtool_remove) then
          ! remove the shape under the mouse
          k = planar_shape_at(r,inp%xm,objedit_hit_px*inp%pxs)
          if (k > 0) then
             call planar_delete(r%planar,k)
             changed = .true.
          end if
       elseif (ikind >= 1 .and. ikind <= planarkind_NUM) then
          ! start a new shape: an ellipse, rectangle or arrow grows from
          ! the press as if its far corner (handle 3) or its tip (handle 2)
          ! were dragged
          w%oe%x0 = inp%xm
          w%oe%moved = .false.
          w%oe%planar%sh = planar_template(r%planar)
          w%oe%planar%sh%kind = ikind
          w%oe%planar%sh%xc = inp%xm
          if (ikind == planarkind_ellipse .or. ikind == planarkind_rect) then
             w%oe%op = objop_drag
             w%oe%ih = 3
          elseif (ikind == planarkind_arrow .or. ikind == planarkind_curve) then
             ! the tip is the last point; a curve keeps its middle point
             ! at the default bend while it is drawn (the heads, like the
             ! rest of the style, come from the selected shape)
             w%oe%op = objop_drag
             call planar_addpoint(w%oe%planar%sh,inp%xm)
             call planar_addpoint(w%oe%planar%sh,inp%xm)
             if (ikind == planarkind_curve) call planar_addpoint(w%oe%planar%sh,inp%xm)
             w%oe%ih = w%oe%planar%sh%npt
          elseif (ikind == planarkind_freehand) then
             w%oe%op = objop_free
             call planar_addpoint(w%oe%planar%sh,inp%xm)
          else
             w%oe%op = objop_clicks
             call planar_addpoint(w%oe%planar%sh,inp%xm)
          end if
          w%oe%planar%sh0 = w%oe%planar%sh
       end if
    end if

    ! the exit bind finishes a polygon or polyline in progress
    if (inp%exitev .and. w%oe%op == objop_clicks) then
       call planar_finish_clicks(w,r)
       inp%exitev = .false.
    end if

    ! polygon/polyline in progress: finish it, or remove its last point
    if (w%oe%op == objop_clicks .and. inp%hover) then
       if (is_bind_event(BIND_OBJEDIT_FINISH,norepeat=.true.,iview=w%id)) then
          call planar_finish_clicks(w,r)
       elseif (is_bind_event(BIND_OBJEDIT_DELPOINT,iview=w%id)) then
          w%oe%planar%sh%npt = w%oe%planar%sh%npt - 1
          if (w%oe%planar%sh%npt == 0) w%oe%op = objop_none
       end if
    end if

    ! drag
    if (w%oe%op == objop_drag .or. w%oe%op == objop_free .or. w%oe%op == objop_handle .or.&
       w%oe%op == objop_move) then
       if (norm2(inp%xm - w%oe%x0) > objedit_drag_px * inp%pxs) w%oe%moved = .true.
       if (inp%down) then
          if (w%oe%op == objop_drag .or. w%oe%op == objop_handle) then
             w%oe%planar%sh = w%oe%planar%sh0
             call planar_drag_handle(w%oe%planar%sh,w%oe%planar%sh0,w%oe%ih,inp%xm,inp%constrain)
             if (w%oe%op == objop_drag .and. w%oe%planar%sh%kind == planarkind_curve) &
                call planar_curve_default_bend(w%oe%planar%sh)
          elseif (w%oe%op == objop_free) then
             if (norm2(inp%xm - w%oe%planar%sh%x(:,w%oe%planar%sh%npt)) > free_px * inp%pxs) &
                call planar_addpoint(w%oe%planar%sh,inp%xm)
          elseif (w%oe%op == objop_move) then
             w%oe%planar%sh = w%oe%planar%sh0
             d = inp%xm - w%oe%x0
             if (inp%constrain) then
                ! along the dominant direction only
                if (abs(d(1)) > abs(d(2))) then
                   d(2) = 0d0
                else
                   d(1) = 0d0
                end if
             end if
             call planar_move(w%oe%planar%sh,d)
          end if
       else
          ! release: commit
          if (w%oe%op == objop_drag) then
             ok = w%oe%moved
             if (ok) then
                if (w%oe%planar%sh%kind == planarkind_arrow .or. w%oe%planar%sh%kind == planarkind_curve) then
                   ok = (norm2(w%oe%planar%sh%x(:,w%oe%planar%sh%npt) - w%oe%planar%sh%x(:,1)) > objedit_drag_px * inp%pxs)
                else
                   ok = all(w%oe%planar%sh%hs > 0d0)
                end if
             end if
             if (ok) then
                call planar_append(r%planar,w%oe%planar%sh)
                changed = .true.
             end if
          elseif (w%oe%op == objop_free) then
             call planar_simplify(w%oe%planar%sh%npt,w%oe%planar%sh%x,simplify_px * inp%pxs)
             if (w%oe%planar%sh%npt >= 2) then
                call planar_append(r%planar,w%oe%planar%sh)
                changed = .true.
             end if
          elseif (w%oe%moved) then
             r%planar%shape(w%oe%iitem) = w%oe%planar%sh
             changed = .true.
          end if
          w%oe%op = objop_none
       end if
    end if

    ! delete the selected shape (draw_view leaves Delete to this mode)
    if (inp%hover .and. inp%itool == objtool_select .and. w%oe%op == objop_none .and.&
       r%planar%isel > 0) then
       if (is_bind_event(BIND_OBJEDIT_DELETE,norepeat=.true.,iview=w%id)) then
          call planar_delete(r%planar,r%planar%isel)
          changed = .true.
       end if
    end if

    ! the shapes are drawn by the scene: rebuild it
    if (changed) w%sc%forcebuildlists = .true.

    ! the preview and the handles
    call planar_view_overlay(w,r,inp%itool,inp%xm,inp%pxs)

  end subroutine planar_events

  !> The object editing view mode (vm_objedit) for the 3D shapes r, with
  !> the input of this frame inp: draw new shapes with the tool (a press
  !> and a drag), select them and edit them (move them, drag their
  !> handles), or remove them. The points snap to the atom under the
  !> mouse unless the no-snap bind is held; otherwise they lie on the
  !> plane facing the camera through the scene center (a new shape) or
  !> through the dragged point. The shape being drawn or dragged is shown
  !> as a wireframe over the view and goes to the object on release.
  module subroutine shapes_events(w,r,inp)
    use representations, only: representation, shapekind_NUM, shapekind_sphere,&
       shapekind_box, shapes_template, shapes_append, shapes_delete
    use keybindings, only: is_bind_event, BIND_OBJEDIT_DELETE
    class(window), intent(inout), target :: w
    type(representation), intent(inout) :: r
    type(objedit_input), intent(inout) :: inp

    integer :: ikind, nh, i, k
    real*8 :: tol, p(3), d(3), xh(3,4)
    real(c_float) :: t(3), zc
    logical :: changed, ok

    ikind = inp%itool - objtool_kind0
    changed = .false.
    tol = objedit_hit_px * inp%pxs * 0.5d0 * w%FBOside ! pick radius in texture pixels
    if (r%shapes%isel > r%shapes%nshape) r%shapes%isel = 0
    if (w%oe%op == objop_handle .or. w%oe%op == objop_move) then
       if (w%oe%iitem < 1 .or. w%oe%iitem > r%shapes%nshape) w%oe%op = objop_none
    end if

    associate (st => w%oe%shapes)
      ! press
      if (inp%press) then
         w%oe%x0 = inp%xm
         w%oe%moved = .false.
         if (inp%itool == objtool_select) then
            w%oe%op = objop_none

            ! a handle of the selected (and shown) shape
            if (r%shapes%isel > 0) then
               if (r%shapes%shape(r%shapes%isel)%shown) then
                  call shape_handles(w,r%shapes%shape(r%shapes%isel),nh,xh)
                  do i = nh, 1, -1
                     t = to_tex(w,xh(:,i))
                     if (norm2(inp%tex - t(1:2)) <= tol) then
                        w%oe%op = objop_handle
                        w%oe%ih = i
                        w%oe%iitem = r%shapes%isel
                        st%zdep = t(3)
                        exit
                     end if
                  end do
               end if
            end if

            ! a shape: the nearest under the mouse
            if (w%oe%op == objop_none) then
               k = shape_at(w,r,inp%tex,tol,st%zdep)
               r%shapes%isel = k
               if (k > 0) then
                  w%oe%op = objop_move
                  w%oe%iitem = k
                  st%xg = from_tex(w,inp%tex,st%zdep)
               end if
            end if
            if (w%oe%op /= objop_none) then
               st%sh0 = r%shapes%shape(w%oe%iitem)
               st%sh = st%sh0
            end if
         elseif (inp%itool == objtool_remove) then
            ! remove the shape under the mouse
            k = shape_at(w,r,inp%tex,tol,zc)
            if (k > 0) then
               call shapes_delete(r%shapes,k)
               changed = .true.
            end if
         elseif (ikind >= 1 .and. ikind <= shapekind_NUM) then
            ! start a new shape at the press: its anchor, which the drag
            ! then gives a size, on the plane through the scene center
            t = to_tex(w,real(w%sc%scenecenter,8) + molshift(w))
            p = snap_point(w,inp,t(3))
            t = to_tex(w,p)
            st%zdep = t(3)
            st%sh = shapes_template(r%shapes,ikind)
            st%sh%x1 = p
            st%sh0 = st%sh
            w%oe%op = objop_drag
         end if
      end if

      ! drag
      if (w%oe%op == objop_drag .or. w%oe%op == objop_handle .or. w%oe%op == objop_move) then
         if (norm2(inp%xm - w%oe%x0) > objedit_drag_px * inp%pxs) w%oe%moved = .true.
         if (inp%down) then
            st%sh = st%sh0
            if (w%oe%op == objop_drag) then
               call drag_new(w,st%sh,snap_point(w,inp,st%zdep),inp%constrain)
            elseif (w%oe%op == objop_handle) then
               call drag_handle(w,st%sh,st%sh0,w%oe%ih,snap_point(w,inp,st%zdep),inp%constrain)
            else
               ! the whole shape, on the plane of the grabbed point
               d = from_tex(w,inp%tex,st%zdep) - st%xg
               if (inp%constrain) d = axis_snap(w,d)
               st%sh%x1 = st%sh0%x1 + d
            end if
         else
            ! release: commit
            if (w%oe%op == objop_drag) then
               ok = w%oe%moved
               if (ok) then
                  if (st%sh%kind == shapekind_sphere) then
                     ok = (st%sh%rad > 1d-6)
                  elseif (st%sh%kind == shapekind_box) then
                     ok = all(norm2(st%sh%v,1) > 1d-6)
                  else
                     ok = (norm2(st%sh%v(:,1)) > 1d-6)
                  end if
               end if
               if (ok) then
                  call shapes_append(r%shapes,st%sh)
                  changed = .true.
               end if
            elseif (w%oe%moved) then
               r%shapes%shape(w%oe%iitem) = st%sh
               changed = .true.
            end if
            w%oe%op = objop_none
         end if
      end if
    end associate

    ! delete the selected shape (draw_view leaves Delete to this mode)
    if (inp%hover .and. inp%itool == objtool_select .and. w%oe%op == objop_none .and.&
       r%shapes%isel > 0) then
       if (is_bind_event(BIND_OBJEDIT_DELETE,norepeat=.true.,iview=w%id)) then
          call shapes_delete(r%shapes,r%shapes%isel)
          changed = .true.
       end if
    end if

    ! the shapes are drawn by the scene: rebuild it, keeping the camera
    ! where it is (the shapes count in the scene extent, which would
    ! otherwise move the view under the mouse)
    if (changed) then
       w%sc%forcebuildlists = .true.
       w%sc%nextbuildlists_fixcam = .true.
    end if

    ! the preview and the handles
    call shapes_overlay(w,r,inp%itool)

  end subroutine shapes_events

  !xx! private procedures

  !> The topmost shown shape of r at position x (vm_objedit), 0 if none:
  !> the outlines (within tol of their edge) first, then the insides of
  !> the closed shapes.
  function planar_shape_at(r,x,tol) result(k)
    use representations, only: representation, planar_hit
    type(representation), intent(in) :: r
    real*8, intent(in) :: x(2), tol
    integer :: k

    integer :: ipass
    logical :: onstroke

    do ipass = 1, 2
       onstroke = (ipass == 1)
       do k = r%planar%nshape, 1, -1
          if (.not.r%planar%shape(k)%shown) cycle
          if (planar_hit(r%planar%shape(k),x,merge(tol,0d0,onstroke),onstroke)) return
       end do
    end do
    k = 0

  end function planar_shape_at

  !> Whether the polygon or polyline being drawn (vm_objedit) has enough points.
  function clicks_complete(w) result(ok)
    use representations, only: planarkind_polygon
    class(window), intent(in) :: w
    logical :: ok

    ok = (w%oe%planar%sh%npt >= merge(3,2,w%oe%planar%sh%kind == planarkind_polygon))

  end function clicks_complete

  !> Finish the polygon or polyline being drawn point by point (vm_objedit),
  !> adding it to the shapes of r if it has enough points.
  subroutine planar_finish_clicks(w,r)
    use representations, only: representation, planar_append
    class(window), intent(inout) :: w
    type(representation), intent(inout) :: r

    if (clicks_complete(w)) then
       call planar_append(r%planar,w%oe%planar%sh)
       w%sc%forcebuildlists = .true.
    end if
    w%oe%op = objop_none

  end subroutine planar_finish_clicks

  !> Append the point x to the shape sh.
  subroutine planar_addpoint(sh,x)
    use representations, only: planar_shape
    use types, only: realloc
    type(planar_shape), intent(inout) :: sh
    real*8, intent(in) :: x(2)

    if (.not.allocated(sh%x)) then
       allocate(sh%x(2,16))
    elseif (sh%npt >= size(sh%x,2)) then
       call realloc(sh%x,2,2*size(sh%x,2))
    end if
    sh%npt = sh%npt + 1
    sh%x(:,sh%npt) = x

  end subroutine planar_addpoint

  !> Draw over the view w the planar shape being drawn or dragged, and
  !> the outline and handles of the selected shape (select tool) of the
  !> planar shapes representation r. itool = the editor's tool, xm =
  !> the mouse position (NDC of the render buffer), pxs = NDC size of a
  !> screen pixel.
  subroutine planar_view_overlay(w,r,itool,xm,pxs)
    use representations, only: representation, planar_shape, planar_handles,&
       planar_rotation_handle, planarkind_polyline
    use gui_main, only: ColorHighlightSelectScene
    class(window), intent(inout), target :: w
    type(representation), intent(in) :: r
    integer, intent(in) :: itool
    real*8, intent(in) :: xm(2), pxs

    type(c_ptr) :: dl
    type(ImVec2) :: p1, p2
    integer(c_int) :: colw, colk, colh

    dl = igGetWindowDrawList()
    call ImDrawList_PushClipRect(dl,w%v_rmin,w%v_rmax,.true._c_bool)
    colw = igGetColorU32_Vec4(ImVec4(1._c_float,1._c_float,1._c_float,1._c_float))
    colk = igGetColorU32_Vec4(ImVec4(0._c_float,0._c_float,0._c_float,1._c_float))
    colh = igGetColorU32_Vec4(ImVec4(ColorHighlightSelectScene(1),ColorHighlightSelectScene(2),&
       ColorHighlightSelectScene(3),1._c_float))

    ! the shape being drawn or dragged, in its own color and width
    if (w%oe%op /= objop_none) then
       call preview_shape(w%oe%planar%sh,max(real(w%oe%planar%sh%width/pxs,c_float),1._c_float),&
          igGetColorU32_Vec4(ImVec4(w%oe%planar%sh%rgb(1),w%oe%planar%sh%rgb(2),w%oe%planar%sh%rgb(3),&
          max(w%oe%planar%sh%alpha,0.3_c_float))))
       ! the segment to the mouse of a polygon/polyline in progress
       if (w%oe%op == objop_clicks) then
          p1 = ndc_to_mouse(w%oe%planar%sh%x(:,w%oe%planar%sh%npt))
          p2 = ndc_to_mouse(xm)
          call ImDrawList_AddLine(dl,p1,p2,colh,1._c_float)
          if (w%oe%planar%sh%kind /= planarkind_polyline .and. w%oe%planar%sh%npt >= 2) then
             p1 = ndc_to_mouse(w%oe%planar%sh%x(:,1))
             call ImDrawList_AddLine(dl,p1,p2,colh,1._c_float)
          end if
       end if
    end if

    ! the outline and the handles of the selected shape, if shown
    if (itool == objtool_select) then
       if (w%oe%op == objop_handle .or. w%oe%op == objop_move) then
          call selected_shape(w%oe%planar%sh)
       elseif (r%planar%isel >= 1 .and. r%planar%isel <= r%planar%nshape) then
          if (r%planar%shape(r%planar%isel)%shown) call selected_shape(r%planar%shape(r%planar%isel))
       end if
    end if
    call ImDrawList_PopClipRect(dl)

  contains
    !> Draw the outline (thin, highlighted) and the handles of the
    !> selected shape sh.
    subroutine selected_shape(sh)
      type(planar_shape), intent(in) :: sh

      integer :: nh, i, irot
      real*8, allocatable :: xh(:,:)
      type(ImVec2) :: q1, q2

      call preview_shape(sh,1._c_float,colh)
      call planar_handles(sh,nh,xh)
      irot = planar_rotation_handle(sh)
      if (irot > 0) then
         ! the rotation handle hangs from the midpoint of the top edge (handle 7)
         q1 = ndc_to_mouse(xh(:,7))
         q2 = ndc_to_mouse(xh(:,irot))
         call ImDrawList_AddLine(dl,q1,q2,colh,1._c_float)
      end if
      do i = 1, nh
         call draw_handle(dl,ndc_to_mouse(xh(:,i)),i == irot,colw,colk)
      end do

    end subroutine selected_shape

    !> Draw the center line of the shape sh with thickness thick (pixels)
    !> and color col, closed if the shape is.
    subroutine preview_shape(sh,thick,col)
      use representations, only: planar_path, planar_isclosed
      type(planar_shape), intent(in) :: sh
      real(c_float), intent(in) :: thick
      integer(c_int), intent(in) :: col

      integer :: n, j, nseg
      real*8, allocatable :: xp(:,:)
      type(ImVec2) :: q1, q2

      call planar_path(sh,n,xp)
      if (n < 2) return
      nseg = n - 1
      if (planar_isclosed(sh)) nseg = n
      do j = 1, nseg
         q1 = ndc_to_mouse(xp(:,j))
         q2 = ndc_to_mouse(xp(:,modulo(j,n)+1))
         call ImDrawList_AddLine(dl,q1,q2,col,thick)
      end do

    end subroutine preview_shape

    !> Mouse (screen) position of the point x in the NDC of the render buffer.
    function ndc_to_mouse(x) result(p)
      real*8, intent(in) :: x(2)
      type(ImVec2) :: p

      p%x = real(0.5d0 * (x(1) + 1d0) * w%FBOside,c_float)
      p%y = real(0.5d0 * (x(2) + 1d0) * w%FBOside,c_float)
      call w%texpos_to_mousepos(p)

    end function ndc_to_mouse

  end subroutine planar_view_overlay

  !> Shift from the scene coordinates of view w to the absolute frame of
  !> the shapes (cartesian bohr): the molecular origin for a molecule,
  !> zero for a crystal.
  function molshift(w) result(x)
    use systems, only: sys
    class(window), intent(in) :: w
    real*8 :: x(3)

    x = 0d0
    if (sys(w%isys)%c%ismolecule) x = sys(w%isys)%c%molx0

  end function molshift

  !> Texture position (pixels, and depth) of the absolute-frame point x
  !> in view w.
  function to_tex(w,x) result(t)
    use utils, only: mult
    class(window), intent(inout), target :: w
    real*8, intent(in) :: x(3)
    real(c_float) :: t(3)

    real(c_float) :: v(3)

    v = real(x - molshift(w),c_float)
    call mult(t,w%sc%world,v)
    call w%world_to_texpos(t)

  end function to_tex

  !> Absolute-frame point at the texture position tex (pixels) and depth
  !> z in view w.
  function from_tex(w,tex,z) result(x)
    use utils, only: invmult
    class(window), intent(inout), target :: w
    real(c_float), intent(in) :: tex(2), z
    real*8 :: x(3)

    real(c_float) :: v(3)

    v = (/tex(1),tex(2),z/)
    call w%texpos_to_world(v)
    call invmult(v,w%sc%world)
    x = real(v,8) + molshift(w)

  end function from_tex

  !> The point under the mouse (input inp) in view w: the atom under it,
  !> unless the no-snap bind is held, or the point at texture depth z.
  function snap_point(w,inp,z) result(p)
    use systems, only: sys
    class(window), intent(inout), target :: w
    type(objedit_input), intent(in) :: inp
    real(c_float), intent(in) :: z
    real*8 :: p(3)

    integer :: i

    i = w%mousepos_idx(1)
    if (.not.inp%nosnap .and. i >= 1 .and. i <= sys(w%isys)%c%ncel) then
       p = sys(w%isys)%c%x2c(sys(w%isys)%c%atcel(i)%x + w%mousepos_idx(2:4)) + molshift(w)
    else
       p = from_tex(w,inp%tex,z)
    end if

  end function snap_point

  !> Unit vector along the screen right direction at the absolute-frame
  !> point x of view w.
  function cam_right(w,x) result(u)
    class(window), intent(inout), target :: w
    real*8, intent(in) :: x(3)
    real*8 :: u(3)

    real(c_float) :: t(3)

    t = to_tex(w,x)
    u = from_tex(w,(/t(1)+10._c_float,t(2)/),t(3)) - x
    if (norm2(u) > 1d-12) then
       u = u / norm2(u)
    else
       u = (/1d0,0d0,0d0/)
    end if

  end function cam_right

  !> The point of the sphere sh on its rim to the right on the screen of
  !> view w (absolute frame).
  function sphere_rim(w,sh) result(x)
    use representations, only: rep_shape
    class(window), intent(inout), target :: w
    type(rep_shape), intent(in) :: sh
    real*8 :: x(3)

    x = sh%x1 + sh%rad * cam_right(w,sh%x1)

  end function sphere_rim

  !> Texture positions (pixels, and depth) of the corners of the box sh
  !> in view w: corner i is the anchor plus the edge vectors v_j for the
  !> bits j-1 set in i-1.
  function box_corners(w,sh) result(c)
    use representations, only: rep_shape
    class(window), intent(inout), target :: w
    type(rep_shape), intent(in) :: sh
    real(c_float) :: c(3,8)

    integer :: i, j
    real*8 :: x(3)

    do i = 1, 8
       x = sh%x1
       do j = 1, 3
          if (btest(i-1,j-1)) x = x + sh%v(:,j)
       end do
       c(:,i) = to_tex(w,x)
    end do

  end function box_corners

  !> The axes the drags snap to in view w (columns, unit vectors): the
  !> cell axes of a crystal, the cartesian axes of a molecule.
  function snap_axes(w) result(a)
    use systems, only: sys
    class(window), intent(in) :: w
    real*8 :: a(3,3)

    integer :: k

    if (sys(w%isys)%c%ismolecule) then
       a = 0d0
       do k = 1, 3
          a(k,k) = 1d0
       end do
    else
       do k = 1, 3
          a(:,k) = sys(w%isys)%c%m_x2c(:,k) / norm2(sys(w%isys)%c%m_x2c(:,k))
       end do
    end if

  end function snap_axes

  !> The component of the vector d along the snap axis of view w closest
  !> to its direction.
  function axis_snap(w,d) result(ds)
    class(window), intent(in) :: w
    real*8, intent(in) :: d(3)
    real*8 :: ds(3)

    real*8 :: a(3,3)
    integer :: k

    a = snap_axes(w)
    k = maxloc(abs(matmul(d,a)),1)
    ds = dot_product(d,a(:,k)) * a(:,k)

  end function axis_snap

  !> The shape sh (a new one, at its anchor sh%x1) drawn out to the
  !> point p: the radius of a sphere, the end of an arrow, cone or
  !> cylinder (constrain: along a snap axis), or the opposite corner of a
  !> box with edges along the snap axes. A box drawn on the plane facing
  !> the camera gets, along the third axis, the shorter of the two drawn
  !> edges (constrain: a cube).
  subroutine drag_new(w,sh,p,constrain)
    use representations, only: rep_shape, shapekind_sphere, shapekind_box
    class(window), intent(in) :: w
    type(rep_shape), intent(inout) :: sh
    real*8, intent(in) :: p(3)
    logical, intent(in) :: constrain

    real*8 :: a(3,3), q(3), e(3), d(3)
    integer :: k

    d = p - sh%x1
    if (sh%kind == shapekind_sphere) then
       sh%rad = norm2(d)
    elseif (sh%kind == shapekind_box) then
       ! the drag in the basis of the snap axes
       a = snap_axes(w)
       q = solve3(a,d)
       e = abs(q)
       if (constrain) then
          e = maxval(e)
       else
          ! the shortest edge (the one along the view) takes the middle one
          e(minloc(e,1)) = sum(e) - maxval(e) - minval(e)
       end if
       do k = 1, 3
          sh%v(:,k) = sign(e(k),q(k)) * a(:,k)
       end do
    else
       if (constrain) d = axis_snap(w,d)
       sh%v(:,1) = d
    end if

  end subroutine drag_new

  !> Drag handle ih (see shape_handles) of the shape to the point p. sh0
  !> is the shape as it was at the press, and the result goes to sh.
  !> Moving the anchor keeps the other ends in place; constrain keeps an
  !> end on a snap axis from the anchor.
  subroutine drag_handle(w,sh,sh0,ih,p,constrain)
    use representations, only: rep_shape, shapekind_sphere, shapekind_box
    class(window), intent(in) :: w
    type(rep_shape), intent(inout) :: sh
    type(rep_shape), intent(in) :: sh0
    integer, intent(in) :: ih
    real*8, intent(in) :: p(3)
    logical, intent(in) :: constrain

    integer :: k
    real*8 :: d(3)

    if (sh0%kind == shapekind_sphere) then
       if (ih == 1) then
          sh%x1 = p
       else
          sh%rad = norm2(p - sh0%x1)
       end if
    elseif (ih == 1) then
       ! the anchor; the ends stay where they are
       if (sh0%kind == shapekind_box) then
          do k = 1, 3
             sh%v(:,k) = sh0%v(:,k) + sh0%x1 - p
          end do
       else
          sh%v(:,1) = sh0%v(:,1) + sh0%x1 - p
       end if
       sh%x1 = p
    else
       d = p - sh0%x1
       if (constrain) d = axis_snap(w,d)
       sh%v(:,ih-1) = d
    end if

  end subroutine drag_handle

  !> The editing handles of the shape sh in view w (absolute frame): the
  !> anchor (1) and, for a sphere, a radius knob on its rim to the right
  !> on the screen (2), for an arrow, cone or cylinder, its other end
  !> (2), and for a box, the ends of its three edges from the anchor
  !> (2-4).
  subroutine shape_handles(w,sh,nh,xh)
    use representations, only: rep_shape, shapekind_sphere, shapekind_box
    class(window), intent(inout), target :: w
    type(rep_shape), intent(in) :: sh
    integer, intent(out) :: nh
    real*8, intent(out) :: xh(3,4)

    integer :: k

    xh(:,1) = sh%x1
    if (sh%kind == shapekind_sphere) then
       nh = 2
       xh(:,2) = sphere_rim(w,sh)
    elseif (sh%kind == shapekind_box) then
       nh = 4
       do k = 1, 3
          xh(:,k+1) = sh%x1 + sh%v(:,k)
       end do
    else
       nh = 2
       xh(:,2) = sh%x1 + sh%v(:,1)
    end if

  end subroutine shape_handles

  !> The shown shape of r nearest to the camera under the texture
  !> position tex in view w, within tol pixels of its outline, and the
  !> texture depth of the hit (z); 0 if none.
  function shape_at(w,r,tex,tol,z) result(ksel)
    use representations, only: representation
    class(window), intent(inout), target :: w
    type(representation), intent(in) :: r
    real(c_float), intent(in) :: tex(2)
    real*8, intent(in) :: tol
    real(c_float), intent(out) :: z
    integer :: ksel

    integer :: k
    real(c_float) :: zk

    ksel = 0
    z = 1._c_float
    do k = 1, r%shapes%nshape
       if (.not.r%shapes%shape(k)%shown) cycle
       if (.not.shape_hit(w,r%shapes%shape(k),tex,tol,zk)) cycle
       if (ksel == 0 .or. zk < z) then
          ksel = k
          z = zk
       end if
    end do

  end function shape_at

  !> Whether the texture position tex hits the shape sh in view w
  !> (within tol pixels): inside the outline of a sphere, near the axis
  !> of an arrow, cone or cylinder (within its projected thickness),
  !> near the edges of a box or inside one of its faces. z is the
  !> texture depth of the hit.
  function shape_hit(w,sh,tex,tol,z) result(ok)
    use representations, only: rep_shape, shapekind_sphere, shapekind_box,&
       shapekind_arrow
    class(window), intent(inout), target :: w
    type(rep_shape), intent(in) :: sh
    real(c_float), intent(in) :: tex(2)
    real*8, intent(in) :: tol
    real(c_float), intent(out) :: z
    logical :: ok

    real(c_float) :: a(3), b(3), c(3,8), m(3), e(3)
    real*8 :: rpx, xm(3), radv
    integer :: i

    integer, parameter :: faces(4,6) = reshape((/1,2,4,3, 5,6,8,7, 1,2,6,5, 3,4,8,7,&
       1,3,7,5, 2,4,8,6/),shape(faces))

    ok = .false.
    z = 1._c_float
    if (sh%kind == shapekind_sphere) then
       a = to_tex(w,sh%x1)
       b = to_tex(w,sphere_rim(w,sh))
       rpx = norm2(b(1:2) - a(1:2))
       ok = (norm2(tex - a(1:2)) <= rpx + tol)
       z = a(3)
    elseif (sh%kind == shapekind_box) then
       c = box_corners(w,sh)
       z = minval(c(3,:))
       do i = 1, 12
          if (segdist2(tex,c(1:2,box_edges(1,i)),c(1:2,box_edges(2,i))) <= tol) ok = .true.
       end do
       do i = 1, 6
          if (in_quad(tex,c(1:2,faces(:,i)))) ok = .true.
       end do
    else
       a = to_tex(w,sh%x1)
       b = to_tex(w,sh%x1 + sh%v(:,1))
       xm = sh%x1 + 0.5d0 * sh%v(:,1)
       radv = sh%rad
       if (sh%kind == shapekind_arrow) radv = sh%rad * max(sh%headr,1d0)
       m = to_tex(w,xm)
       e = to_tex(w,xm + radv * cam_right(w,xm))
       rpx = norm2(e(1:2) - m(1:2))
       ok = (segdist2(tex,a(1:2),b(1:2)) <= rpx + tol)
       z = min(a(3),b(3))
    end if

  end function shape_hit

  !> Distance from the 2D point p to the segment ab.
  function segdist2(p,a,b) result(d)
    real(c_float), intent(in) :: p(2), a(2), b(2)
    real*8 :: d

    real*8 :: ab(2), t, l2

    ab = b - a
    l2 = dot_product(ab,ab)
    if (l2 < 1d-12) then
       d = norm2(real(p - a,8))
    else
       t = min(max(dot_product(real(p - a,8),ab) / l2,0d0),1d0)
       d = norm2(real(p - a,8) - t * ab)
    end if

  end function segdist2

  !> Whether the 2D point p is inside the convex quadrilateral q (either
  !> orientation).
  function in_quad(p,q) result(ok)
    real(c_float), intent(in) :: p(2), q(2,4)
    logical :: ok

    real*8 :: s(4)
    integer :: i, j

    do i = 1, 4
       j = modulo(i,4) + 1
       s(i) = (q(1,j)-q(1,i)) * (p(2)-q(2,i)) - (q(2,j)-q(2,i)) * (p(1)-q(1,i))
    end do
    ok = all(s >= 0d0) .or. all(s <= 0d0)
    if (ok) ok = any(abs(s) > 1d-12)

  end function in_quad

  !> The solution x of a x = d, for the 3x3 matrix a (zero if singular).
  function solve3(a,d) result(x)
    use tools_math, only: cross
    real*8, intent(in) :: a(3,3), d(3)
    real*8 :: x(3)

    real*8 :: det

    det = dot_product(a(:,1),cross(a(:,2),a(:,3)))
    x = 0d0
    if (abs(det) < 1d-14) return
    x(1) = dot_product(d,cross(a(:,2),a(:,3))) / det
    x(2) = dot_product(a(:,1),cross(d,a(:,3))) / det
    x(3) = dot_product(a(:,1),cross(a(:,2),d)) / det

  end function solve3

  !> Draw over the view w the shape being drawn or dragged (in its own
  !> color), and the outline and handles of the selected shape of r
  !> (select tool itool), as wireframes on their screen projection.
  subroutine shapes_overlay(w,r,itool)
    use representations, only: representation, rep_shape, shapekind_sphere
    use gui_main, only: ColorHighlightSelectScene
    class(window), intent(inout), target :: w
    type(representation), intent(in) :: r
    integer, intent(in) :: itool

    type(c_ptr) :: dl
    integer(c_int) :: colw, colk, colh

    dl = igGetWindowDrawList()
    call ImDrawList_PushClipRect(dl,w%v_rmin,w%v_rmax,.true._c_bool)
    colw = igGetColorU32_Vec4(ImVec4(1._c_float,1._c_float,1._c_float,1._c_float))
    colk = igGetColorU32_Vec4(ImVec4(0._c_float,0._c_float,0._c_float,1._c_float))
    colh = igGetColorU32_Vec4(ImVec4(ColorHighlightSelectScene(1),ColorHighlightSelectScene(2),&
       ColorHighlightSelectScene(3),1._c_float))

    ! the shape being drawn or dragged
    if (w%oe%op /= objop_none) then
       associate (sh => w%oe%shapes%sh)
         call wireframe(sh,igGetColorU32_Vec4(ImVec4(sh%rgb(1),sh%rgb(2),sh%rgb(3),1._c_float)),&
            2._c_float)
       end associate
    end if

    ! the outline and the handles of the selected shape, if shown
    if (itool == objtool_select) then
       if (w%oe%op == objop_handle .or. w%oe%op == objop_move) then
          call selected(w%oe%shapes%sh)
       elseif (r%shapes%isel >= 1 .and. r%shapes%isel <= r%shapes%nshape) then
          if (r%shapes%shape(r%shapes%isel)%shown) call selected(r%shapes%shape(r%shapes%isel))
       end if
    end if
    call ImDrawList_PopClipRect(dl)

  contains
    !> The outline and the handles of the selected shape sh.
    subroutine selected(sh)
      type(rep_shape), intent(in) :: sh

      integer :: nh, i
      real*8 :: xh(3,4)

      call wireframe(sh,colh,1._c_float)
      call shape_handles(w,sh,nh,xh)
      do i = 1, nh
         call draw_handle(dl,tex_to_mouse(to_tex(w,xh(:,i))),sh%kind == shapekind_sphere .and. i == 2,&
            colw,colk)
      end do

    end subroutine selected

    !> The wireframe of the shape sh with color col and thickness thick:
    !> the outline of a sphere, the axis of an arrow, cone or cylinder,
    !> the edges of a box.
    subroutine wireframe(sh,col,thick)
      use representations, only: shapekind_box
      type(rep_shape), intent(in) :: sh
      integer(c_int), intent(in) :: col
      real(c_float), intent(in) :: thick

      integer :: j
      real(c_float) :: c(3,8)
      type(ImVec2) :: p1, p2

      if (sh%kind == shapekind_sphere) then
         p1 = tex_to_mouse(to_tex(w,sh%x1))
         p2 = tex_to_mouse(to_tex(w,sphere_rim(w,sh)))
         call ImDrawList_AddCircle(dl,p1,real(hypot(p2%x-p1%x,p2%y-p1%y),c_float),col,48_c_int,thick)
      elseif (sh%kind == shapekind_box) then
         c = box_corners(w,sh)
         do j = 1, 12
            p1 = tex_to_mouse(c(:,box_edges(1,j)))
            p2 = tex_to_mouse(c(:,box_edges(2,j)))
            call ImDrawList_AddLine(dl,p1,p2,col,thick)
         end do
      else
         p1 = tex_to_mouse(to_tex(w,sh%x1))
         p2 = tex_to_mouse(to_tex(w,sh%x1 + sh%v(:,1)))
         call ImDrawList_AddLine(dl,p1,p2,col,thick)
      end if

    end subroutine wireframe

    !> Mouse (screen) position of the texture position t.
    function tex_to_mouse(t) result(p)
      real(c_float), intent(in) :: t(3)
      type(ImVec2) :: p

      p%x = t(1)
      p%y = t(2)
      call w%texpos_to_mousepos(p)

    end function tex_to_mouse
  end subroutine shapes_overlay

  !> Draw a handle at the screen position q on the draw list dl: a
  !> square, or a disk if round, filled with colw and outlined with colk.
  subroutine draw_handle(dl,q,round,colw,colk)
    type(c_ptr), intent(in) :: dl
    type(ImVec2), intent(in) :: q
    logical, intent(in) :: round
    integer(c_int), intent(in) :: colw, colk

    type(ImVec2) :: q1, q2

    if (round) then
       call ImDrawList_AddCircleFilled(dl,q,handle_px,colw,12_c_int)
       call ImDrawList_AddCircle(dl,q,handle_px,colk,12_c_int,1._c_float)
    else
       q1%x = q%x - handle_px
       q1%y = q%y - handle_px
       q2%x = q%x + handle_px
       q2%y = q%y + handle_px
       call ImDrawList_AddRectFilled(dl,q1,q2,colw,0._c_float,0_c_int)
       call ImDrawList_AddRect(dl,q1,q2,colk,0._c_float,0_c_int,1._c_float)
    end if

  end subroutine draw_handle

end submodule objedit
