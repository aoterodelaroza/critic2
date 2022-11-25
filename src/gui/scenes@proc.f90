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

! Scene object and GL rendering utilities
submodule (scenes) proc
  implicit none

  ! some math parameters
  real(c_float), parameter :: zero = 0._c_float
  real(c_float), parameter :: one = 1._c_float
  real(c_float), parameter :: eye4(4,4) = reshape((/&
     one,zero,zero,zero,&
     zero,one,zero,zero,&
     zero,zero,one,zero,&
     zero,zero,zero,one/),shape(eye4))

  !xx! private procedures: low-level draws
  ! subroutine draw_sphere(x0,rad,rgba,ires)
  ! subroutine draw_cylinder(x1,x2,rad,rgba,ires)

contains

  !xx! scene

  !> Initialize a scene object associated with system isys.
  module subroutine scene_init(s,isys)
    use gui_main, only: nsys, sysc, sys_init, sys
    class(scene), intent(inout), target :: s
    integer, intent(in) :: isys

    if (isys < 1 .or. isys > nsys) return
    if (sysc(isys)%status /= sys_init) return

    ! basic variables
    s%id = isys
    s%isinit = .true.
    s%nc = 1
    s%scenerad = 10d0
    s%scenecenter = 0d0
    s%scenexmin = 0d0
    s%scenexmax = 1d0

    ! resolutions
    s%atom_res = 3
    s%bond_res = 1
    s%uc_res = 1

    ! phong default settings
    s%bgcolor = (/0._c_float,0._c_float,0._c_float,1._c_float/)
    s%lightpos = (/20._c_float,20._c_float,0._c_float/)
    s%lightcolor = (/1._c_float,1._c_float,1._c_float/)
    s%ambient = 0.2_c_float
    s%diffuse = 0.4_c_float
    s%specular = 0.6_c_float
    s%shininess = 8_c_int

    ! initialize representations
    if (allocated(s%rep)) deallocate(s%rep)
    allocate(s%rep(20))
    s%nrep = 0
    if (allocated(s%icount)) deallocate(s%icount)
    allocate(s%icount(0:reptype_NUM))
    s%icount = 0
    if (allocated(s%iord)) deallocate(s%iord)
    allocate(s%iord(20))
    s%iord = 0

    ! atoms
    s%nrep = s%nrep + 1
    call s%rep(s%nrep)%init(s%id,s%nrep,reptype_atoms)

    ! unit cell
    if (.not.sys(isys)%c%ismolecule) then
       s%nrep = s%nrep + 1
       call s%rep(s%nrep)%init(s%id,s%nrep,reptype_unitcell)
    end if

    ! reset the camera later
    s%camratio = 2.5_c_float
    s%forceresetcam = .true.

    ! sort the representations next pass
    s%forcesort = .true.

  end subroutine scene_init

  !> Terminate a scene object
  module subroutine scene_end(s)
    class(scene), intent(inout), target :: s

    s%isinit = .false.
    s%id = 0
    if (allocated(s%rep)) deallocate(s%rep)
    if (allocated(s%icount)) deallocate(s%icount)
    if (allocated(s%iord)) deallocate(s%iord)
    s%nrep = 0

  end subroutine scene_end

  !> Reset the camera position and direction. Sets scenerad,
  !> scenecenter, ortho_fov, persp_fov, v_center, v_up, v_pos, view,
  !> world, projection, and znear.
  module subroutine scene_reset(s)
    use gui_main, only: sys
    use param, only: pi
    class(scene), intent(inout), target :: s

    real(c_float) :: pic
    real*8 :: xmin(3), xmax(3), x(3), hside
    integer :: i, i1, i2, i3

    ! default transformation matrices
    pic = real(pi,c_float)
    s%ortho_fov = 45._c_float
    s%persp_fov = 45._c_float

    ! world matrix
    s%world = eye4

    ! camera distance and view matrix
    hside = 1.1_c_float * 0.5_c_float * max(s%scenexmax(1) - s%scenexmin(1),s%scenexmax(2) - s%scenexmin(2))
    hside = hside * s%camratio
    hside = max(hside,3._c_float)
    s%campos = s%scenecenter
    s%campos(3) = s%campos(3) + hside / tan(0.5_c_float * s%ortho_fov * pic / 180._c_float)
    s%camfront = (/zero,zero,-one/)
    s%camup = (/zero,one,zero/)
    call s%update_view_matrix()

    ! projection matrix
    s%znear = 0._c_float
    s%zfar = 100000._c_float
    call s%update_projection_matrix()

  end subroutine scene_reset

  !> Build the draw lists for the current scene.
  module subroutine scene_build_lists(s)
    use utils, only: translate
    class(scene), intent(inout), target :: s

    integer :: i
    real(c_float) :: xmin(3), xmax(3), maxrad, xc(3)

    ! initialize
    s%nsph = 0
    if (allocated(s%drawlist_sph)) deallocate(s%drawlist_sph)
    allocate(s%drawlist_sph(100))
    s%ncyl = 0
    if (allocated(s%drawlist_cyl)) deallocate(s%drawlist_cyl)
    allocate(s%drawlist_cyl(100))
    s%ncylflat = 0
    if (allocated(s%drawlist_cylflat)) deallocate(s%drawlist_cylflat)
    allocate(s%drawlist_cylflat(10))

    ! add the items by representation
    do i = 1, s%nrep
       call s%rep(i)%add_draw_elements(s%nc,s%nsph,s%drawlist_sph,s%ncyl,s%drawlist_cyl,&
          s%ncylflat,s%drawlist_cylflat)
    end do

    ! recalculate scene center and radius
    maxrad = 0._c_float
    do i = 1, s%nrep
       if (s%rep(i)%shown .and. s%rep(i)%type == reptype_atoms) then
          maxrad = max(maxrad,maxval(s%rep(i)%atom_style(1:s%rep(i)%natom_style)%rad))
       end if
    end do
    if (s%nsph + s%ncyl + s%ncylflat > 0) then
       do i = 1, 3
          xmin(i) = huge(1._c_float)
          xmax(i) = -huge(1._c_float)
          xmin(i) = minval(s%drawlist_sph(1:s%nsph)%x(i)) - maxrad
          xmin(i) = min(xmin(i),minval(s%drawlist_cyl(1:s%ncyl)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%drawlist_cyl(1:s%ncyl)%x2(i)))
          xmin(i) = min(xmin(i),minval(s%drawlist_cylflat(1:s%ncylflat)%x1(i)))
          xmin(i) = min(xmin(i),minval(s%drawlist_cylflat(1:s%ncylflat)%x2(i)))
          xmax(i) = maxval(s%drawlist_sph(1:s%nsph)%x(i)) + maxrad
          xmax(i) = max(xmax(i),maxval(s%drawlist_cyl(1:s%ncyl)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%drawlist_cyl(1:s%ncyl)%x2(i)))
          xmax(i) = max(xmax(i),maxval(s%drawlist_cylflat(1:s%ncylflat)%x1(i)))
          xmax(i) = max(xmax(i),maxval(s%drawlist_cylflat(1:s%ncylflat)%x2(i)))
       end do
    else
       xmin = 0._c_float
       xmax = 1._c_float
    end if

    ! new scene center and center shift in xc
    xc = s%scenecenter
    s%scenecenter = 0.5_c_float * (xmin + xmax)
    xc = s%scenecenter - xc

    ! radius and bounding box
    s%scenerad = 0.5_c_float * norm2(xmax - xmin)
    s%scenexmin = xmin
    s%scenexmax = xmax

    if (s%forceresetcam) then
       ! reset the camera if requested
       call s%reset()
       s%forceresetcam = .false.
    else
       ! translate the scene so the center position remains unchanged
       s%world = translate(s%world,-xc)
    end if

  end subroutine scene_build_lists

  !> Draw the scene
  module subroutine scene_render(s)
    use interfaces_opengl3
    use shapes, only: sphVAO, cylVAO
    use gui_main, only: sysc, sys_init
    use shaders, only: shader_phong, useshader, setuniform_int,&
       setuniform_float, setuniform_vec3, setuniform_vec4, setuniform_mat3,&
       setuniform_mat4
    class(scene), intent(inout), target :: s

    integer :: i

    ! check that the scene and system are initialized
    if (.not.s%isinit) return
    if (.not.sysc(s%id)%status == sys_init) return

    ! build draw lists if not done already
    if (.not.allocated(s%drawlist_sph)) call s%build_lists()

    ! set up the shader and the uniforms
    call useshader(shader_phong)
    call setuniform_int("uselighting",1_c_int)
    call setuniform_vec3("lightPos",s%lightpos)
    call setuniform_vec3("lightColor",s%lightcolor)
    call setuniform_float("ambient",s%ambient)
    call setuniform_float("diffuse",s%diffuse)
    call setuniform_float("specular",s%specular)
    call setuniform_int("shininess",s%shininess)
    call setuniform_mat4("world",s%world)
    call setuniform_mat4("view",s%view)
    call setuniform_mat4("projection",s%projection)

    ! draw the atoms
    call glBindVertexArray(sphVAO(s%atom_res))
    do i = 1, s%nsph
       call draw_sphere(s%drawlist_sph(i)%x,s%drawlist_sph(i)%r,s%drawlist_sph(i)%rgba,s%atom_res)
    end do

    ! draw the bonds
    call glBindVertexArray(cylVAO(s%bond_res))
    do i = 1, s%ncyl
       call draw_cylinder(s%drawlist_cyl(i)%x1,s%drawlist_cyl(i)%x2,s%drawlist_cyl(i)%r,&
          s%drawlist_cyl(i)%rgba,s%bond_res)
    end do

    ! draw the flat cylinders (unit cell)
    call setuniform_int("uselighting",0_c_int)
    call glBindVertexArray(cylVAO(s%uc_res))
    do i = 1, s%ncylflat
       call draw_cylinder(s%drawlist_cylflat(i)%x1,s%drawlist_cylflat(i)%x2,&
          s%drawlist_cylflat(i)%r,s%drawlist_cylflat(i)%rgba,s%uc_res)
    end do
    call glBindVertexArray(0)

  end subroutine scene_render

  !> Show the representation menu (called from view). Return .true.
  !> if the scene needs to be rendered again.
  module function representation_menu(s,idcaller) result(changed)
    use interfaces_cimgui
    use utils, only: iw_text, iw_tooltip
    use windows, only: win, stack_create_window, wintype_editrep, update_window_id
    use gui_main, only: ColorDangerButton, g, sysc
    use tools_io, only: string
    use tools, only: mergesort
    class(scene), intent(inout), target :: s
    integer(c_int), intent(in) :: idcaller
    logical :: changed

    integer :: i, ii, id, ll
    character(kind=c_char,len=:), allocatable, target :: str1, str2, str3
    logical :: discol, doerase
    type(ImVec2) :: szero
    character(kind=c_char,len=1024), target :: txtinp
    integer(c_int) :: flags
    integer, allocatable :: idx(:)

    logical, save :: ttshown = .false. ! tooltip flag

    ! coordinate this with draw_view in windows@view module
    integer(c_int), parameter :: ic_closebutton = 0
    integer(c_int), parameter :: ic_viewbutton = 1
    integer(c_int), parameter :: ic_name = 2
    integer(c_int), parameter :: ic_editbutton = 3

    ! initialization
    szero%x = 0
    szero%y = 0
    changed = .false.

    ! sort the representation array
    if (s%forcesort) then
       allocate(idx(s%nrep))
       do i = 1, s%nrep
          s%iord(i) = i
          idx(i) = s%rep(i)%iord
       end do
       call mergesort(idx,s%iord,1,s%nrep)
       deallocate(idx)
       s%forcesort = .false.
    end if

    ! representation rows
    do ii = 1, s%nrep
       i = s%iord(ii)
       if (.not.s%rep(i)%isinit) cycle

       ! update window ID
       call update_window_id(s%rep(i)%idwin)

       ! close button
       doerase = .false.
       call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float)
       if (igTableSetColumnIndex(ic_closebutton)) then
          str1 = "##2ic_closebutton" // string(ic_closebutton) // "," // string(i) // c_null_char
          if (my_CloseButton(c_loc(str1),ColorDangerButton)) doerase = .true.
       end if

       ! view button
       if (igTableSetColumnIndex(ic_viewbutton)) then
          str1 = "##2ic_viewbutton" // string(ic_viewbutton) // "," // string(i) // c_null_char
          if (igCheckbox(c_loc(str1), s%rep(i)%shown)) &
             changed = .true.
       end if

       ! name
       if (igTableSetColumnIndex(ic_name)) then
          discol = .not.s%rep(i)%shown
          if (discol) &
             call igPushStyleColor_Vec4(ImGuiCol_Text,g%Style%Colors(ImGuiCol_TextDisabled+1))
          str1 = trim(s%rep(i)%name) // "##" // string(ic_name) // "," // string(i) // c_null_char
          flags = ImGuiSelectableFlags_SpanAllColumns
          flags = ior(flags,ImGuiSelectableFlags_AllowItemOverlap)
          if (igSelectable_Bool(c_loc(str1),.false._c_bool,flags,szero)) then
             s%rep(i)%shown = .not.s%rep(i)%shown
             changed = .true.
          end if
          if (discol) call igPopStyleColor(1)

          ! name context menu
          if (igBeginPopupContextItem(c_loc(str1),ImGuiPopupFlags_MouseButtonRight)) then
             ! edit
             str2 = "Edit" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                if (s%rep(i)%idwin == 0) then
                   s%rep(i)%idwin = stack_create_window(wintype_editrep,.true.,isys=s%id,irep=i,idcaller=idcaller)
                else
                   call igSetWindowFocus_Str(c_loc(win(s%rep(i)%idwin)%name))
                end if
             end if
             call iw_tooltip("Edit this representation",ttshown)

             ! duplicate
             str2 = "Duplicate" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                id = sysc(s%id)%sc%get_new_representation_id()
                sysc(s%id)%sc%rep(id) = s%rep(i)
                sysc(s%id)%sc%rep(id)%name = trim(s%rep(i)%name) // " (copy)"
                s%icount(s%rep(i)%type) = s%icount(s%rep(i)%type) + 1
                s%icount(0) = s%icount(0) + 1
                sysc(s%id)%sc%rep(id)%iord = s%icount(0)
                s%forcesort = .true.
                changed = .true.
             end if
             call iw_tooltip("Make a copy of this representation",ttshown)

             ! show/hide
             str2 = "Show/Hide" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                s%rep(i)%shown = .not.s%rep(i)%shown
                changed = .true.
             end if
             call iw_tooltip("Toggle hide/show this representation",ttshown)

             ! rename
             str2 = "Rename" // c_null_char
             if (igBeginMenu(c_loc(str2),.true._c_bool)) then
                str3 = "##inputrenamerep" // c_null_char
                txtinp = trim(adjustl(s%rep(i)%name)) // c_null_char
                call igSetKeyboardFocusHere(0_c_int)
                if (igInputText(c_loc(str3),c_loc(txtinp),1023_c_size_t,ImGuiInputTextFlags_EnterReturnsTrue,&
                   c_null_ptr,c_null_ptr)) then
                   ll = index(txtinp,c_null_char)
                   s%rep(i)%name = txtinp(1:ll-1)
                   call igCloseCurrentPopup()
                end if
                call igEndMenu()
             end if
             call iw_tooltip("Rename this representation",ttshown)

             ! delete
             str2 = "Delete" // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) &
                doerase = .true.
             call iw_tooltip("Delete this representation",ttshown)

             call igEndPopup()
          end if
       end if

       ! edit button
       if (igTableSetColumnIndex(ic_editbutton)) then
          str1 = "Edit##2ic_editbutton" // string(ic_editbutton) // "," // string(i) // c_null_char
          if (igSmallButton(c_loc(str1))) then
             if (s%rep(i)%idwin == 0) then
                s%rep(i)%idwin = stack_create_window(wintype_editrep,.true.,isys=s%id,irep=i,idcaller=idcaller)
             else
                call igSetWindowFocus_Str(c_loc(win(s%rep(i)%idwin)%name))
             end if
          end if
       end if

       ! delete the representation if asked
       if (doerase) then
          call s%rep(i)%end()
          changed = .true.
       end if
    end do

  end function representation_menu

  !> Get the ID for a new representation. If necessary, reallocate the
  !> representations array.
  module function get_new_representation_id(s) result(id)
    use types, only: realloc
    class(scene), intent(inout), target :: s
    integer :: id

    integer :: i
    type(representation), allocatable :: aux(:)

    ! try to find an empty spot
    do i = 1, s%nrep
       if (.not.s%rep(i)%isinit) then
          id = i
          return
       end if
    end do

    ! make new representation at the end
    s%nrep = s%nrep + 1
    if (s%nrep > size(s%rep,1)) then
       allocate(aux(2*s%nrep))
       aux(1:size(s%rep,1)) = s%rep
       call move_alloc(aux,s%rep)
       call realloc(s%iord,2*s%nrep)
    end if
    id = s%nrep

  end function get_new_representation_id

  !> Update the projection matrix from the v_pos
  module subroutine update_projection_matrix(s)
    use utils, only: ortho, mult
    use param, only: pi
    class(scene), intent(inout), target :: s

    real(c_float) :: pic, hw2, sc(3)

    pic = real(pi,c_float)

    ! scene center: world to tworld
    sc = mult(s%world,s%scenecenter)

    ! update the projection matrix
    hw2 = tan(0.5_c_float * s%ortho_fov * pic / 180._c_float) * norm2(s%campos - sc)
    s%projection = ortho(-hw2,hw2,-hw2,hw2,s%znear,s%zfar)

  end subroutine update_projection_matrix

  !> Update the view matrix from the v_pos, v_front, and v_up
  module subroutine update_view_matrix(s)
    use utils, only: lookat
    class(scene), intent(inout), target :: s

    real(c_float) :: front(3), up(3)

    s%view = lookat(s%campos,s%campos+s%camfront,s%camup)

  end subroutine update_view_matrix

  !> Align the view with a given scene axis. a,b,c = 1,2,3 and x,y,z =
  !> -1,-2,-3.
  module subroutine align_view_axis(s,iaxis)
    use gui_main, only: sys
    use tools_math, only: cross
    use utils, only: rotate, translate
    class(scene), intent(inout), target :: s
    integer, intent(in) :: iaxis

    real*8 :: xaxis(3), oaxis(3), raxis(3)
    real(c_float) :: raxis_c(3), angle

    ! alignment axis
    if (iaxis == 1) then
       xaxis = sys(s%id)%c%m_x2c(:,1)
    elseif (iaxis == 2) then
       xaxis = sys(s%id)%c%m_x2c(:,2)
    elseif (iaxis == 3) then
       xaxis = sys(s%id)%c%m_x2c(:,3)
    elseif (iaxis == -1) then
       xaxis = (/1d0,0d0,0d0/)
    elseif (iaxis == -2) then
       xaxis = (/0d0,1d0,0d0/)
    elseif (iaxis == -3) then
       xaxis = (/0d0,0d0,1d0/)
    else
       return
    end if
    xaxis = xaxis / norm2(xaxis)

    oaxis = (/0d0,0d0,1d0/)
    raxis = cross(oaxis,xaxis)
    angle = real(asin(norm2(raxis)),c_float)

    ! reset the camera position
    call s%reset()

    ! set the world matrix
    if (angle > 1e-10_c_float) then
       raxis_c = real(raxis / norm2(raxis),c_float)
       s%world = translate(s%world,s%scenecenter)
       s%world = rotate(s%world,-angle,raxis_c)
       s%world = translate(s%world,-s%scenecenter)
    end if

  end subroutine align_view_axis

  !xx! representation

  !> Initialize a representation. If itype is present and not _none,
  !> fill the representation with the default values for the
  !> corresponding type and set isinit = .true. and shown = .true.
  module subroutine representation_init(r,isys,irep,itype)
    use gui_main, only: sys, sysc
    use tools_io, only: string
    class(representation), intent(inout), target :: r
    integer, intent(in) :: isys
    integer, intent(in) :: irep
    integer, intent(in), optional :: itype

    r%isinit = .false.
    r%shown = .false.
    r%type = reptype_none
    r%id = isys
    r%idrep = irep
    r%idwin = 0
    r%name = ""
    r%filter = ""
    r%goodfilter = .true.
    r%pertype = 1
    r%ncell = 1
    r%border = .true.
    r%onemotif = .false.
    r%atoms_display = .true.
    r%bonds_display = .true.
    r%atom_style_type = 0
    r%natom_style = 0
    r%atom_radii_reset_type = 0
    r%atom_radii_reset_scale = 0.7_c_float
    r%atom_color_reset_type = 0
    r%bond_color_style = 1
    r%bond_rgba = (/1._c_float,0._c_float,0._c_float,1._c_float/)
    r%bond_rad = 0.2_c_float
    r%uc_radius = 0.15_c_float
    r%uc_radiusinner = 0.15_c_float
    r%uc_innersteplen = 2d0
    r%uc_innerstipple = .true.
    r%uc_rgba = (/1._c_float,1._c_float,1._c_float,1._c_float/)
    r%uc_inner = .true.
    r%uc_coloraxes = .true.
    if (present(itype)) then
       if (itype == reptype_atoms) then
          r%isinit = .true.
          r%shown = .true.
          r%type = reptype_atoms
          r%name = "Atoms"
          if (sys(isys)%c%ismolecule) then
             r%ncell = 0
             r%border = .false.
             r%onemotif = .false.
          else
             r%border = .true.
             r%onemotif = sys(isys)%c%ismol3d
             r%ncell = 1
          end if
       elseif (itype == reptype_unitcell) then
          r%isinit = .true.
          r%shown = .true.
          r%type = reptype_unitcell
          r%name = "Unit cell"
       end if

       sysc(isys)%sc%icount(itype) = sysc(isys)%sc%icount(itype) + 1
       if (sysc(isys)%sc%icount(itype) > 1) then
          r%name = trim(r%name) // " " // string(sysc(isys)%sc%icount(itype))
       end if
    end if

    ! increment counter and force sort of the parent scene
    sysc(isys)%sc%icount(0) = sysc(isys)%sc%icount(0) + 1
    r%iord = sysc(isys)%sc%icount(0)
    sysc(isys)%sc%forcesort = .true.

    ! initialize the styles
    call r%reset_atom_style()

  end subroutine representation_init

  !> Terminate a representation
  module subroutine representation_end(r)
    use windows, only: win
    class(representation), intent(inout), target :: r

    r%name = ""
    r%filter = ""
    r%goodfilter = .true.
    r%isinit = .false.
    r%shown = .false.
    r%type = reptype_none
    r%id = 0
    r%idrep = 0
    r%iord = 0
    if (r%idwin > 0) &
       nullify(win(r%idwin)%rep)
    r%idwin = 0

  end subroutine representation_end

  !> Reset atom styles.
  module subroutine reset_atom_style(r)
    use gui_main, only: nsys, sysc, sys, sys_init
    use param, only: jmlcol, atmcov
    class(representation), intent(inout), target :: r

    integer :: i, ispc, iz

    ! set the atom style to zero
    r%natom_style = 0
    if (allocated(r%atom_style)) deallocate(r%atom_style)

    ! check the system is correct and initialized
    if (r%id < 1 .or. r%id > nsys) return
    if (sysc(r%id)%status /= sys_init) return

    ! fill the styles
    if (r%atom_style_type == 0) then
       ! species
       r%natom_style = sys(r%id)%c%nspc
       allocate(r%atom_style(r%natom_style))
       do i = 1, sys(r%id)%c%nspc
          r%atom_style(i)%shown = .true.

          iz = sys(r%id)%c%spc(i)%z
          r%atom_style(i)%rgba(1:3) = real(jmlcol(:,iz),c_float) / 255._c_float
          r%atom_style(i)%rgba(4) = 1._c_float

          r%atom_style(i)%rad = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    elseif (r%atom_style_type == 1) then
       ! nneq
       r%natom_style = sys(r%id)%c%nneq
       allocate(r%atom_style(r%natom_style))
       do i = 1, sys(r%id)%c%nneq
          r%atom_style(i)%shown = .true.

          ispc = sys(r%id)%c%at(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          r%atom_style(i)%rgba(1:3) = real(jmlcol(:,iz),c_float) / 255._c_float
          r%atom_style(i)%rgba(4) = 1._c_float

          r%atom_style(i)%rad = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    elseif (r%atom_style_type == 2) then
       ! ncel
       r%natom_style = sys(r%id)%c%ncel
       allocate(r%atom_style(r%natom_style))
       do i = 1, sys(r%id)%c%ncel
          r%atom_style(i)%shown = .true.

          ispc = sys(r%id)%c%atcel(i)%is
          iz = sys(r%id)%c%spc(ispc)%z
          r%atom_style(i)%rgba(1:3) = real(jmlcol(:,iz),c_float) / 255._c_float
          r%atom_style(i)%rgba(4) = 1._c_float

          r%atom_style(i)%rad = 0.7_c_float * real(atmcov(iz),c_float)
       end do
    end if

  end subroutine reset_atom_style

  !xx! private procedures: low-level draws

  !> Draw a sphere with center x0, radius rad and color rgba. Requires
  !> having the sphere VAO bound.
  subroutine draw_sphere(x0,rad,rgba,ires)
    use interfaces_opengl3
    use shaders, only: setuniform_vec4, setuniform_mat4
    use shapes, only: sphnel
    real(c_float), intent(in) :: x0(3)
    real(c_float), intent(in) :: rad
    real(c_float), intent(in) :: rgba(4)
    integer(c_int), intent(in) :: ires

    real(c_float) :: model(4,4)

    ! the model matrix: scale and translate
    model = eye4
    model(1:3,4) = x0
    model(1,1) = rad
    model(2,2) = rad
    model(3,3) = rad

    ! draw the sphere
    call setuniform_vec4("vColor",rgba)
    call setuniform_mat4("model",model)
    call glDrawElements(GL_TRIANGLES, int(3*sphnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)

  end subroutine draw_sphere

  !> Draw a cylinder from x1 to x2 with radius rad and color
  !> rgba. Requires having the cylinder VAO bound.
  subroutine draw_cylinder(x1,x2,rad,rgba,ires)
    use interfaces_opengl3
    use tools_math, only: cross_cfloat
    use shaders, only: setuniform_vec4, setuniform_mat4
    use shapes, only: cylnel
    real(c_float), intent(in) :: x1(3)
    real(c_float), intent(in) :: x2(3)
    real(c_float), intent(in) :: rad
    real(c_float), intent(in) :: rgba(4)
    integer(c_int), intent(in) :: ires

    real(c_float) :: xmid(3), xdif(3), up(3), crs(3), model(4,4), blen
    real(c_float) :: a, ca, sa, axis(3), temp(3)

    xmid = 0.5_c_float * (x1 + x2)
    xdif = x2 - x1
    blen = norm2(xdif)
    xdif = xdif / blen
    up = (/0._c_float,0._c_float,1._c_float/)
    crs = cross_cfloat(up,xdif)

    ! the model matrix
    model = eye4
    model(1:3,4) = xmid ! translate(m_model,xmid);
    if (dot_product(crs,crs) > 1e-14_c_float) then
       ! m_model = m_model * rotate(acos(dot(xdif,up)),crs);
       a = acos(dot_product(xdif,up))
       ca = cos(a)
       sa = sin(a)
       axis = crs / norm2(crs)
       temp = (1._c_float - ca) * axis

       model(1,1) = ca + temp(1) * axis(1)
       model(2,1) = temp(1) * axis(2) + sa * axis(3)
       model(3,1) = temp(1) * axis(3) - sa * axis(2)

       model(1,2) = temp(2) * axis(1) - sa * axis(3)
       model(2,2) = ca + temp(2) * axis(2)
       model(3,2) = temp(2) * axis(3) + sa * axis(1)

       model(1,3) = temp(3) * axis(1) + sa * axis(2)
       model(2,3) = temp(3) * axis(2) - sa * axis(1)
       model(3,3) = ca + temp(3) * axis(3)
    end if
    ! m_model = scale(m_model,glm::vec3(rad,rad,blen));
    model(:,1) = model(:,1) * rad
    model(:,2) = model(:,2) * rad
    model(:,3) = model(:,3) * blen

    ! draw the cylinder
    call setuniform_vec4("vColor",rgba)
    call setuniform_mat4("model",model)
    call glDrawElements(GL_TRIANGLES, int(3*cylnel(ires),c_int), GL_UNSIGNED_INT, c_null_ptr)

  end subroutine draw_cylinder

  !> Add the spheres, cylinder, etc. to the draw lists.
  module subroutine add_draw_elements(r,nc,nsph,drawlist_sph,ncyl,drawlist_cyl,ncylflat,drawlist_cylflat)
    use gui_main, only: sys
    use tools_io, only: string
    use hashmod, only: hash
    class(representation), intent(inout), target :: r
    integer, intent(in) :: nc(3)
    integer, intent(inout) :: nsph
    type(dl_sphere), intent(inout), allocatable :: drawlist_sph(:)
    integer, intent(inout) :: ncyl
    type(dl_cylinder), intent(inout), allocatable :: drawlist_cyl(:)
    integer, intent(inout) :: ncylflat
    type(dl_cylinder), intent(inout), allocatable :: drawlist_cylflat(:)

    type(hash) :: shown_atoms
    logical :: havefilter, step, ok, isedge(3)
    integer :: n(3), i, j, k, imol, lvec(3), id, idaux, n0(3), n1(3), i1, i2, i3, ix(3)
    integer :: ib, ineigh, ixn(3), ix1(3), ix2(3), nstep
    real(c_float) :: rgba(4), rad
    real*8 :: xx(3), x0(3), x1(3), x2(3), res
    type(dl_sphere), allocatable :: auxsph(:)
    type(dl_cylinder), allocatable :: auxcyl(:)
    character(len=:), allocatable :: atcode

    real*8, parameter :: rthr = 0.01d0
    real*8, parameter :: rthr1 = 1-rthr
    integer, parameter :: uc(3,2,12) = reshape((/&
       0,0,0,  1,0,0,&
       0,0,0,  0,1,0,&
       0,0,0,  0,0,1,&
       1,1,0,  1,1,1,&
       1,0,1,  1,1,1,&
       0,1,1,  1,1,1,&
       1,0,0,  1,1,0,&
       0,1,0,  1,1,0,&
       0,1,0,  0,1,1,&
       0,0,1,  0,1,1,&
       0,0,1,  1,0,1,&
       1,0,0,  1,0,1/),shape(uc))
    integer, parameter :: ucdir(12) = (/1, 2, 3, 3, 2, 1, 2, 1, 3, 2, 1, 3/)

    ! return if not initialized
    if (.not.r%isinit) return
    if (.not.r%shown) return

    ! initialize the drawlists if not done already
    if (.not.allocated(drawlist_sph)) then
       allocate(drawlist_sph(100))
       nsph = 0
    end if
    if (.not.allocated(drawlist_cyl)) then
       allocate(drawlist_cyl(100))
       ncyl = 0
    end if
    if (.not.allocated(drawlist_cylflat)) then
       allocate(drawlist_cylflat(100))
       ncylflat = 0
    end if

    if (r%type == reptype_atoms) then
       !!! atoms and bonds representation !!!

       !! first, the atoms
       ! do we have a filter?
       havefilter = (len_trim(r%filter) > 0) .and. r%goodfilter

       ! calculate the periodicity
       n = 1
       if (r%pertype == 1) then
          n = nc
       elseif (r%pertype == 2) then
          n = r%ncell
       end if

       ! run over atoms, either directly or per-molecule
       i = 0
       imol = 0
       do while(i < sys(r%id)%c%ncel)
          if (r%onemotif .and. sys(r%id)%c%ismol3d) then
             step = (imol == 0)
             if (.not.step) step = (k == sys(r%id)%c%mol(imol)%nat)
             if (step) then
                imol = imol + 1
                k = 0
             end if
             if (imol > sys(r%id)%c%nmol) exit
             k = k + 1
             i = sys(r%id)%c%mol(imol)%at(k)%cidx
             lvec = sys(r%id)%c%mol(imol)%at(k)%lvec
          else
             i = i + 1
             lvec = 0
          end if

          ! skip hidden atoms
          if (r%atom_style_type == 0) then
             id = sys(r%id)%c%atcel(i)%is
          elseif (r%atom_style_type == 1) then
             id = sys(r%id)%c%atcel(i)%idx
          else
             id = i
          end if
          if (.not.r%atom_style(id)%shown) cycle

          ! calculate the border
          n0 = 0
          n1 = n-1
          if (r%border.and..not.r%onemotif) then
             xx = sys(r%id)%c%atcel(i)%x
             do j = 1, 3
                if (xx(j) < rthr) then
                   n1(j) = n(j)
                elseif (xx(j) > rthr1) then
                   n0(j) = -1
                end if
             end do
          end if

          ! draw the spheres and cylinders
          rgba = r%atom_style(id)%rgba
          rad = r%atom_style(id)%rad
          do i1 = n0(1), n1(1)
             do i2 = n0(2), n1(2)
                do i3 = n0(3), n1(3)
                   ! atoms
                   ix = (/i1,i2,i3/) + lvec
                   xx = sys(r%id)%c%atcel(i)%x + ix
                   xx = sys(r%id)%c%x2c(xx)

                   ! apply the filter
                   if (havefilter) then
                      res = sys(r%id)%eval(r%filter,.false.,ok,xx)
                      if (ok) then
                         if (res == 0d0) cycle
                      end if
                   end if

                   ! draw the atom, reallocate if necessary
                   if (r%atoms_display) then
                      nsph = nsph + 1
                      if (nsph > size(drawlist_sph,1)) then
                         allocate(auxsph(2*nsph))
                         auxsph(1:size(drawlist_sph,1)) = drawlist_sph
                         call move_alloc(auxsph,drawlist_sph)
                      end if

                      ! write down the sphere
                      drawlist_sph(nsph)%x = real(xx,c_float)
                      drawlist_sph(nsph)%r = rad
                      drawlist_sph(nsph)%rgba = rgba
                   end if

                   ! start with bonds
                   if (.not.r%bonds_display) cycle

                   ! add this atom to the hash
                   atcode = string(i) // "_" // string(ix(1)) // "_" // string(ix(2)) // "_" // string(ix(3))
                   call shown_atoms%put(atcode,1)

                   ! bonds
                   do ib = 1, sys(r%id)%c%nstar(i)%ncon
                      ineigh = sys(r%id)%c%nstar(i)%idcon(ib)
                      ixn = ix + sys(r%id)%c%nstar(i)%lcon(:,ib)

                      ! check if the atom has been represented already
                      atcode = string(ineigh) // "_" // string(ixn(1)) // "_" // string(ixn(2)) // "_" // string(ixn(3))
                      if (.not.shown_atoms%iskey(atcode)) cycle

                      ! draw the bond, reallocate if necessary
                      if (r%bond_color_style == 0) then
                         ncyl = ncyl + 1
                      else
                         ncyl = ncyl + 2
                      end if
                      if (ncyl > size(drawlist_cyl,1)) then
                         allocate(auxcyl(2*ncyl))
                         auxcyl(1:size(drawlist_cyl,1)) = drawlist_cyl
                         call move_alloc(auxcyl,drawlist_cyl)
                      end if

                      x1 = xx
                      x2 = sys(r%id)%c%atcel(ineigh)%x + ixn
                      x2 = sys(r%id)%c%x2c(x2)
                      if (r%bond_color_style == 0) then
                         drawlist_cyl(ncyl)%x1 = real(x1,c_float)
                         drawlist_cyl(ncyl)%x2 = real(x2,c_float)
                         drawlist_cyl(ncyl)%r = r%bond_rad
                         drawlist_cyl(ncyl)%rgba = r%bond_rgba
                      else
                         x0 = 0.5d0 * (x1 + x2)
                         drawlist_cyl(ncyl-1)%x1 = real(x1,c_float)
                         drawlist_cyl(ncyl-1)%x2 = real(x0,c_float)
                         drawlist_cyl(ncyl-1)%r = r%bond_rad
                         drawlist_cyl(ncyl-1)%rgba = r%atom_style(id)%rgba

                         if (r%atom_style_type == 0) then
                            idaux = sys(r%id)%c%atcel(ineigh)%is
                         elseif (r%atom_style_type == 1) then
                            idaux = sys(r%id)%c%atcel(ineigh)%idx
                         else
                            idaux = ineigh
                         end if
                         drawlist_cyl(ncyl)%x1 = real(x0,c_float)
                         drawlist_cyl(ncyl)%x2 = real(x2,c_float)
                         drawlist_cyl(ncyl)%r = r%bond_rad
                         drawlist_cyl(ncyl)%rgba = r%atom_style(idaux)%rgba
                      end if
                   end do ! ncon

                end do ! i3
             end do ! i2
          end do ! i1
       end do ! loop over complete atom list (i)
    elseif (r%type == reptype_unitcell) then
       !!! unit cell representation !!!

       ! number of cells
       n = 1
       if (r%pertype == 1) then
          n = nc
       elseif (r%pertype == 2) then
          n = r%ncell
       end if

       ! external cell
       do i = 1, 12
          x1 = real(uc(:,1,i) * n,8)
          x1 = sys(r%id)%c%x2c(x1)
          x2 = real(uc(:,2,i) * n,8)
          x2 = sys(r%id)%c%x2c(x2)

          call increase_ncylflat()
          drawlist_cylflat(ncylflat)%x1 = real(x1,c_float)
          drawlist_cylflat(ncylflat)%x2 = real(x2,c_float)
          drawlist_cylflat(ncylflat)%r = r%uc_radius
          if (r%uc_coloraxes.and.i==1) then
             drawlist_cylflat(ncylflat)%rgba = (/1._c_float,0._c_float,0._c_float,1._c_float/)
          elseif (r%uc_coloraxes.and.i==2) then
             drawlist_cylflat(ncylflat)%rgba = (/0._c_float,1._c_float,0._c_float,1._c_float/)
          elseif (r%uc_coloraxes.and.i==3) then
             drawlist_cylflat(ncylflat)%rgba = (/0._c_float,0._c_float,1._c_float,1._c_float/)
          else
             drawlist_cylflat(ncylflat)%rgba = r%uc_rgba
          end if
       end do

       ! draw inner cylinders
       if (r%uc_inner) then
          do i1 = 0, n(1)-1
             do i2 = 0, n(2)-1
                do i3 = 0, n(3)-1
                   do i = 1, 12
                      ix1 = (/i1,i2,i3/) + uc(:,1,i)
                      ix2 = (/i1,i2,i3/) + uc(:,2,i)

                      ! skip outer cylinders
                      isedge = (ix1 == 0 .and. ix2 == 0) .or. (ix1 == n .and. ix2 == n)
                      isedge(ucdir(i)) = .true.
                      if (all(isedge)) cycle

                      x1 = real(ix1,8)
                      x1 = sys(r%id)%c%x2c(x1)
                      x2 = real(ix2,8)
                      x2 = sys(r%id)%c%x2c(x2)

                      ! logical :: uc_innerstipple ! stippled lines for the inner lines
                      if (r%uc_innerstipple) then
                         nstep = ceiling(norm2(x2 - x1) / r%uc_innersteplen)
                         do j = 1, nstep
                            call increase_ncylflat()
                            drawlist_cylflat(ncylflat)%x1 = real(x1 + real(2*j-1,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            drawlist_cylflat(ncylflat)%x2 = real(x1 + real(2*j,8)/real(2*nstep,8) * (x2-x1) ,c_float)
                            drawlist_cylflat(ncylflat)%r = r%uc_radiusinner
                            drawlist_cylflat(ncylflat)%rgba = r%uc_rgba
                         end do
                      else
                         call increase_ncylflat()
                         drawlist_cylflat(ncylflat)%x1 = real(x1 ,c_float)
                         drawlist_cylflat(ncylflat)%x2 = real(x2 ,c_float)
                         drawlist_cylflat(ncylflat)%r = r%uc_radiusinner
                         drawlist_cylflat(ncylflat)%rgba = r%uc_rgba
                      end if
                   end do
                end do
             end do
          end do
       end if
    end if ! reptype
  contains
    subroutine increase_ncylflat()

      ncylflat = ncylflat + 1
      if (ncylflat > size(drawlist_cylflat,1)) then
         allocate(auxcyl(2*ncylflat))
         auxcyl(1:size(drawlist_cylflat,1)) = drawlist_cylflat
         call move_alloc(auxcyl,drawlist_cylflat)
      end if

    end subroutine increase_ncylflat
  end subroutine add_draw_elements

end submodule proc
