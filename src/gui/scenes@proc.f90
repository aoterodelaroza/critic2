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

  real(c_float) :: min_scenerad = 10._c_float

  ! object resolution
  integer, parameter :: icylres = 3 ! cylinder

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
  ! subroutine draw_cylinder(x1,x2,rad,rgba)

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

    ! phong default settings
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

    ! bonds
    s%nrep = s%nrep + 1
    call s%rep(s%nrep)%init(s%id,s%nrep,reptype_bonds)

    ! unit cell
    if (.not.sys(isys)%c%ismolecule) then
       s%nrep = s%nrep + 1
       call s%rep(s%nrep)%init(s%id,s%nrep,reptype_unitcell)
    end if

    ! reset the camera
    call s%reset()

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
    real*8 :: xmin(3), xmax(3), x(3)
    integer :: isys
    integer :: i, i1, i2, i3

    ! initialize
    isys = s%id

    ! use xmin and xmax
    ! ! draw all initialized and visible representations
    ! do i = 1, s%nrep
    !    if (s%rep(i)%isinit .and. s%rep(i)%shown) then
    !       call s%rep(i)%draw()
    !    end if
    ! end do

    ! calculate the initial scene radius
    xmin = 0d0
    xmax = 0d0
    if (sys(isys)%c%ncel > 0) then
       xmin = sys(isys)%c%atcel(1)%r
       xmax = sys(isys)%c%atcel(1)%r
       do i = 2, sys(isys)%c%ncel
          xmax = max(sys(isys)%c%atcel(i)%r,xmax)
          xmin = min(sys(isys)%c%atcel(i)%r,xmin)
       end do

       if (.not.sys(isys)%c%ismolecule) then
          do i1 = 0, 1
             do i2 = 0, 1
                do i3 = 0, 1
                   x = real((/i1,i2,i3/),8)
                   x = sys(isys)%c%x2c(x)
                   xmax = max(x,xmax)
                   xmin = min(x,xmin)
                end do
             end do
          end do
       end if
    end if
    s%scenerad = real(norm2(xmax-xmin),c_float)
    s%scenerad = max(s%scenerad,min_scenerad)

    ! calculate the scene center
    x = (/0.5d0,0.5d0,0.5d0/)
    x = sys(isys)%c%x2c(x)
    s%scenecenter = real(x,c_float)

    ! default transformation matrices
    pic = real(pi,c_float)
    s%ortho_fov = 45._c_float
    s%persp_fov = 45._c_float

    ! world matrix
    s%world = eye4

    ! camera distance and view matrix
    s%campos = s%scenecenter
    ! s%campos(3) = s%campos(3) + 1.1_c_float * s%scenerad * &
    !    sqrt(real(maxval(s%ncell),c_float)) / tan(0.5_c_float * s%ortho_fov * pic / 180._c_float)
    s%campos(3) = s%campos(3) + 1.1_c_float * s%scenerad / tan(0.5_c_float * s%ortho_fov * pic / 180._c_float)
    call s%update_view_matrix()

    ! projection matrix
    s%znear = 0.1_c_float
    s%zfar = 100000._c_float
    call s%update_projection_matrix()

  end subroutine scene_reset

  !> Draw the scene
  module subroutine scene_render(s)
    use interfaces_opengl3
    use gui_main, only: sysc, sys_init, sys
    use shaders, only: shader_phong, useshader, setuniform_int,&
       setuniform_float, setuniform_vec3, setuniform_vec4, setuniform_mat3,&
       setuniform_mat4
    class(scene), intent(inout), target :: s

    integer :: i

    ! check that the scene and system are initialized
    if (.not.s%isinit) return
    if (.not.sysc(s%id)%status == sys_init) return

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

    ! draw all initialized and visible representations
    do i = 1, s%nrep
       if (s%rep(i)%isinit .and. s%rep(i)%shown) then
          call s%rep(i)%draw()
       end if
    end do

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
    logical(c_bool) :: ldum
    logical :: discol, doerase
    type(ImVec2) :: szero, sz
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

             ! show/hide
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
          str1 = "E##2ic_editbutton" // string(ic_editbutton) // "," // string(i) // c_null_char
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
    use utils, only: ortho
    use param, only: pi
    class(scene), intent(inout), target :: s

    real(c_float) :: pic, hw2

    pic = real(pi,c_float)
    hw2 = tan(0.5_c_float * s%ortho_fov * pic / 180._c_float) * s%campos(3)
    s%projection = ortho(-hw2,hw2,-hw2,hw2,s%znear,s%zfar)

  end subroutine update_projection_matrix

  !> Update the view matrix from the v_pos, v_front, and v_up
  module subroutine update_view_matrix(s)
    use utils, only: lookat
    class(scene), intent(inout), target :: s

    real(c_float) :: front(3), up(3)

    front = (/zero,zero,-one/)
    up = (/zero,one,zero/)
    s%view = lookat(s%campos,s%campos+front,up)

  end subroutine update_view_matrix

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
    r%atom_style_type = 0
    r%natom_style = 0
    r%atom_radii_reset_type = 0
    r%atom_radii_reset_scale = 0.7_c_float
    r%atom_color_reset_type = 0
    r%atom_res = 3
    r%bond_scale = 1d0
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
       elseif (itype == reptype_bonds) then
          r%isinit = .true.
          r%shown = .true.
          r%type = reptype_bonds
          r%name = "Bonds"
          if (sys(isys)%c%ismolecule) then
             r%ncell = 0
             r%border = .false.
             r%onemotif = .false.
          else
             r%border = .true.
             r%onemotif = sys(isys)%c%ismol3d
             r%ncell = 1
          end if
          r%bond_scale = 1d0
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

  !> Draw a representation. If xmin and xmax are present, calculate
  !> the bounding box instead of drawing.
  module subroutine representation_draw(r,xmin,xmax)
    class(representation), intent(inout), target :: r
    real*8, optional, intent(inout) :: xmin(3)
    real*8, optional, intent(inout) :: xmax(3)

    ! initial checks
    if (.not.r%isinit) return
    if (.not.r%shown) return

    ! draw the objects
    if (r%type == reptype_atoms) then
       call r%draw_atoms(xmin,xmax)
    elseif (r%type == reptype_bonds) then
       call r%draw_bonds(xmin,xmax)
    elseif (r%type == reptype_unitcell) then
       call r%draw_unitcell(xmin,xmax)
    elseif (r%type == reptype_labels) then
       call r%draw_labels(xmin,xmax)
    end if

  end subroutine representation_draw

  !> Draw an atoms representation. If xmin and xmax are present,
  !> calculate the bounding box instead of drawing.
  module subroutine draw_atoms(r,xmin,xmax)
    use interfaces_opengl3
    use gui_main, only: sys
    use shapes, only: sphVAO
    class(representation), intent(inout), target :: r
    real*8, optional, intent(inout) :: xmin(3)
    real*8, optional, intent(inout) :: xmax(3)

    real(c_float) :: rgba(4), x0(3), x1(3), rad
    integer :: i, id
    logical :: havefilter, ok
    real*8 :: res

    havefilter = (len_trim(r%filter) > 0) .and. r%goodfilter

    call glBindVertexArray(sphVAO(r%atom_res))
    do i = 1, sys(r%id)%c%ncel
       ! apply the filter
       if (havefilter) then
          res = sys(r%id)%eval(r%filter,.false.,ok,sys(r%id)%c%atcel(i)%r)
          if (ok) then
             if (res == 0d0) cycle
          end if
       end if
       if (r%atom_style_type == 0) then
          id = sys(r%id)%c%atcel(i)%is
       elseif (r%atom_style_type == 1) then
          id = sys(r%id)%c%atcel(i)%idx
       else
          id = i
       end if

       if (.not.r%atom_style(id)%shown) cycle
       rgba = r%atom_style(id)%rgba
       rad = r%atom_style(id)%rad

       x0 = real(sys(r%id)%c%atcel(i)%r,c_float)
       call draw_sphere(x0,rad,rgba,r%atom_res)
    end do
    call glBindVertexArray(0)

  end subroutine draw_atoms

  !> Draw a bonds representation. If xmin and xmax are present,
  !> calculate the bounding box instead of drawing.
  module subroutine draw_bonds(r,xmin,xmax)
    class(representation), intent(inout), target :: r
    real*8, optional, intent(inout) :: xmin(3)
    real*8, optional, intent(inout) :: xmax(3)


  end subroutine draw_bonds

  !> Draw a unit cell representation. If xmin and xmax are present,
  !> calculate the bounding box instead of drawing.
  module subroutine draw_unitcell(r,xmin,xmax)
    class(representation), intent(inout), target :: r
    real*8, optional, intent(inout) :: xmin(3)
    real*8, optional, intent(inout) :: xmax(3)

    ! use shapes, only: cylVAO

    ! real*8, parameter :: uc(3,2,12) = reshape((/&
    !    0._c_float,0._c_float,0._c_float,  1._c_float,0._c_float,0._c_float,&
    !    0._c_float,0._c_float,0._c_float,  0._c_float,1._c_float,0._c_float,&
    !    0._c_float,0._c_float,0._c_float,  0._c_float,0._c_float,1._c_float,&
    !    1._c_float,1._c_float,1._c_float,  1._c_float,1._c_float,0._c_float,&
    !    1._c_float,1._c_float,1._c_float,  1._c_float,0._c_float,1._c_float,&
    !    1._c_float,1._c_float,1._c_float,  0._c_float,1._c_float,1._c_float,&
    !    1._c_float,0._c_float,0._c_float,  1._c_float,1._c_float,0._c_float,&
    !    0._c_float,1._c_float,0._c_float,  1._c_float,1._c_float,0._c_float,&
    !    0._c_float,1._c_float,0._c_float,  0._c_float,1._c_float,1._c_float,&
    !    0._c_float,0._c_float,1._c_float,  0._c_float,1._c_float,1._c_float,&
    !    0._c_float,0._c_float,1._c_float,  1._c_float,0._c_float,1._c_float,&
    !    1._c_float,0._c_float,0._c_float,  1._c_float,0._c_float,1._c_float/),shape(uc))

    ! ! ! draw the unit cell
    ! ! if (s%showcell) then
    ! !    call glBindVertexArray(cylVAO(icylres))
    ! !    rad = 0.2_c_float
    ! !    rgba = (/1._c_float,0._c_float,0._c_float,1._c_float/)
    ! !    do i = 1, 12
    ! !       x0 = real(sys(s%id)%c%x2c(uc(:,1,i)),c_float)
    ! !       x1 = real(sys(s%id)%c%x2c(uc(:,2,i)),c_float)
    ! !       call s%draw_cylinder(x0,x1,rad,rgba)
    ! !    end do
    ! !    call glBindVertexArray(0)
    ! ! end if

  end subroutine draw_unitcell

  !> Draw a labels representation. If xmin and xmax are present,
  !> calculate the bounding box instead of drawing.
  module subroutine draw_labels(r,xmin,xmax)
    class(representation), intent(inout), target :: r
    real*8, optional, intent(inout) :: xmin(3)
    real*8, optional, intent(inout) :: xmax(3)


  end subroutine draw_labels

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
  subroutine draw_cylinder(x1,x2,rad,rgba)
    use interfaces_opengl3
    use tools_math, only: cross_cfloat
    use shaders, only: setuniform_vec4, setuniform_mat4
    use shapes, only: cylnel
    real(c_float), intent(in) :: x1(3)
    real(c_float), intent(in) :: x2(3)
    real(c_float), intent(in) :: rad
    real(c_float), intent(in) :: rgba(4)

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
    call glDrawElements(GL_TRIANGLES, int(3*cylnel(icylres),c_int), GL_UNSIGNED_INT, c_null_ptr)

  end subroutine draw_cylinder

end submodule proc
