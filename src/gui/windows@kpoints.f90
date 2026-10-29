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


! Routines for the k-point grid window (the KPOINTS keyword).
submodule (windows) kpoints
  use interfaces_cimgui
  use param, only: bohrtoa
  implicit none

contains

  !> Draw the k-point grid window: the dimensions of the uniform
  !> k-point grid that a length parameter (rklength) gives for the
  !> system shown by the anchor view, and the table of rklength ranges
  !> that give the same grid, as KPOINTS calculates them.
  module subroutine draw_kpoints(w)
    use systems, only: sys, sysc, sys_init, ok_system
    use utils, only: iw_text, iw_button, iw_tooltip, iw_dragfloat_real8, iw_close_event,&
       iw_setpos_bottomright, iw_table_column, iw_table_headers_row, iw_calcheight,&
       iw_separatortext, iw_highlight_selectable, iw_cell_right
    use crystalmod, only: kpoints_rkstep, kpoints_rkmax_def
    use tools_io, only: string
    class(window), intent(inout), target :: w

    integer, parameter :: itable_ranges = 1 ! the table, for table_write_*

    logical :: doquit, goodsys, iscrystal, syschanged, ldum, clicked
    integer :: i, j, isys, iview, iwrite, nkcur(3), nrange
    real*8, allocatable :: rk1(:), rk2(:)
    integer, allocatable :: nk(:,:)
    integer(c_int) :: tflags
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
    iscrystal = goodsys
    if (iscrystal) iscrystal = .not.sys(isys)%c%ismolecule

    ! initialize the form on the first pass
    if (w%firstpass) then
       w%kp = kpoints_state()
       w%kp%rkmax = kpoints_rkmax_def
       w%errmsg = ""
    end if

    ! the system
    if (goodsys) then
       call iw_text("System",highlight=.true.)
       call iw_text(string(isys) // ": " // trim(sysc(isys)%seed%name),sameline=.true.)
    end if

    ! a k-point grid needs a periodic cell
    if (goodsys .and. .not.iscrystal) &
       call iw_text("K-point grids can only be calculated for crystals",danger=.true.,wrap=.true.)

    if (iscrystal) then
       ! the grid for the current rklength
       call sys(isys)%c%get_kpoints(w%kp%rk,nkcur)
       call sys(isys)%c%get_kpoints_ranges(w%kp%rkmax,rk1,rk2,nk)
       nrange = size(rk1,1)
       call iw_separatortext("Grid",w%kp%secgrid)
       if (w%kp%secgrid%open) then
          call iw_text("Rk (bohr)",alignframe=.true.)
          ldum = iw_dragfloat_real8("##kpointsrk",x1=w%kp%rk,speed=0.1d0,min=kpoints_rkstep,max=1d3,&
             decimal=1,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.)
          call iw_tooltip("Length parameter: the grid has max(1,nint(Rk·|ai*|)) points along&
             & each reciprocal lattice vector ai*, so the spacing between k-points is about&
             & 1/Rk (no 2π)",ttshown)

          call iw_text("Grid")
          call iw_text(string(nkcur(1)) // " × " // string(nkcur(2)) // " × " // string(nkcur(3)),&
             sameline=.true.)
          call iw_tooltip("Number of k-points along each reciprocal lattice vector",ttshown)
          call iw_text("Total",gap=.true.)
          call iw_text(string(product(nkcur)),sameline=.true.)
          call iw_tooltip("Number of k-points in the grid",ttshown)

          call iw_text("Spacing (1/Å)")
          s = ""
          do i = 1, 3
             s = s // string(sys(isys)%c%ar(i) / bohrtoa / nkcur(i),'f',decimal=4) // "  "
          end do
          call iw_text(trim(s),sameline=.true.)
          call iw_tooltip("Distance between neighboring k-points along each reciprocal lattice&
             & vector, |ai*|/ni (no 2π)",ttshown)

          call iw_text("Supercell (Å)")
          s = ""
          do i = 1, 3
             s = s // string(nkcur(i) * sys(isys)%c%aa(i) * bohrtoa,'f',decimal=3) // "  "
          end do
          call iw_text(trim(s),sameline=.true.)
          call iw_tooltip("Lengths of the supercell that corresponds to the grid (the&
             & Born-von Kármán cell), ni·|ai|. In an orthogonal cell they are close to Rk",ttshown)
       end if

       ! the table of rklength ranges
       call iw_separatortext("Rk ranges (" // string(nrange) // ")",w%kp%secranges,hasmore=.true.)
       if (w%kp%secranges%open) then
          if (w%kp%secranges%more) then
             call iw_text("Rk max (bohr)",alignframe=.true.)
             ldum = iw_dragfloat_real8("##kpointsrkmax",x1=w%kp%rkmax,speed=1d0,min=1d0,max=1d3,&
                decimal=1,flags=ImGuiSliderFlags_AlwaysClamp,sameline=.true.,notlive=.true.)
             call iw_tooltip("The table covers the length parameters up to this value",ttshown)
          end if

          tflags = ImGuiTableFlags_None
          tflags = ior(tflags,ImGuiTableFlags_NoSavedSettings)
          tflags = ior(tflags,ImGuiTableFlags_RowBg)
          tflags = ior(tflags,ImGuiTableFlags_Borders)
          tflags = ior(tflags,ImGuiTableFlags_SizingFixedFit)
          tflags = ior(tflags,ImGuiTableFlags_ScrollY)
          str1 = "##kpointstable" // c_null_char
          ! the rest of the window, but at least a few rows (the stretch
          ! beyond them is elastic)
          call igGetContentRegionAvail(sz0)
          sz0%x = 0._c_float
          sz0%y = max(sz0%y,iw_calcheight(5,0,.false.))
          w%heightslack = sz0%y - iw_calcheight(5,0,.false.)
          if (igBeginTable(c_loc(str1),6,tflags,sz0,0._c_float)) then
             call iw_table_column("Rk from",id=0_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Rk to",id=1_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("n1",id=2_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("n2",id=3_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("n3",id=4_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_column("Total",id=5_c_int,flags=ImGuiTableColumnFlags_WidthFixed)
             call iw_table_headers_row(freezetop=.true.,writemenu=iwrite,ttshown=ttshown)

             clipper = ImGuiListClipper_ImGuiListClipper()
             call ImGuiListClipper_Begin(clipper,nrange,igGetTextLineHeightWithSpacing())
             call w%table_write_begin(itable_ranges,iwrite,clipper)
             do while (ImGuiListClipper_Step(clipper))
                call c_f_pointer(clipper,clipper_f)
                do i = clipper_f%DisplayStart+1, clipper_f%DisplayEnd
                   call igTableNextRow(ImGuiTableRowFlags_None,0._c_float)
                   if (igTableSetColumnIndex(0_c_int)) &
                      call iw_cell_right(string(rk1(i),'f',decimal=1))
                   ! selectable spans the row (outside the column block, so it
                   ! survives column 0 being clipped); clicking it picks its grid
                   ldum = iw_highlight_selectable("##kpointsrow" // string(i),clicked=clicked,&
                      selected=all(nk(:,i) == nkcur))
                   if (clicked) w%kp%rk = rk1(i)
                   if (igTableSetColumnIndex(1_c_int)) &
                      call iw_cell_right(string(rk2(i),'f',decimal=1))
                   do j = 1, 3
                      if (igTableSetColumnIndex(int(j+1,c_int))) &
                         call iw_cell_right(string(nk(j,i)))
                   end do
                   if (igTableSetColumnIndex(5_c_int)) &
                      call iw_cell_right(string(product(nk(:,i))))
                end do
             end do
             call ImGuiListClipper_End(clipper)
             call ImGuiListClipper_destroy(clipper)
             call w%table_write_end(itable_ranges,table_title("Rk ranges (bohr)",isys))
             call igEndTable()
          end if
       end if
    end if

    ! the error message from writing the table
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    ! close button
    call iw_setpos_bottomright(5,1)
    if (iw_button("Close")) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused(),okcloses=.false.)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_kpoints

end submodule kpoints
