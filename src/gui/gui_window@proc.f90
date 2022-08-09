! Copyright (c) 2019 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Some utilities for building the GUI (e.g. wrappers around ImGui routines).
submodule (gui_window) proc
  use gui_interfaces_cimgui
  implicit none

  ! Count ID for keeping track of windows and widgets
  integer :: idcount = 0

contains

  !> Initialize a window of the given type. If isiopen, initialize it
  !> as open.
  module subroutine window_init(w,type,isopen)
    class(window), intent(inout) :: w
    integer, intent(in) :: type
    logical, intent(in) :: isopen

    w%isinit = .true.
    w%isopen = isopen
    w%type = type
    w%id = -1

  end subroutine window_init

  !> Draw an ImGui window.
  module subroutine window_draw(w)
    use tools_io, only: string
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str

    ! First pass: assign ID, name, and flags
    if (w%id < 0) then
       idcount = idcount + 1
       w%id = idcount
       if (w%type == wintype_tree) then
          ! w%name = "Tree###" // string(w%id) // c_null_char
          w%name = "Tree" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_view) then
          ! w%name = "View###" // string(w%id) // c_null_char
          w%name = "View" // c_null_char
          w%flags = ImGuiWindowFlags_None
       elseif (w%type == wintype_console) then
          ! w%name = "Console###" // string(w%id) // c_null_char
          w%name = "Console" // c_null_char
          w%flags = ImGuiWindowFlags_None
       end if
    end if

    if (w%isopen) then
       if (igBegin(c_loc(w%name),w%isopen,w%flags)) then
          ! assign the pointer ID
          w%ptr = igGetCurrentWindow()

          ! draw the window contents, depending on type
          if (w%type == wintype_tree) then
             call w%draw_tree()
          elseif (w%type == wintype_view) then
             str = "Hello View!"
             call igText(c_loc(str))
          elseif (w%type == wintype_console) then
             str = "Hello Console!"
             call igText(c_loc(str))
          end if
       end if
       call igEnd()
    end if

  end subroutine window_draw

  !> Draw the contents of a tree window
  module subroutine draw_tree(w)
    use gui_main, only: nsys, sys, sys_status, sys_empty, sys_init, sys_loaded_not_init,&
       sys_seed, system_initialize
    use tools_io, only: string
    use param, only: bohrtoa
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str
    type(ImVec2) :: zero2
    integer(c_int) :: mycol_id, flags
    integer :: i
    logical(c_bool) :: ldum, selected

    ! initialize the currently selected system
    if (w%table_selected_sys > 0 .and. w%table_selected_sys <= nsys) then
       if (sys_status(w%table_selected_sys) /= sys_init) &
          call system_initialize(w%table_selected_sys)
    end if

    ! two zeros
    zero2%x = 0
    zero2%y = 0

    ! set up the table
    str = "Structures" // c_null_char
    flags = ImGuiTableFlags_Borders
    flags = ior(flags,ImGuiTableFlags_Resizable)
    flags = ior(flags,ImGuiTableFlags_Reorderable)
    flags = ior(flags,ImGuiTableFlags_Hideable)
    flags = ior(flags,ImGuiTableFlags_Sortable)
    flags = ior(flags,ImGuiTableFlags_SortMulti)
    flags = ior(flags,ImGuiTableFlags_SizingFixedFit)
    ! flags = ior(flags,ImGuiTableFlags_NoHostExtendX)
    ! flags = ior(flags,ImGuiTableFlags_ScrollX)
    ! flags = ior(flags,)
    if (igBeginTable(c_loc(str),13,flags,zero2,0._c_float)) then

       ! set up the columns
       ! ID - name - spg - volume - nneq - ncel - nmol - a - b - c - alpha - beta - gamma
       str = "ID" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultSort
       flags = ior(flags,ImGuiTableColumnFlags_PreferSortDescending)
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "Name" // c_null_char
       flags = ImGuiTableColumnFlags_WidthStretch
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "spg" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "V(A^3)" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "nneq" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "ncel" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "nmol" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "a" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "b" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "c" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "alpha" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "beta" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

       str = "gamma" // c_null_char
       flags = ImGuiTableColumnFlags_DefaultHide
       call igTableSetupColumn(c_loc(str),flags,0.0_c_float,0)

    ! // Declare columns
    ! // We use the "user_id" parameter of TableSetupColumn() to specify a user id that will be stored in the sort specifications.
    ! // This is so our sort function can identify a column given our own identifier. We could also identify them based on their index!
    ! ImGui::TableSetupColumn("ID",           ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_NoHide, 0.0f, MyItemColumnID_ID);
    ! ImGui::TableSetupColumn("Name",         ImGuiTableColumnFlags_WidthFixed, 0.0f, MyItemColumnID_Name);
    ! ImGui::TableSetupColumn("Action",       ImGuiTableColumnFlags_NoSort | ImGuiTableColumnFlags_WidthFixed, 0.0f, MyItemColumnID_Action);
    ! ImGui::TableSetupColumn("Quantity",     ImGuiTableColumnFlags_PreferSortDescending, 0.0f, MyItemColumnID_Quantity);
    ! ImGui::TableSetupColumn("Description",  (flags & ImGuiTableFlags_NoHostExtendX) ? 0 : ImGuiTableColumnFlags_WidthStretch, 0.0f, MyItemColumnID_Description);
    ! ImGui::TableSetupColumn("Hidden",       ImGuiTableColumnFlags_DefaultHide | ImGuiTableColumnFlags_NoSort);
    ! ImGui::TableSetupScrollFreeze(freeze_cols, freeze_rows);

    ! // Sort our data if sort specs have been changed!
    ! ImGuiTableSortSpecs* sorts_specs = ImGui::TableGetSortSpecs();
    ! if (sorts_specs && sorts_specs->SpecsDirty)
    !     items_need_sort = true;
    ! if (sorts_specs && items_need_sort && items.Size > 1)
    ! {
    !     MyItem::s_current_sort_specs = sorts_specs; // Store in variable accessible by the sort function.
    !     qsort(&items[0], (size_t)items.Size, sizeof(items[0]), MyItem::CompareWithSortSpecs);
    !     MyItem::s_current_sort_specs = NULL;
    !     sorts_specs->SpecsDirty = false;
    ! }
    ! items_need_sort = false;

    ! // Take note of whether we are currently sorting based on the Quantity field,
    ! // we will use this to trigger sorting when we know the data of this column has been modified.
    ! const bool sorts_specs_using_quantity = (ImGui::TableGetColumnFlags(3) & ImGuiTableColumnFlags_IsSorted) != 0;

       ! draw the header
       call igTableHeadersRow()

       ! draw the rows
       do i = 1, nsys
          if (sys_status(i) == sys_empty) cycle
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float);

          ! ID
          if (igTableSetColumnIndex(0)) then
             str = string(i) // c_null_char
             flags = ImGuiSelectableFlags_SpanAllColumns
             selected = (w%table_selected_sys==i)
             if (igSelectable_Bool(c_loc(str),selected,flags,zero2)) then
                w%table_selected_sys = i
             end if
          end if

          ! name
          if (igTableSetColumnIndex(1)) then
             str = trim(sys_seed(i)%file) // c_null_char
             call igTextWrapped(c_loc(str))
          end if

          if (sys_status(i) == sys_init) then
             if (igTableSetColumnIndex(2)) then ! spg
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                elseif (.not.sys(i)%c%spgavail) then
                   str = "n/a" // c_null_char
                else
                   str = trim(sys(i)%c%spg%international_symbol) // c_null_char
                end if
                call igText(c_loc(str))
             end if

             if (igTableSetColumnIndex(3)) then ! volume
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                else
                   str = string(sys(i)%c%omega,'f',decimal=2) // c_null_char
                end if
                call igText(c_loc(str))
             end if

             if (igTableSetColumnIndex(4)) then ! nneq
                str = string(sys(i)%c%nneq) // c_null_char
                call igText(c_loc(str))
             end if

             if (igTableSetColumnIndex(5)) then ! ncel
                str = string(sys(i)%c%ncel) // c_null_char
                call igText(c_loc(str))
             end if

             if (igTableSetColumnIndex(6)) then ! nmol
                str = string(sys(i)%c%nmol) // c_null_char
                call igText(c_loc(str))
             end if

             if (igTableSetColumnIndex(7)) then ! a
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                else
                   str = string(sys(i)%c%aa(1)*bohrtoa,'f',decimal=4) // c_null_char
                end if
                call igText(c_loc(str))
             end if
             if (igTableSetColumnIndex(8)) then ! b
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                else
                   str = string(sys(i)%c%aa(2)*bohrtoa,'f',decimal=4) // c_null_char
                end if
                call igText(c_loc(str))
             end if
             if (igTableSetColumnIndex(9)) then ! c
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                else
                   str = string(sys(i)%c%aa(3)*bohrtoa,'f',decimal=4) // c_null_char
                end if
                call igText(c_loc(str))
             end if
             if (igTableSetColumnIndex(10)) then ! alpha
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                else
                   str = string(sys(i)%c%bb(1),'f',decimal=2) // c_null_char
                end if
                call igText(c_loc(str))
             end if
             if (igTableSetColumnIndex(11)) then ! beta
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                else
                   str = string(sys(i)%c%bb(2),'f',decimal=2) // c_null_char
                end if
                call igText(c_loc(str))
             end if
             if (igTableSetColumnIndex(12)) then ! gamma
                if (sys(i)%c%ismolecule) then
                   str = "<mol>" // c_null_char
                else
                   str = string(sys(i)%c%bb(3),'f',decimal=2) // c_null_char
                end if
                call igText(c_loc(str))
             end if

          end if

       end do

       call igEndTable()
    end if

  end subroutine draw_tree

  !xx! private procedures

end submodule proc
