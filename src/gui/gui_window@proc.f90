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
       sys_seed
    use tools_io, only: string
    class(window), intent(inout), target :: w

    character(kind=c_char,len=:), allocatable, target :: str
    type(ImVec2) :: outer_size
    integer(c_int) :: mycol_id, flags
    integer :: i
    integer, save :: selected_item = 0
    logical(c_bool) :: ldum

    str = "Structures" // c_null_char
    outer_size%x = 0
    outer_size%y = 0
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
    if (igBeginTable(c_loc(str),8,flags,outer_size,0._c_float)) then

       ! set up the columns
       ! ID - name - spg - Z - nneq - ncel - Volume - nmol
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

       str = "Volume" // c_null_char
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

       str = "Z" // c_null_char
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
          !    const bool item_is_selected = selection.contains(item->ID);
          call igTableNextRow(ImGuiTableRowFlags_None, 0._c_float);

          ! ID
          ldum = igTableSetColumnIndex(0)
          str = string(i) // c_null_char
          call igText(c_loc(str))
          ! name
          ldum = igTableSetColumnIndex(1)
          str = trim(sys_seed(i)%name) // c_null_char
          call igTextWrapped(c_loc(str))
          if (sys_status(i) == sys_init) then
             write (*,*) "write me, in the tree!"
             stop 1
          end if

       end do

       call igEndTable()
    end if


  end subroutine draw_tree

  !xx! private procedures

end submodule proc
