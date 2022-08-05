! -*-F90-*-

!xx! cimgui procedures

! [SECTION/imgui.cpp] MAIN CODE (most of the code! lots of stuff, needs tidying up!)
! ImGuiContext* igGetCurrentContext(void);
function igGetCurrentContext() bind(c,name="igGetCurrentContext")
  import c_ptr
  type(c_ptr) :: igGetCurrentContext
end function igGetCurrentContext
! ImGuiContext* igCreateContext(ImFontAtlas* shared_font_atlas);
function igCreateContext(shared_font_atlas) bind(c,name="igCreateContext")
  import c_ptr
  type(c_ptr), value :: shared_font_atlas
  type(c_ptr) :: igCreateContext
end function igCreateContext
! ImGuiIO* igGetIO(void);
function igGetIO() bind(c,name="igGetIO")
  import c_ptr
  type(c_ptr) :: igGetIO
end function igGetIO
! ImDrawData* igGetDrawData(void);
function igGetDrawData() bind(c,name="igGetDrawData")
  import c_ptr
  type(c_ptr) :: igGetDrawData
end function igGetDrawData
! void igDestroyContext(ImGuiContext* ctx);
subroutine igDestroyContext(ctx) bind(c,name="igDestroyContext")
  import c_ptr
  type(c_ptr), value :: ctx
end subroutine igDestroyContext
! void igNewFrame(void);
subroutine igNewFrame() bind(c,name="igNewFrame")
end subroutine igNewFrame
! void igRender(void);
subroutine igRender() bind(c,name="igRender")
end subroutine igRender

! [SECTION/imgui.cpp] INPUTS
! bool igIsKeyDown(ImGuiKey key);
function igIsKeyDown(key) bind(c,name="igIsKeyDown")
  import c_int, c_bool
  integer(c_int), value :: key
  logical(c_bool) :: igIsKeyDown
end function igIsKeyDown
! bool igIsKeyPressed(ImGuiKey key,bool repeat);
function igIsKeyPressed(key,repeat) bind(c,name="igIsKeyPressed")
  import c_int, c_bool
  integer(c_int), value :: key
  logical(c_bool), value :: repeat
  logical(c_bool) :: igIsKeyPressed
end function igIsKeyPressed
! bool igIsKeyReleased(ImGuiKey key);
function igIsKeyReleased(key) bind(c,name="igIsKeyReleased")
  import c_int, c_bool
  integer(c_int), value :: key
  logical(c_bool) :: igIsKeyReleased
end function igIsKeyReleased
! const char* igGetKeyName(ImGuiKey key);
function igGetKeyName(key) bind(c,name="igGetKeyName")
  import c_int, c_ptr
  integer(c_int), value :: key
  type(c_ptr) :: igGetKeyName
end function igGetKeyName

! [SECTION/imgui.cpp] LAYOUT
! void igSameLine(float offset_from_start_x,float spacing);
subroutine igSameLine(offset_from_start_x, spacing) bind(c,name="igSameLine")
  import c_float
  real(c_float), value :: offset_from_start_x
  real(c_float), value :: spacing
end subroutine igSameLine
! void igGetContentRegionAvail(ImVec2 *pOut);
subroutine igGetContentRegionAvail(pOut) bind(c,name="igGetContentRegionAvail")
  import ImVec2
  type(ImVec2) :: pOut
end subroutine igGetContentRegionAvail
! void igText(const char* fmt,...); (Fortran handles the formatting)
subroutine igText(fmt) bind(c,name="igText")
  import c_ptr
  type(c_ptr), intent(in), value :: fmt
end subroutine igText

! [SECTION/imgui.cpp] USER FACING STRUCTURES (ImGuiStyle, ImGuiIO)
! void igStyleColorsDark(ImGuiStyle* dst);
subroutine igStyleColorsDark(dst) bind(c,name="igStyleColorsDark")
  import c_ptr
  type(c_ptr), value :: dst
end subroutine igStyleColorsDark

! [SECTION/imgui_demo.cpp]
! void igShowDemoWindow(bool* p_open);
subroutine igShowDemoWindow(p_open) bind(c,name="igShowDemoWindow")
  import c_bool
  logical(c_bool) :: p_open
end subroutine igShowDemoWindow

! [SECTION/header] Dear ImGui end-user API functions
! bool igBeginMainMenuBar(void);
function igBeginMainMenuBar() bind(c,name="igBeginMainMenuBar")
  import c_bool
  logical(c_bool) :: igBeginMainMenuBar
end function igBeginMainMenuBar
! void igEndMainMenuBar(void);
subroutine igEndMainMenuBar() bind(c,name="igEndMainMenuBar")
end subroutine igEndMainMenuBar
! bool igBeginMenu(const char* label,bool enabled);
function igBeginMenu(label, enabled) bind(c,name="igBeginMenu")
  import c_ptr, c_bool
  type(c_ptr), intent(in), value :: label
  logical(c_bool), value :: enabled
  logical(c_bool) :: igBeginMenu
end function igBeginMenu
! void igEndMenu(void);
subroutine igEndMenu() bind(c,name="igEndMenu")
end subroutine igEndMenu
! bool igMenuItem_Bool(const char* label,const char* shortcut,bool selected,bool enabled);
function igMenuItem_Bool(label, shortcut, selected, enabled) bind(c,name="igMenuItem_Bool")
  import c_ptr, c_bool
  type(c_ptr), value, intent(in) :: label
  type(c_ptr), value, intent(in) :: shortcut
  logical(c_bool), value :: selected
  logical(c_bool), value :: enabled
  logical(c_bool) :: igMenuItem_Bool
end function igMenuItem_Bool
! bool igMenuItem_BoolPtr(const char* label,const char* shortcut,bool* p_selected,bool enabled);
function igMenuItem_BoolPtr(label, shortcut, p_selected, enabled) bind(c,name="igMenuItem_BoolPtr")
  import c_ptr, c_bool
  type(c_ptr), value, intent(in) :: label
  type(c_ptr), value, intent(in) :: shortcut
  logical(c_bool) :: p_selected
  logical(c_bool), value :: enabled
  logical(c_bool) :: igMenuItem_BoolPtr
end function igMenuItem_BoolPtr
