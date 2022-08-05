! -*-F90-*-

!xx! cimgui procedures

! [SECTION] MAIN CODE (most of the code! lots of stuff, needs tidying up!)
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

! [SECTION] INPUTS
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

! [SECTION] USER FACING STRUCTURES (ImGuiStyle, ImGuiIO)
! void igStyleColorsDark(ImGuiStyle* dst);
subroutine igStyleColorsDark(dst) bind(c,name="igStyleColorsDark")
  import c_ptr
  type(c_ptr), value :: dst
end subroutine igStyleColorsDark

! [SECTION] imgui_demo.cpp
! void igShowDemoWindow(bool* p_open);
subroutine igShowDemoWindow(p_open) bind(c,name="igShowDemoWindow")
  import c_bool
  logical(c_bool) :: p_open
end subroutine igShowDemoWindow

