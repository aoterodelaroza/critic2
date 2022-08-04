! -*-F90-*-

!xx! cimgui procedures
! ImGuiContext* igGetCurrentContext(void);
function igGetCurrentContext() bind(c,name="igGetCurrentContext")
  import c_ptr
  type(c_ptr) :: igGetCurrentContext
end function igGetCurrentContext
! ImGuiIO* igGetIO(void);
function igGetIO() bind(c,name="igGetIO")
  import c_ptr
  type(c_ptr) :: igGetIO
end function igGetIO
! ImGuiContext* igCreateContext(ImFontAtlas* shared_font_atlas);
function igCreateContext(shared_font_atlas) bind(c,name="igCreateContext")
  import c_ptr
  type(c_ptr), value :: shared_font_atlas
  type(c_ptr) :: igCreateContext
end function igCreateContext
! void igDestroyContext(ImGuiContext* ctx);
subroutine igDestroyContext(ctx) bind(c,name="igDestroyContext")
  import c_ptr
  type(c_ptr), value :: ctx
end subroutine igDestroyContext
! void igStyleColorsDark(ImGuiStyle* dst);
subroutine igStyleColorsDark(dst) bind(c,name="igStyleColorsDark")
  import c_ptr
  type(c_ptr), value :: dst
end subroutine igStyleColorsDark
! void igNewFrame(void);
subroutine igNewFrame() bind(c,name="igNewFrame")
end subroutine igNewFrame
! void igShowDemoWindow(bool* p_open);
subroutine igShowDemoWindow(p_open) bind(c,name="igShowDemoWindow")
  import c_bool
  logical(c_bool) :: p_open
end subroutine igShowDemoWindow
! void igRender(void);
subroutine igRender() bind(c,name="igRender")
end subroutine igRender
! ImDrawData* igGetDrawData(void);
function igGetDrawData() bind(c,name="igGetDrawData")
  import c_ptr
  type(c_ptr) :: igGetDrawData
end function igGetDrawData
