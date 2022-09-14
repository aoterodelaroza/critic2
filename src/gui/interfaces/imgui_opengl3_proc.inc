! -*-F90-*-

!xx! ImGui OpenGL3 backend procedures (imgui_impl_opengl3)
! bool ImGui_ImplOpenGL3_Init(const char* glsl_version);
function ImGui_ImplOpenGL3_Init(glsl_version) bind(c,name="ImGui_ImplOpenGL3_Init")
  import c_bool, c_ptr
  type(c_ptr), value, intent(in) :: glsl_version
  logical(c_bool) :: ImGui_ImplOpenGL3_Init
end function ImGui_ImplOpenGL3_Init
! void ImGui_ImplOpenGL3_NewFrame(void);
subroutine ImGui_ImplOpenGL3_NewFrame() bind(c,name="ImGui_ImplOpenGL3_NewFrame")
end subroutine ImGui_ImplOpenGL3_NewFrame
! void ImGui_ImplOpenGL3_RenderDrawData(ImDrawData* draw_data);
subroutine ImGui_ImplOpenGL3_RenderDrawData(draw_data) bind(c,name="ImGui_ImplOpenGL3_RenderDrawData")
  import c_ptr
  type(c_ptr), value :: draw_data
end subroutine ImGui_ImplOpenGL3_RenderDrawData
! void ImGui_ImplOpenGL3_Shutdown(void);
subroutine ImGui_ImplOpenGL3_Shutdown() bind(c,name="ImGui_ImplOpenGL3_Shutdown")
end subroutine ImGui_ImplOpenGL3_Shutdown
