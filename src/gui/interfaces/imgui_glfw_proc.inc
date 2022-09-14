! -*-F90-*-

!xx! ImGui GLFW backend procedures (imgui_impl_glfw)
! bool ImGui_ImplGlfw_InitForOpenGL(GLFWwindow* window, bool install_callbacks);
function ImGui_ImplGlfw_InitForOpenGL(window, install_callbacks) bind(c,name="ImGui_ImplGlfw_InitForOpenGL")
  import c_bool, c_ptr
  type(c_ptr), value :: window
  logical(c_bool), value :: install_callbacks
  logical(c_bool) :: ImGui_ImplGlfw_InitForOpenGL
end function ImGui_ImplGlfw_InitForOpenGL
! void ImGui_ImplGlfw_NewFrame(void);
subroutine ImGui_ImplGlfw_NewFrame() bind(c,name="ImGui_ImplGlfw_NewFrame")
end subroutine ImGui_ImplGlfw_NewFrame
! void ImGui_ImplGlfw_Shutdown(void);
subroutine ImGui_ImplGlfw_Shutdown() bind(c,name="ImGui_ImplGlfw_Shutdown")
end subroutine ImGui_ImplGlfw_Shutdown
