! -*-F90-*-

!xx! cimgui procedures

!x! ImVec2* ImVec2_ImVec2_Nil(void);
!x! void ImVec2_destroy(ImVec2* self);
!x! ImVec2* ImVec2_ImVec2_Float(float _x,float _y);
!x! ImVec4* ImVec4_ImVec4_Nil(void);
!x! void ImVec4_destroy(ImVec4* self);
!x! ImVec4* ImVec4_ImVec4_Float(float _x,float _y,float _z,float _w);
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
! ImGuiContext* igGetCurrentContext(void);
function igGetCurrentContext() bind(c,name="igGetCurrentContext")
  import c_ptr
  type(c_ptr) :: igGetCurrentContext
end function igGetCurrentContext
!x! void igSetCurrentContext(ImGuiContext* ctx);
! ImGuiIO* igGetIO(void);
function igGetIO() bind(c,name="igGetIO")
  import c_ptr
  type(c_ptr) :: igGetIO
end function igGetIO
!x! ImGuiStyle* igGetStyle(void);
! void igNewFrame(void);
subroutine igNewFrame() bind(c,name="igNewFrame")
end subroutine igNewFrame
!x! void igEndFrame(void);
! void igRender(void);
subroutine igRender() bind(c,name="igRender")
end subroutine igRender
! ImDrawData* igGetDrawData(void);
function igGetDrawData() bind(c,name="igGetDrawData")
  import c_ptr
  type(c_ptr) :: igGetDrawData
end function igGetDrawData
! void igShowDemoWindow(bool* p_open);
subroutine igShowDemoWindow(p_open) bind(c,name="igShowDemoWindow")
  import c_bool
  logical(c_bool) :: p_open
end subroutine igShowDemoWindow
!x! void igShowMetricsWindow(bool* p_open);
!x! void igShowDebugLogWindow(bool* p_open);
! void igShowStackToolWindow(bool* p_open);
subroutine igShowStackToolWindow(p_open) bind(c,name="igShowStackToolWindow")
  import c_bool
  logical(c_bool) :: p_open
end subroutine igShowStackToolWindow
!x! void igShowAboutWindow(bool* p_open);
!x! void igShowStyleEditor(ImGuiStyle* ref);
!x! bool igShowStyleSelector(const char* label);
!x! void igShowFontSelector(const char* label);
!x! void igShowUserGuide(void);
!x! const char* igGetVersion(void);
! void igStyleColorsDark(ImGuiStyle* dst);
subroutine igStyleColorsDark(dst) bind(c,name="igStyleColorsDark")
  import c_ptr
  type(c_ptr), value :: dst
end subroutine igStyleColorsDark
!x! void igStyleColorsLight(ImGuiStyle* dst);
!x! void igStyleColorsClassic(ImGuiStyle* dst);
! bool igBegin(const char* name,bool* p_open,ImGuiWindowFlags flags);
function igBegin(name,p_open,flags) bind(c,name="igBegin")
  import c_ptr, c_bool, c_int
  type(c_ptr), intent(in), value :: name
  logical(c_bool) :: p_open
  integer(c_int), value :: flags
  logical(c_bool) :: igBegin
end function igBegin
! void igEnd(void);
subroutine igEnd() bind(c,name="igEnd")
end subroutine igEnd
!x! bool igBeginChild_Str(const char* str_id,const ImVec2 size,bool border,ImGuiWindowFlags flags);
!x! bool igBeginChild_ID(ImGuiID id,const ImVec2 size,bool border,ImGuiWindowFlags flags);
!x! void igEndChild(void);
!x! bool igIsWindowAppearing(void);
!x! bool igIsWindowCollapsed(void);
!x! bool igIsWindowFocused(ImGuiFocusedFlags flags);
!x! bool igIsWindowHovered(ImGuiHoveredFlags flags);
!x! ImDrawList* igGetWindowDrawList(void);
!x! float igGetWindowDpiScale(void);
!x! void igGetWindowPos(ImVec2 *pOut);
!x! void igGetWindowSize(ImVec2 *pOut);
!x! float igGetWindowWidth(void);
!x! float igGetWindowHeight(void);
!x! ImGuiViewport* igGetWindowViewport(void);
!x! void igSetNextWindowPos(const ImVec2 pos,ImGuiCond cond,const ImVec2 pivot);
!x! void igSetNextWindowSize(const ImVec2 size,ImGuiCond cond);
!x! void igSetNextWindowSizeConstraints(const ImVec2 size_min,const ImVec2 size_max,ImGuiSizeCallback custom_callback,void* custom_callback_data);
!x! void igSetNextWindowContentSize(const ImVec2 size);
!x! void igSetNextWindowCollapsed(bool collapsed,ImGuiCond cond);
!x! void igSetNextWindowFocus(void);
!x! void igSetNextWindowBgAlpha(float alpha);
!x! void igSetNextWindowViewport(ImGuiID viewport_id);
!x! void igSetWindowPos_Vec2(const ImVec2 pos,ImGuiCond cond);
!x! void igSetWindowSize_Vec2(const ImVec2 size,ImGuiCond cond);
!x! void igSetWindowCollapsed_Bool(bool collapsed,ImGuiCond cond);
!x! void igSetWindowFocus_Nil(void);
!x! void igSetWindowFontScale(float scale);
!x! void igSetWindowPos_Str(const char* name,const ImVec2 pos,ImGuiCond cond);
!x! void igSetWindowSize_Str(const char* name,const ImVec2 size,ImGuiCond cond);
!x! void igSetWindowCollapsed_Str(const char* name,bool collapsed,ImGuiCond cond);
!x! void igSetWindowFocus_Str(const char* name);
! void igGetContentRegionAvail(ImVec2 *pOut);
subroutine igGetContentRegionAvail(pOut) bind(c,name="igGetContentRegionAvail")
  import ImVec2
  type(ImVec2) :: pOut
end subroutine igGetContentRegionAvail
!x! void igGetContentRegionMax(ImVec2 *pOut);
!x! void igGetWindowContentRegionMin(ImVec2 *pOut);
!x! void igGetWindowContentRegionMax(ImVec2 *pOut);
!x! float igGetScrollX(void);
!x! float igGetScrollY(void);
!x! void igSetScrollX_Float(float scroll_x);
!x! void igSetScrollY_Float(float scroll_y);
!x! float igGetScrollMaxX(void);
!x! float igGetScrollMaxY(void);
!x! void igSetScrollHereX(float center_x_ratio);
!x! void igSetScrollHereY(float center_y_ratio);
!x! void igSetScrollFromPosX_Float(float local_x,float center_x_ratio);
!x! void igSetScrollFromPosY_Float(float local_y,float center_y_ratio);
!x! void igPushFont(ImFont* font);
!x! void igPopFont(void);
!x! void igPushStyleColor_U32(ImGuiCol idx,ImU32 col);
!x! void igPushStyleColor_Vec4(ImGuiCol idx,const ImVec4 col);
!x! void igPopStyleColor(int count);
!x! void igPushStyleVar_Float(ImGuiStyleVar idx,float val);
!x! void igPushStyleVar_Vec2(ImGuiStyleVar idx,const ImVec2 val);
!x! void igPopStyleVar(int count);
!x! void igPushAllowKeyboardFocus(bool allow_keyboard_focus);
!x! void igPopAllowKeyboardFocus(void);
!x! void igPushButtonRepeat(bool repeat);
!x! void igPopButtonRepeat(void);
!x! void igPushItemWidth(float item_width);
!x! void igPopItemWidth(void);
!x! void igSetNextItemWidth(float item_width);
!x! float igCalcItemWidth(void);
!x! void igPushTextWrapPos(float wrap_local_pos_x);
!x! void igPopTextWrapPos(void);
!x! ImFont* igGetFont(void);
!x! float igGetFontSize(void);
!x! void igGetFontTexUvWhitePixel(ImVec2 *pOut);
!x! ImU32 igGetColorU32_Col(ImGuiCol idx,float alpha_mul);
!x! ImU32 igGetColorU32_Vec4(const ImVec4 col);
!x! ImU32 igGetColorU32_U32(ImU32 col);
!x! const ImVec4* igGetStyleColorVec4(ImGuiCol idx);
!x! void igSeparator(void);
! void igSameLine(float offset_from_start_x,float spacing);
subroutine igSameLine(offset_from_start_x, spacing) bind(c,name="igSameLine")
  import c_float
  real(c_float), value :: offset_from_start_x
  real(c_float), value :: spacing
end subroutine igSameLine
!x! void igNewLine(void);
!x! void igSpacing(void);
!x! void igDummy(const ImVec2 size);
!x! void igIndent(float indent_w);
!x! void igUnindent(float indent_w);
!x! void igBeginGroup(void);
!x! void igEndGroup(void);
!x! void igGetCursorPos(ImVec2 *pOut);
!x! float igGetCursorPosX(void);
!x! float igGetCursorPosY(void);
!x! void igSetCursorPos(const ImVec2 local_pos);
!x! void igSetCursorPosX(float local_x);
!x! void igSetCursorPosY(float local_y);
!x! void igGetCursorStartPos(ImVec2 *pOut);
!x! void igGetCursorScreenPos(ImVec2 *pOut);
!x! void igSetCursorScreenPos(const ImVec2 pos);
!x! void igAlignTextToFramePadding(void);
!x! float igGetTextLineHeight(void);
!x! float igGetTextLineHeightWithSpacing(void);
!x! float igGetFrameHeight(void);
!x! float igGetFrameHeightWithSpacing(void);
!x! void igPushID_Str(const char* str_id);
!x! void igPushID_StrStr(const char* str_id_begin,const char* str_id_end);
!x! void igPushID_Ptr(const void* ptr_id);
!x! void igPushID_Int(int int_id);
!x! void igPopID(void);
!x! ImGuiID igGetID_Str(const char* str_id);
!x! ImGuiID igGetID_StrStr(const char* str_id_begin,const char* str_id_end);
!x! ImGuiID igGetID_Ptr(const void* ptr_id);
!x! void igTextUnformatted(const char* text,const char* text_end);
! void igText(const char* fmt,...);
subroutine igText(fmt) bind(c,name="igText")
  import c_ptr
  type(c_ptr), intent(in), value :: fmt
end subroutine igText
!x! void igTextV(const char* fmt,va_list args);
!x! void igTextColored(const ImVec4 col,const char* fmt,...);
!x! void igTextColoredV(const ImVec4 col,const char* fmt,va_list args);
!x! void igTextDisabled(const char* fmt,...);
!x! void igTextDisabledV(const char* fmt,va_list args);
!x! void igTextWrapped(const char* fmt,...);
!x! void igTextWrappedV(const char* fmt,va_list args);
!x! void igLabelText(const char* label,const char* fmt,...);
!x! void igLabelTextV(const char* label,const char* fmt,va_list args);
!x! void igBulletText(const char* fmt,...);
!x! void igBulletTextV(const char* fmt,va_list args);
!x! bool igButton(const char* label,const ImVec2 size);
!x! bool igSmallButton(const char* label);
!x! bool igInvisibleButton(const char* str_id,const ImVec2 size,ImGuiButtonFlags flags);
!x! bool igArrowButton(const char* str_id,ImGuiDir dir);
!x! void igImage(ImTextureID user_texture_id,const ImVec2 size,const ImVec2 uv0,const ImVec2 uv1,const ImVec4 tint_col,const ImVec4 border_col);
!x! bool igImageButton(ImTextureID user_texture_id,const ImVec2 size,const ImVec2 uv0,const ImVec2 uv1,int frame_padding,const ImVec4 bg_col,const ImVec4 tint_col);
!x! bool igCheckbox(const char* label,bool* v);
!x! bool igCheckboxFlags_IntPtr(const char* label,int* flags,int flags_value);
!x! bool igCheckboxFlags_UintPtr(const char* label,unsigned int* flags,unsigned int flags_value);
!x! bool igRadioButton_Bool(const char* label,bool active);
!x! bool igRadioButton_IntPtr(const char* label,int* v,int v_button);
!x! void igProgressBar(float fraction,const ImVec2 size_arg,const char* overlay);
!x! void igBullet(void);
!x! bool igBeginCombo(const char* label,const char* preview_value,ImGuiComboFlags flags);
!x! void igEndCombo(void);
!x! bool igCombo_Str_arr(const char* label,int* current_item,const char* const items[],int items_count,int popup_max_height_in_items);
!x! bool igCombo_Str(const char* label,int* current_item,const char* items_separated_by_zeros,int popup_max_height_in_items);
!x! bool igCombo_FnBoolPtr(const char* label,int* current_item,bool(*items_getter)(void* data,int idx,const char** out_text),void* data,int items_count,int popup_max_height_in_items);
!x! bool igDragFloat(const char* label,float* v,float v_speed,float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragFloat2(const char* label,float v[2],float v_speed,float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragFloat3(const char* label,float v[3],float v_speed,float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragFloat4(const char* label,float v[4],float v_speed,float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragFloatRange2(const char* label,float* v_current_min,float* v_current_max,float v_speed,float v_min,float v_max,const char* format,const char* format_max,ImGuiSliderFlags flags);
!x! bool igDragInt(const char* label,int* v,float v_speed,int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragInt2(const char* label,int v[2],float v_speed,int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragInt3(const char* label,int v[3],float v_speed,int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragInt4(const char* label,int v[4],float v_speed,int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragIntRange2(const char* label,int* v_current_min,int* v_current_max,float v_speed,int v_min,int v_max,const char* format,const char* format_max,ImGuiSliderFlags flags);
!x! bool igDragScalar(const char* label,ImGuiDataType data_type,void* p_data,float v_speed,const void* p_min,const void* p_max,const char* format,ImGuiSliderFlags flags);
!x! bool igDragScalarN(const char* label,ImGuiDataType data_type,void* p_data,int components,float v_speed,const void* p_min,const void* p_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderFloat(const char* label,float* v,float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderFloat2(const char* label,float v[2],float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderFloat3(const char* label,float v[3],float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderFloat4(const char* label,float v[4],float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderAngle(const char* label,float* v_rad,float v_degrees_min,float v_degrees_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderInt(const char* label,int* v,int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderInt2(const char* label,int v[2],int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderInt3(const char* label,int v[3],int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderInt4(const char* label,int v[4],int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderScalar(const char* label,ImGuiDataType data_type,void* p_data,const void* p_min,const void* p_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderScalarN(const char* label,ImGuiDataType data_type,void* p_data,int components,const void* p_min,const void* p_max,const char* format,ImGuiSliderFlags flags);
!x! bool igVSliderFloat(const char* label,const ImVec2 size,float* v,float v_min,float v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igVSliderInt(const char* label,const ImVec2 size,int* v,int v_min,int v_max,const char* format,ImGuiSliderFlags flags);
!x! bool igVSliderScalar(const char* label,const ImVec2 size,ImGuiDataType data_type,void* p_data,const void* p_min,const void* p_max,const char* format,ImGuiSliderFlags flags);
!x! bool igInputText(const char* label,char* buf,size_t buf_size,ImGuiInputTextFlags flags,ImGuiInputTextCallback callback,void* user_data);
!x! bool igInputTextMultiline(const char* label,char* buf,size_t buf_size,const ImVec2 size,ImGuiInputTextFlags flags,ImGuiInputTextCallback callback,void* user_data);
!x! bool igInputTextWithHint(const char* label,const char* hint,char* buf,size_t buf_size,ImGuiInputTextFlags flags,ImGuiInputTextCallback callback,void* user_data);
!x! bool igInputFloat(const char* label,float* v,float step,float step_fast,const char* format,ImGuiInputTextFlags flags);
!x! bool igInputFloat2(const char* label,float v[2],const char* format,ImGuiInputTextFlags flags);
!x! bool igInputFloat3(const char* label,float v[3],const char* format,ImGuiInputTextFlags flags);
!x! bool igInputFloat4(const char* label,float v[4],const char* format,ImGuiInputTextFlags flags);
!x! bool igInputInt(const char* label,int* v,int step,int step_fast,ImGuiInputTextFlags flags);
!x! bool igInputInt2(const char* label,int v[2],ImGuiInputTextFlags flags);
!x! bool igInputInt3(const char* label,int v[3],ImGuiInputTextFlags flags);
!x! bool igInputInt4(const char* label,int v[4],ImGuiInputTextFlags flags);
!x! bool igInputDouble(const char* label,double* v,double step,double step_fast,const char* format,ImGuiInputTextFlags flags);
!x! bool igInputScalar(const char* label,ImGuiDataType data_type,void* p_data,const void* p_step,const void* p_step_fast,const char* format,ImGuiInputTextFlags flags);
!x! bool igInputScalarN(const char* label,ImGuiDataType data_type,void* p_data,int components,const void* p_step,const void* p_step_fast,const char* format,ImGuiInputTextFlags flags);
!x! bool igColorEdit3(const char* label,float col[3],ImGuiColorEditFlags flags);
!x! bool igColorEdit4(const char* label,float col[4],ImGuiColorEditFlags flags);
!x! bool igColorPicker3(const char* label,float col[3],ImGuiColorEditFlags flags);
!x! bool igColorPicker4(const char* label,float col[4],ImGuiColorEditFlags flags,const float* ref_col);
!x! bool igColorButton(const char* desc_id,const ImVec4 col,ImGuiColorEditFlags flags,const ImVec2 size);
!x! void igSetColorEditOptions(ImGuiColorEditFlags flags);
! bool igTreeNode_Str(const char* label);
function igTreeNode_Str(label) bind(c,name="igTreeNode_Str")
  import c_ptr, c_bool
  type(c_ptr), intent(in), value :: label
  logical(c_bool) :: igTreeNode_Str
end function igTreeNode_Str
!x! bool igTreeNode_StrStr(const char* str_id,const char* fmt,...);
!x! bool igTreeNode_Ptr(const void* ptr_id,const char* fmt,...);
!x! bool igTreeNodeV_Str(const char* str_id,const char* fmt,va_list args);
!x! bool igTreeNodeV_Ptr(const void* ptr_id,const char* fmt,va_list args);
!x! bool igTreeNodeEx_Str(const char* label,ImGuiTreeNodeFlags flags);
!x! bool igTreeNodeEx_StrStr(const char* str_id,ImGuiTreeNodeFlags flags,const char* fmt,...);
!x! bool igTreeNodeEx_Ptr(const void* ptr_id,ImGuiTreeNodeFlags flags,const char* fmt,...);
!x! bool igTreeNodeExV_Str(const char* str_id,ImGuiTreeNodeFlags flags,const char* fmt,va_list args);
!x! bool igTreeNodeExV_Ptr(const void* ptr_id,ImGuiTreeNodeFlags flags,const char* fmt,va_list args);
!x! void igTreePush_Str(const char* str_id);
!x! void igTreePush_Ptr(const void* ptr_id);
! void igTreePop(void);
subroutine igTreePop() bind(c,name="igTreePop")
end subroutine igTreePop
!x! float igGetTreeNodeToLabelSpacing(void);
!x! bool igCollapsingHeader_TreeNodeFlags(const char* label,ImGuiTreeNodeFlags flags);
!x! bool igCollapsingHeader_BoolPtr(const char* label,bool* p_visible,ImGuiTreeNodeFlags flags);
!x! void igSetNextItemOpen(bool is_open,ImGuiCond cond);
!x! bool igSelectable_Bool(const char* label,bool selected,ImGuiSelectableFlags flags,const ImVec2 size);
!x! bool igSelectable_BoolPtr(const char* label,bool* p_selected,ImGuiSelectableFlags flags,const ImVec2 size);
!x! bool igBeginListBox(const char* label,const ImVec2 size);
!x! void igEndListBox(void);
!x! bool igListBox_Str_arr(const char* label,int* current_item,const char* const items[],int items_count,int height_in_items);
!x! bool igListBox_FnBoolPtr(const char* label,int* current_item,bool(*items_getter)(void* data,int idx,const char** out_text),void* data,int items_count,int height_in_items);
!x! void igPlotLines_FloatPtr(const char* label,const float* values,int values_count,int values_offset,const char* overlay_text,float scale_min,float scale_max,ImVec2 graph_size,int stride);
!x! void igPlotLines_FnFloatPtr(const char* label,float(*values_getter)(void* data,int idx),void* data,int values_count,int values_offset,const char* overlay_text,float scale_min,float scale_max,ImVec2 graph_size);
!x! void igPlotHistogram_FloatPtr(const char* label,const float* values,int values_count,int values_offset,const char* overlay_text,float scale_min,float scale_max,ImVec2 graph_size,int stride);
!x! void igPlotHistogram_FnFloatPtr(const char* label,float(*values_getter)(void* data,int idx),void* data,int values_count,int values_offset,const char* overlay_text,float scale_min,float scale_max,ImVec2 graph_size);
!x! void igValue_Bool(const char* prefix,bool b);
!x! void igValue_Int(const char* prefix,int v);
!x! void igValue_Uint(const char* prefix,unsigned int v);
!x! void igValue_Float(const char* prefix,float v,const char* float_format);
!x! bool igBeginMenuBar(void);
!x! void igEndMenuBar(void);
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
!x! void igBeginTooltip(void);
!x! void igEndTooltip(void);
! void igSetTooltip(const char* fmt,...);
subroutine igSetTooltip(fmt) bind(c,name="igSetTooltip")
  import c_ptr
  type(c_ptr), intent(in), value :: fmt
end subroutine igSetTooltip
!x! void igSetTooltipV(const char* fmt,va_list args);
!x! bool igBeginPopup(const char* str_id,ImGuiWindowFlags flags);
!x! bool igBeginPopupModal(const char* name,bool* p_open,ImGuiWindowFlags flags);
!x! void igEndPopup(void);
!x! void igOpenPopup_Str(const char* str_id,ImGuiPopupFlags popup_flags);
!x! void igOpenPopup_ID(ImGuiID id,ImGuiPopupFlags popup_flags);
!x! void igOpenPopupOnItemClick(const char* str_id,ImGuiPopupFlags popup_flags);
!x! void igCloseCurrentPopup(void);
!x! bool igBeginPopupContextItem(const char* str_id,ImGuiPopupFlags popup_flags);
!x! bool igBeginPopupContextWindow(const char* str_id,ImGuiPopupFlags popup_flags);
!x! bool igBeginPopupContextVoid(const char* str_id,ImGuiPopupFlags popup_flags);
!x! bool igIsPopupOpen_Str(const char* str_id,ImGuiPopupFlags flags);
!x! bool igBeginTable(const char* str_id,int column,ImGuiTableFlags flags,const ImVec2 outer_size,float inner_width);
!x! void igEndTable(void);
!x! void igTableNextRow(ImGuiTableRowFlags row_flags,float min_row_height);
!x! bool igTableNextColumn(void);
!x! bool igTableSetColumnIndex(int column_n);
!x! void igTableSetupColumn(const char* label,ImGuiTableColumnFlags flags,float init_width_or_weight,ImGuiID user_id);
!x! void igTableSetupScrollFreeze(int cols,int rows);
!x! void igTableHeadersRow(void);
!x! void igTableHeader(const char* label);
!x! ImGuiTableSortSpecs* igTableGetSortSpecs(void);
!x! int igTableGetColumnCount(void);
!x! int igTableGetColumnIndex(void);
!x! int igTableGetRowIndex(void);
!x! const char* igTableGetColumnName_Int(int column_n);
!x! ImGuiTableColumnFlags igTableGetColumnFlags(int column_n);
!x! void igTableSetColumnEnabled(int column_n,bool v);
!x! void igTableSetBgColor(ImGuiTableBgTarget target,ImU32 color,int column_n);
!x! void igColumns(int count,const char* id,bool border);
!x! void igNextColumn(void);
!x! int igGetColumnIndex(void);
!x! float igGetColumnWidth(int column_index);
!x! void igSetColumnWidth(int column_index,float width);
!x! float igGetColumnOffset(int column_index);
!x! void igSetColumnOffset(int column_index,float offset_x);
!x! int igGetColumnsCount(void);
!x! bool igBeginTabBar(const char* str_id,ImGuiTabBarFlags flags);
!x! void igEndTabBar(void);
!x! bool igBeginTabItem(const char* label,bool* p_open,ImGuiTabItemFlags flags);
!x! void igEndTabItem(void);
!x! bool igTabItemButton(const char* label,ImGuiTabItemFlags flags);
!x! void igSetTabItemClosed(const char* tab_or_docked_window_label);
!x! ImGuiID igDockSpace(ImGuiID id,const ImVec2 size,ImGuiDockNodeFlags flags,const ImGuiWindowClass* window_class);
! ImGuiID igDockSpaceOverViewport(const ImGuiViewport* viewport,ImGuiDockNodeFlags flags,const ImGuiWindowClass* window_class);
function igDockSpaceOverViewport(viewport,flags,window_class) bind(c,name="igDockSpaceOverViewport")
  import c_int, c_ptr
  integer(c_int) :: igDockSpaceOverViewport
  type(c_ptr), intent(in), value :: viewport
  integer(c_int), value :: flags
  type(c_ptr), intent(in), value :: window_class
end function igDockSpaceOverViewport
!x! void igSetNextWindowDockID(ImGuiID dock_id,ImGuiCond cond);
!x! void igSetNextWindowClass(const ImGuiWindowClass* window_class);
!x! ImGuiID igGetWindowDockID(void);
!x! bool igIsWindowDocked(void);
!x! void igLogToTTY(int auto_open_depth);
!x! void igLogToFile(int auto_open_depth,const char* filename);
!x! void igLogToClipboard(int auto_open_depth);
!x! void igLogFinish(void);
!x! void igLogButtons(void);
!x! void igLogTextV(const char* fmt,va_list args);
!x! bool igBeginDragDropSource(ImGuiDragDropFlags flags);
!x! bool igSetDragDropPayload(const char* type,const void* data,size_t sz,ImGuiCond cond);
!x! void igEndDragDropSource(void);
!x! bool igBeginDragDropTarget(void);
!x! const ImGuiPayload* igAcceptDragDropPayload(const char* type,ImGuiDragDropFlags flags);
!x! void igEndDragDropTarget(void);
!x! const ImGuiPayload* igGetDragDropPayload(void);
!x! void igBeginDisabled(bool disabled);
!x! void igEndDisabled(void);
!x! void igPushClipRect(const ImVec2 clip_rect_min,const ImVec2 clip_rect_max,bool intersect_with_current_clip_rect);
!x! void igPopClipRect(void);
!x! void igSetItemDefaultFocus(void);
!x! void igSetKeyboardFocusHere(int offset);
! bool igIsItemHovered(ImGuiHoveredFlags flags);
function igIsItemHovered(flags) bind(c,name="igIsItemHovered")
  import c_int, c_bool
  integer(c_int), value :: flags
  logical(c_bool) :: igIsItemHovered
end function igIsItemHovered
!x! bool igIsItemActive(void);
!x! bool igIsItemFocused(void);
!x! bool igIsItemClicked(ImGuiMouseButton mouse_button);
!x! bool igIsItemVisible(void);
!x! bool igIsItemEdited(void);
!x! bool igIsItemActivated(void);
!x! bool igIsItemDeactivated(void);
!x! bool igIsItemDeactivatedAfterEdit(void);
!x! bool igIsItemToggledOpen(void);
!x! bool igIsAnyItemHovered(void);
!x! bool igIsAnyItemActive(void);
!x! bool igIsAnyItemFocused(void);
!x! void igGetItemRectMin(ImVec2 *pOut);
!x! void igGetItemRectMax(ImVec2 *pOut);
!x! void igGetItemRectSize(ImVec2 *pOut);
!x! void igSetItemAllowOverlap(void);
! ImGuiViewport* igGetMainViewport(void);
function igGetMainViewport() bind(c,name="igGetMainViewport")
  import c_ptr
  type(c_ptr) :: igGetMainViewport
end function igGetMainViewport
!x! ImDrawList* igGetBackgroundDrawList_Nil(void);
!x! ImDrawList* igGetForegroundDrawList_Nil(void);
!x! ImDrawList* igGetBackgroundDrawList_ViewportPtr(ImGuiViewport* viewport);
!x! ImDrawList* igGetForegroundDrawList_ViewportPtr(ImGuiViewport* viewport);
!x! bool igIsRectVisible_Nil(const ImVec2 size);
!x! bool igIsRectVisible_Vec2(const ImVec2 rect_min,const ImVec2 rect_max);
!x! double igGetTime(void);
!x! int igGetFrameCount(void);
!x! ImDrawListSharedData* igGetDrawListSharedData(void);
!x! const char* igGetStyleColorName(ImGuiCol idx);
!x! void igSetStateStorage(ImGuiStorage* storage);
!x! ImGuiStorage* igGetStateStorage(void);
!x! bool igBeginChildFrame(ImGuiID id,const ImVec2 size,ImGuiWindowFlags flags);
!x! void igEndChildFrame(void);
!x! void igCalcTextSize(ImVec2 *pOut,const char* text,const char* text_end,bool hide_text_after_double_hash,float wrap_width);
!x! void igColorConvertU32ToFloat4(ImVec4 *pOut,ImU32 in);
!x! ImU32 igColorConvertFloat4ToU32(const ImVec4 in);
!x! void igColorConvertRGBtoHSV(float r,float g,float b,float* out_h,float* out_s,float* out_v);
!x! void igColorConvertHSVtoRGB(float h,float s,float v,float* out_r,float* out_g,float* out_b);
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
!x! int igGetKeyPressedAmount(ImGuiKey key,float repeat_delay,float rate);
! const char* igGetKeyName(ImGuiKey key);
function igGetKeyName(key) bind(c,name="igGetKeyName")
  import c_int, c_ptr
  integer(c_int), value :: key
  type(c_ptr) :: igGetKeyName
end function igGetKeyName
!x! void igSetNextFrameWantCaptureKeyboard(bool want_capture_keyboard);
!x! bool igIsMouseDown(ImGuiMouseButton button);
!x! bool igIsMouseClicked(ImGuiMouseButton button,bool repeat);
!x! bool igIsMouseReleased(ImGuiMouseButton button);
!x! bool igIsMouseDoubleClicked(ImGuiMouseButton button);
!x! int igGetMouseClickedCount(ImGuiMouseButton button);
!x! bool igIsMouseHoveringRect(const ImVec2 r_min,const ImVec2 r_max,bool clip);
!x! bool igIsMousePosValid(const ImVec2* mouse_pos);
!x! bool igIsAnyMouseDown(void);
!x! void igGetMousePos(ImVec2 *pOut);
!x! void igGetMousePosOnOpeningCurrentPopup(ImVec2 *pOut);
!x! bool igIsMouseDragging(ImGuiMouseButton button,float lock_threshold);
!x! void igGetMouseDragDelta(ImVec2 *pOut,ImGuiMouseButton button,float lock_threshold);
!x! void igResetMouseDragDelta(ImGuiMouseButton button);
!x! ImGuiMouseCursor igGetMouseCursor(void);
!x! void igSetMouseCursor(ImGuiMouseCursor cursor_type);
!x! void igSetNextFrameWantCaptureMouse(bool want_capture_mouse);
!x! const char* igGetClipboardText(void);
!x! void igSetClipboardText(const char* text);
!x! void igLoadIniSettingsFromDisk(const char* ini_filename);
!x! void igLoadIniSettingsFromMemory(const char* ini_data,size_t ini_size);
!x! void igSaveIniSettingsToDisk(const char* ini_filename);
!x! const char* igSaveIniSettingsToMemory(size_t* out_ini_size);
!x! void igDebugTextEncoding(const char* text);
!x! bool igDebugCheckVersionAndDataLayout(const char* version_str,size_t sz_io,size_t sz_style,size_t sz_vec2,size_t sz_vec4,size_t sz_drawvert,size_t sz_drawidx);
!x! void igSetAllocatorFunctions(ImGuiMemAllocFunc alloc_func,ImGuiMemFreeFunc free_func,void* user_data);
!x! void igGetAllocatorFunctions(ImGuiMemAllocFunc* p_alloc_func,ImGuiMemFreeFunc* p_free_func,void** p_user_data);
!x! void* igMemAlloc(size_t size);
!x! void igMemFree(void* ptr);
!x! ImGuiPlatformIO* igGetPlatformIO(void);
!x! void igUpdatePlatformWindows(void);
!x! void igRenderPlatformWindowsDefault(void* platform_render_arg,void* renderer_render_arg);
!x! void igDestroyPlatformWindows(void);
!x! ImGuiViewport* igFindViewportByID(ImGuiID id);
!x! ImGuiViewport* igFindViewportByPlatformHandle(void* platform_handle);
!x! ImGuiStyle* ImGuiStyle_ImGuiStyle(void);
!x! void ImGuiStyle_destroy(ImGuiStyle* self);
!x! void ImGuiStyle_ScaleAllSizes(ImGuiStyle* self,float scale_factor);
!x! void ImGuiIO_AddKeyEvent(ImGuiIO* self,ImGuiKey key,bool down);
!x! void ImGuiIO_AddKeyAnalogEvent(ImGuiIO* self,ImGuiKey key,bool down,float v);
!x! void ImGuiIO_AddMousePosEvent(ImGuiIO* self,float x,float y);
!x! void ImGuiIO_AddMouseButtonEvent(ImGuiIO* self,int button,bool down);
!x! void ImGuiIO_AddMouseWheelEvent(ImGuiIO* self,float wh_x,float wh_y);
!x! void ImGuiIO_AddMouseViewportEvent(ImGuiIO* self,ImGuiID id);
!x! void ImGuiIO_AddFocusEvent(ImGuiIO* self,bool focused);
!x! void ImGuiIO_AddInputCharacter(ImGuiIO* self,unsigned int c);
!x! void ImGuiIO_AddInputCharacterUTF16(ImGuiIO* self,ImWchar16 c);
!x! void ImGuiIO_AddInputCharactersUTF8(ImGuiIO* self,const char* str);
!x! void ImGuiIO_SetKeyEventNativeData(ImGuiIO* self,ImGuiKey key,int native_keycode,int native_scancode,int native_legacy_index);
!x! void ImGuiIO_SetAppAcceptingEvents(ImGuiIO* self,bool accepting_events);
!x! void ImGuiIO_ClearInputCharacters(ImGuiIO* self);
!x! void ImGuiIO_ClearInputKeys(ImGuiIO* self);
!x! ImGuiIO* ImGuiIO_ImGuiIO(void);
!x! void ImGuiIO_destroy(ImGuiIO* self);
!x! ImGuiInputTextCallbackData* ImGuiInputTextCallbackData_ImGuiInputTextCallbackData(void);
!x! void ImGuiInputTextCallbackData_destroy(ImGuiInputTextCallbackData* self);
!x! void ImGuiInputTextCallbackData_DeleteChars(ImGuiInputTextCallbackData* self,int pos,int bytes_count);
!x! void ImGuiInputTextCallbackData_InsertChars(ImGuiInputTextCallbackData* self,int pos,const char* text,const char* text_end);
!x! void ImGuiInputTextCallbackData_SelectAll(ImGuiInputTextCallbackData* self);
!x! void ImGuiInputTextCallbackData_ClearSelection(ImGuiInputTextCallbackData* self);
!x! bool ImGuiInputTextCallbackData_HasSelection(ImGuiInputTextCallbackData* self);
!x! ImGuiWindowClass* ImGuiWindowClass_ImGuiWindowClass(void);
!x! void ImGuiWindowClass_destroy(ImGuiWindowClass* self);
!x! ImGuiPayload* ImGuiPayload_ImGuiPayload(void);
!x! void ImGuiPayload_destroy(ImGuiPayload* self);
!x! void ImGuiPayload_Clear(ImGuiPayload* self);
!x! bool ImGuiPayload_IsDataType(ImGuiPayload* self,const char* type);
!x! bool ImGuiPayload_IsPreview(ImGuiPayload* self);
!x! bool ImGuiPayload_IsDelivery(ImGuiPayload* self);
!x! ImGuiTableColumnSortSpecs* ImGuiTableColumnSortSpecs_ImGuiTableColumnSortSpecs(void);
!x! void ImGuiTableColumnSortSpecs_destroy(ImGuiTableColumnSortSpecs* self);
!x! ImGuiTableSortSpecs* ImGuiTableSortSpecs_ImGuiTableSortSpecs(void);
!x! void ImGuiTableSortSpecs_destroy(ImGuiTableSortSpecs* self);
!x! ImGuiOnceUponAFrame* ImGuiOnceUponAFrame_ImGuiOnceUponAFrame(void);
!x! void ImGuiOnceUponAFrame_destroy(ImGuiOnceUponAFrame* self);
!x! ImGuiTextFilter* ImGuiTextFilter_ImGuiTextFilter(const char* default_filter);
!x! void ImGuiTextFilter_destroy(ImGuiTextFilter* self);
!x! bool ImGuiTextFilter_Draw(ImGuiTextFilter* self,const char* label,float width);
!x! bool ImGuiTextFilter_PassFilter(ImGuiTextFilter* self,const char* text,const char* text_end);
!x! void ImGuiTextFilter_Build(ImGuiTextFilter* self);
!x! void ImGuiTextFilter_Clear(ImGuiTextFilter* self);
!x! bool ImGuiTextFilter_IsActive(ImGuiTextFilter* self);
!x! ImGuiTextRange* ImGuiTextRange_ImGuiTextRange_Nil(void);
!x! void ImGuiTextRange_destroy(ImGuiTextRange* self);
!x! ImGuiTextRange* ImGuiTextRange_ImGuiTextRange_Str(const char* _b,const char* _e);
!x! bool ImGuiTextRange_empty(ImGuiTextRange* self);
!x! void ImGuiTextRange_split(ImGuiTextRange* self,char separator,ImVector_ImGuiTextRange* out);
!x! ImGuiTextBuffer* ImGuiTextBuffer_ImGuiTextBuffer(void);
!x! void ImGuiTextBuffer_destroy(ImGuiTextBuffer* self);
!x! const char* ImGuiTextBuffer_begin(ImGuiTextBuffer* self);
!x! const char* ImGuiTextBuffer_end(ImGuiTextBuffer* self);
!x! int ImGuiTextBuffer_size(ImGuiTextBuffer* self);
!x! bool ImGuiTextBuffer_empty(ImGuiTextBuffer* self);
!x! void ImGuiTextBuffer_clear(ImGuiTextBuffer* self);
!x! void ImGuiTextBuffer_reserve(ImGuiTextBuffer* self,int capacity);
!x! const char* ImGuiTextBuffer_c_str(ImGuiTextBuffer* self);
!x! void ImGuiTextBuffer_append(ImGuiTextBuffer* self,const char* str,const char* str_end);
!x! void ImGuiTextBuffer_appendfv(ImGuiTextBuffer* self,const char* fmt,va_list args);
!x! ImGuiStoragePair* ImGuiStoragePair_ImGuiStoragePair_Int(ImGuiID _key,int _val_i);
!x! void ImGuiStoragePair_destroy(ImGuiStoragePair* self);
!x! ImGuiStoragePair* ImGuiStoragePair_ImGuiStoragePair_Float(ImGuiID _key,float _val_f);
!x! ImGuiStoragePair* ImGuiStoragePair_ImGuiStoragePair_Ptr(ImGuiID _key,void* _val_p);
!x! void ImGuiStorage_Clear(ImGuiStorage* self);
!x! int ImGuiStorage_GetInt(ImGuiStorage* self,ImGuiID key,int default_val);
!x! void ImGuiStorage_SetInt(ImGuiStorage* self,ImGuiID key,int val);
!x! bool ImGuiStorage_GetBool(ImGuiStorage* self,ImGuiID key,bool default_val);
!x! void ImGuiStorage_SetBool(ImGuiStorage* self,ImGuiID key,bool val);
!x! float ImGuiStorage_GetFloat(ImGuiStorage* self,ImGuiID key,float default_val);
!x! void ImGuiStorage_SetFloat(ImGuiStorage* self,ImGuiID key,float val);
!x! void* ImGuiStorage_GetVoidPtr(ImGuiStorage* self,ImGuiID key);
!x! void ImGuiStorage_SetVoidPtr(ImGuiStorage* self,ImGuiID key,void* val);
!x! int* ImGuiStorage_GetIntRef(ImGuiStorage* self,ImGuiID key,int default_val);
!x! bool* ImGuiStorage_GetBoolRef(ImGuiStorage* self,ImGuiID key,bool default_val);
!x! float* ImGuiStorage_GetFloatRef(ImGuiStorage* self,ImGuiID key,float default_val);
!x! void** ImGuiStorage_GetVoidPtrRef(ImGuiStorage* self,ImGuiID key,void* default_val);
!x! void ImGuiStorage_SetAllInt(ImGuiStorage* self,int val);
!x! void ImGuiStorage_BuildSortByKey(ImGuiStorage* self);
!x! ImGuiListClipper* ImGuiListClipper_ImGuiListClipper(void);
!x! void ImGuiListClipper_destroy(ImGuiListClipper* self);
!x! void ImGuiListClipper_Begin(ImGuiListClipper* self,int items_count,float items_height);
!x! void ImGuiListClipper_End(ImGuiListClipper* self);
!x! bool ImGuiListClipper_Step(ImGuiListClipper* self);
!x! void ImGuiListClipper_ForceDisplayRangeByIndices(ImGuiListClipper* self,int item_min,int item_max);
!x! ImColor* ImColor_ImColor_Nil(void);
!x! void ImColor_destroy(ImColor* self);
!x! ImColor* ImColor_ImColor_Float(float r,float g,float b,float a);
!x! ImColor* ImColor_ImColor_Vec4(const ImVec4 col);
!x! ImColor* ImColor_ImColor_Int(int r,int g,int b,int a);
!x! ImColor* ImColor_ImColor_U32(ImU32 rgba);
!x! void ImColor_SetHSV(ImColor* self,float h,float s,float v,float a);
!x! void ImColor_HSV(ImColor *pOut,float h,float s,float v,float a);
!x! ImDrawCmd* ImDrawCmd_ImDrawCmd(void);
!x! void ImDrawCmd_destroy(ImDrawCmd* self);
!x! ImTextureID ImDrawCmd_GetTexID(ImDrawCmd* self);
!x! ImDrawListSplitter* ImDrawListSplitter_ImDrawListSplitter(void);
!x! void ImDrawListSplitter_destroy(ImDrawListSplitter* self);
!x! void ImDrawListSplitter_Clear(ImDrawListSplitter* self);
!x! void ImDrawListSplitter_ClearFreeMemory(ImDrawListSplitter* self);
!x! void ImDrawListSplitter_Split(ImDrawListSplitter* self,ImDrawList* draw_list,int count);
!x! void ImDrawListSplitter_Merge(ImDrawListSplitter* self,ImDrawList* draw_list);
!x! void ImDrawListSplitter_SetCurrentChannel(ImDrawListSplitter* self,ImDrawList* draw_list,int channel_idx);
!x! ImDrawList* ImDrawList_ImDrawList(const ImDrawListSharedData* shared_data);
!x! void ImDrawList_destroy(ImDrawList* self);
!x! void ImDrawList_PushClipRect(ImDrawList* self,const ImVec2 clip_rect_min,const ImVec2 clip_rect_max,bool intersect_with_current_clip_rect);
!x! void ImDrawList_PushClipRectFullScreen(ImDrawList* self);
!x! void ImDrawList_PopClipRect(ImDrawList* self);
!x! void ImDrawList_PushTextureID(ImDrawList* self,ImTextureID texture_id);
!x! void ImDrawList_PopTextureID(ImDrawList* self);
!x! void ImDrawList_GetClipRectMin(ImVec2 *pOut,ImDrawList* self);
!x! void ImDrawList_GetClipRectMax(ImVec2 *pOut,ImDrawList* self);
!x! void ImDrawList_AddLine(ImDrawList* self,const ImVec2 p1,const ImVec2 p2,ImU32 col,float thickness);
!x! void ImDrawList_AddRect(ImDrawList* self,const ImVec2 p_min,const ImVec2 p_max,ImU32 col,float rounding,ImDrawFlags flags,float thickness);
!x! void ImDrawList_AddRectFilled(ImDrawList* self,const ImVec2 p_min,const ImVec2 p_max,ImU32 col,float rounding,ImDrawFlags flags);
!x! void ImDrawList_AddRectFilledMultiColor(ImDrawList* self,const ImVec2 p_min,const ImVec2 p_max,ImU32 col_upr_left,ImU32 col_upr_right,ImU32 col_bot_right,ImU32 col_bot_left);
!x! void ImDrawList_AddQuad(ImDrawList* self,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,ImU32 col,float thickness);
!x! void ImDrawList_AddQuadFilled(ImDrawList* self,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,ImU32 col);
!x! void ImDrawList_AddTriangle(ImDrawList* self,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,ImU32 col,float thickness);
!x! void ImDrawList_AddTriangleFilled(ImDrawList* self,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,ImU32 col);
!x! void ImDrawList_AddCircle(ImDrawList* self,const ImVec2 center,float radius,ImU32 col,int num_segments,float thickness);
!x! void ImDrawList_AddCircleFilled(ImDrawList* self,const ImVec2 center,float radius,ImU32 col,int num_segments);
!x! void ImDrawList_AddNgon(ImDrawList* self,const ImVec2 center,float radius,ImU32 col,int num_segments,float thickness);
!x! void ImDrawList_AddNgonFilled(ImDrawList* self,const ImVec2 center,float radius,ImU32 col,int num_segments);
!x! void ImDrawList_AddText_Vec2(ImDrawList* self,const ImVec2 pos,ImU32 col,const char* text_begin,const char* text_end);
!x! void ImDrawList_AddText_FontPtr(ImDrawList* self,const ImFont* font,float font_size,const ImVec2 pos,ImU32 col,const char* text_begin,const char* text_end,float wrap_width,const ImVec4* cpu_fine_clip_rect);
!x! void ImDrawList_AddPolyline(ImDrawList* self,const ImVec2* points,int num_points,ImU32 col,ImDrawFlags flags,float thickness);
!x! void ImDrawList_AddConvexPolyFilled(ImDrawList* self,const ImVec2* points,int num_points,ImU32 col);
!x! void ImDrawList_AddBezierCubic(ImDrawList* self,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,ImU32 col,float thickness,int num_segments);
!x! void ImDrawList_AddBezierQuadratic(ImDrawList* self,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,ImU32 col,float thickness,int num_segments);
!x! void ImDrawList_AddImage(ImDrawList* self,ImTextureID user_texture_id,const ImVec2 p_min,const ImVec2 p_max,const ImVec2 uv_min,const ImVec2 uv_max,ImU32 col);
!x! void ImDrawList_AddImageQuad(ImDrawList* self,ImTextureID user_texture_id,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,const ImVec2 uv1,const ImVec2 uv2,const ImVec2 uv3,const ImVec2 uv4,ImU32 col);
!x! void ImDrawList_AddImageRounded(ImDrawList* self,ImTextureID user_texture_id,const ImVec2 p_min,const ImVec2 p_max,const ImVec2 uv_min,const ImVec2 uv_max,ImU32 col,float rounding,ImDrawFlags flags);
!x! void ImDrawList_PathClear(ImDrawList* self);
!x! void ImDrawList_PathLineTo(ImDrawList* self,const ImVec2 pos);
!x! void ImDrawList_PathLineToMergeDuplicate(ImDrawList* self,const ImVec2 pos);
!x! void ImDrawList_PathFillConvex(ImDrawList* self,ImU32 col);
!x! void ImDrawList_PathStroke(ImDrawList* self,ImU32 col,ImDrawFlags flags,float thickness);
!x! void ImDrawList_PathArcTo(ImDrawList* self,const ImVec2 center,float radius,float a_min,float a_max,int num_segments);
!x! void ImDrawList_PathArcToFast(ImDrawList* self,const ImVec2 center,float radius,int a_min_of_12,int a_max_of_12);
!x! void ImDrawList_PathBezierCubicCurveTo(ImDrawList* self,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,int num_segments);
!x! void ImDrawList_PathBezierQuadraticCurveTo(ImDrawList* self,const ImVec2 p2,const ImVec2 p3,int num_segments);
!x! void ImDrawList_PathRect(ImDrawList* self,const ImVec2 rect_min,const ImVec2 rect_max,float rounding,ImDrawFlags flags);
!x! void ImDrawList_AddCallback(ImDrawList* self,ImDrawCallback callback,void* callback_data);
!x! void ImDrawList_AddDrawCmd(ImDrawList* self);
!x! ImDrawList* ImDrawList_CloneOutput(ImDrawList* self);
!x! void ImDrawList_ChannelsSplit(ImDrawList* self,int count);
!x! void ImDrawList_ChannelsMerge(ImDrawList* self);
!x! void ImDrawList_ChannelsSetCurrent(ImDrawList* self,int n);
!x! void ImDrawList_PrimReserve(ImDrawList* self,int idx_count,int vtx_count);
!x! void ImDrawList_PrimUnreserve(ImDrawList* self,int idx_count,int vtx_count);
!x! void ImDrawList_PrimRect(ImDrawList* self,const ImVec2 a,const ImVec2 b,ImU32 col);
!x! void ImDrawList_PrimRectUV(ImDrawList* self,const ImVec2 a,const ImVec2 b,const ImVec2 uv_a,const ImVec2 uv_b,ImU32 col);
!x! void ImDrawList_PrimQuadUV(ImDrawList* self,const ImVec2 a,const ImVec2 b,const ImVec2 c,const ImVec2 d,const ImVec2 uv_a,const ImVec2 uv_b,const ImVec2 uv_c,const ImVec2 uv_d,ImU32 col);
!x! void ImDrawList_PrimWriteVtx(ImDrawList* self,const ImVec2 pos,const ImVec2 uv,ImU32 col);
!x! void ImDrawList_PrimWriteIdx(ImDrawList* self,ImDrawIdx idx);
!x! void ImDrawList_PrimVtx(ImDrawList* self,const ImVec2 pos,const ImVec2 uv,ImU32 col);
!x! void ImDrawList__ResetForNewFrame(ImDrawList* self);
!x! void ImDrawList__ClearFreeMemory(ImDrawList* self);
!x! void ImDrawList__PopUnusedDrawCmd(ImDrawList* self);
!x! void ImDrawList__TryMergeDrawCmds(ImDrawList* self);
!x! void ImDrawList__OnChangedClipRect(ImDrawList* self);
!x! void ImDrawList__OnChangedTextureID(ImDrawList* self);
!x! void ImDrawList__OnChangedVtxOffset(ImDrawList* self);
!x! int ImDrawList__CalcCircleAutoSegmentCount(ImDrawList* self,float radius);
!x! void ImDrawList__PathArcToFastEx(ImDrawList* self,const ImVec2 center,float radius,int a_min_sample,int a_max_sample,int a_step);
!x! void ImDrawList__PathArcToN(ImDrawList* self,const ImVec2 center,float radius,float a_min,float a_max,int num_segments);
!x! ImDrawData* ImDrawData_ImDrawData(void);
!x! void ImDrawData_destroy(ImDrawData* self);
!x! void ImDrawData_Clear(ImDrawData* self);
!x! void ImDrawData_DeIndexAllBuffers(ImDrawData* self);
!x! void ImDrawData_ScaleClipRects(ImDrawData* self,const ImVec2 fb_scale);
!x! ImFontConfig* ImFontConfig_ImFontConfig(void);
!x! void ImFontConfig_destroy(ImFontConfig* self);
!x! ImFontGlyphRangesBuilder* ImFontGlyphRangesBuilder_ImFontGlyphRangesBuilder(void);
!x! void ImFontGlyphRangesBuilder_destroy(ImFontGlyphRangesBuilder* self);
!x! void ImFontGlyphRangesBuilder_Clear(ImFontGlyphRangesBuilder* self);
!x! bool ImFontGlyphRangesBuilder_GetBit(ImFontGlyphRangesBuilder* self,size_t n);
!x! void ImFontGlyphRangesBuilder_SetBit(ImFontGlyphRangesBuilder* self,size_t n);
!x! void ImFontGlyphRangesBuilder_AddChar(ImFontGlyphRangesBuilder* self,ImWchar c);
!x! void ImFontGlyphRangesBuilder_AddText(ImFontGlyphRangesBuilder* self,const char* text,const char* text_end);
!x! void ImFontGlyphRangesBuilder_AddRanges(ImFontGlyphRangesBuilder* self,const ImWchar* ranges);
!x! void ImFontGlyphRangesBuilder_BuildRanges(ImFontGlyphRangesBuilder* self,ImVector_ImWchar* out_ranges);
!x! ImFontAtlasCustomRect* ImFontAtlasCustomRect_ImFontAtlasCustomRect(void);
!x! void ImFontAtlasCustomRect_destroy(ImFontAtlasCustomRect* self);
!x! bool ImFontAtlasCustomRect_IsPacked(ImFontAtlasCustomRect* self);
!x! ImFontAtlas* ImFontAtlas_ImFontAtlas(void);
!x! void ImFontAtlas_destroy(ImFontAtlas* self);
!x! ImFont* ImFontAtlas_AddFont(ImFontAtlas* self,const ImFontConfig* font_cfg);
!x! ImFont* ImFontAtlas_AddFontDefault(ImFontAtlas* self,const ImFontConfig* font_cfg);
!x! ImFont* ImFontAtlas_AddFontFromFileTTF(ImFontAtlas* self,const char* filename,float size_pixels,const ImFontConfig* font_cfg,const ImWchar* glyph_ranges);
!x! ImFont* ImFontAtlas_AddFontFromMemoryTTF(ImFontAtlas* self,void* font_data,int font_size,float size_pixels,const ImFontConfig* font_cfg,const ImWchar* glyph_ranges);
!x! ImFont* ImFontAtlas_AddFontFromMemoryCompressedTTF(ImFontAtlas* self,const void* compressed_font_data,int compressed_font_size,float size_pixels,const ImFontConfig* font_cfg,const ImWchar* glyph_ranges);
!x! ImFont* ImFontAtlas_AddFontFromMemoryCompressedBase85TTF(ImFontAtlas* self,const char* compressed_font_data_base85,float size_pixels,const ImFontConfig* font_cfg,const ImWchar* glyph_ranges);
!x! void ImFontAtlas_ClearInputData(ImFontAtlas* self);
!x! void ImFontAtlas_ClearTexData(ImFontAtlas* self);
!x! void ImFontAtlas_ClearFonts(ImFontAtlas* self);
!x! void ImFontAtlas_Clear(ImFontAtlas* self);
!x! bool ImFontAtlas_Build(ImFontAtlas* self);
!x! void ImFontAtlas_GetTexDataAsAlpha8(ImFontAtlas* self,unsigned char** out_pixels,int* out_width,int* out_height,int* out_bytes_per_pixel);
!x! void ImFontAtlas_GetTexDataAsRGBA32(ImFontAtlas* self,unsigned char** out_pixels,int* out_width,int* out_height,int* out_bytes_per_pixel);
!x! bool ImFontAtlas_IsBuilt(ImFontAtlas* self);
!x! void ImFontAtlas_SetTexID(ImFontAtlas* self,ImTextureID id);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesDefault(ImFontAtlas* self);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesKorean(ImFontAtlas* self);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesJapanese(ImFontAtlas* self);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesChineseFull(ImFontAtlas* self);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesChineseSimplifiedCommon(ImFontAtlas* self);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesCyrillic(ImFontAtlas* self);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesThai(ImFontAtlas* self);
!x! const ImWchar* ImFontAtlas_GetGlyphRangesVietnamese(ImFontAtlas* self);
!x! int ImFontAtlas_AddCustomRectRegular(ImFontAtlas* self,int width,int height);
!x! int ImFontAtlas_AddCustomRectFontGlyph(ImFontAtlas* self,ImFont* font,ImWchar id,int width,int height,float advance_x,const ImVec2 offset);
!x! ImFontAtlasCustomRect* ImFontAtlas_GetCustomRectByIndex(ImFontAtlas* self,int index);
!x! void ImFontAtlas_CalcCustomRectUV(ImFontAtlas* self,const ImFontAtlasCustomRect* rect,ImVec2* out_uv_min,ImVec2* out_uv_max);
!x! bool ImFontAtlas_GetMouseCursorTexData(ImFontAtlas* self,ImGuiMouseCursor cursor,ImVec2* out_offset,ImVec2* out_size,ImVec2 out_uv_border[2],ImVec2 out_uv_fill[2]);
!x! ImFont* ImFont_ImFont(void);
!x! void ImFont_destroy(ImFont* self);
!x! const ImFontGlyph* ImFont_FindGlyph(ImFont* self,ImWchar c);
!x! const ImFontGlyph* ImFont_FindGlyphNoFallback(ImFont* self,ImWchar c);
!x! float ImFont_GetCharAdvance(ImFont* self,ImWchar c);
!x! bool ImFont_IsLoaded(ImFont* self);
!x! const char* ImFont_GetDebugName(ImFont* self);
!x! void ImFont_CalcTextSizeA(ImVec2 *pOut,ImFont* self,float size,float max_width,float wrap_width,const char* text_begin,const char* text_end,const char** remaining);
!x! const char* ImFont_CalcWordWrapPositionA(ImFont* self,float scale,const char* text,const char* text_end,float wrap_width);
!x! void ImFont_RenderChar(ImFont* self,ImDrawList* draw_list,float size,const ImVec2 pos,ImU32 col,ImWchar c);
!x! void ImFont_RenderText(ImFont* self,ImDrawList* draw_list,float size,const ImVec2 pos,ImU32 col,const ImVec4 clip_rect,const char* text_begin,const char* text_end,float wrap_width,bool cpu_fine_clip);
!x! void ImFont_BuildLookupTable(ImFont* self);
!x! void ImFont_ClearOutputData(ImFont* self);
!x! void ImFont_GrowIndex(ImFont* self,int new_size);
!x! void ImFont_AddGlyph(ImFont* self,const ImFontConfig* src_cfg,ImWchar c,float x0,float y0,float x1,float y1,float u0,float v0,float u1,float v1,float advance_x);
!x! void ImFont_AddRemapChar(ImFont* self,ImWchar dst,ImWchar src,bool overwrite_dst);
!x! void ImFont_SetGlyphVisible(ImFont* self,ImWchar c,bool visible);
!x! bool ImFont_IsGlyphRangeUnused(ImFont* self,unsigned int c_begin,unsigned int c_last);
!x! ImGuiViewport* ImGuiViewport_ImGuiViewport(void);
!x! void ImGuiViewport_destroy(ImGuiViewport* self);
!x! void ImGuiViewport_GetCenter(ImVec2 *pOut,ImGuiViewport* self);
!x! void ImGuiViewport_GetWorkCenter(ImVec2 *pOut,ImGuiViewport* self);
!x! ImGuiPlatformIO* ImGuiPlatformIO_ImGuiPlatformIO(void);
!x! void ImGuiPlatformIO_destroy(ImGuiPlatformIO* self);
!x! ImGuiPlatformMonitor* ImGuiPlatformMonitor_ImGuiPlatformMonitor(void);
!x! void ImGuiPlatformMonitor_destroy(ImGuiPlatformMonitor* self);
!x! ImGuiPlatformImeData* ImGuiPlatformImeData_ImGuiPlatformImeData(void);
!x! void ImGuiPlatformImeData_destroy(ImGuiPlatformImeData* self);
!x! int igGetKeyIndex(ImGuiKey key);
!x! ImGuiID igImHashData(const void* data,size_t data_size,ImU32 seed);
!x! ImGuiID igImHashStr(const char* data,size_t data_size,ImU32 seed);
!x! void igImQsort(void* base,size_t count,size_t size_of_element,int(*compare_func)(void const*,void const*));
!x! ImU32 igImAlphaBlendColors(ImU32 col_a,ImU32 col_b);
!x! bool igImIsPowerOfTwo_Int(int v);
!x! bool igImIsPowerOfTwo_U64(ImU64 v);
!x! int igImUpperPowerOfTwo(int v);
!x! int igImStricmp(const char* str1,const char* str2);
!x! int igImStrnicmp(const char* str1,const char* str2,size_t count);
!x! void igImStrncpy(char* dst,const char* src,size_t count);
!x! char* igImStrdup(const char* str);
!x! char* igImStrdupcpy(char* dst,size_t* p_dst_size,const char* str);
!x! const char* igImStrchrRange(const char* str_begin,const char* str_end,char c);
!x! int igImStrlenW(const ImWchar* str);
!x! const char* igImStreolRange(const char* str,const char* str_end);
!x! const ImWchar* igImStrbolW(const ImWchar* buf_mid_line,const ImWchar* buf_begin);
!x! const char* igImStristr(const char* haystack,const char* haystack_end,const char* needle,const char* needle_end);
!x! void igImStrTrimBlanks(char* str);
!x! const char* igImStrSkipBlank(const char* str);
!x! bool igImCharIsBlankA(char c);
!x! bool igImCharIsBlankW(unsigned int c);
!x! int igImFormatString(char* buf,size_t buf_size,const char* fmt,...);
!x! int igImFormatStringV(char* buf,size_t buf_size,const char* fmt,va_list args);
!x! void igImFormatStringToTempBuffer(const char** out_buf,const char** out_buf_end,const char* fmt,...);
!x! void igImFormatStringToTempBufferV(const char** out_buf,const char** out_buf_end,const char* fmt,va_list args);
!x! const char* igImParseFormatFindStart(const char* format);
!x! const char* igImParseFormatFindEnd(const char* format);
!x! const char* igImParseFormatTrimDecorations(const char* format,char* buf,size_t buf_size);
!x! void igImParseFormatSanitizeForPrinting(const char* fmt_in,char* fmt_out,size_t fmt_out_size);
!x! const char* igImParseFormatSanitizeForScanning(const char* fmt_in,char* fmt_out,size_t fmt_out_size);
!x! int igImParseFormatPrecision(const char* format,int default_value);
!x! const char* igImTextCharToUtf8(char out_buf[5],unsigned int c);
!x! int igImTextStrToUtf8(char* out_buf,int out_buf_size,const ImWchar* in_text,const ImWchar* in_text_end);
!x! int igImTextCharFromUtf8(unsigned int* out_char,const char* in_text,const char* in_text_end);
!x! int igImTextStrFromUtf8(ImWchar* out_buf,int out_buf_size,const char* in_text,const char* in_text_end,const char** in_remaining);
!x! int igImTextCountCharsFromUtf8(const char* in_text,const char* in_text_end);
!x! int igImTextCountUtf8BytesFromChar(const char* in_text,const char* in_text_end);
!x! int igImTextCountUtf8BytesFromStr(const ImWchar* in_text,const ImWchar* in_text_end);
!x! ImFileHandle igImFileOpen(const char* filename,const char* mode);
!x! bool igImFileClose(ImFileHandle file);
!x! ImU64 igImFileGetSize(ImFileHandle file);
!x! ImU64 igImFileRead(void* data,ImU64 size,ImU64 count,ImFileHandle file);
!x! ImU64 igImFileWrite(const void* data,ImU64 size,ImU64 count,ImFileHandle file);
!x! void* igImFileLoadToMemory(const char* filename,const char* mode,size_t* out_file_size,int padding_bytes);
!x! float igImPow_Float(float x,float y);
!x! double igImPow_double(double x,double y);
!x! float igImLog_Float(float x);
!x! double igImLog_double(double x);
!x! int igImAbs_Int(int x);
!x! float igImAbs_Float(float x);
!x! double igImAbs_double(double x);
!x! float igImSign_Float(float x);
!x! double igImSign_double(double x);
!x! float igImRsqrt_Float(float x);
!x! double igImRsqrt_double(double x);
!x! void igImMin(ImVec2 *pOut,const ImVec2 lhs,const ImVec2 rhs);
!x! void igImMax(ImVec2 *pOut,const ImVec2 lhs,const ImVec2 rhs);
!x! void igImClamp(ImVec2 *pOut,const ImVec2 v,const ImVec2 mn,ImVec2 mx);
!x! void igImLerp_Vec2Float(ImVec2 *pOut,const ImVec2 a,const ImVec2 b,float t);
!x! void igImLerp_Vec2Vec2(ImVec2 *pOut,const ImVec2 a,const ImVec2 b,const ImVec2 t);
!x! void igImLerp_Vec4(ImVec4 *pOut,const ImVec4 a,const ImVec4 b,float t);
!x! float igImSaturate(float f);
!x! float igImLengthSqr_Vec2(const ImVec2 lhs);
!x! float igImLengthSqr_Vec4(const ImVec4 lhs);
!x! float igImInvLength(const ImVec2 lhs,float fail_value);
!x! float igImFloor_Float(float f);
!x! float igImFloorSigned_Float(float f);
!x! void igImFloor_Vec2(ImVec2 *pOut,const ImVec2 v);
!x! void igImFloorSigned_Vec2(ImVec2 *pOut,const ImVec2 v);
!x! int igImModPositive(int a,int b);
!x! float igImDot(const ImVec2 a,const ImVec2 b);
!x! void igImRotate(ImVec2 *pOut,const ImVec2 v,float cos_a,float sin_a);
!x! float igImLinearSweep(float current,float target,float speed);
!x! void igImMul(ImVec2 *pOut,const ImVec2 lhs,const ImVec2 rhs);
!x! bool igImIsFloatAboveGuaranteedIntegerPrecision(float f);
!x! void igImBezierCubicCalc(ImVec2 *pOut,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,float t);
!x! void igImBezierCubicClosestPoint(ImVec2 *pOut,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,const ImVec2 p,int num_segments);
!x! void igImBezierCubicClosestPointCasteljau(ImVec2 *pOut,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,const ImVec2 p4,const ImVec2 p,float tess_tol);
!x! void igImBezierQuadraticCalc(ImVec2 *pOut,const ImVec2 p1,const ImVec2 p2,const ImVec2 p3,float t);
!x! void igImLineClosestPoint(ImVec2 *pOut,const ImVec2 a,const ImVec2 b,const ImVec2 p);
!x! bool igImTriangleContainsPoint(const ImVec2 a,const ImVec2 b,const ImVec2 c,const ImVec2 p);
!x! void igImTriangleClosestPoint(ImVec2 *pOut,const ImVec2 a,const ImVec2 b,const ImVec2 c,const ImVec2 p);
!x! void igImTriangleBarycentricCoords(const ImVec2 a,const ImVec2 b,const ImVec2 c,const ImVec2 p,float* out_u,float* out_v,float* out_w);
!x! float igImTriangleArea(const ImVec2 a,const ImVec2 b,const ImVec2 c);
!x! ImGuiDir igImGetDirQuadrantFromDelta(float dx,float dy);
!x! ImVec1* ImVec1_ImVec1_Nil(void);
!x! void ImVec1_destroy(ImVec1* self);
!x! ImVec1* ImVec1_ImVec1_Float(float _x);
!x! ImVec2ih* ImVec2ih_ImVec2ih_Nil(void);
!x! void ImVec2ih_destroy(ImVec2ih* self);
!x! ImVec2ih* ImVec2ih_ImVec2ih_short(short _x,short _y);
!x! ImVec2ih* ImVec2ih_ImVec2ih_Vec2(const ImVec2 rhs);
!x! ImRect* ImRect_ImRect_Nil(void);
!x! void ImRect_destroy(ImRect* self);
!x! ImRect* ImRect_ImRect_Vec2(const ImVec2 min,const ImVec2 max);
!x! ImRect* ImRect_ImRect_Vec4(const ImVec4 v);
!x! ImRect* ImRect_ImRect_Float(float x1,float y1,float x2,float y2);
!x! void ImRect_GetCenter(ImVec2 *pOut,ImRect* self);
!x! void ImRect_GetSize(ImVec2 *pOut,ImRect* self);
!x! float ImRect_GetWidth(ImRect* self);
!x! float ImRect_GetHeight(ImRect* self);
!x! float ImRect_GetArea(ImRect* self);
!x! void ImRect_GetTL(ImVec2 *pOut,ImRect* self);
!x! void ImRect_GetTR(ImVec2 *pOut,ImRect* self);
!x! void ImRect_GetBL(ImVec2 *pOut,ImRect* self);
!x! void ImRect_GetBR(ImVec2 *pOut,ImRect* self);
!x! bool ImRect_Contains_Vec2(ImRect* self,const ImVec2 p);
!x! bool ImRect_Contains_Rect(ImRect* self,const ImRect r);
!x! bool ImRect_Overlaps(ImRect* self,const ImRect r);
!x! void ImRect_Add_Vec2(ImRect* self,const ImVec2 p);
!x! void ImRect_Add_Rect(ImRect* self,const ImRect r);
!x! void ImRect_Expand_Float(ImRect* self,const float amount);
!x! void ImRect_Expand_Vec2(ImRect* self,const ImVec2 amount);
!x! void ImRect_Translate(ImRect* self,const ImVec2 d);
!x! void ImRect_TranslateX(ImRect* self,float dx);
!x! void ImRect_TranslateY(ImRect* self,float dy);
!x! void ImRect_ClipWith(ImRect* self,const ImRect r);
!x! void ImRect_ClipWithFull(ImRect* self,const ImRect r);
!x! void ImRect_Floor(ImRect* self);
!x! bool ImRect_IsInverted(ImRect* self);
!x! void ImRect_ToVec4(ImVec4 *pOut,ImRect* self);
!x! bool igImBitArrayTestBit(const ImU32* arr,int n);
!x! void igImBitArrayClearBit(ImU32* arr,int n);
!x! void igImBitArraySetBit(ImU32* arr,int n);
!x! void igImBitArraySetBitRange(ImU32* arr,int n,int n2);
!x! void ImBitVector_Create(ImBitVector* self,int sz);
!x! void ImBitVector_Clear(ImBitVector* self);
!x! bool ImBitVector_TestBit(ImBitVector* self,int n);
!x! void ImBitVector_SetBit(ImBitVector* self,int n);
!x! void ImBitVector_ClearBit(ImBitVector* self,int n);
!x! ImDrawListSharedData* ImDrawListSharedData_ImDrawListSharedData(void);
!x! void ImDrawListSharedData_destroy(ImDrawListSharedData* self);
!x! void ImDrawListSharedData_SetCircleTessellationMaxError(ImDrawListSharedData* self,float max_error);
!x! void ImDrawDataBuilder_Clear(ImDrawDataBuilder* self);
!x! void ImDrawDataBuilder_ClearFreeMemory(ImDrawDataBuilder* self);
!x! int ImDrawDataBuilder_GetDrawListCount(ImDrawDataBuilder* self);
!x! void ImDrawDataBuilder_FlattenIntoSingleLayer(ImDrawDataBuilder* self);
!x! ImGuiStyleMod* ImGuiStyleMod_ImGuiStyleMod_Int(ImGuiStyleVar idx,int v);
!x! void ImGuiStyleMod_destroy(ImGuiStyleMod* self);
!x! ImGuiStyleMod* ImGuiStyleMod_ImGuiStyleMod_Float(ImGuiStyleVar idx,float v);
!x! ImGuiStyleMod* ImGuiStyleMod_ImGuiStyleMod_Vec2(ImGuiStyleVar idx,ImVec2 v);
!x! ImGuiComboPreviewData* ImGuiComboPreviewData_ImGuiComboPreviewData(void);
!x! void ImGuiComboPreviewData_destroy(ImGuiComboPreviewData* self);
!x! ImGuiMenuColumns* ImGuiMenuColumns_ImGuiMenuColumns(void);
!x! void ImGuiMenuColumns_destroy(ImGuiMenuColumns* self);
!x! void ImGuiMenuColumns_Update(ImGuiMenuColumns* self,float spacing,bool window_reappearing);
!x! float ImGuiMenuColumns_DeclColumns(ImGuiMenuColumns* self,float w_icon,float w_label,float w_shortcut,float w_mark);
!x! void ImGuiMenuColumns_CalcNextTotalWidth(ImGuiMenuColumns* self,bool update_offsets);
!x! ImGuiInputTextState* ImGuiInputTextState_ImGuiInputTextState(void);
!x! void ImGuiInputTextState_destroy(ImGuiInputTextState* self);
!x! void ImGuiInputTextState_ClearText(ImGuiInputTextState* self);
!x! void ImGuiInputTextState_ClearFreeMemory(ImGuiInputTextState* self);
!x! int ImGuiInputTextState_GetUndoAvailCount(ImGuiInputTextState* self);
!x! int ImGuiInputTextState_GetRedoAvailCount(ImGuiInputTextState* self);
!x! void ImGuiInputTextState_OnKeyPressed(ImGuiInputTextState* self,int key);
!x! void ImGuiInputTextState_CursorAnimReset(ImGuiInputTextState* self);
!x! void ImGuiInputTextState_CursorClamp(ImGuiInputTextState* self);
!x! bool ImGuiInputTextState_HasSelection(ImGuiInputTextState* self);
!x! void ImGuiInputTextState_ClearSelection(ImGuiInputTextState* self);
!x! int ImGuiInputTextState_GetCursorPos(ImGuiInputTextState* self);
!x! int ImGuiInputTextState_GetSelectionStart(ImGuiInputTextState* self);
!x! int ImGuiInputTextState_GetSelectionEnd(ImGuiInputTextState* self);
!x! void ImGuiInputTextState_SelectAll(ImGuiInputTextState* self);
!x! ImGuiPopupData* ImGuiPopupData_ImGuiPopupData(void);
!x! void ImGuiPopupData_destroy(ImGuiPopupData* self);
!x! ImGuiNextWindowData* ImGuiNextWindowData_ImGuiNextWindowData(void);
!x! void ImGuiNextWindowData_destroy(ImGuiNextWindowData* self);
!x! void ImGuiNextWindowData_ClearFlags(ImGuiNextWindowData* self);
!x! ImGuiNextItemData* ImGuiNextItemData_ImGuiNextItemData(void);
!x! void ImGuiNextItemData_destroy(ImGuiNextItemData* self);
!x! void ImGuiNextItemData_ClearFlags(ImGuiNextItemData* self);
!x! ImGuiLastItemData* ImGuiLastItemData_ImGuiLastItemData(void);
!x! void ImGuiLastItemData_destroy(ImGuiLastItemData* self);
!x! ImGuiStackSizes* ImGuiStackSizes_ImGuiStackSizes(void);
!x! void ImGuiStackSizes_destroy(ImGuiStackSizes* self);
!x! void ImGuiStackSizes_SetToCurrentState(ImGuiStackSizes* self);
!x! void ImGuiStackSizes_CompareWithCurrentState(ImGuiStackSizes* self);
!x! ImGuiPtrOrIndex* ImGuiPtrOrIndex_ImGuiPtrOrIndex_Ptr(void* ptr);
!x! void ImGuiPtrOrIndex_destroy(ImGuiPtrOrIndex* self);
!x! ImGuiPtrOrIndex* ImGuiPtrOrIndex_ImGuiPtrOrIndex_Int(int index);
!x! ImGuiInputEvent* ImGuiInputEvent_ImGuiInputEvent(void);
!x! void ImGuiInputEvent_destroy(ImGuiInputEvent* self);
!x! ImGuiListClipperRange ImGuiListClipperRange_FromIndices(int min,int max);
!x! ImGuiListClipperRange ImGuiListClipperRange_FromPositions(float y1,float y2,int off_min,int off_max);
!x! ImGuiListClipperData* ImGuiListClipperData_ImGuiListClipperData(void);
!x! void ImGuiListClipperData_destroy(ImGuiListClipperData* self);
!x! void ImGuiListClipperData_Reset(ImGuiListClipperData* self,ImGuiListClipper* clipper);
!x! ImGuiNavItemData* ImGuiNavItemData_ImGuiNavItemData(void);
!x! void ImGuiNavItemData_destroy(ImGuiNavItemData* self);
!x! void ImGuiNavItemData_Clear(ImGuiNavItemData* self);
!x! ImGuiOldColumnData* ImGuiOldColumnData_ImGuiOldColumnData(void);
!x! void ImGuiOldColumnData_destroy(ImGuiOldColumnData* self);
!x! ImGuiOldColumns* ImGuiOldColumns_ImGuiOldColumns(void);
!x! void ImGuiOldColumns_destroy(ImGuiOldColumns* self);
!x! ImGuiDockNode* ImGuiDockNode_ImGuiDockNode(ImGuiID id);
!x! void ImGuiDockNode_destroy(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsRootNode(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsDockSpace(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsFloatingNode(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsCentralNode(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsHiddenTabBar(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsNoTabBar(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsSplitNode(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsLeafNode(ImGuiDockNode* self);
!x! bool ImGuiDockNode_IsEmpty(ImGuiDockNode* self);
!x! void ImGuiDockNode_Rect(ImRect *pOut,ImGuiDockNode* self);
!x! void ImGuiDockNode_SetLocalFlags(ImGuiDockNode* self,ImGuiDockNodeFlags flags);
!x! void ImGuiDockNode_UpdateMergedFlags(ImGuiDockNode* self);
!x! ImGuiDockContext* ImGuiDockContext_ImGuiDockContext(void);
!x! void ImGuiDockContext_destroy(ImGuiDockContext* self);
!x! ImGuiViewportP* ImGuiViewportP_ImGuiViewportP(void);
!x! void ImGuiViewportP_destroy(ImGuiViewportP* self);
!x! void ImGuiViewportP_ClearRequestFlags(ImGuiViewportP* self);
!x! void ImGuiViewportP_CalcWorkRectPos(ImVec2 *pOut,ImGuiViewportP* self,const ImVec2 off_min);
!x! void ImGuiViewportP_CalcWorkRectSize(ImVec2 *pOut,ImGuiViewportP* self,const ImVec2 off_min,const ImVec2 off_max);
!x! void ImGuiViewportP_UpdateWorkRect(ImGuiViewportP* self);
!x! void ImGuiViewportP_GetMainRect(ImRect *pOut,ImGuiViewportP* self);
!x! void ImGuiViewportP_GetWorkRect(ImRect *pOut,ImGuiViewportP* self);
!x! void ImGuiViewportP_GetBuildWorkRect(ImRect *pOut,ImGuiViewportP* self);
!x! ImGuiWindowSettings* ImGuiWindowSettings_ImGuiWindowSettings(void);
!x! void ImGuiWindowSettings_destroy(ImGuiWindowSettings* self);
!x! char* ImGuiWindowSettings_GetName(ImGuiWindowSettings* self);
!x! ImGuiSettingsHandler* ImGuiSettingsHandler_ImGuiSettingsHandler(void);
!x! void ImGuiSettingsHandler_destroy(ImGuiSettingsHandler* self);
!x! ImGuiMetricsConfig* ImGuiMetricsConfig_ImGuiMetricsConfig(void);
!x! void ImGuiMetricsConfig_destroy(ImGuiMetricsConfig* self);
!x! ImGuiStackLevelInfo* ImGuiStackLevelInfo_ImGuiStackLevelInfo(void);
!x! void ImGuiStackLevelInfo_destroy(ImGuiStackLevelInfo* self);
!x! ImGuiStackTool* ImGuiStackTool_ImGuiStackTool(void);
!x! void ImGuiStackTool_destroy(ImGuiStackTool* self);
!x! ImGuiContextHook* ImGuiContextHook_ImGuiContextHook(void);
!x! void ImGuiContextHook_destroy(ImGuiContextHook* self);
!x! ImGuiContext* ImGuiContext_ImGuiContext(ImFontAtlas* shared_font_atlas);
!x! void ImGuiContext_destroy(ImGuiContext* self);
!x! ImGuiWindow* ImGuiWindow_ImGuiWindow(ImGuiContext* context,const char* name);
!x! void ImGuiWindow_destroy(ImGuiWindow* self);
!x! ImGuiID ImGuiWindow_GetID_Str(ImGuiWindow* self,const char* str,const char* str_end);
!x! ImGuiID ImGuiWindow_GetID_Ptr(ImGuiWindow* self,const void* ptr);
!x! ImGuiID ImGuiWindow_GetID_Int(ImGuiWindow* self,int n);
!x! ImGuiID ImGuiWindow_GetIDFromRectangle(ImGuiWindow* self,const ImRect r_abs);
!x! void ImGuiWindow_Rect(ImRect *pOut,ImGuiWindow* self);
!x! float ImGuiWindow_CalcFontSize(ImGuiWindow* self);
!x! float ImGuiWindow_TitleBarHeight(ImGuiWindow* self);
!x! void ImGuiWindow_TitleBarRect(ImRect *pOut,ImGuiWindow* self);
!x! float ImGuiWindow_MenuBarHeight(ImGuiWindow* self);
!x! void ImGuiWindow_MenuBarRect(ImRect *pOut,ImGuiWindow* self);
!x! ImGuiTabItem* ImGuiTabItem_ImGuiTabItem(void);
!x! void ImGuiTabItem_destroy(ImGuiTabItem* self);
!x! ImGuiTabBar* ImGuiTabBar_ImGuiTabBar(void);
!x! void ImGuiTabBar_destroy(ImGuiTabBar* self);
!x! int ImGuiTabBar_GetTabOrder(ImGuiTabBar* self,const ImGuiTabItem* tab);
!x! const char* ImGuiTabBar_GetTabName(ImGuiTabBar* self,const ImGuiTabItem* tab);
!x! ImGuiTableColumn* ImGuiTableColumn_ImGuiTableColumn(void);
!x! void ImGuiTableColumn_destroy(ImGuiTableColumn* self);
!x! ImGuiTableInstanceData* ImGuiTableInstanceData_ImGuiTableInstanceData(void);
!x! void ImGuiTableInstanceData_destroy(ImGuiTableInstanceData* self);
!x! ImGuiTable* ImGuiTable_ImGuiTable(void);
!x! void ImGuiTable_destroy(ImGuiTable* self);
!x! ImGuiTableTempData* ImGuiTableTempData_ImGuiTableTempData(void);
!x! void ImGuiTableTempData_destroy(ImGuiTableTempData* self);
!x! ImGuiTableColumnSettings* ImGuiTableColumnSettings_ImGuiTableColumnSettings(void);
!x! void ImGuiTableColumnSettings_destroy(ImGuiTableColumnSettings* self);
!x! ImGuiTableSettings* ImGuiTableSettings_ImGuiTableSettings(void);
!x! void ImGuiTableSettings_destroy(ImGuiTableSettings* self);
!x! ImGuiTableColumnSettings* ImGuiTableSettings_GetColumnSettings(ImGuiTableSettings* self);
!x! ImGuiWindow* igGetCurrentWindowRead(void);
! ImGuiWindow* igGetCurrentWindow(void);
function igGetCurrentWindow() bind(c,name="igGetCurrentWindow")
  import c_ptr
  type(c_ptr) :: igGetCurrentWindow
end function igGetCurrentWindow
!x! ImGuiWindow* igFindWindowByID(ImGuiID id);
!x! ImGuiWindow* igFindWindowByName(const char* name);
!x! void igUpdateWindowParentAndRootLinks(ImGuiWindow* window,ImGuiWindowFlags flags,ImGuiWindow* parent_window);
!x! void igCalcWindowNextAutoFitSize(ImVec2 *pOut,ImGuiWindow* window);
!x! bool igIsWindowChildOf(ImGuiWindow* window,ImGuiWindow* potential_parent,bool popup_hierarchy,bool dock_hierarchy);
!x! bool igIsWindowWithinBeginStackOf(ImGuiWindow* window,ImGuiWindow* potential_parent);
!x! bool igIsWindowAbove(ImGuiWindow* potential_above,ImGuiWindow* potential_below);
!x! bool igIsWindowNavFocusable(ImGuiWindow* window);
!x! void igSetWindowPos_WindowPtr(ImGuiWindow* window,const ImVec2 pos,ImGuiCond cond);
!x! void igSetWindowSize_WindowPtr(ImGuiWindow* window,const ImVec2 size,ImGuiCond cond);
!x! void igSetWindowCollapsed_WindowPtr(ImGuiWindow* window,bool collapsed,ImGuiCond cond);
!x! void igSetWindowHitTestHole(ImGuiWindow* window,const ImVec2 pos,const ImVec2 size);
!x! void igWindowRectAbsToRel(ImRect *pOut,ImGuiWindow* window,const ImRect r);
!x! void igWindowRectRelToAbs(ImRect *pOut,ImGuiWindow* window,const ImRect r);
!x! void igFocusWindow(ImGuiWindow* window);
!x! void igFocusTopMostWindowUnderOne(ImGuiWindow* under_this_window,ImGuiWindow* ignore_window);
!x! void igBringWindowToFocusFront(ImGuiWindow* window);
!x! void igBringWindowToDisplayFront(ImGuiWindow* window);
!x! void igBringWindowToDisplayBack(ImGuiWindow* window);
!x! void igBringWindowToDisplayBehind(ImGuiWindow* window,ImGuiWindow* above_window);
!x! int igFindWindowDisplayIndex(ImGuiWindow* window);
!x! ImGuiWindow* igFindBottomMostVisibleWindowWithinBeginStack(ImGuiWindow* window);
!x! void igSetCurrentFont(ImFont* font);
!x! ImFont* igGetDefaultFont(void);
!x! ImDrawList* igGetForegroundDrawList_WindowPtr(ImGuiWindow* window);
!x! void igInitialize(void);
!x! void igShutdown(void);
!x! void igUpdateInputEvents(bool trickle_fast_inputs);
!x! void igUpdateHoveredWindowAndCaptureFlags(void);
!x! void igStartMouseMovingWindow(ImGuiWindow* window);
!x! void igStartMouseMovingWindowOrNode(ImGuiWindow* window,ImGuiDockNode* node,bool undock_floating_node);
!x! void igUpdateMouseMovingWindowNewFrame(void);
!x! void igUpdateMouseMovingWindowEndFrame(void);
!x! ImGuiID igAddContextHook(ImGuiContext* context,const ImGuiContextHook* hook);
!x! void igRemoveContextHook(ImGuiContext* context,ImGuiID hook_to_remove);
!x! void igCallContextHooks(ImGuiContext* context,ImGuiContextHookType type);
!x! void igTranslateWindowsInViewport(ImGuiViewportP* viewport,const ImVec2 old_pos,const ImVec2 new_pos);
!x! void igScaleWindowsInViewport(ImGuiViewportP* viewport,float scale);
!x! void igDestroyPlatformWindow(ImGuiViewportP* viewport);
!x! void igSetWindowViewport(ImGuiWindow* window,ImGuiViewportP* viewport);
!x! void igSetCurrentViewport(ImGuiWindow* window,ImGuiViewportP* viewport);
!x! const ImGuiPlatformMonitor* igGetViewportPlatformMonitor(ImGuiViewport* viewport);
!x! ImGuiViewportP* igFindHoveredViewportFromPlatformWindowStack(const ImVec2 mouse_platform_pos);
!x! void igMarkIniSettingsDirty_Nil(void);
!x! void igMarkIniSettingsDirty_WindowPtr(ImGuiWindow* window);
!x! void igClearIniSettings(void);
!x! ImGuiWindowSettings* igCreateNewWindowSettings(const char* name);
!x! ImGuiWindowSettings* igFindWindowSettings(ImGuiID id);
!x! ImGuiWindowSettings* igFindOrCreateWindowSettings(const char* name);
!x! void igAddSettingsHandler(const ImGuiSettingsHandler* handler);
!x! void igRemoveSettingsHandler(const char* type_name);
!x! ImGuiSettingsHandler* igFindSettingsHandler(const char* type_name);
!x! void igSetNextWindowScroll(const ImVec2 scroll);
!x! void igSetScrollX_WindowPtr(ImGuiWindow* window,float scroll_x);
!x! void igSetScrollY_WindowPtr(ImGuiWindow* window,float scroll_y);
!x! void igSetScrollFromPosX_WindowPtr(ImGuiWindow* window,float local_x,float center_x_ratio);
!x! void igSetScrollFromPosY_WindowPtr(ImGuiWindow* window,float local_y,float center_y_ratio);
!x! void igScrollToItem(ImGuiScrollFlags flags);
!x! void igScrollToRect(ImGuiWindow* window,const ImRect rect,ImGuiScrollFlags flags);
!x! void igScrollToRectEx(ImVec2 *pOut,ImGuiWindow* window,const ImRect rect,ImGuiScrollFlags flags);
!x! void igScrollToBringRectIntoView(ImGuiWindow* window,const ImRect rect);
!x! ImGuiID igGetItemID(void);
!x! ImGuiItemStatusFlags igGetItemStatusFlags(void);
!x! ImGuiItemFlags igGetItemFlags(void);
!x! ImGuiID igGetActiveID(void);
!x! ImGuiID igGetFocusID(void);
!x! void igSetActiveID(ImGuiID id,ImGuiWindow* window);
!x! void igSetFocusID(ImGuiID id,ImGuiWindow* window);
!x! void igClearActiveID(void);
!x! ImGuiID igGetHoveredID(void);
!x! void igSetHoveredID(ImGuiID id);
!x! void igKeepAliveID(ImGuiID id);
!x! void igMarkItemEdited(ImGuiID id);
!x! void igPushOverrideID(ImGuiID id);
!x! ImGuiID igGetIDWithSeed(const char* str_id_begin,const char* str_id_end,ImGuiID seed);
!x! void igItemSize_Vec2(const ImVec2 size,float text_baseline_y);
!x! void igItemSize_Rect(const ImRect bb,float text_baseline_y);
!x! bool igItemAdd(const ImRect bb,ImGuiID id,const ImRect* nav_bb,ImGuiItemFlags extra_flags);
!x! bool igItemHoverable(const ImRect bb,ImGuiID id);
!x! bool igIsClippedEx(const ImRect bb,ImGuiID id);
!x! void igSetLastItemData(ImGuiID item_id,ImGuiItemFlags in_flags,ImGuiItemStatusFlags status_flags,const ImRect item_rect);
!x! void igCalcItemSize(ImVec2 *pOut,ImVec2 size,float default_w,float default_h);
!x! float igCalcWrapWidthForPos(const ImVec2 pos,float wrap_pos_x);
!x! void igPushMultiItemsWidths(int components,float width_full);
!x! bool igIsItemToggledSelection(void);
!x! void igGetContentRegionMaxAbs(ImVec2 *pOut);
!x! void igShrinkWidths(ImGuiShrinkWidthItem* items,int count,float width_excess);
!x! void igPushItemFlag(ImGuiItemFlags option,bool enabled);
!x! void igPopItemFlag(void);
!x! void igLogBegin(ImGuiLogType type,int auto_open_depth);
!x! void igLogToBuffer(int auto_open_depth);
!x! void igLogRenderedText(const ImVec2* ref_pos,const char* text,const char* text_end);
!x! void igLogSetNextTextDecoration(const char* prefix,const char* suffix);
!x! bool igBeginChildEx(const char* name,ImGuiID id,const ImVec2 size_arg,bool border,ImGuiWindowFlags flags);
!x! void igOpenPopupEx(ImGuiID id,ImGuiPopupFlags popup_flags);
!x! void igClosePopupToLevel(int remaining,bool restore_focus_to_window_under_popup);
!x! void igClosePopupsOverWindow(ImGuiWindow* ref_window,bool restore_focus_to_window_under_popup);
!x! void igClosePopupsExceptModals(void);
!x! bool igIsPopupOpen_ID(ImGuiID id,ImGuiPopupFlags popup_flags);
!x! bool igBeginPopupEx(ImGuiID id,ImGuiWindowFlags extra_flags);
!x! void igBeginTooltipEx(ImGuiTooltipFlags tooltip_flags,ImGuiWindowFlags extra_window_flags);
!x! void igGetPopupAllowedExtentRect(ImRect *pOut,ImGuiWindow* window);
!x! ImGuiWindow* igGetTopMostPopupModal(void);
!x! ImGuiWindow* igGetTopMostAndVisiblePopupModal(void);
!x! void igFindBestWindowPosForPopup(ImVec2 *pOut,ImGuiWindow* window);
!x! void igFindBestWindowPosForPopupEx(ImVec2 *pOut,const ImVec2 ref_pos,const ImVec2 size,ImGuiDir* last_dir,const ImRect r_outer,const ImRect r_avoid,ImGuiPopupPositionPolicy policy);
!x! bool igBeginViewportSideBar(const char* name,ImGuiViewport* viewport,ImGuiDir dir,float size,ImGuiWindowFlags window_flags);
!x! bool igBeginMenuEx(const char* label,const char* icon,bool enabled);
!x! bool igMenuItemEx(const char* label,const char* icon,const char* shortcut,bool selected,bool enabled);
!x! bool igBeginComboPopup(ImGuiID popup_id,const ImRect bb,ImGuiComboFlags flags);
!x! bool igBeginComboPreview(void);
!x! void igEndComboPreview(void);
!x! void igNavInitWindow(ImGuiWindow* window,bool force_reinit);
!x! void igNavInitRequestApplyResult(void);
!x! bool igNavMoveRequestButNoResultYet(void);
!x! void igNavMoveRequestSubmit(ImGuiDir move_dir,ImGuiDir clip_dir,ImGuiNavMoveFlags move_flags,ImGuiScrollFlags scroll_flags);
!x! void igNavMoveRequestForward(ImGuiDir move_dir,ImGuiDir clip_dir,ImGuiNavMoveFlags move_flags,ImGuiScrollFlags scroll_flags);
!x! void igNavMoveRequestResolveWithLastItem(ImGuiNavItemData* result);
!x! void igNavMoveRequestCancel(void);
!x! void igNavMoveRequestApplyResult(void);
!x! void igNavMoveRequestTryWrapping(ImGuiWindow* window,ImGuiNavMoveFlags move_flags);
!x! const char* igGetNavInputName(ImGuiNavInput n);
!x! float igGetNavInputAmount(ImGuiNavInput n,ImGuiNavReadMode mode);
!x! void igGetNavInputAmount2d(ImVec2 *pOut,ImGuiNavDirSourceFlags dir_sources,ImGuiNavReadMode mode,float slow_factor,float fast_factor);
!x! int igCalcTypematicRepeatAmount(float t0,float t1,float repeat_delay,float repeat_rate);
!x! void igActivateItem(ImGuiID id);
!x! void igSetNavWindow(ImGuiWindow* window);
!x! void igSetNavID(ImGuiID id,ImGuiNavLayer nav_layer,ImGuiID focus_scope_id,const ImRect rect_rel);
!x! void igPushFocusScope(ImGuiID id);
!x! void igPopFocusScope(void);
!x! ImGuiID igGetFocusedFocusScope(void);
!x! ImGuiID igGetFocusScope(void);
!x! bool igIsNamedKey(ImGuiKey key);
!x! bool igIsLegacyKey(ImGuiKey key);
!x! bool igIsGamepadKey(ImGuiKey key);
!x! ImGuiKeyData* igGetKeyData(ImGuiKey key);
!x! void igSetItemUsingMouseWheel(void);
!x! void igSetActiveIdUsingNavAndKeys(void);
!x! bool igIsActiveIdUsingNavDir(ImGuiDir dir);
!x! bool igIsActiveIdUsingNavInput(ImGuiNavInput input);
!x! bool igIsActiveIdUsingKey(ImGuiKey key);
!x! void igSetActiveIdUsingKey(ImGuiKey key);
!x! bool igIsMouseDragPastThreshold(ImGuiMouseButton button,float lock_threshold);
!x! bool igIsNavInputDown(ImGuiNavInput n);
!x! bool igIsNavInputTest(ImGuiNavInput n,ImGuiNavReadMode rm);
!x! ImGuiModFlags igGetMergedModFlags(void);
!x! bool igIsKeyPressedMap(ImGuiKey key,bool repeat);
!x! void igDockContextInitialize(ImGuiContext* ctx);
!x! void igDockContextShutdown(ImGuiContext* ctx);
!x! void igDockContextClearNodes(ImGuiContext* ctx,ImGuiID root_id,bool clear_settings_refs);
!x! void igDockContextRebuildNodes(ImGuiContext* ctx);
!x! void igDockContextNewFrameUpdateUndocking(ImGuiContext* ctx);
!x! void igDockContextNewFrameUpdateDocking(ImGuiContext* ctx);
!x! void igDockContextEndFrame(ImGuiContext* ctx);
!x! ImGuiID igDockContextGenNodeID(ImGuiContext* ctx);
!x! void igDockContextQueueDock(ImGuiContext* ctx,ImGuiWindow* target,ImGuiDockNode* target_node,ImGuiWindow* payload,ImGuiDir split_dir,float split_ratio,bool split_outer);
!x! void igDockContextQueueUndockWindow(ImGuiContext* ctx,ImGuiWindow* window);
!x! void igDockContextQueueUndockNode(ImGuiContext* ctx,ImGuiDockNode* node);
!x! bool igDockContextCalcDropPosForDocking(ImGuiWindow* target,ImGuiDockNode* target_node,ImGuiWindow* payload,ImGuiDir split_dir,bool split_outer,ImVec2* out_pos);
!x! bool igDockNodeBeginAmendTabBar(ImGuiDockNode* node);
!x! void igDockNodeEndAmendTabBar(void);
!x! ImGuiDockNode* igDockNodeGetRootNode(ImGuiDockNode* node);
!x! bool igDockNodeIsInHierarchyOf(ImGuiDockNode* node,ImGuiDockNode* parent);
!x! int igDockNodeGetDepth(const ImGuiDockNode* node);
!x! ImGuiID igDockNodeGetWindowMenuButtonId(const ImGuiDockNode* node);
!x! ImGuiDockNode* igGetWindowDockNode(void);
!x! bool igGetWindowAlwaysWantOwnTabBar(ImGuiWindow* window);
!x! void igBeginDocked(ImGuiWindow* window,bool* p_open);
!x! void igBeginDockableDragDropSource(ImGuiWindow* window);
!x! void igBeginDockableDragDropTarget(ImGuiWindow* window);
!x! void igSetWindowDock(ImGuiWindow* window,ImGuiID dock_id,ImGuiCond cond);
!x! void igDockBuilderDockWindow(const char* window_name,ImGuiID node_id);
!x! ImGuiDockNode* igDockBuilderGetNode(ImGuiID node_id);
!x! ImGuiDockNode* igDockBuilderGetCentralNode(ImGuiID node_id);
!x! ImGuiID igDockBuilderAddNode(ImGuiID node_id,ImGuiDockNodeFlags flags);
!x! void igDockBuilderRemoveNode(ImGuiID node_id);
!x! void igDockBuilderRemoveNodeDockedWindows(ImGuiID node_id,bool clear_settings_refs);
!x! void igDockBuilderRemoveNodeChildNodes(ImGuiID node_id);
!x! void igDockBuilderSetNodePos(ImGuiID node_id,ImVec2 pos);
!x! void igDockBuilderSetNodeSize(ImGuiID node_id,ImVec2 size);
!x! ImGuiID igDockBuilderSplitNode(ImGuiID node_id,ImGuiDir split_dir,float size_ratio_for_node_at_dir,ImGuiID* out_id_at_dir,ImGuiID* out_id_at_opposite_dir);
!x! void igDockBuilderCopyDockSpace(ImGuiID src_dockspace_id,ImGuiID dst_dockspace_id,ImVector_const_charPtr* in_window_remap_pairs);
!x! void igDockBuilderCopyNode(ImGuiID src_node_id,ImGuiID dst_node_id,ImVector_ImGuiID* out_node_remap_pairs);
!x! void igDockBuilderCopyWindowSettings(const char* src_name,const char* dst_name);
!x! void igDockBuilderFinish(ImGuiID node_id);
!x! bool igIsDragDropActive(void);
!x! bool igBeginDragDropTargetCustom(const ImRect bb,ImGuiID id);
!x! void igClearDragDrop(void);
!x! bool igIsDragDropPayloadBeingAccepted(void);
!x! void igSetWindowClipRectBeforeSetChannel(ImGuiWindow* window,const ImRect clip_rect);
!x! void igBeginColumns(const char* str_id,int count,ImGuiOldColumnFlags flags);
!x! void igEndColumns(void);
!x! void igPushColumnClipRect(int column_index);
!x! void igPushColumnsBackground(void);
!x! void igPopColumnsBackground(void);
!x! ImGuiID igGetColumnsID(const char* str_id,int count);
!x! ImGuiOldColumns* igFindOrCreateColumns(ImGuiWindow* window,ImGuiID id);
!x! float igGetColumnOffsetFromNorm(const ImGuiOldColumns* columns,float offset_norm);
!x! float igGetColumnNormFromOffset(const ImGuiOldColumns* columns,float offset);
!x! void igTableOpenContextMenu(int column_n);
!x! void igTableSetColumnWidth(int column_n,float width);
!x! void igTableSetColumnSortDirection(int column_n,ImGuiSortDirection sort_direction,bool append_to_sort_specs);
!x! int igTableGetHoveredColumn(void);
!x! float igTableGetHeaderRowHeight(void);
!x! void igTablePushBackgroundChannel(void);
!x! void igTablePopBackgroundChannel(void);
!x! ImGuiTable* igGetCurrentTable(void);
!x! ImGuiTable* igTableFindByID(ImGuiID id);
!x! bool igBeginTableEx(const char* name,ImGuiID id,int columns_count,ImGuiTableFlags flags,const ImVec2 outer_size,float inner_width);
!x! void igTableBeginInitMemory(ImGuiTable* table,int columns_count);
!x! void igTableBeginApplyRequests(ImGuiTable* table);
!x! void igTableSetupDrawChannels(ImGuiTable* table);
!x! void igTableUpdateLayout(ImGuiTable* table);
!x! void igTableUpdateBorders(ImGuiTable* table);
!x! void igTableUpdateColumnsWeightFromWidth(ImGuiTable* table);
!x! void igTableDrawBorders(ImGuiTable* table);
!x! void igTableDrawContextMenu(ImGuiTable* table);
!x! void igTableMergeDrawChannels(ImGuiTable* table);
!x! ImGuiTableInstanceData* igTableGetInstanceData(ImGuiTable* table,int instance_no);
!x! void igTableSortSpecsSanitize(ImGuiTable* table);
!x! void igTableSortSpecsBuild(ImGuiTable* table);
!x! ImGuiSortDirection igTableGetColumnNextSortDirection(ImGuiTableColumn* column);
!x! void igTableFixColumnSortDirection(ImGuiTable* table,ImGuiTableColumn* column);
!x! float igTableGetColumnWidthAuto(ImGuiTable* table,ImGuiTableColumn* column);
!x! void igTableBeginRow(ImGuiTable* table);
!x! void igTableEndRow(ImGuiTable* table);
!x! void igTableBeginCell(ImGuiTable* table,int column_n);
!x! void igTableEndCell(ImGuiTable* table);
!x! void igTableGetCellBgRect(ImRect *pOut,const ImGuiTable* table,int column_n);
!x! const char* igTableGetColumnName_TablePtr(const ImGuiTable* table,int column_n);
!x! ImGuiID igTableGetColumnResizeID(const ImGuiTable* table,int column_n,int instance_no);
!x! float igTableGetMaxColumnWidth(const ImGuiTable* table,int column_n);
!x! void igTableSetColumnWidthAutoSingle(ImGuiTable* table,int column_n);
!x! void igTableSetColumnWidthAutoAll(ImGuiTable* table);
!x! void igTableRemove(ImGuiTable* table);
!x! void igTableGcCompactTransientBuffers_TablePtr(ImGuiTable* table);
!x! void igTableGcCompactTransientBuffers_TableTempDataPtr(ImGuiTableTempData* table);
!x! void igTableGcCompactSettings(void);
!x! void igTableLoadSettings(ImGuiTable* table);
!x! void igTableSaveSettings(ImGuiTable* table);
!x! void igTableResetSettings(ImGuiTable* table);
!x! ImGuiTableSettings* igTableGetBoundSettings(ImGuiTable* table);
!x! void igTableSettingsAddSettingsHandler(void);
!x! ImGuiTableSettings* igTableSettingsCreate(ImGuiID id,int columns_count);
!x! ImGuiTableSettings* igTableSettingsFindByID(ImGuiID id);
!x! bool igBeginTabBarEx(ImGuiTabBar* tab_bar,const ImRect bb,ImGuiTabBarFlags flags,ImGuiDockNode* dock_node);
!x! ImGuiTabItem* igTabBarFindTabByID(ImGuiTabBar* tab_bar,ImGuiID tab_id);
!x! ImGuiTabItem* igTabBarFindMostRecentlySelectedTabForActiveWindow(ImGuiTabBar* tab_bar);
!x! void igTabBarAddTab(ImGuiTabBar* tab_bar,ImGuiTabItemFlags tab_flags,ImGuiWindow* window);
!x! void igTabBarRemoveTab(ImGuiTabBar* tab_bar,ImGuiID tab_id);
!x! void igTabBarCloseTab(ImGuiTabBar* tab_bar,ImGuiTabItem* tab);
!x! void igTabBarQueueReorder(ImGuiTabBar* tab_bar,const ImGuiTabItem* tab,int offset);
!x! void igTabBarQueueReorderFromMousePos(ImGuiTabBar* tab_bar,const ImGuiTabItem* tab,ImVec2 mouse_pos);
!x! bool igTabBarProcessReorder(ImGuiTabBar* tab_bar);
!x! bool igTabItemEx(ImGuiTabBar* tab_bar,const char* label,bool* p_open,ImGuiTabItemFlags flags,ImGuiWindow* docked_window);
!x! void igTabItemCalcSize(ImVec2 *pOut,const char* label,bool has_close_button);
!x! void igTabItemBackground(ImDrawList* draw_list,const ImRect bb,ImGuiTabItemFlags flags,ImU32 col);
!x! void igTabItemLabelAndCloseButton(ImDrawList* draw_list,const ImRect bb,ImGuiTabItemFlags flags,ImVec2 frame_padding,const char* label,ImGuiID tab_id,ImGuiID close_button_id,bool is_contents_visible,bool* out_just_closed,bool* out_text_clipped);
!x! void igRenderText(ImVec2 pos,const char* text,const char* text_end,bool hide_text_after_hash);
!x! void igRenderTextWrapped(ImVec2 pos,const char* text,const char* text_end,float wrap_width);
!x! void igRenderTextClipped(const ImVec2 pos_min,const ImVec2 pos_max,const char* text,const char* text_end,const ImVec2* text_size_if_known,const ImVec2 align,const ImRect* clip_rect);
!x! void igRenderTextClippedEx(ImDrawList* draw_list,const ImVec2 pos_min,const ImVec2 pos_max,const char* text,const char* text_end,const ImVec2* text_size_if_known,const ImVec2 align,const ImRect* clip_rect);
!x! void igRenderTextEllipsis(ImDrawList* draw_list,const ImVec2 pos_min,const ImVec2 pos_max,float clip_max_x,float ellipsis_max_x,const char* text,const char* text_end,const ImVec2* text_size_if_known);
!x! void igRenderFrame(ImVec2 p_min,ImVec2 p_max,ImU32 fill_col,bool border,float rounding);
!x! void igRenderFrameBorder(ImVec2 p_min,ImVec2 p_max,float rounding);
!x! void igRenderColorRectWithAlphaCheckerboard(ImDrawList* draw_list,ImVec2 p_min,ImVec2 p_max,ImU32 fill_col,float grid_step,ImVec2 grid_off,float rounding,ImDrawFlags flags);
!x! void igRenderNavHighlight(const ImRect bb,ImGuiID id,ImGuiNavHighlightFlags flags);
!x! const char* igFindRenderedTextEnd(const char* text,const char* text_end);
!x! void igRenderMouseCursor(ImVec2 pos,float scale,ImGuiMouseCursor mouse_cursor,ImU32 col_fill,ImU32 col_border,ImU32 col_shadow);
!x! void igRenderArrow(ImDrawList* draw_list,ImVec2 pos,ImU32 col,ImGuiDir dir,float scale);
!x! void igRenderBullet(ImDrawList* draw_list,ImVec2 pos,ImU32 col);
!x! void igRenderCheckMark(ImDrawList* draw_list,ImVec2 pos,ImU32 col,float sz);
!x! void igRenderArrowPointingAt(ImDrawList* draw_list,ImVec2 pos,ImVec2 half_sz,ImGuiDir direction,ImU32 col);
!x! void igRenderArrowDockMenu(ImDrawList* draw_list,ImVec2 p_min,float sz,ImU32 col);
!x! void igRenderRectFilledRangeH(ImDrawList* draw_list,const ImRect rect,ImU32 col,float x_start_norm,float x_end_norm,float rounding);
!x! void igRenderRectFilledWithHole(ImDrawList* draw_list,const ImRect outer,const ImRect inner,ImU32 col,float rounding);
!x! ImDrawFlags igCalcRoundingFlagsForRectInRect(const ImRect r_in,const ImRect r_outer,float threshold);
!x! void igTextEx(const char* text,const char* text_end,ImGuiTextFlags flags);
!x! bool igButtonEx(const char* label,const ImVec2 size_arg,ImGuiButtonFlags flags);
!x! bool igCloseButton(ImGuiID id,const ImVec2 pos);
!x! bool igCollapseButton(ImGuiID id,const ImVec2 pos,ImGuiDockNode* dock_node);
!x! bool igArrowButtonEx(const char* str_id,ImGuiDir dir,ImVec2 size_arg,ImGuiButtonFlags flags);
!x! void igScrollbar(ImGuiAxis axis);
!x! bool igScrollbarEx(const ImRect bb,ImGuiID id,ImGuiAxis axis,ImS64* p_scroll_v,ImS64 avail_v,ImS64 contents_v,ImDrawFlags flags);
!x! bool igImageButtonEx(ImGuiID id,ImTextureID texture_id,const ImVec2 size,const ImVec2 uv0,const ImVec2 uv1,const ImVec2 padding,const ImVec4 bg_col,const ImVec4 tint_col);
!x! void igGetWindowScrollbarRect(ImRect *pOut,ImGuiWindow* window,ImGuiAxis axis);
!x! ImGuiID igGetWindowScrollbarID(ImGuiWindow* window,ImGuiAxis axis);
!x! ImGuiID igGetWindowResizeCornerID(ImGuiWindow* window,int n);
!x! ImGuiID igGetWindowResizeBorderID(ImGuiWindow* window,ImGuiDir dir);
!x! void igSeparatorEx(ImGuiSeparatorFlags flags);
!x! bool igCheckboxFlags_S64Ptr(const char* label,ImS64* flags,ImS64 flags_value);
!x! bool igCheckboxFlags_U64Ptr(const char* label,ImU64* flags,ImU64 flags_value);
!x! bool igButtonBehavior(const ImRect bb,ImGuiID id,bool* out_hovered,bool* out_held,ImGuiButtonFlags flags);
!x! bool igDragBehavior(ImGuiID id,ImGuiDataType data_type,void* p_v,float v_speed,const void* p_min,const void* p_max,const char* format,ImGuiSliderFlags flags);
!x! bool igSliderBehavior(const ImRect bb,ImGuiID id,ImGuiDataType data_type,void* p_v,const void* p_min,const void* p_max,const char* format,ImGuiSliderFlags flags,ImRect* out_grab_bb);
!x! bool igSplitterBehavior(const ImRect bb,ImGuiID id,ImGuiAxis axis,float* size1,float* size2,float min_size1,float min_size2,float hover_extend,float hover_visibility_delay,ImU32 bg_col);
!x! bool igTreeNodeBehavior(ImGuiID id,ImGuiTreeNodeFlags flags,const char* label,const char* label_end);
!x! bool igTreeNodeBehaviorIsOpen(ImGuiID id,ImGuiTreeNodeFlags flags);
!x! void igTreePushOverrideID(ImGuiID id);
!x! const ImGuiDataTypeInfo* igDataTypeGetInfo(ImGuiDataType data_type);
!x! int igDataTypeFormatString(char* buf,int buf_size,ImGuiDataType data_type,const void* p_data,const char* format);
!x! void igDataTypeApplyOp(ImGuiDataType data_type,int op,void* output,const void* arg_1,const void* arg_2);
!x! bool igDataTypeApplyFromText(const char* buf,ImGuiDataType data_type,void* p_data,const char* format);
!x! int igDataTypeCompare(ImGuiDataType data_type,const void* arg_1,const void* arg_2);
!x! bool igDataTypeClamp(ImGuiDataType data_type,void* p_data,const void* p_min,const void* p_max);
!x! bool igInputTextEx(const char* label,const char* hint,char* buf,int buf_size,const ImVec2 size_arg,ImGuiInputTextFlags flags,ImGuiInputTextCallback callback,void* user_data);
!x! bool igTempInputText(const ImRect bb,ImGuiID id,const char* label,char* buf,int buf_size,ImGuiInputTextFlags flags);
!x! bool igTempInputScalar(const ImRect bb,ImGuiID id,const char* label,ImGuiDataType data_type,void* p_data,const char* format,const void* p_clamp_min,const void* p_clamp_max);
!x! bool igTempInputIsActive(ImGuiID id);
!x! ImGuiInputTextState* igGetInputTextState(ImGuiID id);
!x! void igColorTooltip(const char* text,const float* col,ImGuiColorEditFlags flags);
!x! void igColorEditOptionsPopup(const float* col,ImGuiColorEditFlags flags);
!x! void igColorPickerOptionsPopup(const float* ref_col,ImGuiColorEditFlags flags);
!x! int igPlotEx(ImGuiPlotType plot_type,const char* label,float(*values_getter)(void* data,int idx),void* data,int values_count,int values_offset,const char* overlay_text,float scale_min,float scale_max,ImVec2 frame_size);
!x! void igShadeVertsLinearColorGradientKeepAlpha(ImDrawList* draw_list,int vert_start_idx,int vert_end_idx,ImVec2 gradient_p0,ImVec2 gradient_p1,ImU32 col0,ImU32 col1);
!x! void igShadeVertsLinearUV(ImDrawList* draw_list,int vert_start_idx,int vert_end_idx,const ImVec2 a,const ImVec2 b,const ImVec2 uv_a,const ImVec2 uv_b,bool clamp);
!x! void igGcCompactTransientMiscBuffers(void);
!x! void igGcCompactTransientWindowBuffers(ImGuiWindow* window);
!x! void igGcAwakeTransientWindowBuffers(ImGuiWindow* window);
!x! void igDebugLog(const char* fmt,...);
!x! void igDebugLogV(const char* fmt,va_list args);
!x! void igErrorCheckEndFrameRecover(ImGuiErrorLogCallback log_callback,void* user_data);
!x! void igErrorCheckEndWindowRecover(ImGuiErrorLogCallback log_callback,void* user_data);
!x! void igDebugDrawItemRect(ImU32 col);
!x! void igDebugStartItemPicker(void);
!x! void igShowFontAtlas(ImFontAtlas* atlas);
!x! void igDebugHookIdInfo(ImGuiID id,ImGuiDataType data_type,const void* data_id,const void* data_id_end);
!x! void igDebugNodeColumns(ImGuiOldColumns* columns);
!x! void igDebugNodeDockNode(ImGuiDockNode* node,const char* label);
!x! void igDebugNodeDrawList(ImGuiWindow* window,ImGuiViewportP* viewport,const ImDrawList* draw_list,const char* label);
!x! void igDebugNodeDrawCmdShowMeshAndBoundingBox(ImDrawList* out_draw_list,const ImDrawList* draw_list,const ImDrawCmd* draw_cmd,bool show_mesh,bool show_aabb);
!x! void igDebugNodeFont(ImFont* font);
!x! void igDebugNodeFontGlyph(ImFont* font,const ImFontGlyph* glyph);
!x! void igDebugNodeStorage(ImGuiStorage* storage,const char* label);
!x! void igDebugNodeTabBar(ImGuiTabBar* tab_bar,const char* label);
!x! void igDebugNodeTable(ImGuiTable* table);
!x! void igDebugNodeTableSettings(ImGuiTableSettings* settings);
!x! void igDebugNodeInputTextState(ImGuiInputTextState* state);
!x! void igDebugNodeWindow(ImGuiWindow* window,const char* label);
!x! void igDebugNodeWindowSettings(ImGuiWindowSettings* settings);
!x! void igDebugNodeWindowsList(ImVector_ImGuiWindowPtr* windows,const char* label);
!x! void igDebugNodeWindowsListByBeginStackParent(ImGuiWindow** windows,int windows_size,ImGuiWindow* parent_in_begin_stack);
!x! void igDebugNodeViewport(ImGuiViewportP* viewport);
!x! void igDebugRenderViewportThumbnail(ImDrawList* draw_list,ImGuiViewportP* viewport,const ImRect bb);
!x! void igImFontAtlasBuildInit(ImFontAtlas* atlas);
!x! void igImFontAtlasBuildSetupFont(ImFontAtlas* atlas,ImFont* font,ImFontConfig* font_config,float ascent,float descent);
!x! void igImFontAtlasBuildPackCustomRects(ImFontAtlas* atlas,void* stbrp_context_opaque);
!x! void igImFontAtlasBuildFinish(ImFontAtlas* atlas);
!x! void igImFontAtlasBuildRender8bppRectFromString(ImFontAtlas* atlas,int x,int y,int w,int h,const char* in_str,char in_marker_char,unsigned char in_marker_pixel_value);
!x! void igImFontAtlasBuildRender32bppRectFromString(ImFontAtlas* atlas,int x,int y,int w,int h,const char* in_str,char in_marker_char,unsigned int in_marker_pixel_value);
!x! void igImFontAtlasBuildMultiplyCalcLookupTable(unsigned char out_table[256],float in_multiply_factor);
!x! void igImFontAtlasBuildMultiplyRectAlpha8(const unsigned char table[256],unsigned char* pixels,int x,int y,int w,int h,int stride);
!x! const ImFontBuilderIO* ImGuiFreeType_GetBuilderForFreeType(void);
!x! void ImGuiFreeType_SetAllocatorFunctions(void*(*alloc_func)(size_t sz,void* user_data),void(*free_func)(void* ptr,void* user_data),void* user_data);
!x! void igLogText(CONST char *fmt, ...);
!x! void ImGuiTextBuffer_appendf(struct ImGuiTextBuffer *buffer, const char *fmt, ...);
!x! float igGET_FLT_MAX(void);
!x! float igGET_FLT_MIN(void);
!x! ImVector_ImWchar* ImVector_ImWchar_create(void);
!x! void ImVector_ImWchar_destroy(ImVector_ImWchar* self);
!x! void ImVector_ImWchar_Init(ImVector_ImWchar* p);
!x! void ImVector_ImWchar_UnInit(ImVector_ImWchar* p);
