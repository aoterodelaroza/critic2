  ! -*-F90-*-

  !xx! cimgui user-defined types
  type, bind(c) :: ImVec2
     real(c_float) :: x
     real(c_float) :: y
  end type ImVec2

  type, bind(c) :: ImGuiKeyData
     logical(c_bool) :: Down
     real(c_float) :: DownDuration
     real(c_float) :: DownDurationPrev
     real(c_float) :: AnalogValue
  end type ImGuiKeyData

  type, bind(c) :: ImVector_ImWchar
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImWchar

  type, bind(c) :: ImGuiIO
     integer(c_int) :: ConfigFlags  ! ImGuiConfigFlags
     integer(c_int) :: BackendFlags ! ImGuiBackendFlags
     type(ImVec2) :: DisplaySize
     real(c_float) :: DeltaTime
     real(c_float) :: IniSavingRate
     type(c_ptr) :: IniFilename ! const char *
     type(c_ptr) :: LogFilename ! const char *
     real(c_float) :: MouseDoubleClickTime
     real(c_float) :: MouseDoubleClickMaxDist
     real(c_float) :: MouseDragThreshold
     real(c_float) :: KeyRepeatDelay
     real(c_float) :: KeyRepeatRate
     type(c_ptr) :: UserData ! void *
     type(c_ptr) :: Fonts ! ImFontAtlas*
     real(c_float) :: FontGlobalScale
     logical(c_bool) :: FontAllowUserScaling
     type(c_ptr) :: FontDefault ! ImFont*
     type(ImVec2) :: DisplayFramebufferScale
     logical(c_bool) :: ConfigDockingNoSplit
     logical(c_bool) :: ConfigDockingWithShift
     logical(c_bool) :: ConfigDockingAlwaysTabBar
     logical(c_bool) :: ConfigDockingTransparentPayload
     logical(c_bool) :: ConfigViewportsNoAutoMerge
     logical(c_bool) :: ConfigViewportsNoTaskBarIcon
     logical(c_bool) :: ConfigViewportsNoDecoration
     logical(c_bool) :: ConfigViewportsNoDefaultParent
     logical(c_bool) :: MouseDrawCursor
     logical(c_bool) :: ConfigMacOSXBehaviors
     logical(c_bool) :: ConfigInputTrickleEventQueue
     logical(c_bool) :: ConfigInputTextCursorBlink
     logical(c_bool) :: ConfigDragClickToInputText
     logical(c_bool) :: ConfigWindowsResizeFromEdges
     logical(c_bool) :: ConfigWindowsMoveFromTitleBarOnly
     real(c_float) :: ConfigMemoryCompactTimer
     type(c_ptr) :: BackendPlatformName ! void *
     type(c_ptr) :: BackendRendererName ! void *
     type(c_ptr) :: BackendPlatformUserData ! void *
     type(c_ptr) :: BackendRendererUserData ! void *
     type(c_ptr) :: BackendLanguageUserData ! void *
     type(c_funptr) :: GetClipboardTextFn ! const char* (*GetClipboardTextFn)(void* user_data);
     type(c_funptr) :: SetClipboardTextFn ! void (*SetClipboardTextFn)(void* user_data, const char* text)
     type(c_ptr) :: ClipboardUserData ! void *
     type(c_funptr) :: SetPlatformImeDataFn ! void (*SetPlatformImeDataFn)(ImGuiViewport* viewport, ImGuiPlatformImeData* data)
     type(c_ptr) :: UnusedPadding ! void *
     logical(c_bool) :: WantCaptureMouse
     logical(c_bool) :: WantCaptureKeyboard
     logical(c_bool) :: WantTextInput
     logical(c_bool) :: WantSetMousePos
     logical(c_bool) :: WantSaveIniSettings
     logical(c_bool) :: NavActive
     logical(c_bool) :: NavVisible
     real(c_float) :: Framerate
     integer(c_int) :: MetricsRenderVertices
     integer(c_int) :: MetricsRenderIndices
     integer(c_int) :: MetricsRenderWindows
     integer(c_int) :: MetricsActiveWindows
     integer(c_int) :: MetricsActiveAllocations
     type(ImVec2) :: MouseDelta
     integer(c_int) :: KeyMap(ImGuiKey_COUNT)
     logical(c_bool) :: KeysDown(ImGuiKey_COUNT)
     type(ImVec2) :: MousePos
     logical(c_bool) :: MouseDown(5)
     real(c_float) :: MouseWheel
     real(c_float) :: MouseWheelH
     integer(c_int) :: MouseHoveredViewport ! ImGuiID
     logical(c_bool) :: KeyCtrl
     logical(c_bool) :: KeyShift
     logical(c_bool) :: KeyAlt
     logical(c_bool) :: KeySuper
     real(c_float) :: NavInputs(ImGuiNavInput_COUNT)
     integer(c_int) :: KeyMods ! ImGuiModFlags
     type(ImGuiKeyData) :: KeysData(ImGuiKey_KeysData_SIZE)
     logical(c_bool) :: WantCaptureMouseUnlessPopupClose
     type(ImVec2) :: MousePosPrev
     type(ImVec2) :: MouseClickedPos(5)
     real(c_double) :: MouseClickedTime(5)
     logical(c_bool) :: MouseClicked(5)
     logical(c_bool) :: MouseDoubleClicked(5)
     integer(c_short) :: MouseClickedCount(5)
     integer(c_short) :: MouseClickedLastCount(5)
     logical(c_bool) :: MouseReleased(5)
     logical(c_bool) :: MouseDownOwned(5)
     logical(c_bool) :: MouseDownOwnedUnlessPopupClose(5)
     real(c_float) :: MouseDownDuration(5)
     real(c_float) :: MouseDownDurationPrev(5)
     type(ImVec2) :: MouseDragMaxDistanceAbs(5)
     real(c_float) :: MouseDragMaxDistanceSqr(5)
     real(c_float) :: NavInputsDownDuration(ImGuiNavInput_COUNT)
     real(c_float) :: NavInputsDownDurationPrev(ImGuiNavInput_COUNT)
     real(c_float) :: PenPressure
     logical(c_bool) :: AppFocusLost
     logical(c_bool) :: AppAcceptingEvents
     character(c_char) :: BackendUsingLegacyKeyArrays ! ImS8
     logical(c_bool) :: BackendUsingLegacyNavInputArray
     integer(c_short) :: InputQueueSurrogate ! ImWchar16
     type(ImVector_ImWchar) :: InputQueueCharacters
  end type ImGuiIO
