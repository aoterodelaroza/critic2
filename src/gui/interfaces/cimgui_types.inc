  ! -*-F90-*-

  !xx! cimgui user-defined types
  type, bind(c) :: ImVector_const_charPtr
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_const_charPtr
  type, bind(c) :: ImVec2
     real(c_float) :: x
     real(c_float) :: y
  end type ImVec2
  type, bind(c) :: ImVec4
     real(c_float) :: x
     real(c_float) :: y
     real(c_float) :: z
     real(c_float) :: w
  end type ImVec4

  type, bind(c) :: ImGuiStyle
     real(c_float) :: Alpha
     real(c_float) :: DisabledAlpha
     type(ImVec2) :: WindowPadding
     real(c_float) :: WindowRounding
     real(c_float) :: WindowBorderSize
     type(ImVec2) :: WindowMinSize
     type(ImVec2) :: WindowTitleAlign
     integer(c_int) :: WindowMenuButtonPosition
     real(c_float) :: ChildRounding
     real(c_float) :: ChildBorderSize
     real(c_float) :: PopupRounding
     real(c_float) :: PopupBorderSize
     type(ImVec2) :: FramePadding
     real(c_float) :: FrameRounding
     real(c_float) :: FrameBorderSize
     type(ImVec2) :: ItemSpacing
     type(ImVec2) :: ItemInnerSpacing
     type(ImVec2) :: CellPadding
     type(ImVec2) :: TouchExtraPadding
     real(c_float) :: IndentSpacing
     real(c_float) :: ColumnsMinSpacing
     real(c_float) :: ScrollbarSize
     real(c_float) :: ScrollbarRounding
     real(c_float) :: GrabMinSize
     real(c_float) :: GrabRounding
     real(c_float) :: LogSliderDeadzone
     real(c_float) :: TabRounding
     real(c_float) :: TabBorderSize
     real(c_float) :: TabMinWidthForCloseButton
     integer(c_int) :: ColorButtonPosition
     type(ImVec2) :: ButtonTextAlign
     type(ImVec2) :: SelectableTextAlign
     type(ImVec2) :: DisplayWindowPadding
     type(ImVec2) :: DisplaySafeAreaPadding
     real(c_float) :: MouseCursorScale
     logical(c_bool) :: AntiAliasedLines
     logical(c_bool) :: AntiAliasedLinesUseTex
     logical(c_bool) :: AntiAliasedFill
     real(c_float) :: CurveTessellationTol
     real(c_float) :: CircleTessellationMaxError
     type(ImVec4) :: Colors(ImGuiCol_COUNT)
  end type ImGuiStyle

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
     type(c_funptr) :: GetClipboardTextFn ! const char* (*GetClipboardTextFn)(void* user_data)
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
     integer(c_int) :: MouseHoveredViewport
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

  ! struct ImGuiInputTextCallbackData
  ! {
  !     integer(c_int) :: EventFlag
  !     integer(c_int) :: Flags
  !     type(c_ptr) :: UserData
  !     integer(c_short) :: EventChar
  !     ImGuiKey EventKey
  !     char* Buf
  !     integer(c_int) :: BufTextLen
  !     integer(c_int) :: BufSize
  !     logical(c_bool) :: BufDirty
  !     integer(c_int) :: CursorPos
  !     integer(c_int) :: SelectionStart
  !     integer(c_int) :: SelectionEnd
  ! }
  ! struct ImGuiSizeCallbackData
  ! {
  !     type(c_ptr) :: UserData
  !     type(ImVec2) :: Pos
  !     type(ImVec2) :: CurrentSize
  !     type(ImVec2) :: DesiredSize
  ! }

  type, bind(c) :: ImGuiWindowClass
     integer(c_int) :: ClassId
     integer(c_int) :: ParentViewportId
     integer(c_int) :: ViewportFlagsOverrideSet
     integer(c_int) :: ViewportFlagsOverrideClear
     integer(c_int) :: TabItemFlagsOverrideSet
     integer(c_int) :: DockNodeFlagsOverrideSet
     logical(c_bool) :: DockingAlwaysTabBar
     logical(c_bool) :: DockingAllowUnclassed
  end type ImGuiWindowClass

  type, bind(c) :: ImGuiPayload
     type(c_ptr) :: Data
     integer(c_int) :: DataSize
     integer(c_int) :: SourceId
     integer(c_int) :: SourceParentId
     integer(c_int) :: DataFrameCount
     character(c_char) :: DataType(33)
     logical(c_bool) :: Preview
     logical(c_bool) :: Delivery
  end type ImGuiPayload

  ! struct ImGuiTableColumnSortSpecs
  ! {
  !     integer(c_int) :: ColumnUserID
  !     ImS16 ColumnIndex
  !     ImS16 SortOrder
  !     ImGuiSortDirection SortDirection : 8
  ! }
  ! struct ImGuiTableSortSpecs
  ! {
  !     const ImGuiTableColumnSortSpecs* Specs
  !     integer(c_int) :: SpecsCount
  !     logical(c_bool) :: SpecsDirty
  ! }
  ! struct ImGuiOnceUponAFrame
  ! {
  !      integer(c_int) :: RefFrame
  ! }
  ! struct ImGuiTextRange
  ! {
  !         type(c_ptr) :: b
  !         type(c_ptr) :: e
  ! }
  ! typedef struct ImVector_ImGuiTextRange {int Size;int Capacity;ImGuiTextRange* Data;} ImVector_ImGuiTextRange;
  ! struct ImGuiTextFilter
  ! {
  !     character(c_char) :: InputBuf[256]
  !     ImVector_ImGuiTextRange Filters
  !     integer(c_int) :: CountGrep
  ! }

  type, bind(c) :: ImVector_char
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_char

  type, bind(c) :: ImGuiTextBuffer
     type(ImVector_char) :: Buf
  end type ImGuiTextBuffer

  ! struct ImGuiStoragePair
  ! {
  !         integer(c_int) :: key
  !         union { int val_i; float val_f; void* val_p; };
  ! }

  type, bind(c) :: ImVector_ImGuiStoragePair
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiStoragePair

  type, bind(c) :: ImGuiStorage
     type(ImVector_ImGuiStoragePair) :: Data
  end type ImGuiStorage

  ! struct ImGuiListClipper
  ! {
  !     integer(c_int) :: DisplayStart
  !     integer(c_int) :: DisplayEnd
  !     integer(c_int) :: ItemsCount
  !     real(c_float) :: ItemsHeight
  !     real(c_float) :: StartPosY
  !     type(c_ptr) :: TempData
  ! }
  ! struct ImColor
  ! {
  !     type(ImVec4) :: Value
  ! }
  ! typedef void (*ImDrawCallback)(const ImDrawList* parent_list, const ImDrawCmd* cmd)
  ! struct ImDrawCmd
  ! {
  !     type(ImVec4) :: ClipRect
  !     type(c_ptr) :: TextureId
  !     integer(c_int) :: VtxOffset
  !     integer(c_int) :: IdxOffset
  !     integer(c_int) :: ElemCount
  !     ImDrawCallback UserCallback
  !     type(c_ptr) :: UserCallbackData
  ! }
  ! struct ImDrawVert
  ! {
  !     type(ImVec2) :: pos
  !     type(ImVec2) :: uv
  !     integer(c_int) :: col
  ! }

  type, bind(c) :: ImDrawCmdHeader
     type(ImVec4) :: ClipRect
     type(c_ptr) :: TextureId
     integer(c_int) :: VtxOffset
  end type ImDrawCmdHeader

  type, bind(c) :: ImVector_ImDrawCmd
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImDrawCmd

  type, bind(c) :: ImVector_ImDrawIdx
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImDrawIdx

  ! struct ImDrawChannel
  ! {
  !     type(ImVector_ImDrawCmd) :: _CmdBuffer
  !     type(ImVector_ImDrawIdx) :: _IdxBuffer
  ! }

  type, bind(c) :: ImVector_ImDrawChannel
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImDrawChannel

  type, bind(c) :: ImDrawListSplitter
     integer(c_int) :: Current
     integer(c_int) :: Count
     type(ImVector_ImDrawChannel) :: Channels
  end type ImDrawListSplitter

  type, bind(c) :: ImVector_ImDrawVert
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImDrawVert

  type, bind(c) :: ImVector_ImVec4
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImVec4

  type, bind(c) :: ImVector_ImTextureID
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImTextureID

  type, bind(c) :: ImVector_ImVec2
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImVec2

  type, bind(c) :: ImDrawList
     type(ImVector_ImDrawCmd) :: CmdBuffer
     type(ImVector_ImDrawIdx) :: IdxBuffer
     type(ImVector_ImDrawVert) :: VtxBuffer
     integer(c_int) :: Flags
     integer(c_int) :: VtxCurrentIdx
     type(c_ptr) :: Data
     type(c_ptr) :: OwnerName
     type(c_ptr) :: VtxWritePtr
     type(c_ptr) :: IdxWritePtr
     type(ImVector_ImVec4) :: ClipRectStack
     type(ImVector_ImTextureID) :: TextureIdStack
     type(ImVector_ImVec2) :: Path
     type(ImDrawCmdHeader) :: CmdHeader
     type(ImDrawListSplitter) :: Splitter
     real(c_float) :: FringeScale
  end type ImDrawList

  ! struct ImDrawData
  ! {
  !     logical(c_bool) :: Valid
  !     integer(c_int) :: CmdListsCount
  !     integer(c_int) :: TotalIdxCount
  !     integer(c_int) :: TotalVtxCount
  !     ImDrawList** CmdLists
  !     type(ImVec2) :: DisplayPos
  !     type(ImVec2) :: DisplaySize
  !     type(ImVec2) :: FramebufferScale
  !     ImGuiViewport* OwnerViewport
  ! }
  ! struct ImFontConfig
  ! {
  !     type(c_ptr) :: FontData
  !     integer(c_int) :: FontDataSize
  !     logical(c_bool) :: FontDataOwnedByAtlas
  !     integer(c_int) :: FontNo
  !     real(c_float) :: SizePixels
  !     integer(c_int) :: OversampleH
  !     integer(c_int) :: OversampleV
  !     logical(c_bool) :: PixelSnapH
  !     type(ImVec2) :: GlyphExtraSpacing
  !     type(ImVec2) :: GlyphOffset
  !     const ImWchar* GlyphRanges
  !     real(c_float) :: GlyphMinAdvanceX
  !     real(c_float) :: GlyphMaxAdvanceX
  !     logical(c_bool) :: MergeMode
  !     integer(c_int) :: FontBuilderFlags
  !     real(c_float) :: RasterizerMultiply
  !     integer(c_short) :: EllipsisChar
  !     character(c_char) :: Name[40]
  !     type(c_ptr) :: DstFont
  ! }
  ! struct ImFontGlyph
  ! {
  !     integer(c_int) :: Colored : 1
  !     integer(c_int) :: Visible : 1
  !     integer(c_int) :: Codepoint : 30
  !     real(c_float) :: AdvanceX
  !     real(c_float) :: X0, Y0, X1, Y1
  !     real(c_float) :: U0, V0, U1, V1
  ! }
  ! typedef struct ImVector_ImU32 {int Size;int Capacity;ImU32* Data;} ImVector_ImU32;

  ! struct ImFontGlyphRangesBuilder
  ! {
  !     ImVector_ImU32 UsedChars
  ! }
  ! typedef struct ImFontAtlasCustomRect ImFontAtlasCustomRect
  ! struct ImFontAtlasCustomRect
  ! {
  !     integer(c_short) :: Width, Height
  !     integer(c_short) :: X, Y
  !     integer(c_int) :: GlyphID
  !     real(c_float) :: GlyphAdvanceX
  !     type(ImVec2) :: GlyphOffset
  !     type(c_ptr) :: Font
  ! }

  type, bind(c) :: ImVector_ImFontPtr
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImFontPtr

  ! typedef struct ImVector_ImFontAtlasCustomRect {int Size;int Capacity;ImFontAtlasCustomRect* Data;} ImVector_ImFontAtlasCustomRect;

  ! typedef struct ImVector_ImFontConfig {int Size;int Capacity;ImFontConfig* Data;} ImVector_ImFontConfig;

  ! struct ImFontAtlas
  ! {
  !     ImFontAtlasFlags Flags
  !     type(c_ptr) :: TexID
  !     integer(c_int) :: TexDesiredWidth
  !     integer(c_int) :: TexGlyphPadding
  !     logical(c_bool) :: Locked
  !     logical(c_bool) :: TexReady
  !     logical(c_bool) :: TexPixelsUseColors
  !     unsigned char* TexPixelsAlpha8
  !     unsigned int* TexPixelsRGBA32
  !     integer(c_int) :: TexWidth
  !     integer(c_int) :: TexHeight
  !     type(ImVec2) :: TexUvScale
  !     type(ImVec2) :: TexUvWhitePixel
  !     type(ImVector_ImFontPtr) :: Fonts
  !     ImVector_ImFontAtlasCustomRect CustomRects
  !     ImVector_ImFontConfig ConfigData
  !     type(ImVec4) :: TexUvLines[(63) + 1]
  !     const ImFontBuilderIO* FontBuilderIO
  !     integer(c_int) :: FontBuilderFlags
  !     integer(c_int) :: PackIdMouseCursors
  !     integer(c_int) :: PackIdLines
  ! }

  type, bind(c) :: ImVector_float
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_float

  type, bind(c) :: ImVector_ImFontGlyph
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImFontGlyph

  type, bind(c) :: ImFont
     type(ImVector_float) :: IndexAdvanceX
     real(c_float) :: FallbackAdvanceX
     real(c_float) :: FontSize
     type(ImVector_ImWchar) :: IndexLookup
     type(ImVector_ImFontGlyph) :: Glyphs
     type(c_ptr) :: FallbackGlyph
     type(c_ptr) :: ContainerAtlas
     type(c_ptr) :: ConfigData
     integer(c_short) :: ConfigDataCount
     integer(c_short) :: FallbackChar
     integer(c_short) :: EllipsisChar
     integer(c_short) :: DotChar
     logical(c_bool) :: DirtyLookupTables
     real(c_float) :: Scale
     real(c_float) :: Ascent, Descent
     integer(c_int) :: MetricsTotalSurface
     integer(c_signed_char) :: Used4kPagesMap(2) ! (0xFFFF +1)/4096/8
  end type ImFont

  ! struct ImGuiViewport
  ! {
  !     integer(c_int) :: ID
  !     integer(c_int) :: Flags
  !     type(ImVec2) :: Pos
  !     type(ImVec2) :: Size
  !     type(ImVec2) :: WorkPos
  !     type(ImVec2) :: WorkSize
  !     real(c_float) :: DpiScale
  !     integer(c_int) :: ParentViewportId
  !     ImDrawData* DrawData
  !     type(c_ptr) :: RendererUserData
  !     type(c_ptr) :: PlatformUserData
  !     type(c_ptr) :: PlatformHandle
  !     type(c_ptr) :: PlatformHandleRaw
  !     logical(c_bool) :: PlatformRequestMove
  !     logical(c_bool) :: PlatformRequestResize
  !     logical(c_bool) :: PlatformRequestClose
  ! }

  type, bind(c) :: ImVector_ImGuiPlatformMonitor
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiPlatformMonitor

  type, bind(c) :: ImVector_ImGuiViewportPtr
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiViewportPtr

  type, bind(c) :: ImGuiPlatformIO
     type(c_funptr) :: Platform_CreateWindow
     type(c_funptr) :: Platform_DestroyWindow
     type(c_funptr) :: Platform_ShowWindow
     type(c_funptr) :: Platform_SetWindowPos
     type(c_funptr) :: Platform_GetWindowPos
     type(c_funptr) :: Platform_SetWindowSize
     type(c_funptr) :: Platform_GetWindowSize
     type(c_funptr) :: Platform_SetWindowFocus
     type(c_funptr) :: Platform_GetWindowFocus
     type(c_funptr) :: Platform_GetWindowMinimized
     type(c_funptr) :: Platform_SetWindowTitle
     type(c_funptr) :: Platform_SetWindowAlpha
     type(c_funptr) :: Platform_UpdateWindow
     type(c_funptr) :: Platform_RenderWindow
     type(c_funptr) :: Platform_SwapBuffers
     type(c_funptr) :: Platform_GetWindowDpiScale
     type(c_funptr) :: Platform_OnChangedViewport
     type(c_funptr) :: Platform_CreateVkSurface
     type(c_funptr) :: Renderer_CreateWindow
     type(c_funptr) :: Renderer_DestroyWindow
     type(c_funptr) :: Renderer_SetWindowSize
     type(c_funptr) :: Renderer_RenderWindow
     type(c_funptr) :: Renderer_SwapBuffers
     type(ImVector_ImGuiPlatformMonitor) :: Monitors
     type(ImVector_ImGuiViewportPtr) :: Viewports
  end type ImGuiPlatformIO

  type, bind(c) :: ImGuiPlatformMonitor
     type(ImVec2) :: MainPos
     type(ImVec2) :: MainSize
     type(ImVec2) :: WorkPos
     type(ImVec2) :: WorkSize
     real(c_float) :: DpiScale
  end type ImGuiPlatformMonitor

  type, bind(c) :: ImGuiPlatformImeData
     logical(c_bool) :: WantVisible
     type(ImVec2) :: InputPos
     real(c_float) :: InputLineHeight
  end type ImGuiPlatformImeData

  ! typedef void (*ImGuiErrorLogCallback)(void* user_data, const char* fmt, ...)

  type, bind(c) :: StbUndoRecord
     integer(c_int) :: where
     integer(c_int) :: insert_length
     integer(c_int) :: delete_length
     integer(c_int) :: char_storage
  end type StbUndoRecord

  type, bind(c) :: StbUndoState
     type(StbUndoRecord) :: undo_rec(99)
     integer(c_short) :: undo_char(999)
     integer(c_short) :: undo_point, redo_point
     integer(c_int) :: undo_char_point, redo_char_point
  end type StbUndoState

  type, bind(c) :: STB_TexteditState
     integer(c_int) :: cursor
     integer(c_int) :: select_start
     integer(c_int) :: select_end
     integer(c_signed_char) :: insert_mode
     integer(c_int) :: row_count_per_page
     integer(c_signed_char) :: cursor_at_end_of_line
     integer(c_signed_char) :: initialized
     integer(c_signed_char) :: has_preferred_x
     integer(c_signed_char) :: single_line
     integer(c_signed_char) :: padding1, padding2, padding3
     real(c_float) :: preferred_x
     type(StbUndoState) :: undostate
  end type STB_TexteditState

  ! struct StbTexteditRow
  ! {
  !    real(c_float) :: x0,x1
  !    real(c_float) :: baseline_y_delta
  !    real(c_float) :: ymin,ymax
  !    integer(c_int) :: num_chars
  ! }

  type, bind(c) :: ImVec1
     real(c_float) :: x
  end type ImVec1

  type, bind(c) :: ImVec2ih
     integer(c_short) :: x
     integer(c_short) :: y
  end type ImVec2ih

  type, bind(c) :: ImRect
     type(ImVec2) :: Min
     type(ImVec2) :: Max
  end type ImRect

  ! struct ImBitVector
  ! {
  !     ImVector_ImU32 Storage
  ! }

  type, bind(c) :: ImDrawListSharedData
     type(ImVec2) :: TexUvWhitePixel
     type(c_ptr) :: Font
     real(c_float) :: FontSize
     real(c_float) :: CurveTessellationTol
     real(c_float) :: CircleSegmentMaxError
     type(ImVec4) :: ClipRectFullscreen
     integer(c_int) :: InitialFlags
     type(ImVec2) :: ArcFastVtx(48)
     real(c_float) :: ArcFastRadiusCutoff
     integer(c_signed_char) :: CircleSegmentCounts(64)
     type(c_ptr) :: TexUvLines
  end type ImDrawListSharedData

  ! typedef struct ImVector_ImDrawListPtr {int Size;int Capacity;ImDrawList** Data;} ImVector_ImDrawListPtr;

  ! struct ImDrawDataBuilder
  ! {
  !     ImVector_ImDrawListPtr Layers[2]
  ! }
  ! struct ImGuiDataTypeTempStorage
  ! {
  !     integer(c_signed_char) :: Data[8]
  ! }
  ! struct ImGuiDataTypeInfo
  ! {
  !     size_t Size
  !     type(c_ptr) :: Name
  !     type(c_ptr) :: PrintFmt
  !     type(c_ptr) :: ScanFmt
  ! }
  ! struct ImGuiColorMod
  ! {
  !     ImGuiCol Col
  !     type(ImVec4) :: BackupValue
  ! }
  ! struct ImGuiStyleMod
  ! {
  !     ImGuiStyleVar VarIdx
  !     union { int BackupInt[2]; float BackupFloat[2]; };
  ! }

  type, bind(c) :: ImGuiComboPreviewData
     type(ImRect) :: PreviewRect
     type(ImVec2) :: BackupCursorPos
     type(ImVec2) :: BackupCursorMaxPos
     type(ImVec2) :: BackupCursorPosPrevLine
     real(c_float) :: BackupPrevLineTextBaseOffset
     integer(c_int) :: BackupLayout
  end type ImGuiComboPreviewData

  ! struct ImGuiGroupData
  ! {
  !     integer(c_int) :: WindowID
  !     type(ImVec2) :: BackupCursorPos
  !     type(ImVec2) :: BackupCursorMaxPos
  !     type(ImVec1) :: BackupIndent
  !     type(ImVec1) :: BackupGroupOffset
  !     type(ImVec2) :: BackupCurrLineSize
  !     real(c_float) :: BackupCurrLineTextBaseOffset
  !     integer(c_int) :: BackupActiveIdIsAlive
  !     logical(c_bool) :: BackupActiveIdPreviousFrameIsAlive
  !     logical(c_bool) :: BackupHoveredIdIsAlive
  !     logical(c_bool) :: EmitItem
  ! }

  type, bind(c) :: ImGuiMenuColumns
     integer(c_int) :: TotalWidth
     integer(c_int) :: NextTotalWidth
     integer(c_short) :: Spacing
     integer(c_short) :: OffsetIcon
     integer(c_short) :: OffsetLabel
     integer(c_short) :: OffsetShortcut
     integer(c_short) :: OffsetMark
     integer(c_short) :: Widths(4)
  end type ImGuiMenuColumns

  type, bind(c) :: ImGuiInputTextState
     integer(c_int) :: ID
     integer(c_int) :: CurLenW, CurLenA
     type(ImVector_ImWchar) :: TextW
     type(ImVector_char) :: TextA
     type(ImVector_char) :: InitialTextA
     logical(c_bool) :: TextAIsValid
     integer(c_int) :: BufCapacityA
     real(c_float) :: ScrollX
     type(STB_TexteditState) :: Stb
     real(c_float) :: CursorAnim
     logical(c_bool) :: CursorFollow
     logical(c_bool) :: SelectedAllMouseLock
     logical(c_bool) :: Edited
     integer(c_int) :: Flags
  end type ImGuiInputTextState

  ! struct ImGuiPopupData
  ! {
  !     integer(c_int) :: PopupId
  !     type(c_ptr) :: Window
  !     type(c_ptr) :: SourceWindow
  !     integer(c_int) :: ParentNavLayer
  !     integer(c_int) :: OpenFrameCount
  !     integer(c_int) :: OpenParentId
  !     type(ImVec2) :: OpenPopupPos
  !     type(ImVec2) :: OpenMousePos
  ! }

  type, bind(c) :: ImGuiNextWindowData
     integer(c_int) :: Flags
     integer(c_int) :: PosCond
     integer(c_int) :: SizeCond
     integer(c_int) :: CollapsedCond
     integer(c_int) :: DockCond
     type(ImVec2) :: PosVal
     type(ImVec2) :: PosPivotVal
     type(ImVec2) :: SizeVal
     type(ImVec2) :: ContentSizeVal
     type(ImVec2) :: ScrollVal
     logical(c_bool) :: PosUndock
     logical(c_bool) :: CollapsedVal
     type(ImRect) :: SizeConstraintRect
     type(c_funptr) :: SizeCallback
     type(c_ptr) :: SizeCallbackUserData
     real(c_float) :: BgAlphaVal
     integer(c_int) :: ViewportId
     integer(c_int) :: DockId
     type(ImGuiWindowClass) :: WindowClass
     type(ImVec2) :: MenuBarOffsetMinVal
  end type ImGuiNextWindowData

  type, bind(c) :: ImGuiNextItemData
     integer(c_int) :: Flags
     real(c_float) :: Width
     integer(c_int) :: FocusScopeId
     integer(c_int) :: OpenCond
     logical(c_bool) :: OpenVal
  end type ImGuiNextItemData

  type, bind(c) :: ImGuiLastItemData
     integer(c_int) :: ID
     integer(c_int) :: InFlags
     integer(c_int) :: StatusFlags
     type(ImRect) :: Rect
     type(ImRect) :: NavRect
     type(ImRect) :: DisplayRect
  end type ImGuiLastItemData

  ! struct ImGuiStackSizes
  ! {
  !     integer(c_short) :: SizeOfIDStack
  !     integer(c_short) :: SizeOfColorStack
  !     integer(c_short) :: SizeOfStyleVarStack
  !     integer(c_short) :: SizeOfFontStack
  !     integer(c_short) :: SizeOfFocusScopeStack
  !     integer(c_short) :: SizeOfGroupStack
  !     integer(c_short) :: SizeOfItemFlagsStack
  !     integer(c_short) :: SizeOfBeginPopupStack
  !     integer(c_short) :: SizeOfDisabledStack
  ! }
  ! struct ImGuiWindowStackData
  ! {
  !     type(c_ptr) :: Window
  !     type(ImGuiLastItemData) :: ParentLastItemDataBackup
  !     ImGuiStackSizes StackSizesOnBegin
  ! }
  ! struct ImGuiShrinkWidthItem
  ! {
  !     integer(c_int) :: Index
  !     real(c_float) :: Width
  !     real(c_float) :: InitialWidth
  ! }
  ! struct ImGuiPtrOrIndex
  ! {
  !     type(c_ptr) :: Ptr
  !     integer(c_int) :: Index
  ! }

  type, bind(c) :: ImBitArray_ImGuiKey_NamedKey_COUNT__lessImGuiKey_NamedKey_BEGIN
     integer(c_int) :: Storage(ImGuiKey_NamedKey_COUNT_31shift5)
  end type ImBitArray_ImGuiKey_NamedKey_COUNT__lessImGuiKey_NamedKey_BEGIN

  type, bind(c) :: ImBitArrayForNamedKeys
     integer(c_int) :: Storage(ImGuiKey_NamedKey_COUNT_31shift5)
  end type ImBitArrayForNamedKeys

  ! struct ImGuiInputEventMousePos
  ! { real(c_float) :: PosX, PosY
  ! }
  ! struct ImGuiInputEventMouseWheel
  ! { real(c_float) :: WheelX, WheelY
  ! }
  ! struct ImGuiInputEventMouseButton
  ! { integer(c_int) :: Button logical(c_bool) :: Down
  ! }
  ! struct ImGuiInputEventMouseViewport
  ! { integer(c_int) :: HoveredViewportID
  ! }
  ! struct ImGuiInputEventKey
  ! { ImGuiKey Key; logical(c_bool) :: Down; real(c_float) :: AnalogValue;
  ! }
  ! struct ImGuiInputEventText
  ! { integer(c_int) :: Char
  ! }
  ! struct ImGuiInputEventAppFocused
  ! { logical(c_bool) :: Focused
  ! }
  ! struct ImGuiInputEvent
  ! {
  !     ImGuiInputEventType Type
  !     integer(c_int) :: Source
  !     union
  !     {
  !         ImGuiInputEventMousePos MousePos
  !         ImGuiInputEventMouseWheel MouseWheel
  !         ImGuiInputEventMouseButton MouseButton
  !         ImGuiInputEventMouseViewport MouseViewport
  !         ImGuiInputEventKey Key
  !         ImGuiInputEventText Text
  !         ImGuiInputEventAppFocused AppFocused
  !     }
  !     logical(c_bool) :: AddedByTestEngine
  ! }
  ! struct ImGuiListClipperRange
  ! {
  !     integer(c_int) :: Min
  !     integer(c_int) :: Max
  !     logical(c_bool) :: PosToIndexConvert
  !     ImS8 PosToIndexOffsetMin
  !     ImS8 PosToIndexOffsetMax
  ! }
  ! typedef struct ImVector_ImGuiListClipperRange {int Size;int Capacity;ImGuiListClipperRange* Data;} ImVector_ImGuiListClipperRange;
  ! struct ImGuiListClipperData
  ! {
  !     ImGuiListClipper* ListClipper
  !     real(c_float) :: LossynessOffset
  !     integer(c_int) :: StepNo
  !     integer(c_int) :: ItemsFrozen
  !     ImVector_ImGuiListClipperRange Ranges
  ! }

  type, bind(c) :: ImGuiNavItemData
     type(c_ptr) :: Window
     integer(c_int) :: ID
     integer(c_int) :: FocusScopeId
     type(ImRect) :: RectRel
     integer(c_int) :: InFlags
     real(c_float) :: DistBox
     real(c_float) :: DistCenter
     real(c_float) :: DistAxial
  end type ImGuiNavItemData

  ! struct ImGuiOldColumnData
  ! {
  !     real(c_float) :: OffsetNorm
  !     real(c_float) :: OffsetNormBeforeResize
  !     ImGuiOldColumnFlags Flags
  !     type(ImRect) :: ClipRect
  ! }
  ! typedef struct ImVector_ImGuiOldColumnData {int Size;int Capacity;ImGuiOldColumnData* Data;} ImVector_ImGuiOldColumnData;
  ! struct ImGuiOldColumns
  ! {
  !     integer(c_int) :: ID
  !     ImGuiOldColumnFlags Flags
  !     logical(c_bool) :: IsFirstFrame
  !     logical(c_bool) :: IsBeingResized
  !     integer(c_int) :: Current
  !     integer(c_int) :: Count
  !     real(c_float) :: OffMinX, OffMaxX
  !     real(c_float) :: LineMinY, LineMaxY
  !     real(c_float) :: HostCursorPosY
  !     real(c_float) :: HostCursorMaxPosX
  !     type(ImRect) :: HostInitialClipRect
  !     type(ImRect) :: HostBackupClipRect
  !     type(ImRect) :: HostBackupParentWorkRect
  !     ImVector_ImGuiOldColumnData Columns
  !     type(ImDrawListSplitter) :: Splitter
  ! }

  type, bind(c) :: ImVector_ImGuiWindowPtr
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiWindowPtr

  ! struct ImGuiDockNode
  ! {
  !     integer(c_int) :: ID
  !     integer(c_int) :: SharedFlags
  !     integer(c_int) :: LocalFlags
  !     integer(c_int) :: LocalFlagsInWindows
  !     integer(c_int) :: MergedFlags
  !     ImGuiDockNodeState State
  !     type(c_ptr) :: ParentNode
  !     type(c_ptr) :: ChildNodes[2]
  !     type(ImVector_ImGuiWindowPtr) :: Windows
  !     ImGuiTabBar* TabBar
  !     type(ImVec2) :: Pos
  !     type(ImVec2) :: Size
  !     type(ImVec2) :: SizeRef
  !     ImGuiAxis SplitAxis
  !     type(ImGuiWindowClass) :: WindowClass
  !     integer(c_int) :: LastBgColor
  !     type(c_ptr) :: HostWindow
  !     type(c_ptr) :: VisibleWindow
  !     type(c_ptr) :: CentralNode
  !     type(c_ptr) :: OnlyNodeWithWindows
  !     integer(c_int) :: CountNodeWithWindows
  !     integer(c_int) :: LastFrameAlive
  !     integer(c_int) :: LastFrameActive
  !     integer(c_int) :: LastFrameFocused
  !     integer(c_int) :: LastFocusedNodeId
  !     integer(c_int) :: SelectedTabId
  !     integer(c_int) :: WantCloseTabId
  !     ImGuiDataAuthority AuthorityForPos :3
  !     ImGuiDataAuthority AuthorityForSize :3
  !     ImGuiDataAuthority AuthorityForViewport :3
  !     logical(c_bool) :: IsVisible :1
  !     logical(c_bool) :: IsFocused :1
  !     logical(c_bool) :: IsBgDrawnThisFrame :1
  !     logical(c_bool) :: HasCloseButton :1
  !     logical(c_bool) :: HasWindowMenuButton :1
  !     logical(c_bool) :: HasCentralNodeChild :1
  !     logical(c_bool) :: WantCloseAll :1
  !     logical(c_bool) :: WantLockSizeOnce :1
  !     logical(c_bool) :: WantMouseMove :1
  !     logical(c_bool) :: WantHiddenTabBarUpdate :1
  !     logical(c_bool) :: WantHiddenTabBarToggle :1
  ! }
  ! struct ImGuiWindowDockStyle
  ! {
  !     integer(c_int) :: Colors[ImGuiWindowDockStyleCol_COUNT]
  ! }

  type, bind(c) :: ImVector_ImGuiDockRequest
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiDockRequest

  type, bind(c) :: ImVector_ImGuiDockNodeSettings
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiDockNodeSettings

  type, bind(c) :: ImGuiDockContext
     type(ImGuiStorage) :: Nodes
     type(ImVector_ImGuiDockRequest) :: Requests
     type(ImVector_ImGuiDockNodeSettings) :: NodesSettings
     logical(c_bool) :: WantFullRebuild
  end type ImGuiDockContext

  ! struct ImGuiViewportP
  ! {
  !     ImGuiViewport _ImGuiViewport
  !     integer(c_int) :: Idx
  !     integer(c_int) :: LastFrameActive
  !     integer(c_int) :: LastFrontMostStampCount
  !     integer(c_int) :: LastNameHash
  !     type(ImVec2) :: LastPos
  !     real(c_float) :: Alpha
  !     real(c_float) :: LastAlpha
  !     integer(c_short) :: PlatformMonitor
  !     logical(c_bool) :: PlatformWindowCreated
  !     type(c_ptr) :: Window
  !     integer(c_int) :: DrawListsLastFrame[2]
  !     type(c_ptr) :: DrawLists[2]
  !     ImDrawData DrawDataP
  !     ImDrawDataBuilder DrawDataBuilder
  !     type(ImVec2) :: LastPlatformPos
  !     type(ImVec2) :: LastPlatformSize
  !     type(ImVec2) :: LastRendererSize
  !     type(ImVec2) :: WorkOffsetMin
  !     type(ImVec2) :: WorkOffsetMax
  !     type(ImVec2) :: BuildWorkOffsetMin
  !     type(ImVec2) :: BuildWorkOffsetMax
  ! }
  ! struct ImGuiWindowSettings
  ! {
  !     integer(c_int) :: ID
  !     type(ImVec2ih) :: Pos
  !     type(ImVec2ih) :: Size
  !     type(ImVec2ih) :: ViewportPos
  !     integer(c_int) :: ViewportId
  !     integer(c_int) :: DockId
  !     integer(c_int) :: ClassId
  !     integer(c_short) :: DockOrder
  !     logical(c_bool) :: Collapsed
  !     logical(c_bool) :: WantApply
  ! }
  ! struct ImGuiSettingsHandler
  ! {
  !     type(c_ptr) :: TypeName
  !     integer(c_int) :: TypeHash
  !     void (*ClearAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler)
  !     void (*ReadInitFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler)
  !     type(c_ptr) :: (*ReadOpenFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, const char* name)
  !     void (*ReadLineFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, void* entry, const char* line)
  !     void (*ApplyAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler)
  !     void (*WriteAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, ImGuiTextBuffer* out_buf)
  !     type(c_ptr) :: UserData
  ! }

  type, bind(c) :: ImGuiMetricsConfig
     logical(c_bool) :: ShowDebugLog
     logical(c_bool) :: ShowStackTool
     logical(c_bool) :: ShowWindowsRects
     logical(c_bool) :: ShowWindowsBeginOrder
     logical(c_bool) :: ShowTablesRects
     logical(c_bool) :: ShowDrawCmdMesh
     logical(c_bool) :: ShowDrawCmdBoundingBoxes
     logical(c_bool) :: ShowDockingNodes
     integer(c_int) :: ShowWindowsRectsType
     integer(c_int) :: ShowTablesRectsType
  end type ImGuiMetricsConfig

  ! struct ImGuiStackLevelInfo
  ! {
  !     integer(c_int) :: ID
  !     ImS8 QueryFrameCount
  !     logical(c_bool) :: QuerySuccess
  !     ImGuiDataType DataType : 8
  !     character(c_char) :: Desc[57]
  ! }

  type, bind(c) :: ImVector_ImGuiStackLevelInfo
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiStackLevelInfo

  type, bind(c) :: ImGuiStackTool
     integer(c_int) :: LastActiveFrame
     integer(c_int) :: StackLevel
     integer(c_int) :: QueryId
     type(ImVector_ImGuiStackLevelInfo) :: Results
     logical(c_bool) :: CopyToClipboardOnCtrlC
     real(c_float) :: CopyToClipboardLastTime
  end type ImGuiStackTool

  ! typedef void (*ImGuiContextHookCallback)(ImGuiContext* ctx, ImGuiContextHook* hook)
  ! struct ImGuiContextHook
  ! {
  !     integer(c_int) :: HookId
  !     ImGuiContextHookType Type
  !     integer(c_int) :: Owner
  !     ImGuiContextHookCallback Callback
  !     type(c_ptr) :: UserData
  ! }

  type, bind(c) ::ImVector_ImGuiInputEvent
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiInputEvent

  type, bind(c) :: ImVector_ImGuiWindowStackData
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiWindowStackData

  type, bind(c) :: ImVector_ImGuiColorMod
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiColorMod

  type, bind(c) :: ImVector_ImGuiStyleMod
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiStyleMod

  type, bind(c) :: ImVector_ImGuiID
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiID

  type, bind(c) :: ImVector_ImGuiItemFlags
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiItemFlags

  type, bind(c) :: ImVector_ImGuiGroupData
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiGroupData

  type, bind(c) :: ImVector_ImGuiPopupData
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiPopupData

  type, bind(c) :: ImVector_ImGuiViewportPPtr
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiViewportPPtr

  type, bind(c) :: ImVector_unsigned_char
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_unsigned_char

  type, bind(c) :: ImVector_ImGuiListClipperData
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiListClipperData

  type, bind(c) :: ImVector_ImGuiTableTempData
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiTableTempData

  type, bind(c) :: ImVector_ImGuiTable
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiTable

  type, bind(c) :: ImPool_ImGuiTable
     type(ImVector_ImGuiTable) :: Buf
     type(ImGuiStorage) :: Map
     integer(c_int) :: FreeIdx
     integer(c_int) ::  AliveCount
  end type ImPool_ImGuiTable

  type, bind(c) :: ImVector_ImGuiTabBar
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiTabBar

  type, bind(c) :: ImPool_ImGuiTabBar
     type(ImVector_ImGuiTabBar) :: Buf
     type(ImGuiStorage) :: Map
     integer(c_int) :: FreeIdx
     integer(c_int) :: AliveCount
  end type ImPool_ImGuiTabBar

  type, bind(c) :: ImVector_ImGuiPtrOrIndex
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiPtrOrIndex

  type, bind(c) :: ImVector_ImGuiShrinkWidthItem
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiShrinkWidthItem

  type, bind(c) :: ImVector_ImGuiSettingsHandler
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiSettingsHandler

  type, bind(c) :: ImChunkStream_ImGuiWindowSettings
     type(ImVector_char) :: Buf
  end type ImChunkStream_ImGuiWindowSettings

  type, bind(c) :: ImChunkStream_ImGuiTableSettings
     type(ImVector_char) :: Buf
  end type ImChunkStream_ImGuiTableSettings

  type, bind(c) :: ImVector_ImGuiContextHook
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiContextHook

  type, bind(c) :: ImGuiContext
     logical(c_bool) :: Initialized
     logical(c_bool) :: FontAtlasOwnedByContext
     type(ImGuiIO) :: IO
     type(ImGuiPlatformIO) :: PlatformIO
     type(ImVector_ImGuiInputEvent) :: InputEventsQueue
     type(ImVector_ImGuiInputEvent) :: InputEventsTrail
     type(ImGuiStyle) :: Style
     integer(c_int) :: ConfigFlagsCurrFrame
     integer(c_int) :: ConfigFlagsLastFrame
     type(c_ptr) :: Font
     real(c_float) :: FontSize
     real(c_float) :: FontBaseSize
     type(ImDrawListSharedData) :: DrawListSharedData
     real(c_double) :: Time
     integer(c_int) :: FrameCount
     integer(c_int) :: FrameCountEnded
     integer(c_int) :: FrameCountPlatformEnded
     integer(c_int) :: FrameCountRendered
     logical(c_bool) :: WithinFrameScope
     logical(c_bool) :: WithinFrameScopeWithImplicitWindow
     logical(c_bool) :: WithinEndChild
     logical(c_bool) :: GcCompactAll
     logical(c_bool) :: TestEngineHookItems
     type(c_ptr) :: TestEngine
     type(ImVector_ImGuiWindowPtr) :: Windows
     type(ImVector_ImGuiWindowPtr) :: WindowsFocusOrder
     type(ImVector_ImGuiWindowPtr) :: WindowsTempSortBuffer
     type(ImVector_ImGuiWindowStackData) :: CurrentWindowStack
     type(ImGuiStorage) :: WindowsById
     integer(c_int) :: WindowsActiveCount
     type(ImVec2) :: WindowsHoverPadding
     type(c_ptr) :: CurrentWindow
     type(c_ptr) :: HoveredWindow
     type(c_ptr) :: HoveredWindowUnderMovingWindow
     type(c_ptr) :: HoveredDockNode
     type(c_ptr) :: MovingWindow
     type(c_ptr) :: WheelingWindow
     type(ImVec2) :: WheelingWindowRefMousePos
     real(c_float) :: WheelingWindowTimer
     integer(c_int) :: DebugHookIdInfo
     integer(c_int) :: HoveredId
     integer(c_int) :: HoveredIdPreviousFrame
     logical(c_bool) :: HoveredIdAllowOverlap
     logical(c_bool) :: HoveredIdUsingMouseWheel
     logical(c_bool) :: HoveredIdPreviousFrameUsingMouseWheel
     logical(c_bool) :: HoveredIdDisabled
     real(c_float) :: HoveredIdTimer
     real(c_float) :: HoveredIdNotActiveTimer
     integer(c_int) :: ActiveId
     integer(c_int) :: ActiveIdIsAlive
     real(c_float) :: ActiveIdTimer
     logical(c_bool) :: ActiveIdIsJustActivated
     logical(c_bool) :: ActiveIdAllowOverlap
     logical(c_bool) :: ActiveIdNoClearOnFocusLoss
     logical(c_bool) :: ActiveIdHasBeenPressedBefore
     logical(c_bool) :: ActiveIdHasBeenEditedBefore
     logical(c_bool) :: ActiveIdHasBeenEditedThisFrame
     type(ImVec2) :: ActiveIdClickOffset
     type(c_ptr) :: ActiveIdWindow
     integer(c_int) :: ActiveIdSource
     integer(c_int) :: ActiveIdMouseButton
     integer(c_int) :: ActiveIdPreviousFrame
     logical(c_bool) :: ActiveIdPreviousFrameIsAlive
     logical(c_bool) :: ActiveIdPreviousFrameHasBeenEditedBefore
     type(c_ptr) :: ActiveIdPreviousFrameWindow
     integer(c_int) :: LastActiveId
     real(c_float) :: LastActiveIdTimer
     logical(c_bool) :: ActiveIdUsingMouseWheel
     integer(c_int) :: ActiveIdUsingNavDirMask
     integer(c_int) :: ActiveIdUsingNavInputMask
     type(ImBitArrayForNamedKeys) :: ActiveIdUsingKeyInputMask
     integer(c_int) :: CurrentItemFlags
     type(ImGuiNextItemData) :: NextItemData
     type(ImGuiLastItemData) :: LastItemData
     type(ImGuiNextWindowData) :: NextWindowData
     type(ImVector_ImGuiColorMod) :: ColorStack
     type(ImVector_ImGuiStyleMod) :: StyleVarStack
     type(ImVector_ImFontPtr) :: FontStack
     type(ImVector_ImGuiID) :: FocusScopeStack
     type(ImVector_ImGuiItemFlags) :: ItemFlagsStack
     type(ImVector_ImGuiGroupData) :: GroupStack
     type(ImVector_ImGuiPopupData) :: OpenPopupStack
     type(ImVector_ImGuiPopupData) :: BeginPopupStack
     integer(c_int) :: BeginMenuCount
     type(ImVector_ImGuiViewportPPtr) :: Viewports
     real(c_float) :: CurrentDpiScale
     type(c_ptr) :: CurrentViewport
     type(c_ptr) :: MouseViewport
     type(c_ptr) :: MouseLastHoveredViewport
     integer(c_int) :: PlatformLastFocusedViewportId
     type(ImGuiPlatformMonitor) :: FallbackMonitor
     integer(c_int) :: ViewportFrontMostStampCount
     type(c_ptr) :: NavWindow
     integer(c_int) :: NavId
     integer(c_int) :: NavFocusScopeId
     integer(c_int) :: NavActivateId
     integer(c_int) :: NavActivateDownId
     integer(c_int) :: NavActivatePressedId
     integer(c_int) :: NavActivateInputId
     integer(c_int) :: NavActivateFlags
     integer(c_int) :: NavJustMovedToId
     integer(c_int) :: NavJustMovedToFocusScopeId
     integer(c_int) :: NavJustMovedToKeyMods
     integer(c_int) :: NavNextActivateId
     integer(c_int) :: NavNextActivateFlags
     integer(c_int) :: NavInputSource
     integer(c_int) :: NavLayer
     logical(c_bool) :: NavIdIsAlive
     logical(c_bool) :: NavMousePosDirty
     logical(c_bool) :: NavDisableHighlight
     logical(c_bool) :: NavDisableMouseHover
     logical(c_bool) :: NavAnyRequest
     logical(c_bool) :: NavInitRequest
     logical(c_bool) :: NavInitRequestFromMove
     integer(c_int) :: NavInitResultId
     type(ImRect) :: NavInitResultRectRel
     logical(c_bool) :: NavMoveSubmitted
     logical(c_bool) :: NavMoveScoringItems
     logical(c_bool) :: NavMoveForwardToNextFrame
     integer(c_int) :: NavMoveFlags
     integer(c_int) :: NavMoveScrollFlags
     integer(c_int) :: NavMoveKeyMods
     integer(c_int) :: NavMoveDir
     integer(c_int) :: NavMoveDirForDebug
     integer(c_int) :: NavMoveClipDir
     type(ImRect) :: NavScoringRect
     type(ImRect) :: NavScoringNoClipRect
     integer(c_int) :: NavScoringDebugCount
     integer(c_int) :: NavTabbingDir
     integer(c_int) :: NavTabbingCounter
     type(ImGuiNavItemData) :: NavMoveResultLocal
     type(ImGuiNavItemData) :: NavMoveResultLocalVisible
     type(ImGuiNavItemData) :: NavMoveResultOther
     type(ImGuiNavItemData) :: NavTabbingResultFirst
     type(c_ptr) :: NavWindowingTarget
     type(c_ptr) :: NavWindowingTargetAnim
     type(c_ptr) :: NavWindowingListWindow
     real(c_float) :: NavWindowingTimer
     real(c_float) :: NavWindowingHighlightAlpha
     logical(c_bool) :: NavWindowingToggleLayer
     real(c_float) :: DimBgRatio
     integer(c_int) :: MouseCursor
     logical(c_bool) :: DragDropActive
     logical(c_bool) :: DragDropWithinSource
     logical(c_bool) :: DragDropWithinTarget
     integer(c_int) :: DragDropSourceFlags
     integer(c_int) :: DragDropSourceFrameCount
     integer(c_int) :: DragDropMouseButton
     type(ImGuiPayload) :: DragDropPayload
     type(ImRect) :: DragDropTargetRect
     integer(c_int) :: DragDropTargetId
     integer(c_int) :: DragDropAcceptFlags
     real(c_float) :: DragDropAcceptIdCurrRectSurface
     integer(c_int) :: DragDropAcceptIdCurr
     integer(c_int) :: DragDropAcceptIdPrev
     integer(c_int) :: DragDropAcceptFrameCount
     integer(c_int) :: DragDropHoldJustPressedId
     type(ImVector_unsigned_char) :: DragDropPayloadBufHeap
     integer(c_signed_char) :: DragDropPayloadBufLocal(16)
     integer(c_int) :: ClipperTempDataStacked
     type(ImVector_ImGuiListClipperData) :: ClipperTempData
     type(c_ptr) :: CurrentTable
     integer(c_int) :: TablesTempDataStacked
     type(ImVector_ImGuiTableTempData) :: TablesTempData
     type(ImPool_ImGuiTable) :: Tables
     type(ImVector_float) :: TablesLastTimeActive
     type(ImVector_ImDrawChannel) :: DrawChannelsTempMergeBuffer
     type(c_ptr) :: CurrentTabBar
     type(ImPool_ImGuiTabBar) :: TabBars
     type(ImVector_ImGuiPtrOrIndex) :: CurrentTabBarStack
     type(ImVector_ImGuiShrinkWidthItem) :: ShrinkWidthBuffer
     type(ImVec2) :: MouseLastValidPos
     type(ImGuiInputTextState) :: InputTextState
     type(ImFont) :: InputTextPasswordFont
     integer(c_int) :: TempInputId
     integer(c_int) :: ColorEditOptions
     real(c_float) :: ColorEditLastHue
     real(c_float) :: ColorEditLastSat
     integer(c_int) :: ColorEditLastColor
     type(ImVec4) :: ColorPickerRef
     type(ImGuiComboPreviewData) :: ComboPreviewData
     real(c_float) :: SliderGrabClickOffset
     real(c_float) :: SliderCurrentAccum
     logical(c_bool) :: SliderCurrentAccumDirty
     logical(c_bool) :: DragCurrentAccumDirty
     real(c_float) :: DragCurrentAccum
     real(c_float) :: DragSpeedDefaultRatio
     real(c_float) :: ScrollbarClickDeltaToGrabCenter
     real(c_float) :: DisabledAlphaBackup
     integer(c_short) :: DisabledStackSize
     integer(c_short) :: TooltipOverrideCount
     real(c_float) :: TooltipSlowDelay
     type(ImVector_char) :: ClipboardHandlerData
     type(ImVector_ImGuiID) :: MenusIdSubmittedThisFrame
     type(ImGuiPlatformImeData) :: PlatformImeData
     type(ImGuiPlatformImeData) :: PlatformImeDataPrev
     integer(c_int) :: PlatformImeViewport
     character(c_char) :: PlatformLocaleDecimalPoint
     type(ImGuiDockContext) :: DockContext
     logical(c_bool) :: SettingsLoaded
     real(c_float) :: SettingsDirtyTimer
     type(ImGuiTextBuffer) :: SettingsIniData
     type(ImVector_ImGuiSettingsHandler) :: SettingsHandlers
     type(ImChunkStream_ImGuiWindowSettings) :: SettingsWindows
     type(ImChunkStream_ImGuiTableSettings) :: SettingsTables
     type(ImVector_ImGuiContextHook) :: Hooks
     integer(c_int) :: HookIdNext
     logical(c_bool) :: LogEnabled
     integer(c_int) :: LogType
     type(c_ptr) :: LogFile
     type(ImGuiTextBuffer) :: LogBuffer
     type(c_ptr) :: LogNextPrefix
     type(c_ptr) :: LogNextSuffix
     real(c_float) :: LogLinePosY
     logical(c_bool) :: LogLineFirstItem
     integer(c_int) :: LogDepthRef
     integer(c_int) :: LogDepthToExpand
     integer(c_int) :: LogDepthToExpandDefault
     integer(c_int) :: DebugLogFlags
     type(ImGuiTextBuffer) :: DebugLogBuf
     logical(c_bool) :: DebugItemPickerActive
     integer(c_int) :: DebugItemPickerBreakId
     type(ImGuiMetricsConfig) :: DebugMetricsConfig
     type(ImGuiStackTool) :: DebugStackTool
     real(c_float) :: FramerateSecPerFrame(120)
     integer(c_int) :: FramerateSecPerFrameIdx
     integer(c_int) :: FramerateSecPerFrameCount
     real(c_float) :: FramerateSecPerFrameAccum
     integer(c_int) :: WantCaptureMouseNextFrame
     integer(c_int) :: WantCaptureKeyboardNextFrame
     integer(c_int) :: WantTextInputNextFrame
     type(ImVector_char) :: TempBuffer
  end type ImGuiContext

  type, bind(c) :: ImGuiWindowTempData
     type(ImVec2) :: CursorPos
     type(ImVec2) :: CursorPosPrevLine
     type(ImVec2) :: CursorStartPos
     type(ImVec2) :: CursorMaxPos
     type(ImVec2) :: IdealMaxPos
     type(ImVec2) :: CurrLineSize
     type(ImVec2) :: PrevLineSize
     real(c_float) :: CurrLineTextBaseOffset
     real(c_float) :: PrevLineTextBaseOffset
     logical(c_bool) :: IsSameLine
     type(ImVec1) :: Indent
     type(ImVec1) :: ColumnsOffset
     type(ImVec1) :: GroupOffset
     type(ImVec2) :: CursorStartPosLossyness
     integer(c_int) :: NavLayerCurrent
     integer(c_short) :: NavLayersActiveMask
     integer(c_short) :: NavLayersActiveMaskNext
     integer(c_int) :: NavFocusScopeIdCurrent
     logical(c_bool) :: NavHideHighlightOneFrame
     logical(c_bool) :: NavHasScroll
     logical(c_bool) :: MenuBarAppending
     type(ImVec2) :: MenuBarOffset
     type(ImGuiMenuColumns) :: MenuColumns
     integer(c_int) :: TreeDepth
     integer(c_int) :: TreeJumpToParentOnPopMask
     type(ImVector_ImGuiWindowPtr) :: ChildWindows
     type(c_ptr) :: StateStorage
     type(c_ptr) :: CurrentColumns
     integer(c_int) :: CurrentTableIdx
     integer(c_int) :: LayoutType
     integer(c_int) :: ParentLayoutType
     real(c_float) :: ItemWidth
     real(c_float) :: TextWrapPos
     type(ImVector_float) :: ItemWidthStack
     type(ImVector_float) :: TextWrapPosStack
  end type ImGuiWindowTempData

  type, bind(c) :: ImVector_ImGuiOldColumns
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiOldColumns

  ! struct ImGuiWindow
  ! {
  !     char* Name
  !     integer(c_int) :: ID
  !     ImGuiWindowFlags Flags, FlagsPreviousFrame
  !     type(ImGuiWindowClass) :: WindowClass
  !     type(c_ptr) :: Viewport
  !     integer(c_int) :: ViewportId
  !     type(ImVec2) :: ViewportPos
  !     integer(c_int) :: ViewportAllowPlatformMonitorExtend
  !     type(ImVec2) :: Pos
  !     type(ImVec2) :: Size
  !     type(ImVec2) :: SizeFull
  !     type(ImVec2) :: ContentSize
  !     type(ImVec2) :: ContentSizeIdeal
  !     type(ImVec2) :: ContentSizeExplicit
  !     type(ImVec2) :: WindowPadding
  !     real(c_float) :: WindowRounding
  !     real(c_float) :: WindowBorderSize
  !     integer(c_int) :: NameBufLen
  !     integer(c_int) :: MoveId
  !     integer(c_int) :: TabId
  !     integer(c_int) :: ChildId
  !     type(ImVec2) :: Scroll
  !     type(ImVec2) :: ScrollMax
  !     type(ImVec2) :: ScrollTarget
  !     type(ImVec2) :: ScrollTargetCenterRatio
  !     type(ImVec2) :: ScrollTargetEdgeSnapDist
  !     type(ImVec2) :: ScrollbarSizes
  !     logical(c_bool) :: ScrollbarX, ScrollbarY
  !     logical(c_bool) :: ViewportOwned
  !     logical(c_bool) :: Active
  !     logical(c_bool) :: WasActive
  !     logical(c_bool) :: WriteAccessed
  !     logical(c_bool) :: Collapsed
  !     logical(c_bool) :: WantCollapseToggle
  !     logical(c_bool) :: SkipItems
  !     logical(c_bool) :: Appearing
  !     logical(c_bool) :: Hidden
  !     logical(c_bool) :: IsFallbackWindow
  !     logical(c_bool) :: IsExplicitChild
  !     logical(c_bool) :: HasCloseButton
  !     integer(c_signed_char) :: ResizeBorderHeld
  !     integer(c_short) :: BeginCount
  !     integer(c_short) :: BeginOrderWithinParent
  !     integer(c_short) :: BeginOrderWithinContext
  !     integer(c_short) :: FocusOrder
  !     integer(c_int) :: PopupId
  !     ImS8 AutoFitFramesX, AutoFitFramesY
  !     ImS8 AutoFitChildAxises
  !     logical(c_bool) :: AutoFitOnlyGrows
  !     integer(c_int) :: AutoPosLastDirection
  !     ImS8 HiddenFramesCanSkipItems
  !     ImS8 HiddenFramesCannotSkipItems
  !     ImS8 HiddenFramesForRenderOnly
  !     ImS8 DisableInputsFrames
  !     integer(c_int) :: SetWindowPosAllowFlags : 8
  !     integer(c_int) :: SetWindowSizeAllowFlags : 8
  !     integer(c_int) :: SetWindowCollapsedAllowFlags : 8
  !     integer(c_int) :: SetWindowDockAllowFlags : 8
  !     type(ImVec2) :: SetWindowPosVal
  !     type(ImVec2) :: SetWindowPosPivot
  !     type(ImVector_ImGuiID) :: IDStack
  !     type(ImGuiWindowTempData) :: DC
  !     type(ImRect) :: OuterRectClipped
  !     type(ImRect) :: InnerRect
  !     type(ImRect) :: InnerClipRect
  !     type(ImRect) :: WorkRect
  !     type(ImRect) :: ParentWorkRect
  !     type(ImRect) :: ClipRect
  !     type(ImRect) :: ContentRegionRect
  !     type(ImVec2ih) :: HitTestHoleSize
  !     type(ImVec2ih) :: HitTestHoleOffset
  !     integer(c_int) :: LastFrameActive
  !     integer(c_int) :: LastFrameJustFocused
  !     real(c_float) :: LastTimeActive
  !     real(c_float) :: ItemWidthDefault
  !     type(ImGuiStorage) :: StateStorage
  !     type(ImVector_ImGuiOldColumns) :: ColumnsStorage
  !     real(c_float) :: FontWindowScale
  !     real(c_float) :: FontDpiScale
  !     integer(c_int) :: SettingsOffset
  !     type(c_ptr) :: DrawList
  !     type(ImDrawList) :: DrawListInst
  !     type(c_ptr) :: ParentWindow
  !     type(c_ptr) :: ParentWindowInBeginStack
  !     type(c_ptr) :: RootWindow
  !     type(c_ptr) :: RootWindowPopupTree
  !     type(c_ptr) :: RootWindowDockTree
  !     type(c_ptr) :: RootWindowForTitleBarHighlight
  !     type(c_ptr) :: RootWindowForNav
  !     type(c_ptr) :: NavLastChildNavWindow
  !     integer(c_int) :: NavLastIds[ImGuiNavLayer_COUNT]
  !     type(ImRect) :: NavRectRel[ImGuiNavLayer_COUNT]
  !     integer(c_int) :: MemoryDrawListIdxCapacity
  !     integer(c_int) :: MemoryDrawListVtxCapacity
  !     logical(c_bool) :: MemoryCompacted
  !     logical(c_bool) :: DockIsActive :1
  !     logical(c_bool) :: DockNodeIsVisible :1
  !     logical(c_bool) :: DockTabIsVisible :1
  !     logical(c_bool) :: DockTabWantClose :1
  !     integer(c_short) :: DockOrder
  !     ImGuiWindowDockStyle DockStyle
  !     type(c_ptr) :: DockNode
  !     type(c_ptr) :: DockNodeAsHost
  !     integer(c_int) :: DockId
  !     integer(c_int) :: DockTabItemStatusFlags
  !     type(ImRect) :: DockTabItemRect
  ! }
  ! struct ImGuiTabItem
  ! {
  !     integer(c_int) :: ID
  !     integer(c_int) :: Flags
  !     type(c_ptr) :: Window
  !     integer(c_int) :: LastFrameVisible
  !     integer(c_int) :: LastFrameSelected
  !     real(c_float) :: Offset
  !     real(c_float) :: Width
  !     real(c_float) :: ContentWidth
  !     real(c_float) :: RequestedWidth
  !     ImS32 NameOffset
  !     ImS16 BeginOrder
  !     ImS16 IndexDuringLayout
  !     logical(c_bool) :: WantClose
  ! }
  ! typedef struct ImVector_ImGuiTabItem {int Size;int Capacity;ImGuiTabItem* Data;} ImVector_ImGuiTabItem;

  ! struct ImGuiTabBar
  ! {
  !     ImVector_ImGuiTabItem Tabs
  !     ImGuiTabBarFlags Flags
  !     integer(c_int) :: ID
  !     integer(c_int) :: SelectedTabId
  !     integer(c_int) :: NextSelectedTabId
  !     integer(c_int) :: VisibleTabId
  !     integer(c_int) :: CurrFrameVisible
  !     integer(c_int) :: PrevFrameVisible
  !     type(ImRect) :: BarRect
  !     real(c_float) :: CurrTabsContentsHeight
  !     real(c_float) :: PrevTabsContentsHeight
  !     real(c_float) :: WidthAllTabs
  !     real(c_float) :: WidthAllTabsIdeal
  !     real(c_float) :: ScrollingAnim
  !     real(c_float) :: ScrollingTarget
  !     real(c_float) :: ScrollingTargetDistToVisibility
  !     real(c_float) :: ScrollingSpeed
  !     real(c_float) :: ScrollingRectMinX
  !     real(c_float) :: ScrollingRectMaxX
  !     integer(c_int) :: ReorderRequestTabId
  !     ImS16 ReorderRequestOffset
  !     ImS8 BeginCount
  !     logical(c_bool) :: WantLayout
  !     logical(c_bool) :: VisibleTabWasSubmitted
  !     logical(c_bool) :: TabsAddedNew
  !     ImS16 TabsActiveCount
  !     ImS16 LastTabItemIdx
  !     real(c_float) :: ItemSpacingY
  !     type(ImVec2) :: FramePadding
  !     type(ImVec2) :: BackupCursorPos
  !     type(ImGuiTextBuffer) :: TabsNames
  ! }
  ! struct ImGuiTableColumn
  ! {
  !     ImGuiTableColumnFlags Flags
  !     real(c_float) :: WidthGiven
  !     real(c_float) :: MinX
  !     real(c_float) :: MaxX
  !     real(c_float) :: WidthRequest
  !     real(c_float) :: WidthAuto
  !     real(c_float) :: StretchWeight
  !     real(c_float) :: InitStretchWeightOrWidth
  !     type(ImRect) :: ClipRect
  !     integer(c_int) :: UserID
  !     real(c_float) :: WorkMinX
  !     real(c_float) :: WorkMaxX
  !     real(c_float) :: ItemWidth
  !     real(c_float) :: ContentMaxXFrozen
  !     real(c_float) :: ContentMaxXUnfrozen
  !     real(c_float) :: ContentMaxXHeadersUsed
  !     real(c_float) :: ContentMaxXHeadersIdeal
  !     ImS16 NameOffset
  !     ImGuiTableColumnIdx DisplayOrder
  !     ImGuiTableColumnIdx IndexWithinEnabledSet
  !     ImGuiTableColumnIdx PrevEnabledColumn
  !     ImGuiTableColumnIdx NextEnabledColumn
  !     ImGuiTableColumnIdx SortOrder
  !     ImGuiTableDrawChannelIdx DrawChannelCurrent
  !     ImGuiTableDrawChannelIdx DrawChannelFrozen
  !     ImGuiTableDrawChannelIdx DrawChannelUnfrozen
  !     logical(c_bool) :: IsEnabled
  !     logical(c_bool) :: IsUserEnabled
  !     logical(c_bool) :: IsUserEnabledNextFrame
  !     logical(c_bool) :: IsVisibleX
  !     logical(c_bool) :: IsVisibleY
  !     logical(c_bool) :: IsRequestOutput
  !     logical(c_bool) :: IsSkipItems
  !     logical(c_bool) :: IsPreserveWidthAuto
  !     ImS8 NavLayerCurrent
  !     integer(c_signed_char) :: AutoFitQueue
  !     integer(c_signed_char) :: CannotSkipItemsQueue
  !     integer(c_signed_char) :: SortDirection : 2
  !     integer(c_signed_char) :: SortDirectionsAvailCount : 2
  !     integer(c_signed_char) :: SortDirectionsAvailMask : 4
  !     integer(c_signed_char) :: SortDirectionsAvailList
  ! }
  ! struct ImGuiTableCellData
  ! {
  !     integer(c_int) :: BgColor
  !     ImGuiTableColumnIdx Column
  ! }
  ! struct ImGuiTableInstanceData
  ! {
  !     real(c_float) :: LastOuterHeight
  !     real(c_float) :: LastFirstRowHeight
  ! };
  ! typedef struct ImSpan_ImGuiTableColumn {ImGuiTableColumn* Data;ImGuiTableColumn* DataEnd;} ImSpan_ImGuiTableColumn;

  ! typedef struct ImSpan_ImGuiTableColumnIdx {ImGuiTableColumnIdx* Data;ImGuiTableColumnIdx* DataEnd;} ImSpan_ImGuiTableColumnIdx;

  ! typedef struct ImSpan_ImGuiTableCellData {ImGuiTableCellData* Data;ImGuiTableCellData* DataEnd;} ImSpan_ImGuiTableCellData;

  ! typedef struct ImVector_ImGuiTableInstanceData {int Size;int Capacity;ImGuiTableInstanceData* Data;} ImVector_ImGuiTableInstanceData;

  ! typedef struct ImVector_ImGuiTableColumnSortSpecs {int Size;int Capacity;ImGuiTableColumnSortSpecs* Data;} ImVector_ImGuiTableColumnSortSpecs;

  ! struct ImGuiTable
  ! {
  !     integer(c_int) :: ID
  !     ImGuiTableFlags Flags
  !     type(c_ptr) :: RawData
  !     ImGuiTableTempData* TempData
  !     ImSpan_ImGuiTableColumn Columns
  !     ImSpan_ImGuiTableColumnIdx DisplayOrderToIndex
  !     ImSpan_ImGuiTableCellData RowCellData
  !     ImU64 EnabledMaskByDisplayOrder
  !     ImU64 EnabledMaskByIndex
  !     ImU64 VisibleMaskByIndex
  !     ImU64 RequestOutputMaskByIndex
  !     ImGuiTableFlags SettingsLoadedFlags
  !     integer(c_int) :: SettingsOffset
  !     integer(c_int) :: LastFrameActive
  !     integer(c_int) :: ColumnsCount
  !     integer(c_int) :: CurrentRow
  !     integer(c_int) :: CurrentColumn
  !     ImS16 InstanceCurrent
  !     ImS16 InstanceInteracted
  !     real(c_float) :: RowPosY1
  !     real(c_float) :: RowPosY2
  !     real(c_float) :: RowMinHeight
  !     real(c_float) :: RowTextBaseline
  !     real(c_float) :: RowIndentOffsetX
  !     ImGuiTableRowFlags RowFlags : 16
  !     ImGuiTableRowFlags LastRowFlags : 16
  !     integer(c_int) :: RowBgColorCounter
  !     integer(c_int) :: RowBgColor[2]
  !     integer(c_int) :: BorderColorStrong
  !     integer(c_int) :: BorderColorLight
  !     real(c_float) :: BorderX1
  !     real(c_float) :: BorderX2
  !     real(c_float) :: HostIndentX
  !     real(c_float) :: MinColumnWidth
  !     real(c_float) :: OuterPaddingX
  !     real(c_float) :: CellPaddingX
  !     real(c_float) :: CellPaddingY
  !     real(c_float) :: CellSpacingX1
  !     real(c_float) :: CellSpacingX2
  !     real(c_float) :: InnerWidth
  !     real(c_float) :: ColumnsGivenWidth
  !     real(c_float) :: ColumnsAutoFitWidth
  !     real(c_float) :: ColumnsStretchSumWeights
  !     real(c_float) :: ResizedColumnNextWidth
  !     real(c_float) :: ResizeLockMinContentsX2
  !     real(c_float) :: RefScale
  !     type(ImRect) :: OuterRect
  !     type(ImRect) :: InnerRect
  !     type(ImRect) :: WorkRect
  !     type(ImRect) :: InnerClipRect
  !     type(ImRect) :: BgClipRect
  !     type(ImRect) :: Bg0ClipRectForDrawCmd
  !     type(ImRect) :: Bg2ClipRectForDrawCmd
  !     type(ImRect) :: HostClipRect
  !     type(ImRect) :: HostBackupInnerClipRect
  !     type(c_ptr) :: OuterWindow
  !     type(c_ptr) :: InnerWindow
  !     type(ImGuiTextBuffer) :: ColumnsNames
  !     ImDrawListSplitter* DrawSplitter
  !     ImGuiTableInstanceData InstanceDataFirst
  !     ImVector_ImGuiTableInstanceData InstanceDataExtra
  !     ImGuiTableColumnSortSpecs SortSpecsSingle
  !     ImVector_ImGuiTableColumnSortSpecs SortSpecsMulti
  !     ImGuiTableSortSpecs SortSpecs
  !     ImGuiTableColumnIdx SortSpecsCount
  !     ImGuiTableColumnIdx ColumnsEnabledCount
  !     ImGuiTableColumnIdx ColumnsEnabledFixedCount
  !     ImGuiTableColumnIdx DeclColumnsCount
  !     ImGuiTableColumnIdx HoveredColumnBody
  !     ImGuiTableColumnIdx HoveredColumnBorder
  !     ImGuiTableColumnIdx AutoFitSingleColumn
  !     ImGuiTableColumnIdx ResizedColumn
  !     ImGuiTableColumnIdx LastResizedColumn
  !     ImGuiTableColumnIdx HeldHeaderColumn
  !     ImGuiTableColumnIdx ReorderColumn
  !     ImGuiTableColumnIdx ReorderColumnDir
  !     ImGuiTableColumnIdx LeftMostEnabledColumn
  !     ImGuiTableColumnIdx RightMostEnabledColumn
  !     ImGuiTableColumnIdx LeftMostStretchedColumn
  !     ImGuiTableColumnIdx RightMostStretchedColumn
  !     ImGuiTableColumnIdx ContextPopupColumn
  !     ImGuiTableColumnIdx FreezeRowsRequest
  !     ImGuiTableColumnIdx FreezeRowsCount
  !     ImGuiTableColumnIdx FreezeColumnsRequest
  !     ImGuiTableColumnIdx FreezeColumnsCount
  !     ImGuiTableColumnIdx RowCellDataCurrent
  !     ImGuiTableDrawChannelIdx DummyDrawChannel
  !     ImGuiTableDrawChannelIdx Bg2DrawChannelCurrent
  !     ImGuiTableDrawChannelIdx Bg2DrawChannelUnfrozen
  !     logical(c_bool) :: IsLayoutLocked
  !     logical(c_bool) :: IsInsideRow
  !     logical(c_bool) :: IsInitializing
  !     logical(c_bool) :: IsSortSpecsDirty
  !     logical(c_bool) :: IsUsingHeaders
  !     logical(c_bool) :: IsContextPopupOpen
  !     logical(c_bool) :: IsSettingsRequestLoad
  !     logical(c_bool) :: IsSettingsDirty
  !     logical(c_bool) :: IsDefaultDisplayOrder
  !     logical(c_bool) :: IsResetAllRequest
  !     logical(c_bool) :: IsResetDisplayOrderRequest
  !     logical(c_bool) :: IsUnfrozenRows
  !     logical(c_bool) :: IsDefaultSizingPolicy
  !     logical(c_bool) :: MemoryCompacted
  !     logical(c_bool) :: HostSkipItems
  ! }
  ! struct ImGuiTableTempData
  ! {
  !     integer(c_int) :: TableIndex
  !     real(c_float) :: LastTimeActive
  !     type(ImVec2) :: UserOuterSize
  !     type(ImDrawListSplitter) :: DrawSplitter
  !     type(ImRect) :: HostBackupWorkRect
  !     type(ImRect) :: HostBackupParentWorkRect
  !     type(ImVec2) :: HostBackupPrevLineSize
  !     type(ImVec2) :: HostBackupCurrLineSize
  !     type(ImVec2) :: HostBackupCursorMaxPos
  !     type(ImVec1) :: HostBackupColumnsOffset
  !     real(c_float) :: HostBackupItemWidth
  !     integer(c_int) :: HostBackupItemWidthStackSize
  ! }
  ! struct ImGuiTableColumnSettings
  ! {
  !     real(c_float) :: WidthOrWeight
  !     integer(c_int) :: UserID
  !     ImGuiTableColumnIdx Index
  !     ImGuiTableColumnIdx DisplayOrder
  !     ImGuiTableColumnIdx SortOrder
  !     integer(c_signed_char) :: SortDirection : 2
  !     integer(c_signed_char) :: IsEnabled : 1
  !     integer(c_signed_char) :: IsStretch : 1
  ! }
  ! struct ImGuiTableSettings
  ! {
  !     integer(c_int) :: ID
  !     ImGuiTableFlags SaveFlags
  !     real(c_float) :: RefScale
  !     ImGuiTableColumnIdx ColumnsCount
  !     ImGuiTableColumnIdx ColumnsCountMax
  !     logical(c_bool) :: WantApply
  ! }
  ! struct ImFontBuilderIO
  ! {
  !     bool (*FontBuilder_Build)(ImFontAtlas* atlas)
  ! }
