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
     ! type(ImVec4) :: Colors[ImGuiCol_COUNT]
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

  ! struct ImGuiInputTextCallbackData
  ! {
  !     ImGuiInputTextFlags EventFlag;
  !     ImGuiInputTextFlags Flags;
  !     type(c_ptr) :: UserData;
  !     ImWchar EventChar;
  !     ImGuiKey EventKey;
  !     char* Buf;
  !     integer(c_int) :: BufTextLen;
  !     integer(c_int) :: BufSize;
  !     logical(c_bool) :: BufDirty;
  !     integer(c_int) :: CursorPos;
  !     integer(c_int) :: SelectionStart;
  !     integer(c_int) :: SelectionEnd;
  ! };
  ! struct ImGuiSizeCallbackData
  ! {
  !     type(c_ptr) :: UserData;
  !     type(ImVec2) :: Pos;
  !     type(ImVec2) :: CurrentSize;
  !     type(ImVec2) :: DesiredSize;
  ! };
  ! struct ImGuiWindowClass
  ! {
  !     ImGuiID ClassId;
  !     ImGuiID ParentViewportId;
  !     ImGuiViewportFlags ViewportFlagsOverrideSet;
  !     ImGuiViewportFlags ViewportFlagsOverrideClear;
  !     ImGuiTabItemFlags TabItemFlagsOverrideSet;
  !     ImGuiDockNodeFlags DockNodeFlagsOverrideSet;
  !     logical(c_bool) :: DockingAlwaysTabBar;
  !     logical(c_bool) :: DockingAllowUnclassed;
  ! };
  ! struct ImGuiPayload
  ! {
  !     type(c_ptr) :: Data;
  !     integer(c_int) :: DataSize;
  !     ImGuiID SourceId;
  !     ImGuiID SourceParentId;
  !     integer(c_int) :: DataFrameCount;
  !     character(c_char) :: DataType[32 + 1];
  !     logical(c_bool) :: Preview;
  !     logical(c_bool) :: Delivery;
  ! };
  ! struct ImGuiTableColumnSortSpecs
  ! {
  !     ImGuiID ColumnUserID;
  !     ImS16 ColumnIndex;
  !     ImS16 SortOrder;
  !     ImGuiSortDirection SortDirection : 8;
  ! };
  ! struct ImGuiTableSortSpecs
  ! {
  !     const ImGuiTableColumnSortSpecs* Specs;
  !     integer(c_int) :: SpecsCount;
  !     logical(c_bool) :: SpecsDirty;
  ! };
  ! struct ImGuiOnceUponAFrame
  ! {
  !      integer(c_int) :: RefFrame;
  ! };
  ! struct ImGuiTextRange
  ! {
  !         const char* b;
  !         const char* e;
  ! };
  ! typedef struct ImVector_ImGuiTextRange {int Size;int Capacity;ImGuiTextRange* Data;} ImVector_ImGuiTextRange;
  ! struct ImGuiTextFilter
  ! {
  !     character(c_char) :: InputBuf[256];
  !     ImVector_ImGuiTextRange Filters;
  !     integer(c_int) :: CountGrep;
  ! };
  ! typedef struct ImVector_char {int Size;int Capacity;char* Data;} ImVector_char;
  ! struct ImGuiTextBuffer
  ! {
  !     ImVector_char Buf;
  ! };
  ! struct ImGuiStoragePair
  ! {
  !         ImGuiID key;
  !         union { int val_i; float val_f; void* val_p; };
  ! };
  ! typedef struct ImVector_ImGuiStoragePair {int Size;int Capacity;ImGuiStoragePair* Data;} ImVector_ImGuiStoragePair;
  ! struct ImGuiStorage
  ! {
  !     ImVector_ImGuiStoragePair Data;
  ! };
  ! struct ImGuiListClipper
  ! {
  !     integer(c_int) :: DisplayStart;
  !     integer(c_int) :: DisplayEnd;
  !     integer(c_int) :: ItemsCount;
  !     real(c_float) :: ItemsHeight;
  !     real(c_float) :: StartPosY;
  !     type(c_ptr) :: TempData;
  ! };
  ! struct ImColor
  ! {
  !     type(ImVec4) :: Value;
  ! };
  ! typedef void (*ImDrawCallback)(const ImDrawList* parent_list, const ImDrawCmd* cmd);
  ! struct ImDrawCmd
  ! {
  !     type(ImVec4) :: ClipRect;
  !     ImTextureID TextureId;
  !     integer(c_int) :: VtxOffset;
  !     integer(c_int) :: IdxOffset;
  !     integer(c_int) :: ElemCount;
  !     ImDrawCallback UserCallback;
  !     type(c_ptr) :: UserCallbackData;
  ! };
  ! struct ImDrawVert
  ! {
  !     type(ImVec2) :: pos;
  !     type(ImVec2) :: uv;
  !     ImU32 col;
  ! };
  ! struct ImDrawCmdHeader
  ! {
  !     type(ImVec4) :: ClipRect;
  !     ImTextureID TextureId;
  !     integer(c_int) :: VtxOffset;
  ! };
  ! typedef struct ImVector_ImDrawCmd {int Size;int Capacity;ImDrawCmd* Data;} ImVector_ImDrawCmd;

  ! typedef struct ImVector_ImDrawIdx {int Size;int Capacity;ImDrawIdx* Data;} ImVector_ImDrawIdx;
  ! struct ImDrawChannel
  ! {
  !     ImVector_ImDrawCmd _CmdBuffer;
  !     ImVector_ImDrawIdx _IdxBuffer;
  ! };
  ! typedef struct ImVector_ImDrawChannel {int Size;int Capacity;ImDrawChannel* Data;} ImVector_ImDrawChannel;

  ! struct ImDrawListSplitter
  ! {
  !     integer(c_int) :: _Current;
  !     integer(c_int) :: _Count;
  !     ImVector_ImDrawChannel _Channels;
  ! };
  ! typedef struct ImVector_ImDrawVert {int Size;int Capacity;ImDrawVert* Data;} ImVector_ImDrawVert;

  ! typedef struct ImVector_ImVec4 {int Size;int Capacity;ImVec4* Data;} ImVector_ImVec4;

  ! typedef struct ImVector_ImTextureID {int Size;int Capacity;ImTextureID* Data;} ImVector_ImTextureID;

  ! typedef struct ImVector_ImVec2 {int Size;int Capacity;ImVec2* Data;} ImVector_ImVec2;

  ! struct ImDrawList
  ! {
  !     ImVector_ImDrawCmd CmdBuffer;
  !     ImVector_ImDrawIdx IdxBuffer;
  !     ImVector_ImDrawVert VtxBuffer;
  !     ImDrawListFlags Flags;
  !     integer(c_int) :: _VtxCurrentIdx;
  !     const ImDrawListSharedData* _Data;
  !     const char* _OwnerName;
  !     ImDrawVert* _VtxWritePtr;
  !     ImDrawIdx* _IdxWritePtr;
  !     ImVector_ImVec4 _ClipRectStack;
  !     ImVector_ImTextureID _TextureIdStack;
  !     ImVector_ImVec2 _Path;
  !     ImDrawCmdHeader _CmdHeader;
  !     ImDrawListSplitter _Splitter;
  !     real(c_float) :: _FringeScale;
  ! };
  ! struct ImDrawData
  ! {
  !     logical(c_bool) :: Valid;
  !     integer(c_int) :: CmdListsCount;
  !     integer(c_int) :: TotalIdxCount;
  !     integer(c_int) :: TotalVtxCount;
  !     ImDrawList** CmdLists;
  !     type(ImVec2) :: DisplayPos;
  !     type(ImVec2) :: DisplaySize;
  !     type(ImVec2) :: FramebufferScale;
  !     ImGuiViewport* OwnerViewport;
  ! };
  ! struct ImFontConfig
  ! {
  !     type(c_ptr) :: FontData;
  !     integer(c_int) :: FontDataSize;
  !     logical(c_bool) :: FontDataOwnedByAtlas;
  !     integer(c_int) :: FontNo;
  !     real(c_float) :: SizePixels;
  !     integer(c_int) :: OversampleH;
  !     integer(c_int) :: OversampleV;
  !     logical(c_bool) :: PixelSnapH;
  !     type(ImVec2) :: GlyphExtraSpacing;
  !     type(ImVec2) :: GlyphOffset;
  !     const ImWchar* GlyphRanges;
  !     real(c_float) :: GlyphMinAdvanceX;
  !     real(c_float) :: GlyphMaxAdvanceX;
  !     logical(c_bool) :: MergeMode;
  !     integer(c_int) :: FontBuilderFlags;
  !     real(c_float) :: RasterizerMultiply;
  !     ImWchar EllipsisChar;
  !     character(c_char) :: Name[40];
  !     ImFont* DstFont;
  ! };
  ! struct ImFontGlyph
  ! {
  !     integer(c_int) :: Colored : 1;
  !     integer(c_int) :: Visible : 1;
  !     integer(c_int) :: Codepoint : 30;
  !     real(c_float) :: AdvanceX;
  !     real(c_float) :: X0, Y0, X1, Y1;
  !     real(c_float) :: U0, V0, U1, V1;
  ! };
  ! typedef struct ImVector_ImU32 {int Size;int Capacity;ImU32* Data;} ImVector_ImU32;

  ! struct ImFontGlyphRangesBuilder
  ! {
  !     ImVector_ImU32 UsedChars;
  ! };
  ! typedef struct ImFontAtlasCustomRect ImFontAtlasCustomRect;
  ! struct ImFontAtlasCustomRect
  ! {
  !     integer(c_short) :: Width, Height;
  !     integer(c_short) :: X, Y;
  !     integer(c_int) :: GlyphID;
  !     real(c_float) :: GlyphAdvanceX;
  !     type(ImVec2) :: GlyphOffset;
  !     ImFont* Font;
  ! };
  ! typedef struct ImVector_ImFontPtr {int Size;int Capacity;ImFont** Data;} ImVector_ImFontPtr;

  ! typedef struct ImVector_ImFontAtlasCustomRect {int Size;int Capacity;ImFontAtlasCustomRect* Data;} ImVector_ImFontAtlasCustomRect;

  ! typedef struct ImVector_ImFontConfig {int Size;int Capacity;ImFontConfig* Data;} ImVector_ImFontConfig;

  ! struct ImFontAtlas
  ! {
  !     ImFontAtlasFlags Flags;
  !     ImTextureID TexID;
  !     integer(c_int) :: TexDesiredWidth;
  !     integer(c_int) :: TexGlyphPadding;
  !     logical(c_bool) :: Locked;
  !     logical(c_bool) :: TexReady;
  !     logical(c_bool) :: TexPixelsUseColors;
  !     unsigned char* TexPixelsAlpha8;
  !     unsigned int* TexPixelsRGBA32;
  !     integer(c_int) :: TexWidth;
  !     integer(c_int) :: TexHeight;
  !     type(ImVec2) :: TexUvScale;
  !     type(ImVec2) :: TexUvWhitePixel;
  !     ImVector_ImFontPtr Fonts;
  !     ImVector_ImFontAtlasCustomRect CustomRects;
  !     ImVector_ImFontConfig ConfigData;
  !     type(ImVec4) :: TexUvLines[(63) + 1];
  !     const ImFontBuilderIO* FontBuilderIO;
  !     integer(c_int) :: FontBuilderFlags;
  !     integer(c_int) :: PackIdMouseCursors;
  !     integer(c_int) :: PackIdLines;
  ! };
  ! typedef struct ImVector_float {int Size;int Capacity;float* Data;} ImVector_float;

  ! typedef struct ImVector_ImFontGlyph {int Size;int Capacity;ImFontGlyph* Data;} ImVector_ImFontGlyph;

  ! struct ImFont
  ! {
  !     ImVector_float IndexAdvanceX;
  !     real(c_float) :: FallbackAdvanceX;
  !     real(c_float) :: FontSize;
  !     ImVector_ImWchar IndexLookup;
  !     ImVector_ImFontGlyph Glyphs;
  !     const ImFontGlyph* FallbackGlyph;
  !     ImFontAtlas* ContainerAtlas;
  !     const ImFontConfig* ConfigData;
  !     integer(c_short) :: ConfigDataCount;
  !     ImWchar FallbackChar;
  !     ImWchar EllipsisChar;
  !     ImWchar DotChar;
  !     logical(c_bool) :: DirtyLookupTables;
  !     real(c_float) :: Scale;
  !     real(c_float) :: Ascent, Descent;
  !     integer(c_int) :: MetricsTotalSurface;
  !     ImU8 Used4kPagesMap[(0xFFFF +1)/4096/8];
  ! };
  ! struct ImGuiViewport
  ! {
  !     ImGuiID ID;
  !     ImGuiViewportFlags Flags;
  !     type(ImVec2) :: Pos;
  !     type(ImVec2) :: Size;
  !     type(ImVec2) :: WorkPos;
  !     type(ImVec2) :: WorkSize;
  !     real(c_float) :: DpiScale;
  !     ImGuiID ParentViewportId;
  !     ImDrawData* DrawData;
  !     type(c_ptr) :: RendererUserData;
  !     type(c_ptr) :: PlatformUserData;
  !     type(c_ptr) :: PlatformHandle;
  !     type(c_ptr) :: PlatformHandleRaw;
  !     logical(c_bool) :: PlatformRequestMove;
  !     logical(c_bool) :: PlatformRequestResize;
  !     logical(c_bool) :: PlatformRequestClose;
  ! };
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

  ! struct ImGuiPlatformMonitor
  ! {
  !     type(ImVec2) :: MainPos;
  !     type(ImVec2) :: MainSize;
  !     type(ImVec2) :: WorkPos;
  !     type(ImVec2) :: WorkSize;
  !     real(c_float) :: DpiScale;
  ! };
  ! struct ImGuiPlatformImeData
  ! {
  !     logical(c_bool) :: WantVisible;
  !     type(ImVec2) :: InputPos;
  !     real(c_float) :: InputLineHeight;
  ! };
  ! typedef void (*ImGuiErrorLogCallback)(void* user_data, const char* fmt, ...);
  ! struct StbUndoRecord
  ! {
  !    integer(c_int) :: where;
  !    integer(c_int) :: insert_length;
  !    integer(c_int) :: delete_length;
  !    integer(c_int) :: char_storage;
  ! };
  ! struct StbUndoState
  ! {
  !    StbUndoRecord undo_rec [99];
  !    ImWchar undo_char[999];
  !    integer(c_short) :: undo_point, redo_point;
  !    integer(c_int) :: undo_char_point, redo_char_point;
  ! };
  ! struct STB_TexteditState
  ! {
  !    integer(c_int) :: cursor;
  !    integer(c_int) :: select_start;
  !    integer(c_int) :: select_end;
  !    integer(c_signed_char) :: insert_mode;
  !    integer(c_int) :: row_count_per_page;
  !    integer(c_signed_char) :: cursor_at_end_of_line;
  !    integer(c_signed_char) :: initialized;
  !    integer(c_signed_char) :: has_preferred_x;
  !    integer(c_signed_char) :: single_line;
  !    integer(c_signed_char) :: padding1, padding2, padding3;
  !    real(c_float) :: preferred_x;
  !    StbUndoState undostate;
  ! };
  ! struct StbTexteditRow
  ! {
  !    real(c_float) :: x0,x1;
  !    real(c_float) :: baseline_y_delta;
  !    real(c_float) :: ymin,ymax;
  !    integer(c_int) :: num_chars;
  ! };
  ! struct ImVec1
  ! {
  !     real(c_float) :: x;
  ! };
  ! struct ImVec2ih
  ! {
  !     integer(c_short) :: x, y;
  ! };
  ! struct ImRect
  ! {
  !     type(ImVec2) :: Min;
  !     type(ImVec2) :: Max;
  ! };
  ! struct ImBitVector
  ! {
  !     ImVector_ImU32 Storage;
  ! };
  ! struct ImDrawListSharedData
  ! {
  !     type(ImVec2) :: TexUvWhitePixel;
  !     ImFont* Font;
  !     real(c_float) :: FontSize;
  !     real(c_float) :: CurveTessellationTol;
  !     real(c_float) :: CircleSegmentMaxError;
  !     type(ImVec4) :: ClipRectFullscreen;
  !     ImDrawListFlags InitialFlags;
  !     type(ImVec2) :: ArcFastVtx[48];
  !     real(c_float) :: ArcFastRadiusCutoff;
  !     ImU8 CircleSegmentCounts[64];
  !     const ImVec4* TexUvLines;
  ! };
  ! typedef struct ImVector_ImDrawListPtr {int Size;int Capacity;ImDrawList** Data;} ImVector_ImDrawListPtr;

  ! struct ImDrawDataBuilder
  ! {
  !     ImVector_ImDrawListPtr Layers[2];
  ! };
  ! struct ImGuiDataTypeTempStorage
  ! {
  !     ImU8 Data[8];
  ! };
  ! struct ImGuiDataTypeInfo
  ! {
  !     size_t Size;
  !     const char* Name;
  !     const char* PrintFmt;
  !     const char* ScanFmt;
  ! };
  ! struct ImGuiColorMod
  ! {
  !     ImGuiCol Col;
  !     type(ImVec4) :: BackupValue;
  ! };
  ! struct ImGuiStyleMod
  ! {
  !     ImGuiStyleVar VarIdx;
  !     union { int BackupInt[2]; float BackupFloat[2]; };
  ! };
  ! struct ImGuiComboPreviewData
  ! {
  !     ImRect PreviewRect;
  !     type(ImVec2) :: BackupCursorPos;
  !     type(ImVec2) :: BackupCursorMaxPos;
  !     type(ImVec2) :: BackupCursorPosPrevLine;
  !     real(c_float) :: BackupPrevLineTextBaseOffset;
  !     ImGuiLayoutType BackupLayout;
  ! };
  ! struct ImGuiGroupData
  ! {
  !     ImGuiID WindowID;
  !     type(ImVec2) :: BackupCursorPos;
  !     type(ImVec2) :: BackupCursorMaxPos;
  !     ImVec1 BackupIndent;
  !     ImVec1 BackupGroupOffset;
  !     type(ImVec2) :: BackupCurrLineSize;
  !     real(c_float) :: BackupCurrLineTextBaseOffset;
  !     ImGuiID BackupActiveIdIsAlive;
  !     logical(c_bool) :: BackupActiveIdPreviousFrameIsAlive;
  !     logical(c_bool) :: BackupHoveredIdIsAlive;
  !     logical(c_bool) :: EmitItem;
  ! };
  ! struct ImGuiMenuColumns
  ! {
  !     ImU32 TotalWidth;
  !     ImU32 NextTotalWidth;
  !     ImU16 Spacing;
  !     ImU16 OffsetIcon;
  !     ImU16 OffsetLabel;
  !     ImU16 OffsetShortcut;
  !     ImU16 OffsetMark;
  !     ImU16 Widths[4];
  ! };
  ! struct ImGuiInputTextState
  ! {
  !     ImGuiID ID;
  !     integer(c_int) :: CurLenW, CurLenA;
  !     ImVector_ImWchar TextW;
  !     ImVector_char TextA;
  !     ImVector_char InitialTextA;
  !     logical(c_bool) :: TextAIsValid;
  !     integer(c_int) :: BufCapacityA;
  !     real(c_float) :: ScrollX;
  !     STB_TexteditState Stb;
  !     real(c_float) :: CursorAnim;
  !     logical(c_bool) :: CursorFollow;
  !     logical(c_bool) :: SelectedAllMouseLock;
  !     logical(c_bool) :: Edited;
  !     ImGuiInputTextFlags Flags;
  ! };
  ! struct ImGuiPopupData
  ! {
  !     ImGuiID PopupId;
  !     ImGuiWindow* Window;
  !     ImGuiWindow* SourceWindow;
  !     integer(c_int) :: ParentNavLayer;
  !     integer(c_int) :: OpenFrameCount;
  !     ImGuiID OpenParentId;
  !     type(ImVec2) :: OpenPopupPos;
  !     type(ImVec2) :: OpenMousePos;
  ! };
  ! struct ImGuiNextWindowData
  ! {
  !     ImGuiNextWindowDataFlags Flags;
  !     ImGuiCond PosCond;
  !     ImGuiCond SizeCond;
  !     ImGuiCond CollapsedCond;
  !     ImGuiCond DockCond;
  !     type(ImVec2) :: PosVal;
  !     type(ImVec2) :: PosPivotVal;
  !     type(ImVec2) :: SizeVal;
  !     type(ImVec2) :: ContentSizeVal;
  !     type(ImVec2) :: ScrollVal;
  !     logical(c_bool) :: PosUndock;
  !     logical(c_bool) :: CollapsedVal;
  !     ImRect SizeConstraintRect;
  !     ImGuiSizeCallback SizeCallback;
  !     type(c_ptr) :: SizeCallbackUserData;
  !     real(c_float) :: BgAlphaVal;
  !     ImGuiID ViewportId;
  !     ImGuiID DockId;
  !     ImGuiWindowClass WindowClass;
  !     type(ImVec2) :: MenuBarOffsetMinVal;
  ! };
  ! struct ImGuiNextItemData
  ! {
  !     ImGuiNextItemDataFlags Flags;
  !     real(c_float) :: Width;
  !     ImGuiID FocusScopeId;
  !     ImGuiCond OpenCond;
  !     logical(c_bool) :: OpenVal;
  ! };
  ! struct ImGuiLastItemData
  ! {
  !     ImGuiID ID;
  !     ImGuiItemFlags InFlags;
  !     ImGuiItemStatusFlags StatusFlags;
  !     ImRect Rect;
  !     ImRect NavRect;
  !     ImRect DisplayRect;
  ! };
  ! struct ImGuiStackSizes
  ! {
  !     integer(c_short) :: SizeOfIDStack;
  !     integer(c_short) :: SizeOfColorStack;
  !     integer(c_short) :: SizeOfStyleVarStack;
  !     integer(c_short) :: SizeOfFontStack;
  !     integer(c_short) :: SizeOfFocusScopeStack;
  !     integer(c_short) :: SizeOfGroupStack;
  !     integer(c_short) :: SizeOfItemFlagsStack;
  !     integer(c_short) :: SizeOfBeginPopupStack;
  !     integer(c_short) :: SizeOfDisabledStack;
  ! };
  ! struct ImGuiWindowStackData
  ! {
  !     ImGuiWindow* Window;
  !     ImGuiLastItemData ParentLastItemDataBackup;
  !     ImGuiStackSizes StackSizesOnBegin;
  ! };
  ! struct ImGuiShrinkWidthItem
  ! {
  !     integer(c_int) :: Index;
  !     real(c_float) :: Width;
  !     real(c_float) :: InitialWidth;
  ! };
  ! struct ImGuiPtrOrIndex
  ! {
  !     type(c_ptr) :: Ptr;
  !     integer(c_int) :: Index;
  ! };
  ! typedef struct ImBitArray_ImGuiKey_NamedKey_COUNT__lessImGuiKey_NamedKey_BEGIN {ImU32 Storage[(ImGuiKey_NamedKey_COUNT+31)>>5];} ImBitArray_ImGuiKey_NamedKey_COUNT__lessImGuiKey_NamedKey_BEGIN;
  ! struct ImGuiInputEventMousePos
  ! { real(c_float) :: PosX, PosY;
  ! };
  ! struct ImGuiInputEventMouseWheel
  ! { real(c_float) :: WheelX, WheelY;
  ! };
  ! struct ImGuiInputEventMouseButton
  ! { integer(c_int) :: Button; logical(c_bool) :: Down;
  ! };
  ! struct ImGuiInputEventMouseViewport
  ! { ImGuiID HoveredViewportID;
  ! };
  ! struct ImGuiInputEventKey
  ! { ImGuiKey Key; logical(c_bool) :: Down; real(c_float) :: AnalogValue;
  ! };
  ! struct ImGuiInputEventText
  ! { integer(c_int) :: Char;
  ! };
  ! struct ImGuiInputEventAppFocused
  ! { logical(c_bool) :: Focused;
  ! };
  ! struct ImGuiInputEvent
  ! {
  !     ImGuiInputEventType Type;
  !     ImGuiInputSource Source;
  !     union
  !     {
  !         ImGuiInputEventMousePos MousePos;
  !         ImGuiInputEventMouseWheel MouseWheel;
  !         ImGuiInputEventMouseButton MouseButton;
  !         ImGuiInputEventMouseViewport MouseViewport;
  !         ImGuiInputEventKey Key;
  !         ImGuiInputEventText Text;
  !         ImGuiInputEventAppFocused AppFocused;
  !     };
  !     logical(c_bool) :: AddedByTestEngine;
  ! };
  ! struct ImGuiListClipperRange
  ! {
  !     integer(c_int) :: Min;
  !     integer(c_int) :: Max;
  !     logical(c_bool) :: PosToIndexConvert;
  !     ImS8 PosToIndexOffsetMin;
  !     ImS8 PosToIndexOffsetMax;
  ! };
  ! typedef struct ImVector_ImGuiListClipperRange {int Size;int Capacity;ImGuiListClipperRange* Data;} ImVector_ImGuiListClipperRange;
  ! struct ImGuiListClipperData
  ! {
  !     ImGuiListClipper* ListClipper;
  !     real(c_float) :: LossynessOffset;
  !     integer(c_int) :: StepNo;
  !     integer(c_int) :: ItemsFrozen;
  !     ImVector_ImGuiListClipperRange Ranges;
  ! };
  ! struct ImGuiNavItemData
  ! {
  !     ImGuiWindow* Window;
  !     ImGuiID ID;
  !     ImGuiID FocusScopeId;
  !     ImRect RectRel;
  !     ImGuiItemFlags InFlags;
  !     real(c_float) :: DistBox;
  !     real(c_float) :: DistCenter;
  !     real(c_float) :: DistAxial;
  ! };
  ! struct ImGuiOldColumnData
  ! {
  !     real(c_float) :: OffsetNorm;
  !     real(c_float) :: OffsetNormBeforeResize;
  !     ImGuiOldColumnFlags Flags;
  !     ImRect ClipRect;
  ! };
  ! typedef struct ImVector_ImGuiOldColumnData {int Size;int Capacity;ImGuiOldColumnData* Data;} ImVector_ImGuiOldColumnData;
  ! struct ImGuiOldColumns
  ! {
  !     ImGuiID ID;
  !     ImGuiOldColumnFlags Flags;
  !     logical(c_bool) :: IsFirstFrame;
  !     logical(c_bool) :: IsBeingResized;
  !     integer(c_int) :: Current;
  !     integer(c_int) :: Count;
  !     real(c_float) :: OffMinX, OffMaxX;
  !     real(c_float) :: LineMinY, LineMaxY;
  !     real(c_float) :: HostCursorPosY;
  !     real(c_float) :: HostCursorMaxPosX;
  !     ImRect HostInitialClipRect;
  !     ImRect HostBackupClipRect;
  !     ImRect HostBackupParentWorkRect;
  !     ImVector_ImGuiOldColumnData Columns;
  !     ImDrawListSplitter Splitter;
  ! };
  ! typedef struct ImVector_ImGuiWindowPtr {int Size;int Capacity;ImGuiWindow** Data;} ImVector_ImGuiWindowPtr;
  ! struct ImGuiDockNode
  ! {
  !     ImGuiID ID;
  !     ImGuiDockNodeFlags SharedFlags;
  !     ImGuiDockNodeFlags LocalFlags;
  !     ImGuiDockNodeFlags LocalFlagsInWindows;
  !     ImGuiDockNodeFlags MergedFlags;
  !     ImGuiDockNodeState State;
  !     ImGuiDockNode* ParentNode;
  !     ImGuiDockNode* ChildNodes[2];
  !     ImVector_ImGuiWindowPtr Windows;
  !     ImGuiTabBar* TabBar;
  !     type(ImVec2) :: Pos;
  !     type(ImVec2) :: Size;
  !     type(ImVec2) :: SizeRef;
  !     ImGuiAxis SplitAxis;
  !     ImGuiWindowClass WindowClass;
  !     ImU32 LastBgColor;
  !     ImGuiWindow* HostWindow;
  !     ImGuiWindow* VisibleWindow;
  !     ImGuiDockNode* CentralNode;
  !     ImGuiDockNode* OnlyNodeWithWindows;
  !     integer(c_int) :: CountNodeWithWindows;
  !     integer(c_int) :: LastFrameAlive;
  !     integer(c_int) :: LastFrameActive;
  !     integer(c_int) :: LastFrameFocused;
  !     ImGuiID LastFocusedNodeId;
  !     ImGuiID SelectedTabId;
  !     ImGuiID WantCloseTabId;
  !     ImGuiDataAuthority AuthorityForPos :3;
  !     ImGuiDataAuthority AuthorityForSize :3;
  !     ImGuiDataAuthority AuthorityForViewport :3;
  !     logical(c_bool) :: IsVisible :1;
  !     logical(c_bool) :: IsFocused :1;
  !     logical(c_bool) :: IsBgDrawnThisFrame :1;
  !     logical(c_bool) :: HasCloseButton :1;
  !     logical(c_bool) :: HasWindowMenuButton :1;
  !     logical(c_bool) :: HasCentralNodeChild :1;
  !     logical(c_bool) :: WantCloseAll :1;
  !     logical(c_bool) :: WantLockSizeOnce :1;
  !     logical(c_bool) :: WantMouseMove :1;
  !     logical(c_bool) :: WantHiddenTabBarUpdate :1;
  !     logical(c_bool) :: WantHiddenTabBarToggle :1;
  ! };
  ! struct ImGuiWindowDockStyle
  ! {
  !     ImU32 Colors[ImGuiWindowDockStyleCol_COUNT];
  ! };
  ! typedef struct ImVector_ImGuiDockRequest {int Size;int Capacity;ImGuiDockRequest* Data;} ImVector_ImGuiDockRequest;
  ! typedef struct ImVector_ImGuiDockNodeSettings {int Size;int Capacity;ImGuiDockNodeSettings* Data;} ImVector_ImGuiDockNodeSettings;
  ! struct ImGuiDockContext
  ! {
  !     ImGuiStorage Nodes;
  !     ImVector_ImGuiDockRequest Requests;
  !     ImVector_ImGuiDockNodeSettings NodesSettings;
  !     logical(c_bool) :: WantFullRebuild;
  ! };
  ! struct ImGuiViewportP
  ! {
  !     ImGuiViewport _ImGuiViewport;
  !     integer(c_int) :: Idx;
  !     integer(c_int) :: LastFrameActive;
  !     integer(c_int) :: LastFrontMostStampCount;
  !     ImGuiID LastNameHash;
  !     type(ImVec2) :: LastPos;
  !     real(c_float) :: Alpha;
  !     real(c_float) :: LastAlpha;
  !     integer(c_short) :: PlatformMonitor;
  !     logical(c_bool) :: PlatformWindowCreated;
  !     ImGuiWindow* Window;
  !     integer(c_int) :: DrawListsLastFrame[2];
  !     ImDrawList* DrawLists[2];
  !     ImDrawData DrawDataP;
  !     ImDrawDataBuilder DrawDataBuilder;
  !     type(ImVec2) :: LastPlatformPos;
  !     type(ImVec2) :: LastPlatformSize;
  !     type(ImVec2) :: LastRendererSize;
  !     type(ImVec2) :: WorkOffsetMin;
  !     type(ImVec2) :: WorkOffsetMax;
  !     type(ImVec2) :: BuildWorkOffsetMin;
  !     type(ImVec2) :: BuildWorkOffsetMax;
  ! };
  ! struct ImGuiWindowSettings
  ! {
  !     ImGuiID ID;
  !     ImVec2ih Pos;
  !     ImVec2ih Size;
  !     ImVec2ih ViewportPos;
  !     ImGuiID ViewportId;
  !     ImGuiID DockId;
  !     ImGuiID ClassId;
  !     integer(c_short) :: DockOrder;
  !     logical(c_bool) :: Collapsed;
  !     logical(c_bool) :: WantApply;
  ! };
  ! struct ImGuiSettingsHandler
  ! {
  !     const char* TypeName;
  !     ImGuiID TypeHash;
  !     void (*ClearAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler);
  !     void (*ReadInitFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler);
  !     type(c_ptr) :: (*ReadOpenFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, const char* name);
  !     void (*ReadLineFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, void* entry, const char* line);
  !     void (*ApplyAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler);
  !     void (*WriteAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, ImGuiTextBuffer* out_buf);
  !     type(c_ptr) :: UserData;
  ! };
  ! struct ImGuiMetricsConfig
  ! {
  !     logical(c_bool) :: ShowDebugLog;
  !     logical(c_bool) :: ShowStackTool;
  !     logical(c_bool) :: ShowWindowsRects;
  !     logical(c_bool) :: ShowWindowsBeginOrder;
  !     logical(c_bool) :: ShowTablesRects;
  !     logical(c_bool) :: ShowDrawCmdMesh;
  !     logical(c_bool) :: ShowDrawCmdBoundingBoxes;
  !     logical(c_bool) :: ShowDockingNodes;
  !     integer(c_int) :: ShowWindowsRectsType;
  !     integer(c_int) :: ShowTablesRectsType;
  ! };
  ! struct ImGuiStackLevelInfo
  ! {
  !     ImGuiID ID;
  !     ImS8 QueryFrameCount;
  !     logical(c_bool) :: QuerySuccess;
  !     ImGuiDataType DataType : 8;
  !     character(c_char) :: Desc[57];
  ! };
  ! typedef struct ImVector_ImGuiStackLevelInfo {int Size;int Capacity;ImGuiStackLevelInfo* Data;} ImVector_ImGuiStackLevelInfo;
  ! struct ImGuiStackTool
  ! {
  !     integer(c_int) :: LastActiveFrame;
  !     integer(c_int) :: StackLevel;
  !     ImGuiID QueryId;
  !     ImVector_ImGuiStackLevelInfo Results;
  !     logical(c_bool) :: CopyToClipboardOnCtrlC;
  !     real(c_float) :: CopyToClipboardLastTime;
  ! };
  ! typedef void (*ImGuiContextHookCallback)(ImGuiContext* ctx, ImGuiContextHook* hook);
  ! struct ImGuiContextHook
  ! {
  !     ImGuiID HookId;
  !     ImGuiContextHookType Type;
  !     ImGuiID Owner;
  !     ImGuiContextHookCallback Callback;
  !     type(c_ptr) :: UserData;
  ! };

  type, bind(c) ::ImVector_ImGuiInputEvent
     integer(c_int) :: Size
     integer(c_int) :: Capacity
     type(c_ptr) :: Data
  end type ImVector_ImGuiInputEvent

  ! typedef struct ImVector_ImGuiWindowStackData {int Size;int Capacity;ImGuiWindowStackData* Data;} ImVector_ImGuiWindowStackData;

  ! typedef struct ImVector_ImGuiColorMod {int Size;int Capacity;ImGuiColorMod* Data;} ImVector_ImGuiColorMod;

  ! typedef struct ImVector_ImGuiStyleMod {int Size;int Capacity;ImGuiStyleMod* Data;} ImVector_ImGuiStyleMod;

  ! typedef struct ImVector_ImGuiID {int Size;int Capacity;ImGuiID* Data;} ImVector_ImGuiID;

  ! typedef struct ImVector_ImGuiItemFlags {int Size;int Capacity;ImGuiItemFlags* Data;} ImVector_ImGuiItemFlags;

  ! typedef struct ImVector_ImGuiGroupData {int Size;int Capacity;ImGuiGroupData* Data;} ImVector_ImGuiGroupData;

  ! typedef struct ImVector_ImGuiPopupData {int Size;int Capacity;ImGuiPopupData* Data;} ImVector_ImGuiPopupData;

  ! typedef struct ImVector_ImGuiViewportPPtr {int Size;int Capacity;ImGuiViewportP** Data;} ImVector_ImGuiViewportPPtr;

  ! typedef struct ImVector_unsigned_char {int Size;int Capacity;unsigned char* Data;} ImVector_unsigned_char;

  ! typedef struct ImVector_ImGuiListClipperData {int Size;int Capacity;ImGuiListClipperData* Data;} ImVector_ImGuiListClipperData;

  ! typedef struct ImVector_ImGuiTableTempData {int Size;int Capacity;ImGuiTableTempData* Data;} ImVector_ImGuiTableTempData;

  ! typedef struct ImVector_ImGuiTable {int Size;int Capacity;ImGuiTable* Data;} ImVector_ImGuiTable;

  ! typedef struct ImPool_ImGuiTable {ImVector_ImGuiTable Buf;ImGuiStorage Map;ImPoolIdx FreeIdx;ImPoolIdx AliveCount;} ImPool_ImGuiTable;

  ! typedef struct ImVector_ImGuiTabBar {int Size;int Capacity;ImGuiTabBar* Data;} ImVector_ImGuiTabBar;

  ! typedef struct ImPool_ImGuiTabBar {ImVector_ImGuiTabBar Buf;ImGuiStorage Map;ImPoolIdx FreeIdx;ImPoolIdx AliveCount;} ImPool_ImGuiTabBar;

  ! typedef struct ImVector_ImGuiPtrOrIndex {int Size;int Capacity;ImGuiPtrOrIndex* Data;} ImVector_ImGuiPtrOrIndex;

  ! typedef struct ImVector_ImGuiShrinkWidthItem {int Size;int Capacity;ImGuiShrinkWidthItem* Data;} ImVector_ImGuiShrinkWidthItem;

  ! typedef struct ImVector_ImGuiSettingsHandler {int Size;int Capacity;ImGuiSettingsHandler* Data;} ImVector_ImGuiSettingsHandler;

  ! typedef struct ImChunkStream_ImGuiWindowSettings {ImVector_char Buf;} ImChunkStream_ImGuiWindowSettings;

  ! typedef struct ImChunkStream_ImGuiTableSettings {ImVector_char Buf;} ImChunkStream_ImGuiTableSettings;

  ! typedef struct ImVector_ImGuiContextHook {int Size;int Capacity;ImGuiContextHook* Data;} ImVector_ImGuiContextHook;

  ! struct ImGuiContext
  ! {
  !     logical(c_bool) :: Initialized;
  !     logical(c_bool) :: FontAtlasOwnedByContext;
  !     type(ImGuiIO) :: IO;
  !     type(ImGuiPlatformIO) :: PlatformIO;
  !     type(ImVector_ImGuiInputEvent) :: InputEventsQueue;
  !     type(ImVector_ImGuiInputEvent) :: InputEventsTrail;
  !     ImGuiStyle Style;
  !     ImGuiConfigFlags ConfigFlagsCurrFrame;
  !     ImGuiConfigFlags ConfigFlagsLastFrame;
  !     ImFont* Font;
  !     real(c_float) :: FontSize;
  !     real(c_float) :: FontBaseSize;
  !     ImDrawListSharedData DrawListSharedData;
  !     real(c_double) :: Time;
  !     integer(c_int) :: FrameCount;
  !     integer(c_int) :: FrameCountEnded;
  !     integer(c_int) :: FrameCountPlatformEnded;
  !     integer(c_int) :: FrameCountRendered;
  !     logical(c_bool) :: WithinFrameScope;
  !     logical(c_bool) :: WithinFrameScopeWithImplicitWindow;
  !     logical(c_bool) :: WithinEndChild;
  !     logical(c_bool) :: GcCompactAll;
  !     logical(c_bool) :: TestEngineHookItems;
  !     type(c_ptr) :: TestEngine;
  !     ImVector_ImGuiWindowPtr Windows;
  !     ImVector_ImGuiWindowPtr WindowsFocusOrder;
  !     ImVector_ImGuiWindowPtr WindowsTempSortBuffer;
  !     ImVector_ImGuiWindowStackData CurrentWindowStack;
  !     ImGuiStorage WindowsById;
  !     integer(c_int) :: WindowsActiveCount;
  !     type(ImVec2) :: WindowsHoverPadding;
  !     ImGuiWindow* CurrentWindow;
  !     ImGuiWindow* HoveredWindow;
  !     ImGuiWindow* HoveredWindowUnderMovingWindow;
  !     ImGuiDockNode* HoveredDockNode;
  !     ImGuiWindow* MovingWindow;
  !     ImGuiWindow* WheelingWindow;
  !     type(ImVec2) :: WheelingWindowRefMousePos;
  !     real(c_float) :: WheelingWindowTimer;
  !     ImGuiID DebugHookIdInfo;
  !     ImGuiID HoveredId;
  !     ImGuiID HoveredIdPreviousFrame;
  !     logical(c_bool) :: HoveredIdAllowOverlap;
  !     logical(c_bool) :: HoveredIdUsingMouseWheel;
  !     logical(c_bool) :: HoveredIdPreviousFrameUsingMouseWheel;
  !     logical(c_bool) :: HoveredIdDisabled;
  !     real(c_float) :: HoveredIdTimer;
  !     real(c_float) :: HoveredIdNotActiveTimer;
  !     ImGuiID ActiveId;
  !     ImGuiID ActiveIdIsAlive;
  !     real(c_float) :: ActiveIdTimer;
  !     logical(c_bool) :: ActiveIdIsJustActivated;
  !     logical(c_bool) :: ActiveIdAllowOverlap;
  !     logical(c_bool) :: ActiveIdNoClearOnFocusLoss;
  !     logical(c_bool) :: ActiveIdHasBeenPressedBefore;
  !     logical(c_bool) :: ActiveIdHasBeenEditedBefore;
  !     logical(c_bool) :: ActiveIdHasBeenEditedThisFrame;
  !     type(ImVec2) :: ActiveIdClickOffset;
  !     ImGuiWindow* ActiveIdWindow;
  !     ImGuiInputSource ActiveIdSource;
  !     integer(c_int) :: ActiveIdMouseButton;
  !     ImGuiID ActiveIdPreviousFrame;
  !     logical(c_bool) :: ActiveIdPreviousFrameIsAlive;
  !     logical(c_bool) :: ActiveIdPreviousFrameHasBeenEditedBefore;
  !     ImGuiWindow* ActiveIdPreviousFrameWindow;
  !     ImGuiID LastActiveId;
  !     real(c_float) :: LastActiveIdTimer;
  !     logical(c_bool) :: ActiveIdUsingMouseWheel;
  !     ImU32 ActiveIdUsingNavDirMask;
  !     ImU32 ActiveIdUsingNavInputMask;
  !     ImBitArrayForNamedKeys ActiveIdUsingKeyInputMask;
  !     ImGuiItemFlags CurrentItemFlags;
  !     ImGuiNextItemData NextItemData;
  !     ImGuiLastItemData LastItemData;
  !     ImGuiNextWindowData NextWindowData;
  !     ImVector_ImGuiColorMod ColorStack;
  !     ImVector_ImGuiStyleMod StyleVarStack;
  !     ImVector_ImFontPtr FontStack;
  !     ImVector_ImGuiID FocusScopeStack;
  !     ImVector_ImGuiItemFlags ItemFlagsStack;
  !     ImVector_ImGuiGroupData GroupStack;
  !     ImVector_ImGuiPopupData OpenPopupStack;
  !     ImVector_ImGuiPopupData BeginPopupStack;
  !     integer(c_int) :: BeginMenuCount;
  !     ImVector_ImGuiViewportPPtr Viewports;
  !     real(c_float) :: CurrentDpiScale;
  !     ImGuiViewportP* CurrentViewport;
  !     ImGuiViewportP* MouseViewport;
  !     ImGuiViewportP* MouseLastHoveredViewport;
  !     ImGuiID PlatformLastFocusedViewportId;
  !     ImGuiPlatformMonitor FallbackMonitor;
  !     integer(c_int) :: ViewportFrontMostStampCount;
  !     ImGuiWindow* NavWindow;
  !     ImGuiID NavId;
  !     ImGuiID NavFocusScopeId;
  !     ImGuiID NavActivateId;
  !     ImGuiID NavActivateDownId;
  !     ImGuiID NavActivatePressedId;
  !     ImGuiID NavActivateInputId;
  !     ImGuiActivateFlags NavActivateFlags;
  !     ImGuiID NavJustMovedToId;
  !     ImGuiID NavJustMovedToFocusScopeId;
  !     ImGuiModFlags NavJustMovedToKeyMods;
  !     ImGuiID NavNextActivateId;
  !     ImGuiActivateFlags NavNextActivateFlags;
  !     ImGuiInputSource NavInputSource;
  !     ImGuiNavLayer NavLayer;
  !     logical(c_bool) :: NavIdIsAlive;
  !     logical(c_bool) :: NavMousePosDirty;
  !     logical(c_bool) :: NavDisableHighlight;
  !     logical(c_bool) :: NavDisableMouseHover;
  !     logical(c_bool) :: NavAnyRequest;
  !     logical(c_bool) :: NavInitRequest;
  !     logical(c_bool) :: NavInitRequestFromMove;
  !     ImGuiID NavInitResultId;
  !     ImRect NavInitResultRectRel;
  !     logical(c_bool) :: NavMoveSubmitted;
  !     logical(c_bool) :: NavMoveScoringItems;
  !     logical(c_bool) :: NavMoveForwardToNextFrame;
  !     ImGuiNavMoveFlags NavMoveFlags;
  !     ImGuiScrollFlags NavMoveScrollFlags;
  !     ImGuiModFlags NavMoveKeyMods;
  !     integer(c_int) :: NavMoveDir;
  !     integer(c_int) :: NavMoveDirForDebug;
  !     integer(c_int) :: NavMoveClipDir;
  !     ImRect NavScoringRect;
  !     ImRect NavScoringNoClipRect;
  !     integer(c_int) :: NavScoringDebugCount;
  !     integer(c_int) :: NavTabbingDir;
  !     integer(c_int) :: NavTabbingCounter;
  !     ImGuiNavItemData NavMoveResultLocal;
  !     ImGuiNavItemData NavMoveResultLocalVisible;
  !     ImGuiNavItemData NavMoveResultOther;
  !     ImGuiNavItemData NavTabbingResultFirst;
  !     ImGuiWindow* NavWindowingTarget;
  !     ImGuiWindow* NavWindowingTargetAnim;
  !     ImGuiWindow* NavWindowingListWindow;
  !     real(c_float) :: NavWindowingTimer;
  !     real(c_float) :: NavWindowingHighlightAlpha;
  !     logical(c_bool) :: NavWindowingToggleLayer;
  !     real(c_float) :: DimBgRatio;
  !     ImGuiMouseCursor MouseCursor;
  !     logical(c_bool) :: DragDropActive;
  !     logical(c_bool) :: DragDropWithinSource;
  !     logical(c_bool) :: DragDropWithinTarget;
  !     ImGuiDragDropFlags DragDropSourceFlags;
  !     integer(c_int) :: DragDropSourceFrameCount;
  !     integer(c_int) :: DragDropMouseButton;
  !     ImGuiPayload DragDropPayload;
  !     ImRect DragDropTargetRect;
  !     ImGuiID DragDropTargetId;
  !     ImGuiDragDropFlags DragDropAcceptFlags;
  !     real(c_float) :: DragDropAcceptIdCurrRectSurface;
  !     ImGuiID DragDropAcceptIdCurr;
  !     ImGuiID DragDropAcceptIdPrev;
  !     integer(c_int) :: DragDropAcceptFrameCount;
  !     ImGuiID DragDropHoldJustPressedId;
  !     ImVector_unsigned_char DragDropPayloadBufHeap;
  !     integer(c_signed_char) :: DragDropPayloadBufLocal[16];
  !     integer(c_int) :: ClipperTempDataStacked;
  !     ImVector_ImGuiListClipperData ClipperTempData;
  !     ImGuiTable* CurrentTable;
  !     integer(c_int) :: TablesTempDataStacked;
  !     ImVector_ImGuiTableTempData TablesTempData;
  !     ImPool_ImGuiTable Tables;
  !     ImVector_float TablesLastTimeActive;
  !     ImVector_ImDrawChannel DrawChannelsTempMergeBuffer;
  !     ImGuiTabBar* CurrentTabBar;
  !     ImPool_ImGuiTabBar TabBars;
  !     ImVector_ImGuiPtrOrIndex CurrentTabBarStack;
  !     ImVector_ImGuiShrinkWidthItem ShrinkWidthBuffer;
  !     type(ImVec2) :: MouseLastValidPos;
  !     ImGuiInputTextState InputTextState;
  !     ImFont InputTextPasswordFont;
  !     ImGuiID TempInputId;
  !     ImGuiColorEditFlags ColorEditOptions;
  !     real(c_float) :: ColorEditLastHue;
  !     real(c_float) :: ColorEditLastSat;
  !     ImU32 ColorEditLastColor;
  !     type(ImVec4) :: ColorPickerRef;
  !     ImGuiComboPreviewData ComboPreviewData;
  !     real(c_float) :: SliderGrabClickOffset;
  !     real(c_float) :: SliderCurrentAccum;
  !     logical(c_bool) :: SliderCurrentAccumDirty;
  !     logical(c_bool) :: DragCurrentAccumDirty;
  !     real(c_float) :: DragCurrentAccum;
  !     real(c_float) :: DragSpeedDefaultRatio;
  !     real(c_float) :: ScrollbarClickDeltaToGrabCenter;
  !     real(c_float) :: DisabledAlphaBackup;
  !     integer(c_short) :: DisabledStackSize;
  !     integer(c_short) :: TooltipOverrideCount;
  !     real(c_float) :: TooltipSlowDelay;
  !     ImVector_char ClipboardHandlerData;
  !     ImVector_ImGuiID MenusIdSubmittedThisFrame;
  !     ImGuiPlatformImeData PlatformImeData;
  !     ImGuiPlatformImeData PlatformImeDataPrev;
  !     ImGuiID PlatformImeViewport;
  !     character(c_char) :: PlatformLocaleDecimalPoint;
  !     ImGuiDockContext DockContext;
  !     logical(c_bool) :: SettingsLoaded;
  !     real(c_float) :: SettingsDirtyTimer;
  !     ImGuiTextBuffer SettingsIniData;
  !     ImVector_ImGuiSettingsHandler SettingsHandlers;
  !     ImChunkStream_ImGuiWindowSettings SettingsWindows;
  !     ImChunkStream_ImGuiTableSettings SettingsTables;
  !     ImVector_ImGuiContextHook Hooks;
  !     ImGuiID HookIdNext;
  !     logical(c_bool) :: LogEnabled;
  !     ImGuiLogType LogType;
  !     ImFileHandle LogFile;
  !     ImGuiTextBuffer LogBuffer;
  !     const char* LogNextPrefix;
  !     const char* LogNextSuffix;
  !     real(c_float) :: LogLinePosY;
  !     logical(c_bool) :: LogLineFirstItem;
  !     integer(c_int) :: LogDepthRef;
  !     integer(c_int) :: LogDepthToExpand;
  !     integer(c_int) :: LogDepthToExpandDefault;
  !     ImGuiDebugLogFlags DebugLogFlags;
  !     ImGuiTextBuffer DebugLogBuf;
  !     logical(c_bool) :: DebugItemPickerActive;
  !     ImGuiID DebugItemPickerBreakId;
  !     ImGuiMetricsConfig DebugMetricsConfig;
  !     ImGuiStackTool DebugStackTool;
  !     real(c_float) :: FramerateSecPerFrame[120];
  !     integer(c_int) :: FramerateSecPerFrameIdx;
  !     integer(c_int) :: FramerateSecPerFrameCount;
  !     real(c_float) :: FramerateSecPerFrameAccum;
  !     integer(c_int) :: WantCaptureMouseNextFrame;
  !     integer(c_int) :: WantCaptureKeyboardNextFrame;
  !     integer(c_int) :: WantTextInputNextFrame;
  !     ImVector_char TempBuffer;
  ! };
  ! struct ImGuiWindowTempData
  ! {
  !     type(ImVec2) :: CursorPos;
  !     type(ImVec2) :: CursorPosPrevLine;
  !     type(ImVec2) :: CursorStartPos;
  !     type(ImVec2) :: CursorMaxPos;
  !     type(ImVec2) :: IdealMaxPos;
  !     type(ImVec2) :: CurrLineSize;
  !     type(ImVec2) :: PrevLineSize;
  !     real(c_float) :: CurrLineTextBaseOffset;
  !     real(c_float) :: PrevLineTextBaseOffset;
  !     logical(c_bool) :: IsSameLine;
  !     ImVec1 Indent;
  !     ImVec1 ColumnsOffset;
  !     ImVec1 GroupOffset;
  !     type(ImVec2) :: CursorStartPosLossyness;
  !     ImGuiNavLayer NavLayerCurrent;
  !     integer(c_short) :: NavLayersActiveMask;
  !     integer(c_short) :: NavLayersActiveMaskNext;
  !     ImGuiID NavFocusScopeIdCurrent;
  !     logical(c_bool) :: NavHideHighlightOneFrame;
  !     logical(c_bool) :: NavHasScroll;
  !     logical(c_bool) :: MenuBarAppending;
  !     type(ImVec2) :: MenuBarOffset;
  !     ImGuiMenuColumns MenuColumns;
  !     integer(c_int) :: TreeDepth;
  !     ImU32 TreeJumpToParentOnPopMask;
  !     ImVector_ImGuiWindowPtr ChildWindows;
  !     ImGuiStorage* StateStorage;
  !     ImGuiOldColumns* CurrentColumns;
  !     integer(c_int) :: CurrentTableIdx;
  !     ImGuiLayoutType LayoutType;
  !     ImGuiLayoutType ParentLayoutType;
  !     real(c_float) :: ItemWidth;
  !     real(c_float) :: TextWrapPos;
  !     ImVector_float ItemWidthStack;
  !     ImVector_float TextWrapPosStack;
  ! };
  ! typedef struct ImVector_ImGuiOldColumns {int Size;int Capacity;ImGuiOldColumns* Data;} ImVector_ImGuiOldColumns;

  ! struct ImGuiWindow
  ! {
  !     char* Name;
  !     ImGuiID ID;
  !     ImGuiWindowFlags Flags, FlagsPreviousFrame;
  !     ImGuiWindowClass WindowClass;
  !     ImGuiViewportP* Viewport;
  !     ImGuiID ViewportId;
  !     type(ImVec2) :: ViewportPos;
  !     integer(c_int) :: ViewportAllowPlatformMonitorExtend;
  !     type(ImVec2) :: Pos;
  !     type(ImVec2) :: Size;
  !     type(ImVec2) :: SizeFull;
  !     type(ImVec2) :: ContentSize;
  !     type(ImVec2) :: ContentSizeIdeal;
  !     type(ImVec2) :: ContentSizeExplicit;
  !     type(ImVec2) :: WindowPadding;
  !     real(c_float) :: WindowRounding;
  !     real(c_float) :: WindowBorderSize;
  !     integer(c_int) :: NameBufLen;
  !     ImGuiID MoveId;
  !     ImGuiID TabId;
  !     ImGuiID ChildId;
  !     type(ImVec2) :: Scroll;
  !     type(ImVec2) :: ScrollMax;
  !     type(ImVec2) :: ScrollTarget;
  !     type(ImVec2) :: ScrollTargetCenterRatio;
  !     type(ImVec2) :: ScrollTargetEdgeSnapDist;
  !     type(ImVec2) :: ScrollbarSizes;
  !     logical(c_bool) :: ScrollbarX, ScrollbarY;
  !     logical(c_bool) :: ViewportOwned;
  !     logical(c_bool) :: Active;
  !     logical(c_bool) :: WasActive;
  !     logical(c_bool) :: WriteAccessed;
  !     logical(c_bool) :: Collapsed;
  !     logical(c_bool) :: WantCollapseToggle;
  !     logical(c_bool) :: SkipItems;
  !     logical(c_bool) :: Appearing;
  !     logical(c_bool) :: Hidden;
  !     logical(c_bool) :: IsFallbackWindow;
  !     logical(c_bool) :: IsExplicitChild;
  !     logical(c_bool) :: HasCloseButton;
  !     integer(c_signed_char) :: ResizeBorderHeld;
  !     integer(c_short) :: BeginCount;
  !     integer(c_short) :: BeginOrderWithinParent;
  !     integer(c_short) :: BeginOrderWithinContext;
  !     integer(c_short) :: FocusOrder;
  !     ImGuiID PopupId;
  !     ImS8 AutoFitFramesX, AutoFitFramesY;
  !     ImS8 AutoFitChildAxises;
  !     logical(c_bool) :: AutoFitOnlyGrows;
  !     integer(c_int) :: AutoPosLastDirection;
  !     ImS8 HiddenFramesCanSkipItems;
  !     ImS8 HiddenFramesCannotSkipItems;
  !     ImS8 HiddenFramesForRenderOnly;
  !     ImS8 DisableInputsFrames;
  !     ImGuiCond SetWindowPosAllowFlags : 8;
  !     ImGuiCond SetWindowSizeAllowFlags : 8;
  !     ImGuiCond SetWindowCollapsedAllowFlags : 8;
  !     ImGuiCond SetWindowDockAllowFlags : 8;
  !     type(ImVec2) :: SetWindowPosVal;
  !     type(ImVec2) :: SetWindowPosPivot;
  !     ImVector_ImGuiID IDStack;
  !     ImGuiWindowTempData DC;
  !     ImRect OuterRectClipped;
  !     ImRect InnerRect;
  !     ImRect InnerClipRect;
  !     ImRect WorkRect;
  !     ImRect ParentWorkRect;
  !     ImRect ClipRect;
  !     ImRect ContentRegionRect;
  !     ImVec2ih HitTestHoleSize;
  !     ImVec2ih HitTestHoleOffset;
  !     integer(c_int) :: LastFrameActive;
  !     integer(c_int) :: LastFrameJustFocused;
  !     real(c_float) :: LastTimeActive;
  !     real(c_float) :: ItemWidthDefault;
  !     ImGuiStorage StateStorage;
  !     ImVector_ImGuiOldColumns ColumnsStorage;
  !     real(c_float) :: FontWindowScale;
  !     real(c_float) :: FontDpiScale;
  !     integer(c_int) :: SettingsOffset;
  !     ImDrawList* DrawList;
  !     ImDrawList DrawListInst;
  !     ImGuiWindow* ParentWindow;
  !     ImGuiWindow* ParentWindowInBeginStack;
  !     ImGuiWindow* RootWindow;
  !     ImGuiWindow* RootWindowPopupTree;
  !     ImGuiWindow* RootWindowDockTree;
  !     ImGuiWindow* RootWindowForTitleBarHighlight;
  !     ImGuiWindow* RootWindowForNav;
  !     ImGuiWindow* NavLastChildNavWindow;
  !     ImGuiID NavLastIds[ImGuiNavLayer_COUNT];
  !     ImRect NavRectRel[ImGuiNavLayer_COUNT];
  !     integer(c_int) :: MemoryDrawListIdxCapacity;
  !     integer(c_int) :: MemoryDrawListVtxCapacity;
  !     logical(c_bool) :: MemoryCompacted;
  !     logical(c_bool) :: DockIsActive :1;
  !     logical(c_bool) :: DockNodeIsVisible :1;
  !     logical(c_bool) :: DockTabIsVisible :1;
  !     logical(c_bool) :: DockTabWantClose :1;
  !     integer(c_short) :: DockOrder;
  !     ImGuiWindowDockStyle DockStyle;
  !     ImGuiDockNode* DockNode;
  !     ImGuiDockNode* DockNodeAsHost;
  !     ImGuiID DockId;
  !     ImGuiItemStatusFlags DockTabItemStatusFlags;
  !     ImRect DockTabItemRect;
  ! };
  ! struct ImGuiTabItem
  ! {
  !     ImGuiID ID;
  !     ImGuiTabItemFlags Flags;
  !     ImGuiWindow* Window;
  !     integer(c_int) :: LastFrameVisible;
  !     integer(c_int) :: LastFrameSelected;
  !     real(c_float) :: Offset;
  !     real(c_float) :: Width;
  !     real(c_float) :: ContentWidth;
  !     real(c_float) :: RequestedWidth;
  !     ImS32 NameOffset;
  !     ImS16 BeginOrder;
  !     ImS16 IndexDuringLayout;
  !     logical(c_bool) :: WantClose;
  ! };
  ! typedef struct ImVector_ImGuiTabItem {int Size;int Capacity;ImGuiTabItem* Data;} ImVector_ImGuiTabItem;

  ! struct ImGuiTabBar
  ! {
  !     ImVector_ImGuiTabItem Tabs;
  !     ImGuiTabBarFlags Flags;
  !     ImGuiID ID;
  !     ImGuiID SelectedTabId;
  !     ImGuiID NextSelectedTabId;
  !     ImGuiID VisibleTabId;
  !     integer(c_int) :: CurrFrameVisible;
  !     integer(c_int) :: PrevFrameVisible;
  !     ImRect BarRect;
  !     real(c_float) :: CurrTabsContentsHeight;
  !     real(c_float) :: PrevTabsContentsHeight;
  !     real(c_float) :: WidthAllTabs;
  !     real(c_float) :: WidthAllTabsIdeal;
  !     real(c_float) :: ScrollingAnim;
  !     real(c_float) :: ScrollingTarget;
  !     real(c_float) :: ScrollingTargetDistToVisibility;
  !     real(c_float) :: ScrollingSpeed;
  !     real(c_float) :: ScrollingRectMinX;
  !     real(c_float) :: ScrollingRectMaxX;
  !     ImGuiID ReorderRequestTabId;
  !     ImS16 ReorderRequestOffset;
  !     ImS8 BeginCount;
  !     logical(c_bool) :: WantLayout;
  !     logical(c_bool) :: VisibleTabWasSubmitted;
  !     logical(c_bool) :: TabsAddedNew;
  !     ImS16 TabsActiveCount;
  !     ImS16 LastTabItemIdx;
  !     real(c_float) :: ItemSpacingY;
  !     type(ImVec2) :: FramePadding;
  !     type(ImVec2) :: BackupCursorPos;
  !     ImGuiTextBuffer TabsNames;
  ! };
  ! struct ImGuiTableColumn
  ! {
  !     ImGuiTableColumnFlags Flags;
  !     real(c_float) :: WidthGiven;
  !     real(c_float) :: MinX;
  !     real(c_float) :: MaxX;
  !     real(c_float) :: WidthRequest;
  !     real(c_float) :: WidthAuto;
  !     real(c_float) :: StretchWeight;
  !     real(c_float) :: InitStretchWeightOrWidth;
  !     ImRect ClipRect;
  !     ImGuiID UserID;
  !     real(c_float) :: WorkMinX;
  !     real(c_float) :: WorkMaxX;
  !     real(c_float) :: ItemWidth;
  !     real(c_float) :: ContentMaxXFrozen;
  !     real(c_float) :: ContentMaxXUnfrozen;
  !     real(c_float) :: ContentMaxXHeadersUsed;
  !     real(c_float) :: ContentMaxXHeadersIdeal;
  !     ImS16 NameOffset;
  !     ImGuiTableColumnIdx DisplayOrder;
  !     ImGuiTableColumnIdx IndexWithinEnabledSet;
  !     ImGuiTableColumnIdx PrevEnabledColumn;
  !     ImGuiTableColumnIdx NextEnabledColumn;
  !     ImGuiTableColumnIdx SortOrder;
  !     ImGuiTableDrawChannelIdx DrawChannelCurrent;
  !     ImGuiTableDrawChannelIdx DrawChannelFrozen;
  !     ImGuiTableDrawChannelIdx DrawChannelUnfrozen;
  !     logical(c_bool) :: IsEnabled;
  !     logical(c_bool) :: IsUserEnabled;
  !     logical(c_bool) :: IsUserEnabledNextFrame;
  !     logical(c_bool) :: IsVisibleX;
  !     logical(c_bool) :: IsVisibleY;
  !     logical(c_bool) :: IsRequestOutput;
  !     logical(c_bool) :: IsSkipItems;
  !     logical(c_bool) :: IsPreserveWidthAuto;
  !     ImS8 NavLayerCurrent;
  !     ImU8 AutoFitQueue;
  !     ImU8 CannotSkipItemsQueue;
  !     ImU8 SortDirection : 2;
  !     ImU8 SortDirectionsAvailCount : 2;
  !     ImU8 SortDirectionsAvailMask : 4;
  !     ImU8 SortDirectionsAvailList;
  ! };
  ! struct ImGuiTableCellData
  ! {
  !     ImU32 BgColor;
  !     ImGuiTableColumnIdx Column;
  ! };
  ! struct ImGuiTableInstanceData
  ! {
  !     real(c_float) :: LastOuterHeight;
  !     real(c_float) :: LastFirstRowHeight;
  ! };
  ! typedef struct ImSpan_ImGuiTableColumn {ImGuiTableColumn* Data;ImGuiTableColumn* DataEnd;} ImSpan_ImGuiTableColumn;

  ! typedef struct ImSpan_ImGuiTableColumnIdx {ImGuiTableColumnIdx* Data;ImGuiTableColumnIdx* DataEnd;} ImSpan_ImGuiTableColumnIdx;

  ! typedef struct ImSpan_ImGuiTableCellData {ImGuiTableCellData* Data;ImGuiTableCellData* DataEnd;} ImSpan_ImGuiTableCellData;

  ! typedef struct ImVector_ImGuiTableInstanceData {int Size;int Capacity;ImGuiTableInstanceData* Data;} ImVector_ImGuiTableInstanceData;

  ! typedef struct ImVector_ImGuiTableColumnSortSpecs {int Size;int Capacity;ImGuiTableColumnSortSpecs* Data;} ImVector_ImGuiTableColumnSortSpecs;

  ! struct ImGuiTable
  ! {
  !     ImGuiID ID;
  !     ImGuiTableFlags Flags;
  !     type(c_ptr) :: RawData;
  !     ImGuiTableTempData* TempData;
  !     ImSpan_ImGuiTableColumn Columns;
  !     ImSpan_ImGuiTableColumnIdx DisplayOrderToIndex;
  !     ImSpan_ImGuiTableCellData RowCellData;
  !     ImU64 EnabledMaskByDisplayOrder;
  !     ImU64 EnabledMaskByIndex;
  !     ImU64 VisibleMaskByIndex;
  !     ImU64 RequestOutputMaskByIndex;
  !     ImGuiTableFlags SettingsLoadedFlags;
  !     integer(c_int) :: SettingsOffset;
  !     integer(c_int) :: LastFrameActive;
  !     integer(c_int) :: ColumnsCount;
  !     integer(c_int) :: CurrentRow;
  !     integer(c_int) :: CurrentColumn;
  !     ImS16 InstanceCurrent;
  !     ImS16 InstanceInteracted;
  !     real(c_float) :: RowPosY1;
  !     real(c_float) :: RowPosY2;
  !     real(c_float) :: RowMinHeight;
  !     real(c_float) :: RowTextBaseline;
  !     real(c_float) :: RowIndentOffsetX;
  !     ImGuiTableRowFlags RowFlags : 16;
  !     ImGuiTableRowFlags LastRowFlags : 16;
  !     integer(c_int) :: RowBgColorCounter;
  !     ImU32 RowBgColor[2];
  !     ImU32 BorderColorStrong;
  !     ImU32 BorderColorLight;
  !     real(c_float) :: BorderX1;
  !     real(c_float) :: BorderX2;
  !     real(c_float) :: HostIndentX;
  !     real(c_float) :: MinColumnWidth;
  !     real(c_float) :: OuterPaddingX;
  !     real(c_float) :: CellPaddingX;
  !     real(c_float) :: CellPaddingY;
  !     real(c_float) :: CellSpacingX1;
  !     real(c_float) :: CellSpacingX2;
  !     real(c_float) :: InnerWidth;
  !     real(c_float) :: ColumnsGivenWidth;
  !     real(c_float) :: ColumnsAutoFitWidth;
  !     real(c_float) :: ColumnsStretchSumWeights;
  !     real(c_float) :: ResizedColumnNextWidth;
  !     real(c_float) :: ResizeLockMinContentsX2;
  !     real(c_float) :: RefScale;
  !     ImRect OuterRect;
  !     ImRect InnerRect;
  !     ImRect WorkRect;
  !     ImRect InnerClipRect;
  !     ImRect BgClipRect;
  !     ImRect Bg0ClipRectForDrawCmd;
  !     ImRect Bg2ClipRectForDrawCmd;
  !     ImRect HostClipRect;
  !     ImRect HostBackupInnerClipRect;
  !     ImGuiWindow* OuterWindow;
  !     ImGuiWindow* InnerWindow;
  !     ImGuiTextBuffer ColumnsNames;
  !     ImDrawListSplitter* DrawSplitter;
  !     ImGuiTableInstanceData InstanceDataFirst;
  !     ImVector_ImGuiTableInstanceData InstanceDataExtra;
  !     ImGuiTableColumnSortSpecs SortSpecsSingle;
  !     ImVector_ImGuiTableColumnSortSpecs SortSpecsMulti;
  !     ImGuiTableSortSpecs SortSpecs;
  !     ImGuiTableColumnIdx SortSpecsCount;
  !     ImGuiTableColumnIdx ColumnsEnabledCount;
  !     ImGuiTableColumnIdx ColumnsEnabledFixedCount;
  !     ImGuiTableColumnIdx DeclColumnsCount;
  !     ImGuiTableColumnIdx HoveredColumnBody;
  !     ImGuiTableColumnIdx HoveredColumnBorder;
  !     ImGuiTableColumnIdx AutoFitSingleColumn;
  !     ImGuiTableColumnIdx ResizedColumn;
  !     ImGuiTableColumnIdx LastResizedColumn;
  !     ImGuiTableColumnIdx HeldHeaderColumn;
  !     ImGuiTableColumnIdx ReorderColumn;
  !     ImGuiTableColumnIdx ReorderColumnDir;
  !     ImGuiTableColumnIdx LeftMostEnabledColumn;
  !     ImGuiTableColumnIdx RightMostEnabledColumn;
  !     ImGuiTableColumnIdx LeftMostStretchedColumn;
  !     ImGuiTableColumnIdx RightMostStretchedColumn;
  !     ImGuiTableColumnIdx ContextPopupColumn;
  !     ImGuiTableColumnIdx FreezeRowsRequest;
  !     ImGuiTableColumnIdx FreezeRowsCount;
  !     ImGuiTableColumnIdx FreezeColumnsRequest;
  !     ImGuiTableColumnIdx FreezeColumnsCount;
  !     ImGuiTableColumnIdx RowCellDataCurrent;
  !     ImGuiTableDrawChannelIdx DummyDrawChannel;
  !     ImGuiTableDrawChannelIdx Bg2DrawChannelCurrent;
  !     ImGuiTableDrawChannelIdx Bg2DrawChannelUnfrozen;
  !     logical(c_bool) :: IsLayoutLocked;
  !     logical(c_bool) :: IsInsideRow;
  !     logical(c_bool) :: IsInitializing;
  !     logical(c_bool) :: IsSortSpecsDirty;
  !     logical(c_bool) :: IsUsingHeaders;
  !     logical(c_bool) :: IsContextPopupOpen;
  !     logical(c_bool) :: IsSettingsRequestLoad;
  !     logical(c_bool) :: IsSettingsDirty;
  !     logical(c_bool) :: IsDefaultDisplayOrder;
  !     logical(c_bool) :: IsResetAllRequest;
  !     logical(c_bool) :: IsResetDisplayOrderRequest;
  !     logical(c_bool) :: IsUnfrozenRows;
  !     logical(c_bool) :: IsDefaultSizingPolicy;
  !     logical(c_bool) :: MemoryCompacted;
  !     logical(c_bool) :: HostSkipItems;
  ! };
  ! struct ImGuiTableTempData
  ! {
  !     integer(c_int) :: TableIndex;
  !     real(c_float) :: LastTimeActive;
  !     type(ImVec2) :: UserOuterSize;
  !     ImDrawListSplitter DrawSplitter;
  !     ImRect HostBackupWorkRect;
  !     ImRect HostBackupParentWorkRect;
  !     type(ImVec2) :: HostBackupPrevLineSize;
  !     type(ImVec2) :: HostBackupCurrLineSize;
  !     type(ImVec2) :: HostBackupCursorMaxPos;
  !     ImVec1 HostBackupColumnsOffset;
  !     real(c_float) :: HostBackupItemWidth;
  !     integer(c_int) :: HostBackupItemWidthStackSize;
  ! };
  ! struct ImGuiTableColumnSettings
  ! {
  !     real(c_float) :: WidthOrWeight;
  !     ImGuiID UserID;
  !     ImGuiTableColumnIdx Index;
  !     ImGuiTableColumnIdx DisplayOrder;
  !     ImGuiTableColumnIdx SortOrder;
  !     ImU8 SortDirection : 2;
  !     ImU8 IsEnabled : 1;
  !     ImU8 IsStretch : 1;
  ! };
  ! struct ImGuiTableSettings
  ! {
  !     ImGuiID ID;
  !     ImGuiTableFlags SaveFlags;
  !     real(c_float) :: RefScale;
  !     ImGuiTableColumnIdx ColumnsCount;
  !     ImGuiTableColumnIdx ColumnsCountMax;
  !     logical(c_bool) :: WantApply;
  ! };
  ! struct ImFontBuilderIO
  ! {
  !     bool (*FontBuilder_Build)(ImFontAtlas* atlas);
  ! };
