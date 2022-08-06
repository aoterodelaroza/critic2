  ! -*-F90-*-
  
  !xx! cimgui user-defined types
  ! typedef struct ImVector_const_charPtr {int Size;int Capacity;const char** Data;} ImVector_const_charPtr;
  type, bind(c) :: ImVec2
     real(c_float) :: x
     real(c_float) :: y
  end type ImVec2

  ! struct ImVec4
  ! {
  !     float x, y, z, w;
  ! };
  ! struct ImGuiStyle
  ! {
  !     float Alpha;
  !     float DisabledAlpha;
  !     ImVec2 WindowPadding;
  !     float WindowRounding;
  !     float WindowBorderSize;
  !     ImVec2 WindowMinSize;
  !     ImVec2 WindowTitleAlign;
  !     ImGuiDir WindowMenuButtonPosition;
  !     float ChildRounding;
  !     float ChildBorderSize;
  !     float PopupRounding;
  !     float PopupBorderSize;
  !     ImVec2 FramePadding;
  !     float FrameRounding;
  !     float FrameBorderSize;
  !     ImVec2 ItemSpacing;
  !     ImVec2 ItemInnerSpacing;
  !     ImVec2 CellPadding;
  !     ImVec2 TouchExtraPadding;
  !     float IndentSpacing;
  !     float ColumnsMinSpacing;
  !     float ScrollbarSize;
  !     float ScrollbarRounding;
  !     float GrabMinSize;
  !     float GrabRounding;
  !     float LogSliderDeadzone;
  !     float TabRounding;
  !     float TabBorderSize;
  !     float TabMinWidthForCloseButton;
  !     ImGuiDir ColorButtonPosition;
  !     ImVec2 ButtonTextAlign;
  !     ImVec2 SelectableTextAlign;
  !     ImVec2 DisplayWindowPadding;
  !     ImVec2 DisplaySafeAreaPadding;
  !     float MouseCursorScale;
  !     bool AntiAliasedLines;
  !     bool AntiAliasedLinesUseTex;
  !     bool AntiAliasedFill;
  !     float CurveTessellationTol;
  !     float CircleTessellationMaxError;
  !     ImVec4 Colors[ImGuiCol_COUNT];
  ! };

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
  !     void* UserData;
  !     ImWchar EventChar;
  !     ImGuiKey EventKey;
  !     char* Buf;
  !     int BufTextLen;
  !     int BufSize;
  !     bool BufDirty;
  !     int CursorPos;
  !     int SelectionStart;
  !     int SelectionEnd;
  ! };
  ! struct ImGuiSizeCallbackData
  ! {
  !     void* UserData;
  !     ImVec2 Pos;
  !     ImVec2 CurrentSize;
  !     ImVec2 DesiredSize;
  ! };
  ! struct ImGuiWindowClass
  ! {
  !     ImGuiID ClassId;
  !     ImGuiID ParentViewportId;
  !     ImGuiViewportFlags ViewportFlagsOverrideSet;
  !     ImGuiViewportFlags ViewportFlagsOverrideClear;
  !     ImGuiTabItemFlags TabItemFlagsOverrideSet;
  !     ImGuiDockNodeFlags DockNodeFlagsOverrideSet;
  !     bool DockingAlwaysTabBar;
  !     bool DockingAllowUnclassed;
  ! };
  ! struct ImGuiPayload
  ! {
  !     void* Data;
  !     int DataSize;
  !     ImGuiID SourceId;
  !     ImGuiID SourceParentId;
  !     int DataFrameCount;
  !     char DataType[32 + 1];
  !     bool Preview;
  !     bool Delivery;
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
  !     int SpecsCount;
  !     bool SpecsDirty;
  ! };
  ! struct ImGuiOnceUponAFrame
  ! {
  !      int RefFrame;
  ! };
  ! struct ImGuiTextRange
  ! {
  !         const char* b;
  !         const char* e;
  ! };
  ! typedef struct ImVector_ImGuiTextRange {int Size;int Capacity;ImGuiTextRange* Data;} ImVector_ImGuiTextRange;
  ! struct ImGuiTextFilter
  ! {
  !     char InputBuf[256];
  !     ImVector_ImGuiTextRange Filters;
  !     int CountGrep;
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
  !     int DisplayStart;
  !     int DisplayEnd;
  !     int ItemsCount;
  !     float ItemsHeight;
  !     float StartPosY;
  !     void* TempData;
  ! };
  ! struct ImColor
  ! {
  !     ImVec4 Value;
  ! };
  ! typedef void (*ImDrawCallback)(const ImDrawList* parent_list, const ImDrawCmd* cmd);
  ! struct ImDrawCmd
  ! {
  !     ImVec4 ClipRect;
  !     ImTextureID TextureId;
  !     unsigned int VtxOffset;
  !     unsigned int IdxOffset;
  !     unsigned int ElemCount;
  !     ImDrawCallback UserCallback;
  !     void* UserCallbackData;
  ! };
  ! struct ImDrawVert
  ! {
  !     ImVec2 pos;
  !     ImVec2 uv;
  !     ImU32 col;
  ! };
  ! struct ImDrawCmdHeader
  ! {
  !     ImVec4 ClipRect;
  !     ImTextureID TextureId;
  !     unsigned int VtxOffset;
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
  !     int _Current;
  !     int _Count;
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
  !     unsigned int _VtxCurrentIdx;
  !     const ImDrawListSharedData* _Data;
  !     const char* _OwnerName;
  !     ImDrawVert* _VtxWritePtr;
  !     ImDrawIdx* _IdxWritePtr;
  !     ImVector_ImVec4 _ClipRectStack;
  !     ImVector_ImTextureID _TextureIdStack;
  !     ImVector_ImVec2 _Path;
  !     ImDrawCmdHeader _CmdHeader;
  !     ImDrawListSplitter _Splitter;
  !     float _FringeScale;
  ! };
  ! struct ImDrawData
  ! {
  !     bool Valid;
  !     int CmdListsCount;
  !     int TotalIdxCount;
  !     int TotalVtxCount;
  !     ImDrawList** CmdLists;
  !     ImVec2 DisplayPos;
  !     ImVec2 DisplaySize;
  !     ImVec2 FramebufferScale;
  !     ImGuiViewport* OwnerViewport;
  ! };
  ! struct ImFontConfig
  ! {
  !     void* FontData;
  !     int FontDataSize;
  !     bool FontDataOwnedByAtlas;
  !     int FontNo;
  !     float SizePixels;
  !     int OversampleH;
  !     int OversampleV;
  !     bool PixelSnapH;
  !     ImVec2 GlyphExtraSpacing;
  !     ImVec2 GlyphOffset;
  !     const ImWchar* GlyphRanges;
  !     float GlyphMinAdvanceX;
  !     float GlyphMaxAdvanceX;
  !     bool MergeMode;
  !     unsigned int FontBuilderFlags;
  !     float RasterizerMultiply;
  !     ImWchar EllipsisChar;
  !     char Name[40];
  !     ImFont* DstFont;
  ! };
  ! struct ImFontGlyph
  ! {
  !     unsigned int Colored : 1;
  !     unsigned int Visible : 1;
  !     unsigned int Codepoint : 30;
  !     float AdvanceX;
  !     float X0, Y0, X1, Y1;
  !     float U0, V0, U1, V1;
  ! };
  ! typedef struct ImVector_ImU32 {int Size;int Capacity;ImU32* Data;} ImVector_ImU32;

  ! struct ImFontGlyphRangesBuilder
  ! {
  !     ImVector_ImU32 UsedChars;
  ! };
  ! typedef struct ImFontAtlasCustomRect ImFontAtlasCustomRect;
  ! struct ImFontAtlasCustomRect
  ! {
  !     unsigned short Width, Height;
  !     unsigned short X, Y;
  !     unsigned int GlyphID;
  !     float GlyphAdvanceX;
  !     ImVec2 GlyphOffset;
  !     ImFont* Font;
  ! };
  ! typedef struct ImVector_ImFontPtr {int Size;int Capacity;ImFont** Data;} ImVector_ImFontPtr;

  ! typedef struct ImVector_ImFontAtlasCustomRect {int Size;int Capacity;ImFontAtlasCustomRect* Data;} ImVector_ImFontAtlasCustomRect;

  ! typedef struct ImVector_ImFontConfig {int Size;int Capacity;ImFontConfig* Data;} ImVector_ImFontConfig;

  ! struct ImFontAtlas
  ! {
  !     ImFontAtlasFlags Flags;
  !     ImTextureID TexID;
  !     int TexDesiredWidth;
  !     int TexGlyphPadding;
  !     bool Locked;
  !     bool TexReady;
  !     bool TexPixelsUseColors;
  !     unsigned char* TexPixelsAlpha8;
  !     unsigned int* TexPixelsRGBA32;
  !     int TexWidth;
  !     int TexHeight;
  !     ImVec2 TexUvScale;
  !     ImVec2 TexUvWhitePixel;
  !     ImVector_ImFontPtr Fonts;
  !     ImVector_ImFontAtlasCustomRect CustomRects;
  !     ImVector_ImFontConfig ConfigData;
  !     ImVec4 TexUvLines[(63) + 1];
  !     const ImFontBuilderIO* FontBuilderIO;
  !     unsigned int FontBuilderFlags;
  !     int PackIdMouseCursors;
  !     int PackIdLines;
  ! };
  ! typedef struct ImVector_float {int Size;int Capacity;float* Data;} ImVector_float;

  ! typedef struct ImVector_ImFontGlyph {int Size;int Capacity;ImFontGlyph* Data;} ImVector_ImFontGlyph;

  ! struct ImFont
  ! {
  !     ImVector_float IndexAdvanceX;
  !     float FallbackAdvanceX;
  !     float FontSize;
  !     ImVector_ImWchar IndexLookup;
  !     ImVector_ImFontGlyph Glyphs;
  !     const ImFontGlyph* FallbackGlyph;
  !     ImFontAtlas* ContainerAtlas;
  !     const ImFontConfig* ConfigData;
  !     short ConfigDataCount;
  !     ImWchar FallbackChar;
  !     ImWchar EllipsisChar;
  !     ImWchar DotChar;
  !     bool DirtyLookupTables;
  !     float Scale;
  !     float Ascent, Descent;
  !     int MetricsTotalSurface;
  !     ImU8 Used4kPagesMap[(0xFFFF +1)/4096/8];
  ! };
  ! struct ImGuiViewport
  ! {
  !     ImGuiID ID;
  !     ImGuiViewportFlags Flags;
  !     ImVec2 Pos;
  !     ImVec2 Size;
  !     ImVec2 WorkPos;
  !     ImVec2 WorkSize;
  !     float DpiScale;
  !     ImGuiID ParentViewportId;
  !     ImDrawData* DrawData;
  !     void* RendererUserData;
  !     void* PlatformUserData;
  !     void* PlatformHandle;
  !     void* PlatformHandleRaw;
  !     bool PlatformRequestMove;
  !     bool PlatformRequestResize;
  !     bool PlatformRequestClose;
  ! };
  ! typedef struct ImVector_ImGuiPlatformMonitor {int Size;int Capacity;ImGuiPlatformMonitor* Data;} ImVector_ImGuiPlatformMonitor;

  ! typedef struct ImVector_ImGuiViewportPtr {int Size;int Capacity;ImGuiViewport** Data;} ImVector_ImGuiViewportPtr;

  ! struct ImGuiPlatformIO
  ! {
  !     void (*Platform_CreateWindow)(ImGuiViewport* vp);
  !     void (*Platform_DestroyWindow)(ImGuiViewport* vp);
  !     void (*Platform_ShowWindow)(ImGuiViewport* vp);
  !     void (*Platform_SetWindowPos)(ImGuiViewport* vp, ImVec2 pos);
  !     ImVec2 (*Platform_GetWindowPos)(ImGuiViewport* vp);
  !     void (*Platform_SetWindowSize)(ImGuiViewport* vp, ImVec2 size);
  !     ImVec2 (*Platform_GetWindowSize)(ImGuiViewport* vp);
  !     void (*Platform_SetWindowFocus)(ImGuiViewport* vp);
  !     bool (*Platform_GetWindowFocus)(ImGuiViewport* vp);
  !     bool (*Platform_GetWindowMinimized)(ImGuiViewport* vp);
  !     void (*Platform_SetWindowTitle)(ImGuiViewport* vp, const char* str);
  !     void (*Platform_SetWindowAlpha)(ImGuiViewport* vp, float alpha);
  !     void (*Platform_UpdateWindow)(ImGuiViewport* vp);
  !     void (*Platform_RenderWindow)(ImGuiViewport* vp, void* render_arg);
  !     void (*Platform_SwapBuffers)(ImGuiViewport* vp, void* render_arg);
  !     float (*Platform_GetWindowDpiScale)(ImGuiViewport* vp);
  !     void (*Platform_OnChangedViewport)(ImGuiViewport* vp);
  !     int (*Platform_CreateVkSurface)(ImGuiViewport* vp, ImU64 vk_inst, const void* vk_allocators, ImU64* out_vk_surface);
  !     void (*Renderer_CreateWindow)(ImGuiViewport* vp);
  !     void (*Renderer_DestroyWindow)(ImGuiViewport* vp);
  !     void (*Renderer_SetWindowSize)(ImGuiViewport* vp, ImVec2 size);
  !     void (*Renderer_RenderWindow)(ImGuiViewport* vp, void* render_arg);
  !     void (*Renderer_SwapBuffers)(ImGuiViewport* vp, void* render_arg);
  !     ImVector_ImGuiPlatformMonitor Monitors;
  !     ImVector_ImGuiViewportPtr Viewports;
  ! };
  ! struct ImGuiPlatformMonitor
  ! {
  !     ImVec2 MainPos, MainSize;
  !     ImVec2 WorkPos, WorkSize;
  !     float DpiScale;
  ! };
  ! struct ImGuiPlatformImeData
  ! {
  !     bool WantVisible;
  !     ImVec2 InputPos;
  !     float InputLineHeight;
  ! };
  ! typedef void (*ImGuiErrorLogCallback)(void* user_data, const char* fmt, ...);
  ! struct StbUndoRecord
  ! {
  !    int where;
  !    int insert_length;
  !    int delete_length;
  !    int char_storage;
  ! };
  ! struct StbUndoState
  ! {
  !    StbUndoRecord undo_rec [99];
  !    ImWchar undo_char[999];
  !    short undo_point, redo_point;
  !    int undo_char_point, redo_char_point;
  ! };
  ! struct STB_TexteditState
  ! {
  !    int cursor;
  !    int select_start;
  !    int select_end;
  !    unsigned char insert_mode;
  !    int row_count_per_page;
  !    unsigned char cursor_at_end_of_line;
  !    unsigned char initialized;
  !    unsigned char has_preferred_x;
  !    unsigned char single_line;
  !    unsigned char padding1, padding2, padding3;
  !    float preferred_x;
  !    StbUndoState undostate;
  ! };
  ! struct StbTexteditRow
  ! {
  !    float x0,x1;
  !    float baseline_y_delta;
  !    float ymin,ymax;
  !    int num_chars;
  ! };
  ! struct ImVec1
  ! {
  !     float x;
  ! };
  ! struct ImVec2ih
  ! {
  !     short x, y;
  ! };
  ! struct ImRect
  ! {
  !     ImVec2 Min;
  !     ImVec2 Max;
  ! };
  ! struct ImBitVector
  ! {
  !     ImVector_ImU32 Storage;
  ! };
  ! struct ImDrawListSharedData
  ! {
  !     ImVec2 TexUvWhitePixel;
  !     ImFont* Font;
  !     float FontSize;
  !     float CurveTessellationTol;
  !     float CircleSegmentMaxError;
  !     ImVec4 ClipRectFullscreen;
  !     ImDrawListFlags InitialFlags;
  !     ImVec2 ArcFastVtx[48];
  !     float ArcFastRadiusCutoff;
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
  !     ImVec4 BackupValue;
  ! };
  ! struct ImGuiStyleMod
  ! {
  !     ImGuiStyleVar VarIdx;
  !     union { int BackupInt[2]; float BackupFloat[2]; };
  ! };
  ! struct ImGuiComboPreviewData
  ! {
  !     ImRect PreviewRect;
  !     ImVec2 BackupCursorPos;
  !     ImVec2 BackupCursorMaxPos;
  !     ImVec2 BackupCursorPosPrevLine;
  !     float BackupPrevLineTextBaseOffset;
  !     ImGuiLayoutType BackupLayout;
  ! };
  ! struct ImGuiGroupData
  ! {
  !     ImGuiID WindowID;
  !     ImVec2 BackupCursorPos;
  !     ImVec2 BackupCursorMaxPos;
  !     ImVec1 BackupIndent;
  !     ImVec1 BackupGroupOffset;
  !     ImVec2 BackupCurrLineSize;
  !     float BackupCurrLineTextBaseOffset;
  !     ImGuiID BackupActiveIdIsAlive;
  !     bool BackupActiveIdPreviousFrameIsAlive;
  !     bool BackupHoveredIdIsAlive;
  !     bool EmitItem;
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
  !     int CurLenW, CurLenA;
  !     ImVector_ImWchar TextW;
  !     ImVector_char TextA;
  !     ImVector_char InitialTextA;
  !     bool TextAIsValid;
  !     int BufCapacityA;
  !     float ScrollX;
  !     STB_TexteditState Stb;
  !     float CursorAnim;
  !     bool CursorFollow;
  !     bool SelectedAllMouseLock;
  !     bool Edited;
  !     ImGuiInputTextFlags Flags;
  ! };
  ! struct ImGuiPopupData
  ! {
  !     ImGuiID PopupId;
  !     ImGuiWindow* Window;
  !     ImGuiWindow* SourceWindow;
  !     int ParentNavLayer;
  !     int OpenFrameCount;
  !     ImGuiID OpenParentId;
  !     ImVec2 OpenPopupPos;
  !     ImVec2 OpenMousePos;
  ! };
  ! struct ImGuiNextWindowData
  ! {
  !     ImGuiNextWindowDataFlags Flags;
  !     ImGuiCond PosCond;
  !     ImGuiCond SizeCond;
  !     ImGuiCond CollapsedCond;
  !     ImGuiCond DockCond;
  !     ImVec2 PosVal;
  !     ImVec2 PosPivotVal;
  !     ImVec2 SizeVal;
  !     ImVec2 ContentSizeVal;
  !     ImVec2 ScrollVal;
  !     bool PosUndock;
  !     bool CollapsedVal;
  !     ImRect SizeConstraintRect;
  !     ImGuiSizeCallback SizeCallback;
  !     void* SizeCallbackUserData;
  !     float BgAlphaVal;
  !     ImGuiID ViewportId;
  !     ImGuiID DockId;
  !     ImGuiWindowClass WindowClass;
  !     ImVec2 MenuBarOffsetMinVal;
  ! };
  ! struct ImGuiNextItemData
  ! {
  !     ImGuiNextItemDataFlags Flags;
  !     float Width;
  !     ImGuiID FocusScopeId;
  !     ImGuiCond OpenCond;
  !     bool OpenVal;
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
  !     short SizeOfIDStack;
  !     short SizeOfColorStack;
  !     short SizeOfStyleVarStack;
  !     short SizeOfFontStack;
  !     short SizeOfFocusScopeStack;
  !     short SizeOfGroupStack;
  !     short SizeOfItemFlagsStack;
  !     short SizeOfBeginPopupStack;
  !     short SizeOfDisabledStack;
  ! };
  ! struct ImGuiWindowStackData
  ! {
  !     ImGuiWindow* Window;
  !     ImGuiLastItemData ParentLastItemDataBackup;
  !     ImGuiStackSizes StackSizesOnBegin;
  ! };
  ! struct ImGuiShrinkWidthItem
  ! {
  !     int Index;
  !     float Width;
  !     float InitialWidth;
  ! };
  ! struct ImGuiPtrOrIndex
  ! {
  !     void* Ptr;
  !     int Index;
  ! };
  ! typedef struct ImBitArray_ImGuiKey_NamedKey_COUNT__lessImGuiKey_NamedKey_BEGIN {ImU32 Storage[(ImGuiKey_NamedKey_COUNT+31)>>5];} ImBitArray_ImGuiKey_NamedKey_COUNT__lessImGuiKey_NamedKey_BEGIN;
  ! struct ImGuiInputEventMousePos
  ! { float PosX, PosY;
  ! };
  ! struct ImGuiInputEventMouseWheel
  ! { float WheelX, WheelY;
  ! };
  ! struct ImGuiInputEventMouseButton
  ! { int Button; bool Down;
  ! };
  ! struct ImGuiInputEventMouseViewport
  ! { ImGuiID HoveredViewportID;
  ! };
  ! struct ImGuiInputEventKey
  ! { ImGuiKey Key; bool Down; float AnalogValue;
  ! };
  ! struct ImGuiInputEventText
  ! { unsigned int Char;
  ! };
  ! struct ImGuiInputEventAppFocused
  ! { bool Focused;
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
  !     bool AddedByTestEngine;
  ! };
  ! struct ImGuiListClipperRange
  ! {
  !     int Min;
  !     int Max;
  !     bool PosToIndexConvert;
  !     ImS8 PosToIndexOffsetMin;
  !     ImS8 PosToIndexOffsetMax;
  ! };
  ! typedef struct ImVector_ImGuiListClipperRange {int Size;int Capacity;ImGuiListClipperRange* Data;} ImVector_ImGuiListClipperRange;
  ! struct ImGuiListClipperData
  ! {
  !     ImGuiListClipper* ListClipper;
  !     float LossynessOffset;
  !     int StepNo;
  !     int ItemsFrozen;
  !     ImVector_ImGuiListClipperRange Ranges;
  ! };
  ! struct ImGuiNavItemData
  ! {
  !     ImGuiWindow* Window;
  !     ImGuiID ID;
  !     ImGuiID FocusScopeId;
  !     ImRect RectRel;
  !     ImGuiItemFlags InFlags;
  !     float DistBox;
  !     float DistCenter;
  !     float DistAxial;
  ! };
  ! struct ImGuiOldColumnData
  ! {
  !     float OffsetNorm;
  !     float OffsetNormBeforeResize;
  !     ImGuiOldColumnFlags Flags;
  !     ImRect ClipRect;
  ! };
  ! typedef struct ImVector_ImGuiOldColumnData {int Size;int Capacity;ImGuiOldColumnData* Data;} ImVector_ImGuiOldColumnData;
  ! struct ImGuiOldColumns
  ! {
  !     ImGuiID ID;
  !     ImGuiOldColumnFlags Flags;
  !     bool IsFirstFrame;
  !     bool IsBeingResized;
  !     int Current;
  !     int Count;
  !     float OffMinX, OffMaxX;
  !     float LineMinY, LineMaxY;
  !     float HostCursorPosY;
  !     float HostCursorMaxPosX;
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
  !     ImVec2 Pos;
  !     ImVec2 Size;
  !     ImVec2 SizeRef;
  !     ImGuiAxis SplitAxis;
  !     ImGuiWindowClass WindowClass;
  !     ImU32 LastBgColor;
  !     ImGuiWindow* HostWindow;
  !     ImGuiWindow* VisibleWindow;
  !     ImGuiDockNode* CentralNode;
  !     ImGuiDockNode* OnlyNodeWithWindows;
  !     int CountNodeWithWindows;
  !     int LastFrameAlive;
  !     int LastFrameActive;
  !     int LastFrameFocused;
  !     ImGuiID LastFocusedNodeId;
  !     ImGuiID SelectedTabId;
  !     ImGuiID WantCloseTabId;
  !     ImGuiDataAuthority AuthorityForPos :3;
  !     ImGuiDataAuthority AuthorityForSize :3;
  !     ImGuiDataAuthority AuthorityForViewport :3;
  !     bool IsVisible :1;
  !     bool IsFocused :1;
  !     bool IsBgDrawnThisFrame :1;
  !     bool HasCloseButton :1;
  !     bool HasWindowMenuButton :1;
  !     bool HasCentralNodeChild :1;
  !     bool WantCloseAll :1;
  !     bool WantLockSizeOnce :1;
  !     bool WantMouseMove :1;
  !     bool WantHiddenTabBarUpdate :1;
  !     bool WantHiddenTabBarToggle :1;
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
  !     bool WantFullRebuild;
  ! };
  ! struct ImGuiViewportP
  ! {
  !     ImGuiViewport _ImGuiViewport;
  !     int Idx;
  !     int LastFrameActive;
  !     int LastFrontMostStampCount;
  !     ImGuiID LastNameHash;
  !     ImVec2 LastPos;
  !     float Alpha;
  !     float LastAlpha;
  !     short PlatformMonitor;
  !     bool PlatformWindowCreated;
  !     ImGuiWindow* Window;
  !     int DrawListsLastFrame[2];
  !     ImDrawList* DrawLists[2];
  !     ImDrawData DrawDataP;
  !     ImDrawDataBuilder DrawDataBuilder;
  !     ImVec2 LastPlatformPos;
  !     ImVec2 LastPlatformSize;
  !     ImVec2 LastRendererSize;
  !     ImVec2 WorkOffsetMin;
  !     ImVec2 WorkOffsetMax;
  !     ImVec2 BuildWorkOffsetMin;
  !     ImVec2 BuildWorkOffsetMax;
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
  !     short DockOrder;
  !     bool Collapsed;
  !     bool WantApply;
  ! };
  ! struct ImGuiSettingsHandler
  ! {
  !     const char* TypeName;
  !     ImGuiID TypeHash;
  !     void (*ClearAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler);
  !     void (*ReadInitFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler);
  !     void* (*ReadOpenFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, const char* name);
  !     void (*ReadLineFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, void* entry, const char* line);
  !     void (*ApplyAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler);
  !     void (*WriteAllFn)(ImGuiContext* ctx, ImGuiSettingsHandler* handler, ImGuiTextBuffer* out_buf);
  !     void* UserData;
  ! };
  ! struct ImGuiMetricsConfig
  ! {
  !     bool ShowDebugLog;
  !     bool ShowStackTool;
  !     bool ShowWindowsRects;
  !     bool ShowWindowsBeginOrder;
  !     bool ShowTablesRects;
  !     bool ShowDrawCmdMesh;
  !     bool ShowDrawCmdBoundingBoxes;
  !     bool ShowDockingNodes;
  !     int ShowWindowsRectsType;
  !     int ShowTablesRectsType;
  ! };
  ! struct ImGuiStackLevelInfo
  ! {
  !     ImGuiID ID;
  !     ImS8 QueryFrameCount;
  !     bool QuerySuccess;
  !     ImGuiDataType DataType : 8;
  !     char Desc[57];
  ! };
  ! typedef struct ImVector_ImGuiStackLevelInfo {int Size;int Capacity;ImGuiStackLevelInfo* Data;} ImVector_ImGuiStackLevelInfo;
  ! struct ImGuiStackTool
  ! {
  !     int LastActiveFrame;
  !     int StackLevel;
  !     ImGuiID QueryId;
  !     ImVector_ImGuiStackLevelInfo Results;
  !     bool CopyToClipboardOnCtrlC;
  !     float CopyToClipboardLastTime;
  ! };
  ! typedef void (*ImGuiContextHookCallback)(ImGuiContext* ctx, ImGuiContextHook* hook);
  ! struct ImGuiContextHook
  ! {
  !     ImGuiID HookId;
  !     ImGuiContextHookType Type;
  !     ImGuiID Owner;
  !     ImGuiContextHookCallback Callback;
  !     void* UserData;
  ! };
  ! typedef struct ImVector_ImGuiInputEvent {int Size;int Capacity;ImGuiInputEvent* Data;} ImVector_ImGuiInputEvent;

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
  !     bool Initialized;
  !     bool FontAtlasOwnedByContext;
  !     ImGuiIO IO;
  !     ImGuiPlatformIO PlatformIO;
  !     ImVector_ImGuiInputEvent InputEventsQueue;
  !     ImVector_ImGuiInputEvent InputEventsTrail;
  !     ImGuiStyle Style;
  !     ImGuiConfigFlags ConfigFlagsCurrFrame;
  !     ImGuiConfigFlags ConfigFlagsLastFrame;
  !     ImFont* Font;
  !     float FontSize;
  !     float FontBaseSize;
  !     ImDrawListSharedData DrawListSharedData;
  !     double Time;
  !     int FrameCount;
  !     int FrameCountEnded;
  !     int FrameCountPlatformEnded;
  !     int FrameCountRendered;
  !     bool WithinFrameScope;
  !     bool WithinFrameScopeWithImplicitWindow;
  !     bool WithinEndChild;
  !     bool GcCompactAll;
  !     bool TestEngineHookItems;
  !     void* TestEngine;
  !     ImVector_ImGuiWindowPtr Windows;
  !     ImVector_ImGuiWindowPtr WindowsFocusOrder;
  !     ImVector_ImGuiWindowPtr WindowsTempSortBuffer;
  !     ImVector_ImGuiWindowStackData CurrentWindowStack;
  !     ImGuiStorage WindowsById;
  !     int WindowsActiveCount;
  !     ImVec2 WindowsHoverPadding;
  !     ImGuiWindow* CurrentWindow;
  !     ImGuiWindow* HoveredWindow;
  !     ImGuiWindow* HoveredWindowUnderMovingWindow;
  !     ImGuiDockNode* HoveredDockNode;
  !     ImGuiWindow* MovingWindow;
  !     ImGuiWindow* WheelingWindow;
  !     ImVec2 WheelingWindowRefMousePos;
  !     float WheelingWindowTimer;
  !     ImGuiID DebugHookIdInfo;
  !     ImGuiID HoveredId;
  !     ImGuiID HoveredIdPreviousFrame;
  !     bool HoveredIdAllowOverlap;
  !     bool HoveredIdUsingMouseWheel;
  !     bool HoveredIdPreviousFrameUsingMouseWheel;
  !     bool HoveredIdDisabled;
  !     float HoveredIdTimer;
  !     float HoveredIdNotActiveTimer;
  !     ImGuiID ActiveId;
  !     ImGuiID ActiveIdIsAlive;
  !     float ActiveIdTimer;
  !     bool ActiveIdIsJustActivated;
  !     bool ActiveIdAllowOverlap;
  !     bool ActiveIdNoClearOnFocusLoss;
  !     bool ActiveIdHasBeenPressedBefore;
  !     bool ActiveIdHasBeenEditedBefore;
  !     bool ActiveIdHasBeenEditedThisFrame;
  !     ImVec2 ActiveIdClickOffset;
  !     ImGuiWindow* ActiveIdWindow;
  !     ImGuiInputSource ActiveIdSource;
  !     int ActiveIdMouseButton;
  !     ImGuiID ActiveIdPreviousFrame;
  !     bool ActiveIdPreviousFrameIsAlive;
  !     bool ActiveIdPreviousFrameHasBeenEditedBefore;
  !     ImGuiWindow* ActiveIdPreviousFrameWindow;
  !     ImGuiID LastActiveId;
  !     float LastActiveIdTimer;
  !     bool ActiveIdUsingMouseWheel;
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
  !     int BeginMenuCount;
  !     ImVector_ImGuiViewportPPtr Viewports;
  !     float CurrentDpiScale;
  !     ImGuiViewportP* CurrentViewport;
  !     ImGuiViewportP* MouseViewport;
  !     ImGuiViewportP* MouseLastHoveredViewport;
  !     ImGuiID PlatformLastFocusedViewportId;
  !     ImGuiPlatformMonitor FallbackMonitor;
  !     int ViewportFrontMostStampCount;
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
  !     bool NavIdIsAlive;
  !     bool NavMousePosDirty;
  !     bool NavDisableHighlight;
  !     bool NavDisableMouseHover;
  !     bool NavAnyRequest;
  !     bool NavInitRequest;
  !     bool NavInitRequestFromMove;
  !     ImGuiID NavInitResultId;
  !     ImRect NavInitResultRectRel;
  !     bool NavMoveSubmitted;
  !     bool NavMoveScoringItems;
  !     bool NavMoveForwardToNextFrame;
  !     ImGuiNavMoveFlags NavMoveFlags;
  !     ImGuiScrollFlags NavMoveScrollFlags;
  !     ImGuiModFlags NavMoveKeyMods;
  !     ImGuiDir NavMoveDir;
  !     ImGuiDir NavMoveDirForDebug;
  !     ImGuiDir NavMoveClipDir;
  !     ImRect NavScoringRect;
  !     ImRect NavScoringNoClipRect;
  !     int NavScoringDebugCount;
  !     int NavTabbingDir;
  !     int NavTabbingCounter;
  !     ImGuiNavItemData NavMoveResultLocal;
  !     ImGuiNavItemData NavMoveResultLocalVisible;
  !     ImGuiNavItemData NavMoveResultOther;
  !     ImGuiNavItemData NavTabbingResultFirst;
  !     ImGuiWindow* NavWindowingTarget;
  !     ImGuiWindow* NavWindowingTargetAnim;
  !     ImGuiWindow* NavWindowingListWindow;
  !     float NavWindowingTimer;
  !     float NavWindowingHighlightAlpha;
  !     bool NavWindowingToggleLayer;
  !     float DimBgRatio;
  !     ImGuiMouseCursor MouseCursor;
  !     bool DragDropActive;
  !     bool DragDropWithinSource;
  !     bool DragDropWithinTarget;
  !     ImGuiDragDropFlags DragDropSourceFlags;
  !     int DragDropSourceFrameCount;
  !     int DragDropMouseButton;
  !     ImGuiPayload DragDropPayload;
  !     ImRect DragDropTargetRect;
  !     ImGuiID DragDropTargetId;
  !     ImGuiDragDropFlags DragDropAcceptFlags;
  !     float DragDropAcceptIdCurrRectSurface;
  !     ImGuiID DragDropAcceptIdCurr;
  !     ImGuiID DragDropAcceptIdPrev;
  !     int DragDropAcceptFrameCount;
  !     ImGuiID DragDropHoldJustPressedId;
  !     ImVector_unsigned_char DragDropPayloadBufHeap;
  !     unsigned char DragDropPayloadBufLocal[16];
  !     int ClipperTempDataStacked;
  !     ImVector_ImGuiListClipperData ClipperTempData;
  !     ImGuiTable* CurrentTable;
  !     int TablesTempDataStacked;
  !     ImVector_ImGuiTableTempData TablesTempData;
  !     ImPool_ImGuiTable Tables;
  !     ImVector_float TablesLastTimeActive;
  !     ImVector_ImDrawChannel DrawChannelsTempMergeBuffer;
  !     ImGuiTabBar* CurrentTabBar;
  !     ImPool_ImGuiTabBar TabBars;
  !     ImVector_ImGuiPtrOrIndex CurrentTabBarStack;
  !     ImVector_ImGuiShrinkWidthItem ShrinkWidthBuffer;
  !     ImVec2 MouseLastValidPos;
  !     ImGuiInputTextState InputTextState;
  !     ImFont InputTextPasswordFont;
  !     ImGuiID TempInputId;
  !     ImGuiColorEditFlags ColorEditOptions;
  !     float ColorEditLastHue;
  !     float ColorEditLastSat;
  !     ImU32 ColorEditLastColor;
  !     ImVec4 ColorPickerRef;
  !     ImGuiComboPreviewData ComboPreviewData;
  !     float SliderGrabClickOffset;
  !     float SliderCurrentAccum;
  !     bool SliderCurrentAccumDirty;
  !     bool DragCurrentAccumDirty;
  !     float DragCurrentAccum;
  !     float DragSpeedDefaultRatio;
  !     float ScrollbarClickDeltaToGrabCenter;
  !     float DisabledAlphaBackup;
  !     short DisabledStackSize;
  !     short TooltipOverrideCount;
  !     float TooltipSlowDelay;
  !     ImVector_char ClipboardHandlerData;
  !     ImVector_ImGuiID MenusIdSubmittedThisFrame;
  !     ImGuiPlatformImeData PlatformImeData;
  !     ImGuiPlatformImeData PlatformImeDataPrev;
  !     ImGuiID PlatformImeViewport;
  !     char PlatformLocaleDecimalPoint;
  !     ImGuiDockContext DockContext;
  !     bool SettingsLoaded;
  !     float SettingsDirtyTimer;
  !     ImGuiTextBuffer SettingsIniData;
  !     ImVector_ImGuiSettingsHandler SettingsHandlers;
  !     ImChunkStream_ImGuiWindowSettings SettingsWindows;
  !     ImChunkStream_ImGuiTableSettings SettingsTables;
  !     ImVector_ImGuiContextHook Hooks;
  !     ImGuiID HookIdNext;
  !     bool LogEnabled;
  !     ImGuiLogType LogType;
  !     ImFileHandle LogFile;
  !     ImGuiTextBuffer LogBuffer;
  !     const char* LogNextPrefix;
  !     const char* LogNextSuffix;
  !     float LogLinePosY;
  !     bool LogLineFirstItem;
  !     int LogDepthRef;
  !     int LogDepthToExpand;
  !     int LogDepthToExpandDefault;
  !     ImGuiDebugLogFlags DebugLogFlags;
  !     ImGuiTextBuffer DebugLogBuf;
  !     bool DebugItemPickerActive;
  !     ImGuiID DebugItemPickerBreakId;
  !     ImGuiMetricsConfig DebugMetricsConfig;
  !     ImGuiStackTool DebugStackTool;
  !     float FramerateSecPerFrame[120];
  !     int FramerateSecPerFrameIdx;
  !     int FramerateSecPerFrameCount;
  !     float FramerateSecPerFrameAccum;
  !     int WantCaptureMouseNextFrame;
  !     int WantCaptureKeyboardNextFrame;
  !     int WantTextInputNextFrame;
  !     ImVector_char TempBuffer;
  ! };
  ! struct ImGuiWindowTempData
  ! {
  !     ImVec2 CursorPos;
  !     ImVec2 CursorPosPrevLine;
  !     ImVec2 CursorStartPos;
  !     ImVec2 CursorMaxPos;
  !     ImVec2 IdealMaxPos;
  !     ImVec2 CurrLineSize;
  !     ImVec2 PrevLineSize;
  !     float CurrLineTextBaseOffset;
  !     float PrevLineTextBaseOffset;
  !     bool IsSameLine;
  !     ImVec1 Indent;
  !     ImVec1 ColumnsOffset;
  !     ImVec1 GroupOffset;
  !     ImVec2 CursorStartPosLossyness;
  !     ImGuiNavLayer NavLayerCurrent;
  !     short NavLayersActiveMask;
  !     short NavLayersActiveMaskNext;
  !     ImGuiID NavFocusScopeIdCurrent;
  !     bool NavHideHighlightOneFrame;
  !     bool NavHasScroll;
  !     bool MenuBarAppending;
  !     ImVec2 MenuBarOffset;
  !     ImGuiMenuColumns MenuColumns;
  !     int TreeDepth;
  !     ImU32 TreeJumpToParentOnPopMask;
  !     ImVector_ImGuiWindowPtr ChildWindows;
  !     ImGuiStorage* StateStorage;
  !     ImGuiOldColumns* CurrentColumns;
  !     int CurrentTableIdx;
  !     ImGuiLayoutType LayoutType;
  !     ImGuiLayoutType ParentLayoutType;
  !     float ItemWidth;
  !     float TextWrapPos;
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
  !     ImVec2 ViewportPos;
  !     int ViewportAllowPlatformMonitorExtend;
  !     ImVec2 Pos;
  !     ImVec2 Size;
  !     ImVec2 SizeFull;
  !     ImVec2 ContentSize;
  !     ImVec2 ContentSizeIdeal;
  !     ImVec2 ContentSizeExplicit;
  !     ImVec2 WindowPadding;
  !     float WindowRounding;
  !     float WindowBorderSize;
  !     int NameBufLen;
  !     ImGuiID MoveId;
  !     ImGuiID TabId;
  !     ImGuiID ChildId;
  !     ImVec2 Scroll;
  !     ImVec2 ScrollMax;
  !     ImVec2 ScrollTarget;
  !     ImVec2 ScrollTargetCenterRatio;
  !     ImVec2 ScrollTargetEdgeSnapDist;
  !     ImVec2 ScrollbarSizes;
  !     bool ScrollbarX, ScrollbarY;
  !     bool ViewportOwned;
  !     bool Active;
  !     bool WasActive;
  !     bool WriteAccessed;
  !     bool Collapsed;
  !     bool WantCollapseToggle;
  !     bool SkipItems;
  !     bool Appearing;
  !     bool Hidden;
  !     bool IsFallbackWindow;
  !     bool IsExplicitChild;
  !     bool HasCloseButton;
  !     signed char ResizeBorderHeld;
  !     short BeginCount;
  !     short BeginOrderWithinParent;
  !     short BeginOrderWithinContext;
  !     short FocusOrder;
  !     ImGuiID PopupId;
  !     ImS8 AutoFitFramesX, AutoFitFramesY;
  !     ImS8 AutoFitChildAxises;
  !     bool AutoFitOnlyGrows;
  !     ImGuiDir AutoPosLastDirection;
  !     ImS8 HiddenFramesCanSkipItems;
  !     ImS8 HiddenFramesCannotSkipItems;
  !     ImS8 HiddenFramesForRenderOnly;
  !     ImS8 DisableInputsFrames;
  !     ImGuiCond SetWindowPosAllowFlags : 8;
  !     ImGuiCond SetWindowSizeAllowFlags : 8;
  !     ImGuiCond SetWindowCollapsedAllowFlags : 8;
  !     ImGuiCond SetWindowDockAllowFlags : 8;
  !     ImVec2 SetWindowPosVal;
  !     ImVec2 SetWindowPosPivot;
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
  !     int LastFrameActive;
  !     int LastFrameJustFocused;
  !     float LastTimeActive;
  !     float ItemWidthDefault;
  !     ImGuiStorage StateStorage;
  !     ImVector_ImGuiOldColumns ColumnsStorage;
  !     float FontWindowScale;
  !     float FontDpiScale;
  !     int SettingsOffset;
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
  !     int MemoryDrawListIdxCapacity;
  !     int MemoryDrawListVtxCapacity;
  !     bool MemoryCompacted;
  !     bool DockIsActive :1;
  !     bool DockNodeIsVisible :1;
  !     bool DockTabIsVisible :1;
  !     bool DockTabWantClose :1;
  !     short DockOrder;
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
  !     int LastFrameVisible;
  !     int LastFrameSelected;
  !     float Offset;
  !     float Width;
  !     float ContentWidth;
  !     float RequestedWidth;
  !     ImS32 NameOffset;
  !     ImS16 BeginOrder;
  !     ImS16 IndexDuringLayout;
  !     bool WantClose;
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
  !     int CurrFrameVisible;
  !     int PrevFrameVisible;
  !     ImRect BarRect;
  !     float CurrTabsContentsHeight;
  !     float PrevTabsContentsHeight;
  !     float WidthAllTabs;
  !     float WidthAllTabsIdeal;
  !     float ScrollingAnim;
  !     float ScrollingTarget;
  !     float ScrollingTargetDistToVisibility;
  !     float ScrollingSpeed;
  !     float ScrollingRectMinX;
  !     float ScrollingRectMaxX;
  !     ImGuiID ReorderRequestTabId;
  !     ImS16 ReorderRequestOffset;
  !     ImS8 BeginCount;
  !     bool WantLayout;
  !     bool VisibleTabWasSubmitted;
  !     bool TabsAddedNew;
  !     ImS16 TabsActiveCount;
  !     ImS16 LastTabItemIdx;
  !     float ItemSpacingY;
  !     ImVec2 FramePadding;
  !     ImVec2 BackupCursorPos;
  !     ImGuiTextBuffer TabsNames;
  ! };
  ! struct ImGuiTableColumn
  ! {
  !     ImGuiTableColumnFlags Flags;
  !     float WidthGiven;
  !     float MinX;
  !     float MaxX;
  !     float WidthRequest;
  !     float WidthAuto;
  !     float StretchWeight;
  !     float InitStretchWeightOrWidth;
  !     ImRect ClipRect;
  !     ImGuiID UserID;
  !     float WorkMinX;
  !     float WorkMaxX;
  !     float ItemWidth;
  !     float ContentMaxXFrozen;
  !     float ContentMaxXUnfrozen;
  !     float ContentMaxXHeadersUsed;
  !     float ContentMaxXHeadersIdeal;
  !     ImS16 NameOffset;
  !     ImGuiTableColumnIdx DisplayOrder;
  !     ImGuiTableColumnIdx IndexWithinEnabledSet;
  !     ImGuiTableColumnIdx PrevEnabledColumn;
  !     ImGuiTableColumnIdx NextEnabledColumn;
  !     ImGuiTableColumnIdx SortOrder;
  !     ImGuiTableDrawChannelIdx DrawChannelCurrent;
  !     ImGuiTableDrawChannelIdx DrawChannelFrozen;
  !     ImGuiTableDrawChannelIdx DrawChannelUnfrozen;
  !     bool IsEnabled;
  !     bool IsUserEnabled;
  !     bool IsUserEnabledNextFrame;
  !     bool IsVisibleX;
  !     bool IsVisibleY;
  !     bool IsRequestOutput;
  !     bool IsSkipItems;
  !     bool IsPreserveWidthAuto;
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
  !     float LastOuterHeight;
  !     float LastFirstRowHeight;
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
  !     void* RawData;
  !     ImGuiTableTempData* TempData;
  !     ImSpan_ImGuiTableColumn Columns;
  !     ImSpan_ImGuiTableColumnIdx DisplayOrderToIndex;
  !     ImSpan_ImGuiTableCellData RowCellData;
  !     ImU64 EnabledMaskByDisplayOrder;
  !     ImU64 EnabledMaskByIndex;
  !     ImU64 VisibleMaskByIndex;
  !     ImU64 RequestOutputMaskByIndex;
  !     ImGuiTableFlags SettingsLoadedFlags;
  !     int SettingsOffset;
  !     int LastFrameActive;
  !     int ColumnsCount;
  !     int CurrentRow;
  !     int CurrentColumn;
  !     ImS16 InstanceCurrent;
  !     ImS16 InstanceInteracted;
  !     float RowPosY1;
  !     float RowPosY2;
  !     float RowMinHeight;
  !     float RowTextBaseline;
  !     float RowIndentOffsetX;
  !     ImGuiTableRowFlags RowFlags : 16;
  !     ImGuiTableRowFlags LastRowFlags : 16;
  !     int RowBgColorCounter;
  !     ImU32 RowBgColor[2];
  !     ImU32 BorderColorStrong;
  !     ImU32 BorderColorLight;
  !     float BorderX1;
  !     float BorderX2;
  !     float HostIndentX;
  !     float MinColumnWidth;
  !     float OuterPaddingX;
  !     float CellPaddingX;
  !     float CellPaddingY;
  !     float CellSpacingX1;
  !     float CellSpacingX2;
  !     float InnerWidth;
  !     float ColumnsGivenWidth;
  !     float ColumnsAutoFitWidth;
  !     float ColumnsStretchSumWeights;
  !     float ResizedColumnNextWidth;
  !     float ResizeLockMinContentsX2;
  !     float RefScale;
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
  !     bool IsLayoutLocked;
  !     bool IsInsideRow;
  !     bool IsInitializing;
  !     bool IsSortSpecsDirty;
  !     bool IsUsingHeaders;
  !     bool IsContextPopupOpen;
  !     bool IsSettingsRequestLoad;
  !     bool IsSettingsDirty;
  !     bool IsDefaultDisplayOrder;
  !     bool IsResetAllRequest;
  !     bool IsResetDisplayOrderRequest;
  !     bool IsUnfrozenRows;
  !     bool IsDefaultSizingPolicy;
  !     bool MemoryCompacted;
  !     bool HostSkipItems;
  ! };
  ! struct ImGuiTableTempData
  ! {
  !     int TableIndex;
  !     float LastTimeActive;
  !     ImVec2 UserOuterSize;
  !     ImDrawListSplitter DrawSplitter;
  !     ImRect HostBackupWorkRect;
  !     ImRect HostBackupParentWorkRect;
  !     ImVec2 HostBackupPrevLineSize;
  !     ImVec2 HostBackupCurrLineSize;
  !     ImVec2 HostBackupCursorMaxPos;
  !     ImVec1 HostBackupColumnsOffset;
  !     float HostBackupItemWidth;
  !     int HostBackupItemWidthStackSize;
  ! };
  ! struct ImGuiTableColumnSettings
  ! {
  !     float WidthOrWeight;
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
  !     float RefScale;
  !     ImGuiTableColumnIdx ColumnsCount;
  !     ImGuiTableColumnIdx ColumnsCountMax;
  !     bool WantApply;
  ! };
  ! struct ImFontBuilderIO
  ! {
  !     bool (*FontBuilder_Build)(ImFontAtlas* atlas);
  ! };
