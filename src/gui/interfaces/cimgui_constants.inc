  ! -*-F90-*-

  !xx! cimgui constants

  !! sizes !!!
  !! These need to be the same as in cimgui.h. They are here because
  !! they are sizes for components in C-interoperable user-defined
  !! types, which need to be known at compile-time. Since bind(c) and
  !! parameter are incompatible, it is not possible to use the
  !! "const_" trick below to do this. Hence, copied here.
  !! Note to self: remember to copy them over if imgui/cimgui is updated.
  ! enum ImGuiKey_
  integer(c_int), parameter :: ImGuiKey_COUNT = 645
  integer(c_int), parameter :: ImGuiKey_KeysData_SIZE = 645
  ! enum ImGuiNavInput_
  integer(c_int), parameter :: ImGuiNavInput_COUNT = 20
  !! end of sizes !!

  ! typedef enum {
  !     ImGuiWindowFlags_None = 0,
  !     ImGuiWindowFlags_NoTitleBar = 1 << 0,
  !     ImGuiWindowFlags_NoResize = 1 << 1,
  !     ImGuiWindowFlags_NoMove = 1 << 2,
  !     ImGuiWindowFlags_NoScrollbar = 1 << 3,
  !     ImGuiWindowFlags_NoScrollWithMouse = 1 << 4,
  !     ImGuiWindowFlags_NoCollapse = 1 << 5,
  !     ImGuiWindowFlags_AlwaysAutoResize = 1 << 6,
  !     ImGuiWindowFlags_NoBackground = 1 << 7,
  !     ImGuiWindowFlags_NoSavedSettings = 1 << 8,
  !     ImGuiWindowFlags_NoMouseInputs = 1 << 9,
  !     ImGuiWindowFlags_MenuBar = 1 << 10,
  !     ImGuiWindowFlags_HorizontalScrollbar = 1 << 11,
  !     ImGuiWindowFlags_NoFocusOnAppearing = 1 << 12,
  !     ImGuiWindowFlags_NoBringToFrontOnFocus = 1 << 13,
  !     ImGuiWindowFlags_AlwaysVerticalScrollbar= 1 << 14,
  !     ImGuiWindowFlags_AlwaysHorizontalScrollbar=1<< 15,
  !     ImGuiWindowFlags_AlwaysUseWindowPadding = 1 << 16,
  !     ImGuiWindowFlags_NoNavInputs = 1 << 18,
  !     ImGuiWindowFlags_NoNavFocus = 1 << 19,
  !     ImGuiWindowFlags_UnsavedDocument = 1 << 20,
  !     ImGuiWindowFlags_NoDocking = 1 << 21,
  !     ImGuiWindowFlags_NoNav = ImGuiWindowFlags_NoNavInputs | ImGuiWindowFlags_NoNavFocus,
  !     ImGuiWindowFlags_NoDecoration = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoCollapse,
  !     ImGuiWindowFlags_NoInputs = ImGuiWindowFlags_NoMouseInputs | ImGuiWindowFlags_NoNavInputs | ImGuiWindowFlags_NoNavFocus,
  !     ImGuiWindowFlags_NavFlattened = 1 << 23,
  !     ImGuiWindowFlags_ChildWindow = 1 << 24,
  !     ImGuiWindowFlags_Tooltip = 1 << 25,
  !     ImGuiWindowFlags_Popup = 1 << 26,
  !     ImGuiWindowFlags_Modal = 1 << 27,
  !     ImGuiWindowFlags_ChildMenu = 1 << 28,
  !     ImGuiWindowFlags_DockNodeHost = 1 << 29
  ! }ImGuiWindowFlags_;
  ! typedef enum {
  !     ImGuiInputTextFlags_None = 0,
  !     ImGuiInputTextFlags_CharsDecimal = 1 << 0,
  !     ImGuiInputTextFlags_CharsHexadecimal = 1 << 1,
  !     ImGuiInputTextFlags_CharsUppercase = 1 << 2,
  !     ImGuiInputTextFlags_CharsNoBlank = 1 << 3,
  !     ImGuiInputTextFlags_AutoSelectAll = 1 << 4,
  !     ImGuiInputTextFlags_EnterReturnsTrue = 1 << 5,
  !     ImGuiInputTextFlags_CallbackCompletion = 1 << 6,
  !     ImGuiInputTextFlags_CallbackHistory = 1 << 7,
  !     ImGuiInputTextFlags_CallbackAlways = 1 << 8,
  !     ImGuiInputTextFlags_CallbackCharFilter = 1 << 9,
  !     ImGuiInputTextFlags_AllowTabInput = 1 << 10,
  !     ImGuiInputTextFlags_CtrlEnterForNewLine = 1 << 11,
  !     ImGuiInputTextFlags_NoHorizontalScroll = 1 << 12,
  !     ImGuiInputTextFlags_AlwaysOverwrite = 1 << 13,
  !     ImGuiInputTextFlags_ReadOnly = 1 << 14,
  !     ImGuiInputTextFlags_Password = 1 << 15,
  !     ImGuiInputTextFlags_NoUndoRedo = 1 << 16,
  !     ImGuiInputTextFlags_CharsScientific = 1 << 17,
  !     ImGuiInputTextFlags_CallbackResize = 1 << 18,
  !     ImGuiInputTextFlags_CallbackEdit = 1 << 19
  ! }ImGuiInputTextFlags_;
  ! typedef enum {
  !     ImGuiTreeNodeFlags_None = 0,
  !     ImGuiTreeNodeFlags_Selected = 1 << 0,
  !     ImGuiTreeNodeFlags_Framed = 1 << 1,
  !     ImGuiTreeNodeFlags_AllowItemOverlap = 1 << 2,
  !     ImGuiTreeNodeFlags_NoTreePushOnOpen = 1 << 3,
  !     ImGuiTreeNodeFlags_NoAutoOpenOnLog = 1 << 4,
  !     ImGuiTreeNodeFlags_DefaultOpen = 1 << 5,
  !     ImGuiTreeNodeFlags_OpenOnDoubleClick = 1 << 6,
  !     ImGuiTreeNodeFlags_OpenOnArrow = 1 << 7,
  !     ImGuiTreeNodeFlags_Leaf = 1 << 8,
  !     ImGuiTreeNodeFlags_Bullet = 1 << 9,
  !     ImGuiTreeNodeFlags_FramePadding = 1 << 10,
  !     ImGuiTreeNodeFlags_SpanAvailWidth = 1 << 11,
  !     ImGuiTreeNodeFlags_SpanFullWidth = 1 << 12,
  !     ImGuiTreeNodeFlags_NavLeftJumpsBackHere = 1 << 13,
  !     ImGuiTreeNodeFlags_CollapsingHeader = ImGuiTreeNodeFlags_Framed | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_NoAutoOpenOnLog
  ! }ImGuiTreeNodeFlags_;
  ! typedef enum {
  !     ImGuiPopupFlags_None = 0,
  !     ImGuiPopupFlags_MouseButtonLeft = 0,
  !     ImGuiPopupFlags_MouseButtonRight = 1,
  !     ImGuiPopupFlags_MouseButtonMiddle = 2,
  !     ImGuiPopupFlags_MouseButtonMask_ = 0x1F,
  !     ImGuiPopupFlags_MouseButtonDefault_ = 1,
  !     ImGuiPopupFlags_NoOpenOverExistingPopup = 1 << 5,
  !     ImGuiPopupFlags_NoOpenOverItems = 1 << 6,
  !     ImGuiPopupFlags_AnyPopupId = 1 << 7,
  !     ImGuiPopupFlags_AnyPopupLevel = 1 << 8,
  !     ImGuiPopupFlags_AnyPopup = ImGuiPopupFlags_AnyPopupId | ImGuiPopupFlags_AnyPopupLevel
  ! }ImGuiPopupFlags_;
  ! typedef enum {
  !     ImGuiSelectableFlags_None = 0,
  !     ImGuiSelectableFlags_DontClosePopups = 1 << 0,
  !     ImGuiSelectableFlags_SpanAllColumns = 1 << 1,
  !     ImGuiSelectableFlags_AllowDoubleClick = 1 << 2,
  !     ImGuiSelectableFlags_Disabled = 1 << 3,
  !     ImGuiSelectableFlags_AllowItemOverlap = 1 << 4
  ! }ImGuiSelectableFlags_;
  ! typedef enum {
  !     ImGuiComboFlags_None = 0,
  !     ImGuiComboFlags_PopupAlignLeft = 1 << 0,
  !     ImGuiComboFlags_HeightSmall = 1 << 1,
  !     ImGuiComboFlags_HeightRegular = 1 << 2,
  !     ImGuiComboFlags_HeightLarge = 1 << 3,
  !     ImGuiComboFlags_HeightLargest = 1 << 4,
  !     ImGuiComboFlags_NoArrowButton = 1 << 5,
  !     ImGuiComboFlags_NoPreview = 1 << 6,
  !     ImGuiComboFlags_HeightMask_ = ImGuiComboFlags_HeightSmall | ImGuiComboFlags_HeightRegular | ImGuiComboFlags_HeightLarge | ImGuiComboFlags_HeightLargest
  ! }ImGuiComboFlags_;
  ! typedef enum {
  !     ImGuiTabBarFlags_None = 0,
  !     ImGuiTabBarFlags_Reorderable = 1 << 0,
  !     ImGuiTabBarFlags_AutoSelectNewTabs = 1 << 1,
  !     ImGuiTabBarFlags_TabListPopupButton = 1 << 2,
  !     ImGuiTabBarFlags_NoCloseWithMiddleMouseButton = 1 << 3,
  !     ImGuiTabBarFlags_NoTabListScrollingButtons = 1 << 4,
  !     ImGuiTabBarFlags_NoTooltip = 1 << 5,
  !     ImGuiTabBarFlags_FittingPolicyResizeDown = 1 << 6,
  !     ImGuiTabBarFlags_FittingPolicyScroll = 1 << 7,
  !     ImGuiTabBarFlags_FittingPolicyMask_ = ImGuiTabBarFlags_FittingPolicyResizeDown | ImGuiTabBarFlags_FittingPolicyScroll,
  !     ImGuiTabBarFlags_FittingPolicyDefault_ = ImGuiTabBarFlags_FittingPolicyResizeDown
  ! }ImGuiTabBarFlags_;
  ! typedef enum {
  !     ImGuiTabItemFlags_None = 0,
  !     ImGuiTabItemFlags_UnsavedDocument = 1 << 0,
  !     ImGuiTabItemFlags_SetSelected = 1 << 1,
  !     ImGuiTabItemFlags_NoCloseWithMiddleMouseButton = 1 << 2,
  !     ImGuiTabItemFlags_NoPushId = 1 << 3,
  !     ImGuiTabItemFlags_NoTooltip = 1 << 4,
  !     ImGuiTabItemFlags_NoReorder = 1 << 5,
  !     ImGuiTabItemFlags_Leading = 1 << 6,
  !     ImGuiTabItemFlags_Trailing = 1 << 7
  ! }ImGuiTabItemFlags_;
  ! typedef enum {
  !     ImGuiTableFlags_None = 0,
  !     ImGuiTableFlags_Resizable = 1 << 0,
  !     ImGuiTableFlags_Reorderable = 1 << 1,
  !     ImGuiTableFlags_Hideable = 1 << 2,
  !     ImGuiTableFlags_Sortable = 1 << 3,
  !     ImGuiTableFlags_NoSavedSettings = 1 << 4,
  !     ImGuiTableFlags_ContextMenuInBody = 1 << 5,
  !     ImGuiTableFlags_RowBg = 1 << 6,
  !     ImGuiTableFlags_BordersInnerH = 1 << 7,
  !     ImGuiTableFlags_BordersOuterH = 1 << 8,
  !     ImGuiTableFlags_BordersInnerV = 1 << 9,
  !     ImGuiTableFlags_BordersOuterV = 1 << 10,
  !     ImGuiTableFlags_BordersH = ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_BordersOuterH,
  !     ImGuiTableFlags_BordersV = ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_BordersOuterV,
  !     ImGuiTableFlags_BordersInner = ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_BordersInnerH,
  !     ImGuiTableFlags_BordersOuter = ImGuiTableFlags_BordersOuterV | ImGuiTableFlags_BordersOuterH,
  !     ImGuiTableFlags_Borders = ImGuiTableFlags_BordersInner | ImGuiTableFlags_BordersOuter,
  !     ImGuiTableFlags_NoBordersInBody = 1 << 11,
  !     ImGuiTableFlags_NoBordersInBodyUntilResize = 1 << 12,
  !     ImGuiTableFlags_SizingFixedFit = 1 << 13,
  !     ImGuiTableFlags_SizingFixedSame = 2 << 13,
  !     ImGuiTableFlags_SizingStretchProp = 3 << 13,
  !     ImGuiTableFlags_SizingStretchSame = 4 << 13,
  !     ImGuiTableFlags_NoHostExtendX = 1 << 16,
  !     ImGuiTableFlags_NoHostExtendY = 1 << 17,
  !     ImGuiTableFlags_NoKeepColumnsVisible = 1 << 18,
  !     ImGuiTableFlags_PreciseWidths = 1 << 19,
  !     ImGuiTableFlags_NoClip = 1 << 20,
  !     ImGuiTableFlags_PadOuterX = 1 << 21,
  !     ImGuiTableFlags_NoPadOuterX = 1 << 22,
  !     ImGuiTableFlags_NoPadInnerX = 1 << 23,
  !     ImGuiTableFlags_ScrollX = 1 << 24,
  !     ImGuiTableFlags_ScrollY = 1 << 25,
  !     ImGuiTableFlags_SortMulti = 1 << 26,
  !     ImGuiTableFlags_SortTristate = 1 << 27,
  !     ImGuiTableFlags_SizingMask_ = ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_SizingFixedSame | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_SizingStretchSame
  ! }ImGuiTableFlags_;
  ! typedef enum {
  !     ImGuiTableColumnFlags_None = 0,
  !     ImGuiTableColumnFlags_Disabled = 1 << 0,
  !     ImGuiTableColumnFlags_DefaultHide = 1 << 1,
  !     ImGuiTableColumnFlags_DefaultSort = 1 << 2,
  !     ImGuiTableColumnFlags_WidthStretch = 1 << 3,
  !     ImGuiTableColumnFlags_WidthFixed = 1 << 4,
  !     ImGuiTableColumnFlags_NoResize = 1 << 5,
  !     ImGuiTableColumnFlags_NoReorder = 1 << 6,
  !     ImGuiTableColumnFlags_NoHide = 1 << 7,
  !     ImGuiTableColumnFlags_NoClip = 1 << 8,
  !     ImGuiTableColumnFlags_NoSort = 1 << 9,
  !     ImGuiTableColumnFlags_NoSortAscending = 1 << 10,
  !     ImGuiTableColumnFlags_NoSortDescending = 1 << 11,
  !     ImGuiTableColumnFlags_NoHeaderLabel = 1 << 12,
  !     ImGuiTableColumnFlags_NoHeaderWidth = 1 << 13,
  !     ImGuiTableColumnFlags_PreferSortAscending = 1 << 14,
  !     ImGuiTableColumnFlags_PreferSortDescending = 1 << 15,
  !     ImGuiTableColumnFlags_IndentEnable = 1 << 16,
  !     ImGuiTableColumnFlags_IndentDisable = 1 << 17,
  !     ImGuiTableColumnFlags_IsEnabled = 1 << 24,
  !     ImGuiTableColumnFlags_IsVisible = 1 << 25,
  !     ImGuiTableColumnFlags_IsSorted = 1 << 26,
  !     ImGuiTableColumnFlags_IsHovered = 1 << 27,
  !     ImGuiTableColumnFlags_WidthMask_ = ImGuiTableColumnFlags_WidthStretch | ImGuiTableColumnFlags_WidthFixed,
  !     ImGuiTableColumnFlags_IndentMask_ = ImGuiTableColumnFlags_IndentEnable | ImGuiTableColumnFlags_IndentDisable,
  !     ImGuiTableColumnFlags_StatusMask_ = ImGuiTableColumnFlags_IsEnabled | ImGuiTableColumnFlags_IsVisible | ImGuiTableColumnFlags_IsSorted | ImGuiTableColumnFlags_IsHovered,
  !     ImGuiTableColumnFlags_NoDirectResize_ = 1 << 30
  ! }ImGuiTableColumnFlags_;
  ! typedef enum {
  !     ImGuiTableRowFlags_None = 0,
  !     ImGuiTableRowFlags_Headers = 1 << 0
  ! }ImGuiTableRowFlags_;
  ! typedef enum {
  !     ImGuiTableBgTarget_None = 0,
  !     ImGuiTableBgTarget_RowBg0 = 1,
  !     ImGuiTableBgTarget_RowBg1 = 2,
  !     ImGuiTableBgTarget_CellBg = 3
  ! }ImGuiTableBgTarget_;
  ! typedef enum {
  !     ImGuiFocusedFlags_None = 0,
  !     ImGuiFocusedFlags_ChildWindows = 1 << 0,
  !     ImGuiFocusedFlags_RootWindow = 1 << 1,
  !     ImGuiFocusedFlags_AnyWindow = 1 << 2,
  !     ImGuiFocusedFlags_NoPopupHierarchy = 1 << 3,
  !     ImGuiFocusedFlags_DockHierarchy = 1 << 4,
  !     ImGuiFocusedFlags_RootAndChildWindows = ImGuiFocusedFlags_RootWindow | ImGuiFocusedFlags_ChildWindows
  ! }ImGuiFocusedFlags_;
  ! enum ImGuiHoveredFlags_
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_None") :: ImGuiHoveredFlags_None
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_ChildWindows") :: ImGuiHoveredFlags_ChildWindows
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_RootWindow") :: ImGuiHoveredFlags_RootWindow
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AnyWindow") :: ImGuiHoveredFlags_AnyWindow
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_NoPopupHierarchy") :: ImGuiHoveredFlags_NoPopupHierarchy
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_DockHierarchy") :: ImGuiHoveredFlags_DockHierarchy
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenBlockedByPopup") :: &
     ImGuiHoveredFlags_AllowWhenBlockedByPopup
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenBlockedByActiveItem") :: &
     ImGuiHoveredFlags_AllowWhenBlockedByActiveItem
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenOverlapped") :: ImGuiHoveredFlags_AllowWhenOverlapped
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenDisabled") :: ImGuiHoveredFlags_AllowWhenDisabled
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_NoNavOverride") :: ImGuiHoveredFlags_NoNavOverride
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_RectOnly") :: ImGuiHoveredFlags_RectOnly
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_RootAndChildWindows") :: ImGuiHoveredFlags_RootAndChildWindows
  ! typedef enum {
  !     ImGuiDockNodeFlags_None = 0,
  !     ImGuiDockNodeFlags_KeepAliveOnly = 1 << 0,
  !     ImGuiDockNodeFlags_NoDockingInCentralNode = 1 << 2,
  !     ImGuiDockNodeFlags_PassthruCentralNode = 1 << 3,
  !     ImGuiDockNodeFlags_NoSplit = 1 << 4,
  !     ImGuiDockNodeFlags_NoResize = 1 << 5,
  !     ImGuiDockNodeFlags_AutoHideTabBar = 1 << 6
  ! }ImGuiDockNodeFlags_;
  ! typedef enum {
  !     ImGuiDragDropFlags_None = 0,
  !     ImGuiDragDropFlags_SourceNoPreviewTooltip = 1 << 0,
  !     ImGuiDragDropFlags_SourceNoDisableHover = 1 << 1,
  !     ImGuiDragDropFlags_SourceNoHoldToOpenOthers = 1 << 2,
  !     ImGuiDragDropFlags_SourceAllowNullID = 1 << 3,
  !     ImGuiDragDropFlags_SourceExtern = 1 << 4,
  !     ImGuiDragDropFlags_SourceAutoExpirePayload = 1 << 5,
  !     ImGuiDragDropFlags_AcceptBeforeDelivery = 1 << 10,
  !     ImGuiDragDropFlags_AcceptNoDrawDefaultRect = 1 << 11,
  !     ImGuiDragDropFlags_AcceptNoPreviewTooltip = 1 << 12,
  !     ImGuiDragDropFlags_AcceptPeekOnly = ImGuiDragDropFlags_AcceptBeforeDelivery | ImGuiDragDropFlags_AcceptNoDrawDefaultRect
  ! }ImGuiDragDropFlags_;
  ! typedef enum {
  !     ImGuiDataType_S8,
  !     ImGuiDataType_U8,
  !     ImGuiDataType_S16,
  !     ImGuiDataType_U16,
  !     ImGuiDataType_S32,
  !     ImGuiDataType_U32,
  !     ImGuiDataType_S64,
  !     ImGuiDataType_U64,
  !     ImGuiDataType_Float,
  !     ImGuiDataType_Double,
  !     ImGuiDataType_COUNT
  ! }ImGuiDataType_;
  ! typedef enum {
  !     ImGuiDir_None = -1,
  !     ImGuiDir_Left = 0,
  !     ImGuiDir_Right = 1,
  !     ImGuiDir_Up = 2,
  !     ImGuiDir_Down = 3,
  !     ImGuiDir_COUNT
  ! }ImGuiDir_;
  ! typedef enum {
  !     ImGuiSortDirection_None = 0,
  !     ImGuiSortDirection_Ascending = 1,
  !     ImGuiSortDirection_Descending = 2
  ! }ImGuiSortDirection_;
  ! enum ImGuiKey_
  integer(c_int), bind(C,name="const_ImGuiKey_None") :: ImGuiKey_None
  integer(c_int), bind(C,name="const_ImGuiKey_Tab") :: ImGuiKey_Tab
  integer(c_int), bind(C,name="const_ImGuiKey_LeftArrow") :: ImGuiKey_LeftArrow
  integer(c_int), bind(C,name="const_ImGuiKey_RightArrow") :: ImGuiKey_RightArrow
  integer(c_int), bind(C,name="const_ImGuiKey_UpArrow") :: ImGuiKey_UpArrow
  integer(c_int), bind(C,name="const_ImGuiKey_DownArrow") :: ImGuiKey_DownArrow
  integer(c_int), bind(C,name="const_ImGuiKey_PageUp") :: ImGuiKey_PageUp
  integer(c_int), bind(C,name="const_ImGuiKey_PageDown") :: ImGuiKey_PageDown
  integer(c_int), bind(C,name="const_ImGuiKey_Home") :: ImGuiKey_Home
  integer(c_int), bind(C,name="const_ImGuiKey_End") :: ImGuiKey_End
  integer(c_int), bind(C,name="const_ImGuiKey_Insert") :: ImGuiKey_Insert
  integer(c_int), bind(C,name="const_ImGuiKey_Delete") :: ImGuiKey_Delete
  integer(c_int), bind(C,name="const_ImGuiKey_Backspace") :: ImGuiKey_Backspace
  integer(c_int), bind(C,name="const_ImGuiKey_Space") :: ImGuiKey_Space
  integer(c_int), bind(C,name="const_ImGuiKey_Enter") :: ImGuiKey_Enter
  integer(c_int), bind(C,name="const_ImGuiKey_Escape") :: ImGuiKey_Escape
  integer(c_int), bind(C,name="const_ImGuiKey_LeftCtrl") :: ImGuiKey_LeftCtrl
  integer(c_int), bind(C,name="const_ImGuiKey_LeftShift") :: ImGuiKey_LeftShift
  integer(c_int), bind(C,name="const_ImGuiKey_LeftAlt") :: ImGuiKey_LeftAlt
  integer(c_int), bind(C,name="const_ImGuiKey_LeftSuper") :: ImGuiKey_LeftSuper
  integer(c_int), bind(C,name="const_ImGuiKey_RightCtrl") :: ImGuiKey_RightCtrl
  integer(c_int), bind(C,name="const_ImGuiKey_RightShift") :: ImGuiKey_RightShift
  integer(c_int), bind(C,name="const_ImGuiKey_RightAlt") :: ImGuiKey_RightAlt
  integer(c_int), bind(C,name="const_ImGuiKey_RightSuper") :: ImGuiKey_RightSuper
  integer(c_int), bind(C,name="const_ImGuiKey_Menu") :: ImGuiKey_Menu
  integer(c_int), bind(C,name="const_ImGuiKey_0") :: ImGuiKey_0
  integer(c_int), bind(C,name="const_ImGuiKey_1") :: ImGuiKey_1
  integer(c_int), bind(C,name="const_ImGuiKey_2") :: ImGuiKey_2
  integer(c_int), bind(C,name="const_ImGuiKey_3") :: ImGuiKey_3
  integer(c_int), bind(C,name="const_ImGuiKey_4") :: ImGuiKey_4
  integer(c_int), bind(C,name="const_ImGuiKey_5") :: ImGuiKey_5
  integer(c_int), bind(C,name="const_ImGuiKey_6") :: ImGuiKey_6
  integer(c_int), bind(C,name="const_ImGuiKey_7") :: ImGuiKey_7
  integer(c_int), bind(C,name="const_ImGuiKey_8") :: ImGuiKey_8
  integer(c_int), bind(C,name="const_ImGuiKey_9") :: ImGuiKey_9
  integer(c_int), bind(C,name="const_ImGuiKey_A") :: ImGuiKey_A
  integer(c_int), bind(C,name="const_ImGuiKey_B") :: ImGuiKey_B
  integer(c_int), bind(C,name="const_ImGuiKey_C") :: ImGuiKey_C
  integer(c_int), bind(C,name="const_ImGuiKey_D") :: ImGuiKey_D
  integer(c_int), bind(C,name="const_ImGuiKey_E") :: ImGuiKey_E
  integer(c_int), bind(C,name="const_ImGuiKey_F") :: ImGuiKey_F
  integer(c_int), bind(C,name="const_ImGuiKey_G") :: ImGuiKey_G
  integer(c_int), bind(C,name="const_ImGuiKey_H") :: ImGuiKey_H
  integer(c_int), bind(C,name="const_ImGuiKey_I") :: ImGuiKey_I
  integer(c_int), bind(C,name="const_ImGuiKey_J") :: ImGuiKey_J
  integer(c_int), bind(C,name="const_ImGuiKey_K") :: ImGuiKey_K
  integer(c_int), bind(C,name="const_ImGuiKey_L") :: ImGuiKey_L
  integer(c_int), bind(C,name="const_ImGuiKey_M") :: ImGuiKey_M
  integer(c_int), bind(C,name="const_ImGuiKey_N") :: ImGuiKey_N
  integer(c_int), bind(C,name="const_ImGuiKey_O") :: ImGuiKey_O
  integer(c_int), bind(C,name="const_ImGuiKey_P") :: ImGuiKey_P
  integer(c_int), bind(C,name="const_ImGuiKey_Q") :: ImGuiKey_Q
  integer(c_int), bind(C,name="const_ImGuiKey_R") :: ImGuiKey_R
  integer(c_int), bind(C,name="const_ImGuiKey_S") :: ImGuiKey_S
  integer(c_int), bind(C,name="const_ImGuiKey_T") :: ImGuiKey_T
  integer(c_int), bind(C,name="const_ImGuiKey_U") :: ImGuiKey_U
  integer(c_int), bind(C,name="const_ImGuiKey_V") :: ImGuiKey_V
  integer(c_int), bind(C,name="const_ImGuiKey_W") :: ImGuiKey_W
  integer(c_int), bind(C,name="const_ImGuiKey_X") :: ImGuiKey_X
  integer(c_int), bind(C,name="const_ImGuiKey_Y") :: ImGuiKey_Y
  integer(c_int), bind(C,name="const_ImGuiKey_Z") :: ImGuiKey_Z
  integer(c_int), bind(C,name="const_ImGuiKey_F1") :: ImGuiKey_F1
  integer(c_int), bind(C,name="const_ImGuiKey_F2") :: ImGuiKey_F2
  integer(c_int), bind(C,name="const_ImGuiKey_F3") :: ImGuiKey_F3
  integer(c_int), bind(C,name="const_ImGuiKey_F4") :: ImGuiKey_F4
  integer(c_int), bind(C,name="const_ImGuiKey_F5") :: ImGuiKey_F5
  integer(c_int), bind(C,name="const_ImGuiKey_F6") :: ImGuiKey_F6
  integer(c_int), bind(C,name="const_ImGuiKey_F7") :: ImGuiKey_F7
  integer(c_int), bind(C,name="const_ImGuiKey_F8") :: ImGuiKey_F8
  integer(c_int), bind(C,name="const_ImGuiKey_F9") :: ImGuiKey_F9
  integer(c_int), bind(C,name="const_ImGuiKey_F10") :: ImGuiKey_F10
  integer(c_int), bind(C,name="const_ImGuiKey_F11") :: ImGuiKey_F11
  integer(c_int), bind(C,name="const_ImGuiKey_F12") :: ImGuiKey_F12
  integer(c_int), bind(C,name="const_ImGuiKey_Apostrophe") :: ImGuiKey_Apostrophe
  integer(c_int), bind(C,name="const_ImGuiKey_Comma") :: ImGuiKey_Comma
  integer(c_int), bind(C,name="const_ImGuiKey_Minus") :: ImGuiKey_Minus
  integer(c_int), bind(C,name="const_ImGuiKey_Period") :: ImGuiKey_Period
  integer(c_int), bind(C,name="const_ImGuiKey_Slash") :: ImGuiKey_Slash
  integer(c_int), bind(C,name="const_ImGuiKey_Semicolon") :: ImGuiKey_Semicolon
  integer(c_int), bind(C,name="const_ImGuiKey_Equal") :: ImGuiKey_Equal
  integer(c_int), bind(C,name="const_ImGuiKey_LeftBracket") :: ImGuiKey_LeftBracket
  integer(c_int), bind(C,name="const_ImGuiKey_Backslash") :: ImGuiKey_Backslash
  integer(c_int), bind(C,name="const_ImGuiKey_RightBracket") :: ImGuiKey_RightBracket
  integer(c_int), bind(C,name="const_ImGuiKey_GraveAccent") :: ImGuiKey_GraveAccent
  integer(c_int), bind(C,name="const_ImGuiKey_CapsLock") :: ImGuiKey_CapsLock
  integer(c_int), bind(C,name="const_ImGuiKey_ScrollLock") :: ImGuiKey_ScrollLock
  integer(c_int), bind(C,name="const_ImGuiKey_NumLock") :: ImGuiKey_NumLock
  integer(c_int), bind(C,name="const_ImGuiKey_PrintScreen") :: ImGuiKey_PrintScreen
  integer(c_int), bind(C,name="const_ImGuiKey_Pause") :: ImGuiKey_Pause
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad0") :: ImGuiKey_Keypad0
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad1") :: ImGuiKey_Keypad1
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad2") :: ImGuiKey_Keypad2
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad3") :: ImGuiKey_Keypad3
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad4") :: ImGuiKey_Keypad4
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad5") :: ImGuiKey_Keypad5
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad6") :: ImGuiKey_Keypad6
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad7") :: ImGuiKey_Keypad7
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad8") :: ImGuiKey_Keypad8
  integer(c_int), bind(C,name="const_ImGuiKey_Keypad9") :: ImGuiKey_Keypad9
  integer(c_int), bind(C,name="const_ImGuiKey_KeypadDecimal") :: ImGuiKey_KeypadDecimal
  integer(c_int), bind(C,name="const_ImGuiKey_KeypadDivide") :: ImGuiKey_KeypadDivide
  integer(c_int), bind(C,name="const_ImGuiKey_KeypadMultiply") :: ImGuiKey_KeypadMultiply
  integer(c_int), bind(C,name="const_ImGuiKey_KeypadSubtract") :: ImGuiKey_KeypadSubtract
  integer(c_int), bind(C,name="const_ImGuiKey_KeypadAdd") :: ImGuiKey_KeypadAdd
  integer(c_int), bind(C,name="const_ImGuiKey_KeypadEnter") :: ImGuiKey_KeypadEnter
  integer(c_int), bind(C,name="const_ImGuiKey_KeypadEqual") :: ImGuiKey_KeypadEqual
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadStart") :: ImGuiKey_GamepadStart
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadBack") :: ImGuiKey_GamepadBack
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadFaceUp") :: ImGuiKey_GamepadFaceUp
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadFaceDown") :: ImGuiKey_GamepadFaceDown
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadFaceLeft") :: ImGuiKey_GamepadFaceLeft
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadFaceRight") :: ImGuiKey_GamepadFaceRight
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadDpadUp") :: ImGuiKey_GamepadDpadUp
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadDpadDown") :: ImGuiKey_GamepadDpadDown
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadDpadLeft") :: ImGuiKey_GamepadDpadLeft
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadDpadRight") :: ImGuiKey_GamepadDpadRight
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadL1") :: ImGuiKey_GamepadL1
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadR1") :: ImGuiKey_GamepadR1
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadL2") :: ImGuiKey_GamepadL2
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadR2") :: ImGuiKey_GamepadR2
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadL3") :: ImGuiKey_GamepadL3
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadR3") :: ImGuiKey_GamepadR3
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadLStickUp") :: ImGuiKey_GamepadLStickUp
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadLStickDown") :: ImGuiKey_GamepadLStickDown
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadLStickLeft") :: ImGuiKey_GamepadLStickLeft
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadLStickRight") :: ImGuiKey_GamepadLStickRight
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadRStickUp") :: ImGuiKey_GamepadRStickUp
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadRStickDown") :: ImGuiKey_GamepadRStickDown
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadRStickLeft") :: ImGuiKey_GamepadRStickLeft
  integer(c_int), bind(C,name="const_ImGuiKey_GamepadRStickRight") :: ImGuiKey_GamepadRStickRight
  integer(c_int), bind(C,name="const_ImGuiKey_ModCtrl") :: ImGuiKey_ModCtrl
  integer(c_int), bind(C,name="const_ImGuiKey_ModShift") :: ImGuiKey_ModShift
  integer(c_int), bind(C,name="const_ImGuiKey_ModAlt") :: ImGuiKey_ModAlt
  integer(c_int), bind(C,name="const_ImGuiKey_ModSuper") :: ImGuiKey_ModSuper
  ! integer(c_int), bind(C,name="const_ImGuiKey_COUNT") :: ImGuiKey_COUNT
  integer(c_int), bind(C,name="const_ImGuiKey_NamedKey_BEGIN") :: ImGuiKey_NamedKey_BEGIN
  integer(c_int), bind(C,name="const_ImGuiKey_NamedKey_END") :: ImGuiKey_NamedKey_END
  integer(c_int), bind(C,name="const_ImGuiKey_NamedKey_COUNT") :: ImGuiKey_NamedKey_COUNT
  ! integer(c_int), bind(C,name="const_ImGuiKey_KeysData_SIZE") :: ImGuiKey_KeysData_SIZE
  integer(c_int), bind(C,name="const_ImGuiKey_KeysData_OFFSET") :: ImGuiKey_KeysData_OFFSET
  ! enum ImGuiModFlags_
  integer(c_int), bind(C,name="const_ImGuiModFlags_None") :: ImGuiModFlags_None
  integer(c_int), bind(C,name="const_ImGuiModFlags_Ctrl") :: ImGuiModFlags_Ctrl
  integer(c_int), bind(C,name="const_ImGuiModFlags_Shift") :: ImGuiModFlags_Shift
  integer(c_int), bind(C,name="const_ImGuiModFlags_Alt") :: ImGuiModFlags_Alt
  integer(c_int), bind(C,name="const_ImGuiModFlags_Super") :: ImGuiModFlags_Super
  ! enum ImGuiNavInput_
  integer(c_int), bind(C,name="const_ImGuiNavInput_Activate") :: ImGuiNavInput_Activate
  integer(c_int), bind(C,name="const_ImGuiNavInput_Cancel") :: ImGuiNavInput_Cancel
  integer(c_int), bind(C,name="const_ImGuiNavInput_Input") :: ImGuiNavInput_Input
  integer(c_int), bind(C,name="const_ImGuiNavInput_Menu") :: ImGuiNavInput_Menu
  integer(c_int), bind(C,name="const_ImGuiNavInput_DpadLeft") :: ImGuiNavInput_DpadLeft
  integer(c_int), bind(C,name="const_ImGuiNavInput_DpadRight") :: ImGuiNavInput_DpadRight
  integer(c_int), bind(C,name="const_ImGuiNavInput_DpadUp") :: ImGuiNavInput_DpadUp
  integer(c_int), bind(C,name="const_ImGuiNavInput_DpadDown") :: ImGuiNavInput_DpadDown
  integer(c_int), bind(C,name="const_ImGuiNavInput_LStickLeft") :: ImGuiNavInput_LStickLeft
  integer(c_int), bind(C,name="const_ImGuiNavInput_LStickRight") :: ImGuiNavInput_LStickRight
  integer(c_int), bind(C,name="const_ImGuiNavInput_LStickUp") :: ImGuiNavInput_LStickUp
  integer(c_int), bind(C,name="const_ImGuiNavInput_LStickDown") :: ImGuiNavInput_LStickDown
  integer(c_int), bind(C,name="const_ImGuiNavInput_FocusPrev") :: ImGuiNavInput_FocusPrev
  integer(c_int), bind(C,name="const_ImGuiNavInput_FocusNext") :: ImGuiNavInput_FocusNext
  integer(c_int), bind(C,name="const_ImGuiNavInput_TweakSlow") :: ImGuiNavInput_TweakSlow
  integer(c_int), bind(C,name="const_ImGuiNavInput_TweakFast") :: ImGuiNavInput_TweakFast
  integer(c_int), bind(C,name="const_ImGuiNavInput_KeyLeft_") :: ImGuiNavInput_KeyLeft_
  integer(c_int), bind(C,name="const_ImGuiNavInput_KeyRight_") :: ImGuiNavInput_KeyRight_
  integer(c_int), bind(C,name="const_ImGuiNavInput_KeyUp_") :: ImGuiNavInput_KeyUp_
  integer(c_int), bind(C,name="const_ImGuiNavInput_KeyDown_") :: ImGuiNavInput_KeyDown_
  !integer(c_int), bind(C,name="const_ImGuiNavInput_COUNT") :: ImGuiNavInput_COUNT
  ! enum ImGuiConfigFlags_
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_None") :: ImGuiConfigFlags_None
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_NavEnableKeyboard") ::ImGuiConfigFlags_NavEnableKeyboard
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_NavEnableGamepad") :: ImGuiConfigFlags_NavEnableGamepad
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_NavEnableSetMousePos") :: ImGuiConfigFlags_NavEnableSetMousePos
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_NavNoCaptureKeyboard") :: ImGuiConfigFlags_NavNoCaptureKeyboard
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_NoMouse") :: ImGuiConfigFlags_NoMouse
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_NoMouseCursorChange") :: ImGuiConfigFlags_NoMouseCursorChange
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_DockingEnable") :: ImGuiConfigFlags_DockingEnable
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_ViewportsEnable") :: ImGuiConfigFlags_ViewportsEnable
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_DpiEnableScaleViewport") :: ImGuiConfigFlags_DpiEnableScaleViewport
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_DpiEnableScaleFonts") :: ImGuiConfigFlags_DpiEnableScaleFonts
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_IsSRGB") :: ImGuiConfigFlags_IsSRGB
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_IsTouchScreen") :: ImGuiConfigFlags_IsTouchScreen
  ! typedef enum {
  !     ImGuiBackendFlags_None = 0,
  !     ImGuiBackendFlags_HasGamepad = 1 << 0,
  !     ImGuiBackendFlags_HasMouseCursors = 1 << 1,
  !     ImGuiBackendFlags_HasSetMousePos = 1 << 2,
  !     ImGuiBackendFlags_RendererHasVtxOffset = 1 << 3,
  !     ImGuiBackendFlags_PlatformHasViewports = 1 << 10,
  !     ImGuiBackendFlags_HasMouseHoveredViewport=1 << 11,
  !     ImGuiBackendFlags_RendererHasViewports = 1 << 12
  ! }ImGuiBackendFlags_;
  ! typedef enum {
  !     ImGuiCol_Text,
  !     ImGuiCol_TextDisabled,
  !     ImGuiCol_WindowBg,
  !     ImGuiCol_ChildBg,
  !     ImGuiCol_PopupBg,
  !     ImGuiCol_Border,
  !     ImGuiCol_BorderShadow,
  !     ImGuiCol_FrameBg,
  !     ImGuiCol_FrameBgHovered,
  !     ImGuiCol_FrameBgActive,
  !     ImGuiCol_TitleBg,
  !     ImGuiCol_TitleBgActive,
  !     ImGuiCol_TitleBgCollapsed,
  !     ImGuiCol_MenuBarBg,
  !     ImGuiCol_ScrollbarBg,
  !     ImGuiCol_ScrollbarGrab,
  !     ImGuiCol_ScrollbarGrabHovered,
  !     ImGuiCol_ScrollbarGrabActive,
  !     ImGuiCol_CheckMark,
  !     ImGuiCol_SliderGrab,
  !     ImGuiCol_SliderGrabActive,
  !     ImGuiCol_Button,
  !     ImGuiCol_ButtonHovered,
  !     ImGuiCol_ButtonActive,
  !     ImGuiCol_Header,
  !     ImGuiCol_HeaderHovered,
  !     ImGuiCol_HeaderActive,
  !     ImGuiCol_Separator,
  !     ImGuiCol_SeparatorHovered,
  !     ImGuiCol_SeparatorActive,
  !     ImGuiCol_ResizeGrip,
  !     ImGuiCol_ResizeGripHovered,
  !     ImGuiCol_ResizeGripActive,
  !     ImGuiCol_Tab,
  !     ImGuiCol_TabHovered,
  !     ImGuiCol_TabActive,
  !     ImGuiCol_TabUnfocused,
  !     ImGuiCol_TabUnfocusedActive,
  !     ImGuiCol_DockingPreview,
  !     ImGuiCol_DockingEmptyBg,
  !     ImGuiCol_PlotLines,
  !     ImGuiCol_PlotLinesHovered,
  !     ImGuiCol_PlotHistogram,
  !     ImGuiCol_PlotHistogramHovered,
  !     ImGuiCol_TableHeaderBg,
  !     ImGuiCol_TableBorderStrong,
  !     ImGuiCol_TableBorderLight,
  !     ImGuiCol_TableRowBg,
  !     ImGuiCol_TableRowBgAlt,
  !     ImGuiCol_TextSelectedBg,
  !     ImGuiCol_DragDropTarget,
  !     ImGuiCol_NavHighlight,
  !     ImGuiCol_NavWindowingHighlight,
  !     ImGuiCol_NavWindowingDimBg,
  !     ImGuiCol_ModalWindowDimBg,
  !     ImGuiCol_COUNT
  ! }ImGuiCol_;
  ! typedef enum {
  !     ImGuiStyleVar_Alpha,
  !     ImGuiStyleVar_DisabledAlpha,
  !     ImGuiStyleVar_WindowPadding,
  !     ImGuiStyleVar_WindowRounding,
  !     ImGuiStyleVar_WindowBorderSize,
  !     ImGuiStyleVar_WindowMinSize,
  !     ImGuiStyleVar_WindowTitleAlign,
  !     ImGuiStyleVar_ChildRounding,
  !     ImGuiStyleVar_ChildBorderSize,
  !     ImGuiStyleVar_PopupRounding,
  !     ImGuiStyleVar_PopupBorderSize,
  !     ImGuiStyleVar_FramePadding,
  !     ImGuiStyleVar_FrameRounding,
  !     ImGuiStyleVar_FrameBorderSize,
  !     ImGuiStyleVar_ItemSpacing,
  !     ImGuiStyleVar_ItemInnerSpacing,
  !     ImGuiStyleVar_IndentSpacing,
  !     ImGuiStyleVar_CellPadding,
  !     ImGuiStyleVar_ScrollbarSize,
  !     ImGuiStyleVar_ScrollbarRounding,
  !     ImGuiStyleVar_GrabMinSize,
  !     ImGuiStyleVar_GrabRounding,
  !     ImGuiStyleVar_TabRounding,
  !     ImGuiStyleVar_ButtonTextAlign,
  !     ImGuiStyleVar_SelectableTextAlign,
  !     ImGuiStyleVar_COUNT
  ! }ImGuiStyleVar_;
  ! typedef enum {
  !     ImGuiButtonFlags_None = 0,
  !     ImGuiButtonFlags_MouseButtonLeft = 1 << 0,
  !     ImGuiButtonFlags_MouseButtonRight = 1 << 1,
  !     ImGuiButtonFlags_MouseButtonMiddle = 1 << 2,
  !     ImGuiButtonFlags_MouseButtonMask_ = ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle,
  !     ImGuiButtonFlags_MouseButtonDefault_ = ImGuiButtonFlags_MouseButtonLeft
  ! }ImGuiButtonFlags_;
  ! typedef enum {
  !     ImGuiColorEditFlags_None = 0,
  !     ImGuiColorEditFlags_NoAlpha = 1 << 1,
  !     ImGuiColorEditFlags_NoPicker = 1 << 2,
  !     ImGuiColorEditFlags_NoOptions = 1 << 3,
  !     ImGuiColorEditFlags_NoSmallPreview = 1 << 4,
  !     ImGuiColorEditFlags_NoInputs = 1 << 5,
  !     ImGuiColorEditFlags_NoTooltip = 1 << 6,
  !     ImGuiColorEditFlags_NoLabel = 1 << 7,
  !     ImGuiColorEditFlags_NoSidePreview = 1 << 8,
  !     ImGuiColorEditFlags_NoDragDrop = 1 << 9,
  !     ImGuiColorEditFlags_NoBorder = 1 << 10,
  !     ImGuiColorEditFlags_AlphaBar = 1 << 16,
  !     ImGuiColorEditFlags_AlphaPreview = 1 << 17,
  !     ImGuiColorEditFlags_AlphaPreviewHalf= 1 << 18,
  !     ImGuiColorEditFlags_HDR = 1 << 19,
  !     ImGuiColorEditFlags_DisplayRGB = 1 << 20,
  !     ImGuiColorEditFlags_DisplayHSV = 1 << 21,
  !     ImGuiColorEditFlags_DisplayHex = 1 << 22,
  !     ImGuiColorEditFlags_Uint8 = 1 << 23,
  !     ImGuiColorEditFlags_Float = 1 << 24,
  !     ImGuiColorEditFlags_PickerHueBar = 1 << 25,
  !     ImGuiColorEditFlags_PickerHueWheel = 1 << 26,
  !     ImGuiColorEditFlags_InputRGB = 1 << 27,
  !     ImGuiColorEditFlags_InputHSV = 1 << 28,
  !     ImGuiColorEditFlags_DefaultOptions_ = ImGuiColorEditFlags_Uint8 | ImGuiColorEditFlags_DisplayRGB | ImGuiColorEditFlags_InputRGB | ImGuiColorEditFlags_PickerHueBar,
  !     ImGuiColorEditFlags_DisplayMask_ = ImGuiColorEditFlags_DisplayRGB | ImGuiColorEditFlags_DisplayHSV | ImGuiColorEditFlags_DisplayHex,
  !     ImGuiColorEditFlags_DataTypeMask_ = ImGuiColorEditFlags_Uint8 | ImGuiColorEditFlags_Float,
  !     ImGuiColorEditFlags_PickerMask_ = ImGuiColorEditFlags_PickerHueWheel | ImGuiColorEditFlags_PickerHueBar,
  !     ImGuiColorEditFlags_InputMask_ = ImGuiColorEditFlags_InputRGB | ImGuiColorEditFlags_InputHSV
  ! }ImGuiColorEditFlags_;
  ! typedef enum {
  !     ImGuiSliderFlags_None = 0,
  !     ImGuiSliderFlags_AlwaysClamp = 1 << 4,
  !     ImGuiSliderFlags_Logarithmic = 1 << 5,
  !     ImGuiSliderFlags_NoRoundToFormat = 1 << 6,
  !     ImGuiSliderFlags_NoInput = 1 << 7,
  !     ImGuiSliderFlags_InvalidMask_ = 0x7000000F
  ! }ImGuiSliderFlags_;
  ! typedef enum {
  !     ImGuiMouseButton_Left = 0,
  !     ImGuiMouseButton_Right = 1,
  !     ImGuiMouseButton_Middle = 2,
  !     ImGuiMouseButton_COUNT = 5
  ! }ImGuiMouseButton_;
  ! typedef enum {
  !     ImGuiMouseCursor_None = -1,
  !     ImGuiMouseCursor_Arrow = 0,
  !     ImGuiMouseCursor_TextInput,
  !     ImGuiMouseCursor_ResizeAll,
  !     ImGuiMouseCursor_ResizeNS,
  !     ImGuiMouseCursor_ResizeEW,
  !     ImGuiMouseCursor_ResizeNESW,
  !     ImGuiMouseCursor_ResizeNWSE,
  !     ImGuiMouseCursor_Hand,
  !     ImGuiMouseCursor_NotAllowed,
  !     ImGuiMouseCursor_COUNT
  ! }ImGuiMouseCursor_;
  ! typedef enum {
  !     ImGuiCond_None = 0,
  !     ImGuiCond_Always = 1 << 0,
  !     ImGuiCond_Once = 1 << 1,
  !     ImGuiCond_FirstUseEver = 1 << 2,
  !     ImGuiCond_Appearing = 1 << 3
  ! }ImGuiCond_;
  ! typedef enum {
  !     ImDrawFlags_None = 0,
  !     ImDrawFlags_Closed = 1 << 0,
  !     ImDrawFlags_RoundCornersTopLeft = 1 << 4,
  !     ImDrawFlags_RoundCornersTopRight = 1 << 5,
  !     ImDrawFlags_RoundCornersBottomLeft = 1 << 6,
  !     ImDrawFlags_RoundCornersBottomRight = 1 << 7,
  !     ImDrawFlags_RoundCornersNone = 1 << 8,
  !     ImDrawFlags_RoundCornersTop = ImDrawFlags_RoundCornersTopLeft | ImDrawFlags_RoundCornersTopRight,
  !     ImDrawFlags_RoundCornersBottom = ImDrawFlags_RoundCornersBottomLeft | ImDrawFlags_RoundCornersBottomRight,
  !     ImDrawFlags_RoundCornersLeft = ImDrawFlags_RoundCornersBottomLeft | ImDrawFlags_RoundCornersTopLeft,
  !     ImDrawFlags_RoundCornersRight = ImDrawFlags_RoundCornersBottomRight | ImDrawFlags_RoundCornersTopRight,
  !     ImDrawFlags_RoundCornersAll = ImDrawFlags_RoundCornersTopLeft | ImDrawFlags_RoundCornersTopRight | ImDrawFlags_RoundCornersBottomLeft | ImDrawFlags_RoundCornersBottomRight,
  !     ImDrawFlags_RoundCornersDefault_ = ImDrawFlags_RoundCornersAll,
  !     ImDrawFlags_RoundCornersMask_ = ImDrawFlags_RoundCornersAll | ImDrawFlags_RoundCornersNone
  ! }ImDrawFlags_;
  ! typedef enum {
  !     ImDrawListFlags_None = 0,
  !     ImDrawListFlags_AntiAliasedLines = 1 << 0,
  !     ImDrawListFlags_AntiAliasedLinesUseTex = 1 << 1,
  !     ImDrawListFlags_AntiAliasedFill = 1 << 2,
  !     ImDrawListFlags_AllowVtxOffset = 1 << 3
  ! }ImDrawListFlags_;
  ! typedef enum {
  !     ImFontAtlasFlags_None = 0,
  !     ImFontAtlasFlags_NoPowerOfTwoHeight = 1 << 0,
  !     ImFontAtlasFlags_NoMouseCursors = 1 << 1,
  !     ImFontAtlasFlags_NoBakedLines = 1 << 2
  ! }ImFontAtlasFlags_;
  ! typedef enum {
  !     ImGuiViewportFlags_None = 0,
  !     ImGuiViewportFlags_IsPlatformWindow = 1 << 0,
  !     ImGuiViewportFlags_IsPlatformMonitor = 1 << 1,
  !     ImGuiViewportFlags_OwnedByApp = 1 << 2,
  !     ImGuiViewportFlags_NoDecoration = 1 << 3,
  !     ImGuiViewportFlags_NoTaskBarIcon = 1 << 4,
  !     ImGuiViewportFlags_NoFocusOnAppearing = 1 << 5,
  !     ImGuiViewportFlags_NoFocusOnClick = 1 << 6,
  !     ImGuiViewportFlags_NoInputs = 1 << 7,
  !     ImGuiViewportFlags_NoRendererClear = 1 << 8,
  !     ImGuiViewportFlags_TopMost = 1 << 9,
  !     ImGuiViewportFlags_Minimized = 1 << 10,
  !     ImGuiViewportFlags_NoAutoMerge = 1 << 11,
  !     ImGuiViewportFlags_CanHostOtherWindows = 1 << 12
  ! }ImGuiViewportFlags_;
  ! typedef enum {
  !     ImGuiItemFlags_None = 0,
  !     ImGuiItemFlags_NoTabStop = 1 << 0,
  !     ImGuiItemFlags_ButtonRepeat = 1 << 1,
  !     ImGuiItemFlags_Disabled = 1 << 2,
  !     ImGuiItemFlags_NoNav = 1 << 3,
  !     ImGuiItemFlags_NoNavDefaultFocus = 1 << 4,
  !     ImGuiItemFlags_SelectableDontClosePopup = 1 << 5,
  !     ImGuiItemFlags_MixedValue = 1 << 6,
  !     ImGuiItemFlags_ReadOnly = 1 << 7,
  !     ImGuiItemFlags_Inputable = 1 << 8
  ! }ImGuiItemFlags_;
  ! typedef enum {
  !     ImGuiItemStatusFlags_None = 0,
  !     ImGuiItemStatusFlags_HoveredRect = 1 << 0,
  !     ImGuiItemStatusFlags_HasDisplayRect = 1 << 1,
  !     ImGuiItemStatusFlags_Edited = 1 << 2,
  !     ImGuiItemStatusFlags_ToggledSelection = 1 << 3,
  !     ImGuiItemStatusFlags_ToggledOpen = 1 << 4,
  !     ImGuiItemStatusFlags_HasDeactivated = 1 << 5,
  !     ImGuiItemStatusFlags_Deactivated = 1 << 6,
  !     ImGuiItemStatusFlags_HoveredWindow = 1 << 7,
  !     ImGuiItemStatusFlags_FocusedByTabbing = 1 << 8
  ! }ImGuiItemStatusFlags_;
  ! typedef enum {
  !     ImGuiInputTextFlags_Multiline = 1 << 26,
  !     ImGuiInputTextFlags_NoMarkEdited = 1 << 27,
  !     ImGuiInputTextFlags_MergedItem = 1 << 28
  ! }ImGuiInputTextFlagsPrivate_;
  ! typedef enum {
  !     ImGuiButtonFlags_PressedOnClick = 1 << 4,
  !     ImGuiButtonFlags_PressedOnClickRelease = 1 << 5,
  !     ImGuiButtonFlags_PressedOnClickReleaseAnywhere = 1 << 6,
  !     ImGuiButtonFlags_PressedOnRelease = 1 << 7,
  !     ImGuiButtonFlags_PressedOnDoubleClick = 1 << 8,
  !     ImGuiButtonFlags_PressedOnDragDropHold = 1 << 9,
  !     ImGuiButtonFlags_Repeat = 1 << 10,
  !     ImGuiButtonFlags_FlattenChildren = 1 << 11,
  !     ImGuiButtonFlags_AllowItemOverlap = 1 << 12,
  !     ImGuiButtonFlags_DontClosePopups = 1 << 13,
  !     ImGuiButtonFlags_AlignTextBaseLine = 1 << 15,
  !     ImGuiButtonFlags_NoKeyModifiers = 1 << 16,
  !     ImGuiButtonFlags_NoHoldingActiveId = 1 << 17,
  !     ImGuiButtonFlags_NoNavFocus = 1 << 18,
  !     ImGuiButtonFlags_NoHoveredOnFocus = 1 << 19,
  !     ImGuiButtonFlags_PressedOnMask_ = ImGuiButtonFlags_PressedOnClick | ImGuiButtonFlags_PressedOnClickRelease | ImGuiButtonFlags_PressedOnClickReleaseAnywhere | ImGuiButtonFlags_PressedOnRelease | ImGuiButtonFlags_PressedOnDoubleClick | ImGuiButtonFlags_PressedOnDragDropHold,
  !     ImGuiButtonFlags_PressedOnDefault_ = ImGuiButtonFlags_PressedOnClickRelease
  ! }ImGuiButtonFlagsPrivate_;
  ! typedef enum {
  !     ImGuiComboFlags_CustomPreview = 1 << 20
  ! }ImGuiComboFlagsPrivate_;
  ! typedef enum {
  !     ImGuiSliderFlags_Vertical = 1 << 20,
  !     ImGuiSliderFlags_ReadOnly = 1 << 21
  ! }ImGuiSliderFlagsPrivate_;
  ! typedef enum {
  !     ImGuiSelectableFlags_NoHoldingActiveID = 1 << 20,
  !     ImGuiSelectableFlags_SelectOnNav = 1 << 21,
  !     ImGuiSelectableFlags_SelectOnClick = 1 << 22,
  !     ImGuiSelectableFlags_SelectOnRelease = 1 << 23,
  !     ImGuiSelectableFlags_SpanAvailWidth = 1 << 24,
  !     ImGuiSelectableFlags_DrawHoveredWhenHeld = 1 << 25,
  !     ImGuiSelectableFlags_SetNavIdOnHover = 1 << 26,
  !     ImGuiSelectableFlags_NoPadWithHalfSpacing = 1 << 27
  ! }ImGuiSelectableFlagsPrivate_;
  ! typedef enum {
  !     ImGuiTreeNodeFlags_ClipLabelForTrailingButton = 1 << 20
  ! }ImGuiTreeNodeFlagsPrivate_;
  ! typedef enum {
  !     ImGuiSeparatorFlags_None = 0,
  !     ImGuiSeparatorFlags_Horizontal = 1 << 0,
  !     ImGuiSeparatorFlags_Vertical = 1 << 1,
  !     ImGuiSeparatorFlags_SpanAllColumns = 1 << 2
  ! }ImGuiSeparatorFlags_;
  ! typedef enum {
  !     ImGuiTextFlags_None = 0,
  !     ImGuiTextFlags_NoWidthForLargeClippedText = 1 << 0
  ! }ImGuiTextFlags_;
  ! typedef enum {
  !     ImGuiTooltipFlags_None = 0,
  !     ImGuiTooltipFlags_OverridePreviousTooltip = 1 << 0
  ! }ImGuiTooltipFlags_;
  ! typedef enum {
  !     ImGuiLayoutType_Horizontal = 0,
  !     ImGuiLayoutType_Vertical = 1
  ! }ImGuiLayoutType_;
  ! typedef enum {
  !     ImGuiLogType_None = 0,
  !     ImGuiLogType_TTY,
  !     ImGuiLogType_File,
  !     ImGuiLogType_Buffer,
  !     ImGuiLogType_Clipboard
  ! }ImGuiLogType;
  ! typedef enum {
  !     ImGuiAxis_None = -1,
  !     ImGuiAxis_X = 0,
  !     ImGuiAxis_Y = 1
  ! }ImGuiAxis;
  ! typedef enum {
  !     ImGuiPlotType_Lines,
  !     ImGuiPlotType_Histogram
  ! }ImGuiPlotType;
  ! typedef enum {
  !     ImGuiPopupPositionPolicy_Default,
  !     ImGuiPopupPositionPolicy_ComboBox,
  !     ImGuiPopupPositionPolicy_Tooltip
  ! }ImGuiPopupPositionPolicy;
  ! typedef enum {
  !     ImGuiDataType_String = ImGuiDataType_COUNT + 1,
  !     ImGuiDataType_Pointer,
  !     ImGuiDataType_ID
  ! }ImGuiDataTypePrivate_;
  ! typedef enum {
  !     ImGuiNextWindowDataFlags_None = 0,
  !     ImGuiNextWindowDataFlags_HasPos = 1 << 0,
  !     ImGuiNextWindowDataFlags_HasSize = 1 << 1,
  !     ImGuiNextWindowDataFlags_HasContentSize = 1 << 2,
  !     ImGuiNextWindowDataFlags_HasCollapsed = 1 << 3,
  !     ImGuiNextWindowDataFlags_HasSizeConstraint = 1 << 4,
  !     ImGuiNextWindowDataFlags_HasFocus = 1 << 5,
  !     ImGuiNextWindowDataFlags_HasBgAlpha = 1 << 6,
  !     ImGuiNextWindowDataFlags_HasScroll = 1 << 7,
  !     ImGuiNextWindowDataFlags_HasViewport = 1 << 8,
  !     ImGuiNextWindowDataFlags_HasDock = 1 << 9,
  !     ImGuiNextWindowDataFlags_HasWindowClass = 1 << 10
  ! }ImGuiNextWindowDataFlags_;
  ! typedef enum {
  !     ImGuiNextItemDataFlags_None = 0,
  !     ImGuiNextItemDataFlags_HasWidth = 1 << 0,
  !     ImGuiNextItemDataFlags_HasOpen = 1 << 1
  ! }ImGuiNextItemDataFlags_;
  ! typedef enum {
  !     ImGuiKey_LegacyNativeKey_BEGIN = 0,
  !     ImGuiKey_LegacyNativeKey_END = 512,
  !     ImGuiKey_Gamepad_BEGIN = ImGuiKey_GamepadStart,
  !     ImGuiKey_Gamepad_END = ImGuiKey_GamepadRStickRight + 1
  ! }ImGuiKeyPrivate_;
  ! typedef enum {
  !     ImGuiInputEventType_None = 0,
  !     ImGuiInputEventType_MousePos,
  !     ImGuiInputEventType_MouseWheel,
  !     ImGuiInputEventType_MouseButton,
  !     ImGuiInputEventType_MouseViewport,
  !     ImGuiInputEventType_Key,
  !     ImGuiInputEventType_Text,
  !     ImGuiInputEventType_Focus,
  !     ImGuiInputEventType_COUNT
  ! }ImGuiInputEventType;
  ! typedef enum {
  !     ImGuiInputSource_None = 0,
  !     ImGuiInputSource_Mouse,
  !     ImGuiInputSource_Keyboard,
  !     ImGuiInputSource_Gamepad,
  !     ImGuiInputSource_Clipboard,
  !     ImGuiInputSource_Nav,
  !     ImGuiInputSource_COUNT
  ! }ImGuiInputSource;
  ! typedef enum {
  !     ImGuiNavReadMode_Down,
  !     ImGuiNavReadMode_Pressed,
  !     ImGuiNavReadMode_Released,
  !     ImGuiNavReadMode_Repeat,
  !     ImGuiNavReadMode_RepeatSlow,
  !     ImGuiNavReadMode_RepeatFast
  ! }ImGuiNavReadMode;
  ! typedef enum {
  !     ImGuiActivateFlags_None = 0,
  !     ImGuiActivateFlags_PreferInput = 1 << 0,
  !     ImGuiActivateFlags_PreferTweak = 1 << 1,
  !     ImGuiActivateFlags_TryToPreserveState = 1 << 2
  ! }ImGuiActivateFlags_;
  ! typedef enum {
  !     ImGuiScrollFlags_None = 0,
  !     ImGuiScrollFlags_KeepVisibleEdgeX = 1 << 0,
  !     ImGuiScrollFlags_KeepVisibleEdgeY = 1 << 1,
  !     ImGuiScrollFlags_KeepVisibleCenterX = 1 << 2,
  !     ImGuiScrollFlags_KeepVisibleCenterY = 1 << 3,
  !     ImGuiScrollFlags_AlwaysCenterX = 1 << 4,
  !     ImGuiScrollFlags_AlwaysCenterY = 1 << 5,
  !     ImGuiScrollFlags_NoScrollParent = 1 << 6,
  !     ImGuiScrollFlags_MaskX_ = ImGuiScrollFlags_KeepVisibleEdgeX | ImGuiScrollFlags_KeepVisibleCenterX | ImGuiScrollFlags_AlwaysCenterX,
  !     ImGuiScrollFlags_MaskY_ = ImGuiScrollFlags_KeepVisibleEdgeY | ImGuiScrollFlags_KeepVisibleCenterY | ImGuiScrollFlags_AlwaysCenterY
  ! }ImGuiScrollFlags_;
  ! typedef enum {
  !     ImGuiNavHighlightFlags_None = 0,
  !     ImGuiNavHighlightFlags_TypeDefault = 1 << 0,
  !     ImGuiNavHighlightFlags_TypeThin = 1 << 1,
  !     ImGuiNavHighlightFlags_AlwaysDraw = 1 << 2,
  !     ImGuiNavHighlightFlags_NoRounding = 1 << 3
  ! }ImGuiNavHighlightFlags_;
  ! typedef enum {
  !     ImGuiNavDirSourceFlags_None = 0,
  !     ImGuiNavDirSourceFlags_RawKeyboard = 1 << 0,
  !     ImGuiNavDirSourceFlags_Keyboard = 1 << 1,
  !     ImGuiNavDirSourceFlags_PadDPad = 1 << 2,
  !     ImGuiNavDirSourceFlags_PadLStick = 1 << 3
  ! }ImGuiNavDirSourceFlags_;
  ! typedef enum {
  !     ImGuiNavMoveFlags_None = 0,
  !     ImGuiNavMoveFlags_LoopX = 1 << 0,
  !     ImGuiNavMoveFlags_LoopY = 1 << 1,
  !     ImGuiNavMoveFlags_WrapX = 1 << 2,
  !     ImGuiNavMoveFlags_WrapY = 1 << 3,
  !     ImGuiNavMoveFlags_AllowCurrentNavId = 1 << 4,
  !     ImGuiNavMoveFlags_AlsoScoreVisibleSet = 1 << 5,
  !     ImGuiNavMoveFlags_ScrollToEdgeY = 1 << 6,
  !     ImGuiNavMoveFlags_Forwarded = 1 << 7,
  !     ImGuiNavMoveFlags_DebugNoResult = 1 << 8,
  !     ImGuiNavMoveFlags_FocusApi = 1 << 9,
  !     ImGuiNavMoveFlags_Tabbing = 1 << 10,
  !     ImGuiNavMoveFlags_Activate = 1 << 11,
  !     ImGuiNavMoveFlags_DontSetNavHighlight = 1 << 12
  ! }ImGuiNavMoveFlags_;
  ! typedef enum {
  !     ImGuiNavLayer_Main = 0,
  !     ImGuiNavLayer_Menu = 1,
  !     ImGuiNavLayer_COUNT
  ! }ImGuiNavLayer;
  ! typedef enum {
  !     ImGuiOldColumnFlags_None = 0,
  !     ImGuiOldColumnFlags_NoBorder = 1 << 0,
  !     ImGuiOldColumnFlags_NoResize = 1 << 1,
  !     ImGuiOldColumnFlags_NoPreserveWidths = 1 << 2,
  !     ImGuiOldColumnFlags_NoForceWithinWindow = 1 << 3,
  !     ImGuiOldColumnFlags_GrowParentContentsSize = 1 << 4
  ! }ImGuiOldColumnFlags_;
  ! typedef enum {
  !     ImGuiDockNodeFlags_DockSpace = 1 << 10,
  !     ImGuiDockNodeFlags_CentralNode = 1 << 11,
  !     ImGuiDockNodeFlags_NoTabBar = 1 << 12,
  !     ImGuiDockNodeFlags_HiddenTabBar = 1 << 13,
  !     ImGuiDockNodeFlags_NoWindowMenuButton = 1 << 14,
  !     ImGuiDockNodeFlags_NoCloseButton = 1 << 15,
  !     ImGuiDockNodeFlags_NoDocking = 1 << 16,
  !     ImGuiDockNodeFlags_NoDockingSplitMe = 1 << 17,
  !     ImGuiDockNodeFlags_NoDockingSplitOther = 1 << 18,
  !     ImGuiDockNodeFlags_NoDockingOverMe = 1 << 19,
  !     ImGuiDockNodeFlags_NoDockingOverOther = 1 << 20,
  !     ImGuiDockNodeFlags_NoDockingOverEmpty = 1 << 21,
  !     ImGuiDockNodeFlags_NoResizeX = 1 << 22,
  !     ImGuiDockNodeFlags_NoResizeY = 1 << 23,
  !     ImGuiDockNodeFlags_SharedFlagsInheritMask_ = ~0,
  !     ImGuiDockNodeFlags_NoResizeFlagsMask_ = ImGuiDockNodeFlags_NoResize | ImGuiDockNodeFlags_NoResizeX | ImGuiDockNodeFlags_NoResizeY,
  !     ImGuiDockNodeFlags_LocalFlagsMask_ = ImGuiDockNodeFlags_NoSplit | ImGuiDockNodeFlags_NoResizeFlagsMask_ | ImGuiDockNodeFlags_AutoHideTabBar | ImGuiDockNodeFlags_DockSpace | ImGuiDockNodeFlags_CentralNode | ImGuiDockNodeFlags_NoTabBar | ImGuiDockNodeFlags_HiddenTabBar | ImGuiDockNodeFlags_NoWindowMenuButton | ImGuiDockNodeFlags_NoCloseButton | ImGuiDockNodeFlags_NoDocking,
  !     ImGuiDockNodeFlags_LocalFlagsTransferMask_ = ImGuiDockNodeFlags_LocalFlagsMask_ & ~ImGuiDockNodeFlags_DockSpace,
  !     ImGuiDockNodeFlags_SavedFlagsMask_ = ImGuiDockNodeFlags_NoResizeFlagsMask_ | ImGuiDockNodeFlags_DockSpace | ImGuiDockNodeFlags_CentralNode | ImGuiDockNodeFlags_NoTabBar | ImGuiDockNodeFlags_HiddenTabBar | ImGuiDockNodeFlags_NoWindowMenuButton | ImGuiDockNodeFlags_NoCloseButton | ImGuiDockNodeFlags_NoDocking
  ! }ImGuiDockNodeFlagsPrivate_;
  ! typedef enum {
  !     ImGuiDataAuthority_Auto,
  !     ImGuiDataAuthority_DockNode,
  !     ImGuiDataAuthority_Window
  ! }ImGuiDataAuthority_;
  ! typedef enum {
  !     ImGuiDockNodeState_Unknown,
  !     ImGuiDockNodeState_HostWindowHiddenBecauseSingleWindow,
  !     ImGuiDockNodeState_HostWindowHiddenBecauseWindowsAreResizing,
  !     ImGuiDockNodeState_HostWindowVisible
  ! }ImGuiDockNodeState;
  ! typedef enum {
  !     ImGuiWindowDockStyleCol_Text,
  !     ImGuiWindowDockStyleCol_Tab,
  !     ImGuiWindowDockStyleCol_TabHovered,
  !     ImGuiWindowDockStyleCol_TabActive,
  !     ImGuiWindowDockStyleCol_TabUnfocused,
  !     ImGuiWindowDockStyleCol_TabUnfocusedActive,
  !     ImGuiWindowDockStyleCol_COUNT
  ! }ImGuiWindowDockStyleCol;
  ! typedef enum {
  !     ImGuiDebugLogFlags_None = 0,
  !     ImGuiDebugLogFlags_EventActiveId = 1 << 0,
  !     ImGuiDebugLogFlags_EventFocus = 1 << 1,
  !     ImGuiDebugLogFlags_EventPopup = 1 << 2,
  !     ImGuiDebugLogFlags_EventNav = 1 << 3,
  !     ImGuiDebugLogFlags_EventIO = 1 << 4,
  !     ImGuiDebugLogFlags_EventDocking = 1 << 5,
  !     ImGuiDebugLogFlags_EventViewport = 1 << 6,
  !     ImGuiDebugLogFlags_EventMask_ = ImGuiDebugLogFlags_EventActiveId | ImGuiDebugLogFlags_EventFocus | ImGuiDebugLogFlags_EventPopup | ImGuiDebugLogFlags_EventNav | ImGuiDebugLogFlags_EventIO | ImGuiDebugLogFlags_EventDocking | ImGuiDebugLogFlags_EventViewport,
  !     ImGuiDebugLogFlags_OutputToTTY = 1 << 10
  ! }ImGuiDebugLogFlags_;
  ! typedef enum { ImGuiContextHookType_NewFramePre, ImGuiContextHookType_NewFramePost, ImGuiContextHookType_EndFramePre, ImGuiContextHookType_EndFramePost, ImGuiContextHookType_RenderPre, ImGuiContextHookType_RenderPost, ImGuiContextHookType_Shutdown, ImGuiContextHookType_PendingRemoval_ }ImGuiContextHookType;
  ! typedef enum {
  !     ImGuiTabBarFlags_DockNode = 1 << 20,
  !     ImGuiTabBarFlags_IsFocused = 1 << 21,
  !     ImGuiTabBarFlags_SaveSettings = 1 << 22
  ! }ImGuiTabBarFlagsPrivate_;
  ! typedef enum {
  !     ImGuiTabItemFlags_SectionMask_ = ImGuiTabItemFlags_Leading | ImGuiTabItemFlags_Trailing,
  !     ImGuiTabItemFlags_NoCloseButton = 1 << 20,
  !     ImGuiTabItemFlags_Button = 1 << 21,
  !     ImGuiTabItemFlags_Unsorted = 1 << 22,
  !     ImGuiTabItemFlags_Preview = 1 << 23
  ! }ImGuiTabItemFlagsPrivate_;
  ! typedef enum {
  !     ImGuiFreeTypeBuilderFlags_NoHinting = 1 << 0,
  !     ImGuiFreeTypeBuilderFlags_NoAutoHint = 1 << 1,
  !     ImGuiFreeTypeBuilderFlags_ForceAutoHint = 1 << 2,
  !     ImGuiFreeTypeBuilderFlags_LightHinting = 1 << 3,
  !     ImGuiFreeTypeBuilderFlags_MonoHinting = 1 << 4,
  !     ImGuiFreeTypeBuilderFlags_Bold = 1 << 5,
  !     ImGuiFreeTypeBuilderFlags_Oblique = 1 << 6,
  !     ImGuiFreeTypeBuilderFlags_Monochrome = 1 << 7,
  !     ImGuiFreeTypeBuilderFlags_LoadColor = 1 << 8,
  !     ImGuiFreeTypeBuilderFlags_Bitmap = 1 << 9
  ! }ImGuiFreeTypeBuilderFlags;
