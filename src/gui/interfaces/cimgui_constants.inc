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
  integer(c_int), parameter :: ImGuiKey_NamedKey_COUNT = 133
  integer(c_int), parameter :: ImGuiKey_NamedKey_COUNT_31shift5 = 5 ! (ImGuiKey_NamedKey_COUNT+31)>>5
  integer(c_int), parameter :: ImGuiKey_KeysData_SIZE = 645
  ! enum ImGuiNavInput_
  integer(c_int), parameter :: ImGuiNavInput_COUNT = 20
  ! enum ImGuiStyle
  integer(c_int), parameter :: ImGuiCol_COUNT = 55
  !! end of sizes !!

  ! enum ImGuiWindowFlags_
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_None") :: ImGuiWindowFlags_None
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoTitleBar") :: ImGuiWindowFlags_NoTitleBar
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoResize") :: ImGuiWindowFlags_NoResize
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoMove") :: ImGuiWindowFlags_NoMove
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoScrollbar") :: ImGuiWindowFlags_NoScrollbar
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoScrollWithMouse") :: ImGuiWindowFlags_NoScrollWithMouse
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoCollapse") :: ImGuiWindowFlags_NoCollapse
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_AlwaysAutoResize") :: ImGuiWindowFlags_AlwaysAutoResize
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoBackground") :: ImGuiWindowFlags_NoBackground
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoSavedSettings") :: ImGuiWindowFlags_NoSavedSettings
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoMouseInputs") :: ImGuiWindowFlags_NoMouseInputs
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_MenuBar") :: ImGuiWindowFlags_MenuBar
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_HorizontalScrollbar") :: ImGuiWindowFlags_HorizontalScrollbar
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoFocusOnAppearing") :: ImGuiWindowFlags_NoFocusOnAppearing
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoBringToFrontOnFocus") :: ImGuiWindowFlags_NoBringToFrontOnFocus
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_AlwaysVerticalScrollbar") :: ImGuiWindowFlags_AlwaysVerticalScrollbar
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_AlwaysHorizontalScrollbar") :: ImGuiWindowFlags_AlwaysHorizontalScrollbar
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_AlwaysUseWindowPadding") :: ImGuiWindowFlags_AlwaysUseWindowPadding
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoNavInputs") :: ImGuiWindowFlags_NoNavInputs
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoNavFocus") :: ImGuiWindowFlags_NoNavFocus
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_UnsavedDocument") :: ImGuiWindowFlags_UnsavedDocument
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoDocking") :: ImGuiWindowFlags_NoDocking
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoNav") :: ImGuiWindowFlags_NoNav
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoDecoration") :: ImGuiWindowFlags_NoDecoration
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NoInputs") :: ImGuiWindowFlags_NoInputs
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_NavFlattened") :: ImGuiWindowFlags_NavFlattened
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_ChildWindow") :: ImGuiWindowFlags_ChildWindow
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_Tooltip") :: ImGuiWindowFlags_Tooltip
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_Popup") :: ImGuiWindowFlags_Popup
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_Modal") :: ImGuiWindowFlags_Modal
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_ChildMenu") :: ImGuiWindowFlags_ChildMenu
  integer(c_int), bind(C,name="const_ImGuiWindowFlags_DockNodeHost") :: ImGuiWindowFlags_DockNodeHost
  ! enum ImGuiInputTextFlags_
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_None") :: ImGuiInputTextFlags_None
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CharsDecimal") :: ImGuiInputTextFlags_CharsDecimal
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CharsHexadecimal") :: ImGuiInputTextFlags_CharsHexadecimal
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CharsUppercase") :: ImGuiInputTextFlags_CharsUppercase
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CharsNoBlank") :: ImGuiInputTextFlags_CharsNoBlank
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_AutoSelectAll") :: ImGuiInputTextFlags_AutoSelectAll
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_EnterReturnsTrue") :: ImGuiInputTextFlags_EnterReturnsTrue
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CallbackCompletion") :: ImGuiInputTextFlags_CallbackCompletion
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CallbackHistory") :: ImGuiInputTextFlags_CallbackHistory
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CallbackAlways") :: ImGuiInputTextFlags_CallbackAlways
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CallbackCharFilter") :: ImGuiInputTextFlags_CallbackCharFilter
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_AllowTabInput") :: ImGuiInputTextFlags_AllowTabInput
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CtrlEnterForNewLine") :: ImGuiInputTextFlags_CtrlEnterForNewLine
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_NoHorizontalScroll") :: ImGuiInputTextFlags_NoHorizontalScroll
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_AlwaysOverwrite") :: ImGuiInputTextFlags_AlwaysOverwrite
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_ReadOnly") :: ImGuiInputTextFlags_ReadOnly
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_Password") :: ImGuiInputTextFlags_Password
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_NoUndoRedo") :: ImGuiInputTextFlags_NoUndoRedo
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CharsScientific") :: ImGuiInputTextFlags_CharsScientific
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CallbackResize") :: ImGuiInputTextFlags_CallbackResize
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_CallbackEdit") :: ImGuiInputTextFlags_CallbackEdit
  ! enum ImGuiTreeNodeFlags_;
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_None") :: ImGuiTreeNodeFlags_None
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_Selected") :: ImGuiTreeNodeFlags_Selected
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_Framed") :: ImGuiTreeNodeFlags_Framed
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_AllowItemOverlap") :: ImGuiTreeNodeFlags_AllowItemOverlap
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_NoTreePushOnOpen") :: ImGuiTreeNodeFlags_NoTreePushOnOpen
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_NoAutoOpenOnLog") :: ImGuiTreeNodeFlags_NoAutoOpenOnLog
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_DefaultOpen") :: ImGuiTreeNodeFlags_DefaultOpen
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_OpenOnDoubleClick") :: ImGuiTreeNodeFlags_OpenOnDoubleClick
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_OpenOnArrow") :: ImGuiTreeNodeFlags_OpenOnArrow
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_Leaf") :: ImGuiTreeNodeFlags_Leaf
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_Bullet") :: ImGuiTreeNodeFlags_Bullet
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_FramePadding") :: ImGuiTreeNodeFlags_FramePadding
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_SpanAvailWidth") :: ImGuiTreeNodeFlags_SpanAvailWidth
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_SpanFullWidth") :: ImGuiTreeNodeFlags_SpanFullWidth
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_NavLeftJumpsBackHere") :: ImGuiTreeNodeFlags_NavLeftJumpsBackHere
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_CollapsingHeader") :: ImGuiTreeNodeFlags_CollapsingHeader
  ! enum ImGuiPopupFlags_
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_None") :: ImGuiPopupFlags_None
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_MouseButtonLeft") :: ImGuiPopupFlags_MouseButtonLeft
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_MouseButtonRight") :: ImGuiPopupFlags_MouseButtonRight
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_MouseButtonMiddle") :: ImGuiPopupFlags_MouseButtonMiddle
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_MouseButtonMask_") :: ImGuiPopupFlags_MouseButtonMask_
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_MouseButtonDefault_") :: ImGuiPopupFlags_MouseButtonDefault_
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_NoOpenOverExistingPopup") :: ImGuiPopupFlags_NoOpenOverExistingPopup
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_NoOpenOverItems") :: ImGuiPopupFlags_NoOpenOverItems
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_AnyPopupId") :: ImGuiPopupFlags_AnyPopupId
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_AnyPopupLevel") :: ImGuiPopupFlags_AnyPopupLevel
  integer(c_int), bind(C,name="const_ImGuiPopupFlags_AnyPopup") :: ImGuiPopupFlags_AnyPopup
  ! enum ImGuiSelectableFlags_
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_None") :: ImGuiSelectableFlags_None
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_DontClosePopups") :: ImGuiSelectableFlags_DontClosePopups
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_SpanAllColumns") :: ImGuiSelectableFlags_SpanAllColumns
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_AllowDoubleClick") :: ImGuiSelectableFlags_AllowDoubleClick
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_Disabled") :: ImGuiSelectableFlags_Disabled
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_AllowItemOverlap") :: ImGuiSelectableFlags_AllowItemOverlap
  ! enum ImGuiComboFlags_
  integer(c_int), bind(C,name="const_ImGuiComboFlags_None") :: ImGuiComboFlags_None
  integer(c_int), bind(C,name="const_ImGuiComboFlags_PopupAlignLeft") :: ImGuiComboFlags_PopupAlignLeft
  integer(c_int), bind(C,name="const_ImGuiComboFlags_HeightSmall") :: ImGuiComboFlags_HeightSmall
  integer(c_int), bind(C,name="const_ImGuiComboFlags_HeightRegular") :: ImGuiComboFlags_HeightRegular
  integer(c_int), bind(C,name="const_ImGuiComboFlags_HeightLarge") :: ImGuiComboFlags_HeightLarge
  integer(c_int), bind(C,name="const_ImGuiComboFlags_HeightLargest") :: ImGuiComboFlags_HeightLargest
  integer(c_int), bind(C,name="const_ImGuiComboFlags_NoArrowButton") :: ImGuiComboFlags_NoArrowButton
  integer(c_int), bind(C,name="const_ImGuiComboFlags_NoPreview") :: ImGuiComboFlags_NoPreview
  integer(c_int), bind(C,name="const_ImGuiComboFlags_HeightMask_") :: ImGuiComboFlags_HeightMask_
  ! enum ImGuiTabBarFlags_
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_None") :: ImGuiTabBarFlags_None
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_Reorderable") :: ImGuiTabBarFlags_Reorderable
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_AutoSelectNewTabs") :: ImGuiTabBarFlags_AutoSelectNewTabs
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_TabListPopupButton") :: ImGuiTabBarFlags_TabListPopupButton
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_NoCloseWithMiddleMouseButton") :: &
     ImGuiTabBarFlags_NoCloseWithMiddleMouseButton
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_NoTabListScrollingButtons") :: ImGuiTabBarFlags_NoTabListScrollingButtons
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_NoTooltip") :: ImGuiTabBarFlags_NoTooltip
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_FittingPolicyResizeDown") :: ImGuiTabBarFlags_FittingPolicyResizeDown
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_FittingPolicyScroll") :: ImGuiTabBarFlags_FittingPolicyScroll
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_FittingPolicyMask_") :: ImGuiTabBarFlags_FittingPolicyMask_
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_FittingPolicyDefault_") :: ImGuiTabBarFlags_FittingPolicyDefault_
  ! enum ImGuiTabItemFlags_
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_None") :: ImGuiTabItemFlags_None
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_UnsavedDocument") :: ImGuiTabItemFlags_UnsavedDocument
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_SetSelected") :: ImGuiTabItemFlags_SetSelected
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_NoCloseWithMiddleMouseButton") :: &
     ImGuiTabItemFlags_NoCloseWithMiddleMouseButton
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_NoPushId") :: ImGuiTabItemFlags_NoPushId
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_NoTooltip") :: ImGuiTabItemFlags_NoTooltip
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_NoReorder") :: ImGuiTabItemFlags_NoReorder
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_Leading") :: ImGuiTabItemFlags_Leading
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_Trailing") :: ImGuiTabItemFlags_Trailing
  ! enum ImGuiTableFlags_
  integer(c_int), bind(C,name="const_ImGuiTableFlags_None") :: ImGuiTableFlags_None
  integer(c_int), bind(C,name="const_ImGuiTableFlags_Resizable") :: ImGuiTableFlags_Resizable
  integer(c_int), bind(C,name="const_ImGuiTableFlags_Reorderable") :: ImGuiTableFlags_Reorderable
  integer(c_int), bind(C,name="const_ImGuiTableFlags_Hideable") :: ImGuiTableFlags_Hideable
  integer(c_int), bind(C,name="const_ImGuiTableFlags_Sortable") :: ImGuiTableFlags_Sortable
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoSavedSettings") :: ImGuiTableFlags_NoSavedSettings
  integer(c_int), bind(C,name="const_ImGuiTableFlags_ContextMenuInBody") :: ImGuiTableFlags_ContextMenuInBody
  integer(c_int), bind(C,name="const_ImGuiTableFlags_RowBg") :: ImGuiTableFlags_RowBg
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersInnerH") :: ImGuiTableFlags_BordersInnerH
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersOuterH") :: ImGuiTableFlags_BordersOuterH
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersInnerV") :: ImGuiTableFlags_BordersInnerV
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersOuterV") :: ImGuiTableFlags_BordersOuterV
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersH") :: ImGuiTableFlags_BordersH
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersV") :: ImGuiTableFlags_BordersV
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersInner") :: ImGuiTableFlags_BordersInner
  integer(c_int), bind(C,name="const_ImGuiTableFlags_BordersOuter") :: ImGuiTableFlags_BordersOuter
  integer(c_int), bind(C,name="const_ImGuiTableFlags_Borders") :: ImGuiTableFlags_Borders
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoBordersInBody") :: ImGuiTableFlags_NoBordersInBody
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoBordersInBodyUntilResize") :: ImGuiTableFlags_NoBordersInBodyUntilResize
  integer(c_int), bind(C,name="const_ImGuiTableFlags_SizingFixedFit") :: ImGuiTableFlags_SizingFixedFit
  integer(c_int), bind(C,name="const_ImGuiTableFlags_SizingFixedSame") :: ImGuiTableFlags_SizingFixedSame
  integer(c_int), bind(C,name="const_ImGuiTableFlags_SizingStretchProp") :: ImGuiTableFlags_SizingStretchProp
  integer(c_int), bind(C,name="const_ImGuiTableFlags_SizingStretchSame") :: ImGuiTableFlags_SizingStretchSame
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoHostExtendX") :: ImGuiTableFlags_NoHostExtendX
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoHostExtendY") :: ImGuiTableFlags_NoHostExtendY
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoKeepColumnsVisible") :: ImGuiTableFlags_NoKeepColumnsVisible
  integer(c_int), bind(C,name="const_ImGuiTableFlags_PreciseWidths") :: ImGuiTableFlags_PreciseWidths
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoClip") :: ImGuiTableFlags_NoClip
  integer(c_int), bind(C,name="const_ImGuiTableFlags_PadOuterX") :: ImGuiTableFlags_PadOuterX
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoPadOuterX") :: ImGuiTableFlags_NoPadOuterX
  integer(c_int), bind(C,name="const_ImGuiTableFlags_NoPadInnerX") :: ImGuiTableFlags_NoPadInnerX
  integer(c_int), bind(C,name="const_ImGuiTableFlags_ScrollX") :: ImGuiTableFlags_ScrollX
  integer(c_int), bind(C,name="const_ImGuiTableFlags_ScrollY") :: ImGuiTableFlags_ScrollY
  integer(c_int), bind(C,name="const_ImGuiTableFlags_SortMulti") :: ImGuiTableFlags_SortMulti
  integer(c_int), bind(C,name="const_ImGuiTableFlags_SortTristate") :: ImGuiTableFlags_SortTristate
  integer(c_int), bind(C,name="const_ImGuiTableFlags_SizingMask_") :: ImGuiTableFlags_SizingMask_
  ! enum ImGuiTableColumnFlags_
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_None") :: ImGuiTableColumnFlags_None
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_Disabled") :: ImGuiTableColumnFlags_Disabled
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_DefaultHide") :: ImGuiTableColumnFlags_DefaultHide
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_DefaultSort") :: ImGuiTableColumnFlags_DefaultSort
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_WidthStretch") :: ImGuiTableColumnFlags_WidthStretch
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_WidthFixed") :: ImGuiTableColumnFlags_WidthFixed
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoResize") :: ImGuiTableColumnFlags_NoResize
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoReorder") :: ImGuiTableColumnFlags_NoReorder
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoHide") :: ImGuiTableColumnFlags_NoHide
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoClip") :: ImGuiTableColumnFlags_NoClip
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoSort") :: ImGuiTableColumnFlags_NoSort
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoSortAscending") :: ImGuiTableColumnFlags_NoSortAscending
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoSortDescending") :: ImGuiTableColumnFlags_NoSortDescending
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoHeaderLabel") :: ImGuiTableColumnFlags_NoHeaderLabel
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoHeaderWidth") :: ImGuiTableColumnFlags_NoHeaderWidth
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_PreferSortAscending") :: ImGuiTableColumnFlags_PreferSortAscending
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_PreferSortDescending") :: ImGuiTableColumnFlags_PreferSortDescending
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_IndentEnable") :: ImGuiTableColumnFlags_IndentEnable
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_IndentDisable") :: ImGuiTableColumnFlags_IndentDisable
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_IsEnabled") :: ImGuiTableColumnFlags_IsEnabled
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_IsVisible") :: ImGuiTableColumnFlags_IsVisible
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_IsSorted") :: ImGuiTableColumnFlags_IsSorted
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_IsHovered") :: ImGuiTableColumnFlags_IsHovered
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_WidthMask_") :: ImGuiTableColumnFlags_WidthMask_
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_IndentMask_") :: ImGuiTableColumnFlags_IndentMask_
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_StatusMask_") :: ImGuiTableColumnFlags_StatusMask_
  integer(c_int), bind(C,name="const_ImGuiTableColumnFlags_NoDirectResize_") :: ImGuiTableColumnFlags_NoDirectResize_
  ! enum ImGuiTableRowFlags_
  integer(c_int), bind(C,name="const_ImGuiTableRowFlags_None") :: ImGuiTableRowFlags_None
  integer(c_int), bind(C,name="const_ImGuiTableRowFlags_Headers") :: ImGuiTableRowFlags_Headers
  ! enum ImGuiTableBgTarget_
  integer(c_int), bind(C,name="const_ImGuiTableBgTarget_None") :: ImGuiTableBgTarget_None
  integer(c_int), bind(C,name="const_ImGuiTableBgTarget_RowBg0") :: ImGuiTableBgTarget_RowBg0
  integer(c_int), bind(C,name="const_ImGuiTableBgTarget_RowBg1") :: ImGuiTableBgTarget_RowBg1
  integer(c_int), bind(C,name="const_ImGuiTableBgTarget_CellBg") :: ImGuiTableBgTarget_CellBg
  ! enum ImGuiFocusedFlags_
  integer(c_int), bind(C,name="const_ImGuiFocusedFlags_None") :: ImGuiFocusedFlags_None
  integer(c_int), bind(C,name="const_ImGuiFocusedFlags_ChildWindows") :: ImGuiFocusedFlags_ChildWindows
  integer(c_int), bind(C,name="const_ImGuiFocusedFlags_RootWindow") :: ImGuiFocusedFlags_RootWindow
  integer(c_int), bind(C,name="const_ImGuiFocusedFlags_AnyWindow") :: ImGuiFocusedFlags_AnyWindow
  integer(c_int), bind(C,name="const_ImGuiFocusedFlags_NoPopupHierarchy") :: ImGuiFocusedFlags_NoPopupHierarchy
  integer(c_int), bind(C,name="const_ImGuiFocusedFlags_DockHierarchy") :: ImGuiFocusedFlags_DockHierarchy
  integer(c_int), bind(C,name="const_ImGuiFocusedFlags_RootAndChildWindows") :: ImGuiFocusedFlags_RootAndChildWindows
  ! enum ImGuiHoveredFlags_
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_None") :: ImGuiHoveredFlags_None
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_ChildWindows") :: ImGuiHoveredFlags_ChildWindows
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_RootWindow") :: ImGuiHoveredFlags_RootWindow
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AnyWindow") :: ImGuiHoveredFlags_AnyWindow
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_NoPopupHierarchy") :: ImGuiHoveredFlags_NoPopupHierarchy
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_DockHierarchy") :: ImGuiHoveredFlags_DockHierarchy
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenBlockedByPopup") :: ImGuiHoveredFlags_AllowWhenBlockedByPopup
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenBlockedByActiveItem") :: &
     ImGuiHoveredFlags_AllowWhenBlockedByActiveItem
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenOverlapped") :: ImGuiHoveredFlags_AllowWhenOverlapped
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_AllowWhenDisabled") :: ImGuiHoveredFlags_AllowWhenDisabled
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_NoNavOverride") :: ImGuiHoveredFlags_NoNavOverride
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_RectOnly") :: ImGuiHoveredFlags_RectOnly
  integer(c_int), bind(C,name="const_ImGuiHoveredFlags_RootAndChildWindows") :: ImGuiHoveredFlags_RootAndChildWindows
  ! enum ImGuiDockNodeFlags_
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_None") :: ImGuiDockNodeFlags_None
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_KeepAliveOnly") :: ImGuiDockNodeFlags_KeepAliveOnly
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoDockingInCentralNode") :: ImGuiDockNodeFlags_NoDockingInCentralNode
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_PassthruCentralNode") :: ImGuiDockNodeFlags_PassthruCentralNode
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoSplit") :: ImGuiDockNodeFlags_NoSplit
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoResize") :: ImGuiDockNodeFlags_NoResize
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_AutoHideTabBar") :: ImGuiDockNodeFlags_AutoHideTabBar
  ! enum ImGuiDragDropFlags_
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_None") :: ImGuiDragDropFlags_None
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_SourceNoPreviewTooltip") :: ImGuiDragDropFlags_SourceNoPreviewTooltip
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_SourceNoDisableHover") :: ImGuiDragDropFlags_SourceNoDisableHover
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_SourceNoHoldToOpenOthers") :: ImGuiDragDropFlags_SourceNoHoldToOpenOthers
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_SourceAllowNullID") :: ImGuiDragDropFlags_SourceAllowNullID
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_SourceExtern") :: ImGuiDragDropFlags_SourceExtern
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_SourceAutoExpirePayload") :: ImGuiDragDropFlags_SourceAutoExpirePayload
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_AcceptBeforeDelivery") :: ImGuiDragDropFlags_AcceptBeforeDelivery
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_AcceptNoDrawDefaultRect") :: ImGuiDragDropFlags_AcceptNoDrawDefaultRect
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_AcceptNoPreviewTooltip") :: ImGuiDragDropFlags_AcceptNoPreviewTooltip
  integer(c_int), bind(C,name="const_ImGuiDragDropFlags_AcceptPeekOnly") :: ImGuiDragDropFlags_AcceptPeekOnly
  ! enum ImGuiDataType_
  integer(c_int), bind(C,name="const_ImGuiDataType_S8") :: ImGuiDataType_S8
  integer(c_int), bind(C,name="const_ImGuiDataType_U8") :: ImGuiDataType_U8
  integer(c_int), bind(C,name="const_ImGuiDataType_S16") :: ImGuiDataType_S16
  integer(c_int), bind(C,name="const_ImGuiDataType_U16") :: ImGuiDataType_U16
  integer(c_int), bind(C,name="const_ImGuiDataType_S32") :: ImGuiDataType_S32
  integer(c_int), bind(C,name="const_ImGuiDataType_U32") :: ImGuiDataType_U32
  integer(c_int), bind(C,name="const_ImGuiDataType_S64") :: ImGuiDataType_S64
  integer(c_int), bind(C,name="const_ImGuiDataType_U64") :: ImGuiDataType_U64
  integer(c_int), bind(C,name="const_ImGuiDataType_Float") :: ImGuiDataType_Float
  integer(c_int), bind(C,name="const_ImGuiDataType_Double") :: ImGuiDataType_Double
  integer(c_int), bind(C,name="const_ImGuiDataType_COUNT") :: ImGuiDataType_COUNT
  ! enum ImGuiDir_
  integer(c_int), bind(C,name="const_ImGuiDir_None") :: ImGuiDir_None
  integer(c_int), bind(C,name="const_ImGuiDir_Left") :: ImGuiDir_Left
  integer(c_int), bind(C,name="const_ImGuiDir_Right") :: ImGuiDir_Right
  integer(c_int), bind(C,name="const_ImGuiDir_Up") :: ImGuiDir_Up
  integer(c_int), bind(C,name="const_ImGuiDir_Down") :: ImGuiDir_Down
  integer(c_int), bind(C,name="const_ImGuiDir_COUNT") :: ImGuiDir_COUNT
  ! enum ImGuiSortDirection_
  integer(c_int), bind(C,name="const_ImGuiSortDirection_None") :: ImGuiSortDirection_None
  integer(c_int), bind(C,name="const_ImGuiSortDirection_Ascending") :: ImGuiSortDirection_Ascending
  integer(c_int), bind(C,name="const_ImGuiSortDirection_Descending") :: ImGuiSortDirection_Descending
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
  !! integer(c_int), bind(C,name="const_ImGuiKey_COUNT") :: ImGuiKey_COUNT
  integer(c_int), bind(C,name="const_ImGuiKey_NamedKey_BEGIN") :: ImGuiKey_NamedKey_BEGIN
  integer(c_int), bind(C,name="const_ImGuiKey_NamedKey_END") :: ImGuiKey_NamedKey_END
  !! integer(c_int), bind(C,name="const_ImGuiKey_NamedKey_COUNT") :: ImGuiKey_NamedKey_COUNT
  !! integer(c_int), bind(C,name="const_ImGuiKey_KeysData_SIZE") :: ImGuiKey_KeysData_SIZE
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
  !! integer(c_int), bind(C,name="const_ImGuiNavInput_COUNT") :: ImGuiNavInput_COUNT
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
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_DpiEnableScaleViewports") :: ImGuiConfigFlags_DpiEnableScaleViewports
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_DpiEnableScaleFonts") :: ImGuiConfigFlags_DpiEnableScaleFonts
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_IsSRGB") :: ImGuiConfigFlags_IsSRGB
  integer(c_int), bind(C,name="const_ImGuiConfigFlags_IsTouchScreen") :: ImGuiConfigFlags_IsTouchScreen
  ! enum ImGuiBackendFlags_
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_None") :: ImGuiBackendFlags_None
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_HasGamepad") :: ImGuiBackendFlags_HasGamepad
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_HasMouseCursors") :: ImGuiBackendFlags_HasMouseCursors
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_HasSetMousePos") :: ImGuiBackendFlags_HasSetMousePos
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_RendererHasVtxOffset") :: ImGuiBackendFlags_RendererHasVtxOffset
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_PlatformHasViewports") :: ImGuiBackendFlags_PlatformHasViewports
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_HasMouseHoveredViewport") :: ImGuiBackendFlags_HasMouseHoveredViewport
  integer(c_int), bind(C,name="const_ImGuiBackendFlags_RendererHasViewports") :: ImGuiBackendFlags_RendererHasViewports
  ! enum ImGuiCol_
  integer(c_int), bind(C,name="const_ImGuiCol_Text") :: ImGuiCol_Text
  integer(c_int), bind(C,name="const_ImGuiCol_TextDisabled") :: ImGuiCol_TextDisabled
  integer(c_int), bind(C,name="const_ImGuiCol_WindowBg") :: ImGuiCol_WindowBg
  integer(c_int), bind(C,name="const_ImGuiCol_ChildBg") :: ImGuiCol_ChildBg
  integer(c_int), bind(C,name="const_ImGuiCol_PopupBg") :: ImGuiCol_PopupBg
  integer(c_int), bind(C,name="const_ImGuiCol_Border") :: ImGuiCol_Border
  integer(c_int), bind(C,name="const_ImGuiCol_BorderShadow") :: ImGuiCol_BorderShadow
  integer(c_int), bind(C,name="const_ImGuiCol_FrameBg") :: ImGuiCol_FrameBg
  integer(c_int), bind(C,name="const_ImGuiCol_FrameBgHovered") :: ImGuiCol_FrameBgHovered
  integer(c_int), bind(C,name="const_ImGuiCol_FrameBgActive") :: ImGuiCol_FrameBgActive
  integer(c_int), bind(C,name="const_ImGuiCol_TitleBg") :: ImGuiCol_TitleBg
  integer(c_int), bind(C,name="const_ImGuiCol_TitleBgActive") :: ImGuiCol_TitleBgActive
  integer(c_int), bind(C,name="const_ImGuiCol_TitleBgCollapsed") :: ImGuiCol_TitleBgCollapsed
  integer(c_int), bind(C,name="const_ImGuiCol_MenuBarBg") :: ImGuiCol_MenuBarBg
  integer(c_int), bind(C,name="const_ImGuiCol_ScrollbarBg") :: ImGuiCol_ScrollbarBg
  integer(c_int), bind(C,name="const_ImGuiCol_ScrollbarGrab") :: ImGuiCol_ScrollbarGrab
  integer(c_int), bind(C,name="const_ImGuiCol_ScrollbarGrabHovered") :: ImGuiCol_ScrollbarGrabHovered
  integer(c_int), bind(C,name="const_ImGuiCol_ScrollbarGrabActive") :: ImGuiCol_ScrollbarGrabActive
  integer(c_int), bind(C,name="const_ImGuiCol_CheckMark") :: ImGuiCol_CheckMark
  integer(c_int), bind(C,name="const_ImGuiCol_SliderGrab") :: ImGuiCol_SliderGrab
  integer(c_int), bind(C,name="const_ImGuiCol_SliderGrabActive") :: ImGuiCol_SliderGrabActive
  integer(c_int), bind(C,name="const_ImGuiCol_Button") :: ImGuiCol_Button
  integer(c_int), bind(C,name="const_ImGuiCol_ButtonHovered") :: ImGuiCol_ButtonHovered
  integer(c_int), bind(C,name="const_ImGuiCol_ButtonActive") :: ImGuiCol_ButtonActive
  integer(c_int), bind(C,name="const_ImGuiCol_Header") :: ImGuiCol_Header
  integer(c_int), bind(C,name="const_ImGuiCol_HeaderHovered") :: ImGuiCol_HeaderHovered
  integer(c_int), bind(C,name="const_ImGuiCol_HeaderActive") :: ImGuiCol_HeaderActive
  integer(c_int), bind(C,name="const_ImGuiCol_Separator") :: ImGuiCol_Separator
  integer(c_int), bind(C,name="const_ImGuiCol_SeparatorHovered") :: ImGuiCol_SeparatorHovered
  integer(c_int), bind(C,name="const_ImGuiCol_SeparatorActive") :: ImGuiCol_SeparatorActive
  integer(c_int), bind(C,name="const_ImGuiCol_ResizeGrip") :: ImGuiCol_ResizeGrip
  integer(c_int), bind(C,name="const_ImGuiCol_ResizeGripHovered") :: ImGuiCol_ResizeGripHovered
  integer(c_int), bind(C,name="const_ImGuiCol_ResizeGripActive") :: ImGuiCol_ResizeGripActive
  integer(c_int), bind(C,name="const_ImGuiCol_Tab") :: ImGuiCol_Tab
  integer(c_int), bind(C,name="const_ImGuiCol_TabHovered") :: ImGuiCol_TabHovered
  integer(c_int), bind(C,name="const_ImGuiCol_TabActive") :: ImGuiCol_TabActive
  integer(c_int), bind(C,name="const_ImGuiCol_TabUnfocused") :: ImGuiCol_TabUnfocused
  integer(c_int), bind(C,name="const_ImGuiCol_TabUnfocusedActive") :: ImGuiCol_TabUnfocusedActive
  integer(c_int), bind(C,name="const_ImGuiCol_DockingPreview") :: ImGuiCol_DockingPreview
  integer(c_int), bind(C,name="const_ImGuiCol_DockingEmptyBg") :: ImGuiCol_DockingEmptyBg
  integer(c_int), bind(C,name="const_ImGuiCol_PlotLines") :: ImGuiCol_PlotLines
  integer(c_int), bind(C,name="const_ImGuiCol_PlotLinesHovered") :: ImGuiCol_PlotLinesHovered
  integer(c_int), bind(C,name="const_ImGuiCol_PlotHistogram") :: ImGuiCol_PlotHistogram
  integer(c_int), bind(C,name="const_ImGuiCol_PlotHistogramHovered") :: ImGuiCol_PlotHistogramHovered
  integer(c_int), bind(C,name="const_ImGuiCol_TableHeaderBg") :: ImGuiCol_TableHeaderBg
  integer(c_int), bind(C,name="const_ImGuiCol_TableBorderStrong") :: ImGuiCol_TableBorderStrong
  integer(c_int), bind(C,name="const_ImGuiCol_TableBorderLight") :: ImGuiCol_TableBorderLight
  integer(c_int), bind(C,name="const_ImGuiCol_TableRowBg") :: ImGuiCol_TableRowBg
  integer(c_int), bind(C,name="const_ImGuiCol_TableRowBgAlt") :: ImGuiCol_TableRowBgAlt
  integer(c_int), bind(C,name="const_ImGuiCol_TextSelectedBg") :: ImGuiCol_TextSelectedBg
  integer(c_int), bind(C,name="const_ImGuiCol_DragDropTarget") :: ImGuiCol_DragDropTarget
  integer(c_int), bind(C,name="const_ImGuiCol_NavHighlight") :: ImGuiCol_NavHighlight
  integer(c_int), bind(C,name="const_ImGuiCol_NavWindowingHighlight") :: ImGuiCol_NavWindowingHighlight
  integer(c_int), bind(C,name="const_ImGuiCol_NavWindowingDimBg") :: ImGuiCol_NavWindowingDimBg
  integer(c_int), bind(C,name="const_ImGuiCol_ModalWindowDimBg") :: ImGuiCol_ModalWindowDimBg
  !!integer(c_int), bind(C,name="const_ImGuiCol_COUNT") :: ImGuiCol_COUNT
  ! enum ImGuiStyleVar_
  integer(c_int), bind(C,name="const_ImGuiStyleVar_Alpha") :: ImGuiStyleVar_Alpha
  integer(c_int), bind(C,name="const_ImGuiStyleVar_DisabledAlpha") :: ImGuiStyleVar_DisabledAlpha
  integer(c_int), bind(C,name="const_ImGuiStyleVar_WindowPadding") :: ImGuiStyleVar_WindowPadding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_WindowRounding") :: ImGuiStyleVar_WindowRounding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_WindowBorderSize") :: ImGuiStyleVar_WindowBorderSize
  integer(c_int), bind(C,name="const_ImGuiStyleVar_WindowMinSize") :: ImGuiStyleVar_WindowMinSize
  integer(c_int), bind(C,name="const_ImGuiStyleVar_WindowTitleAlign") :: ImGuiStyleVar_WindowTitleAlign
  integer(c_int), bind(C,name="const_ImGuiStyleVar_ChildRounding") :: ImGuiStyleVar_ChildRounding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_ChildBorderSize") :: ImGuiStyleVar_ChildBorderSize
  integer(c_int), bind(C,name="const_ImGuiStyleVar_PopupRounding") :: ImGuiStyleVar_PopupRounding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_PopupBorderSize") :: ImGuiStyleVar_PopupBorderSize
  integer(c_int), bind(C,name="const_ImGuiStyleVar_FramePadding") :: ImGuiStyleVar_FramePadding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_FrameRounding") :: ImGuiStyleVar_FrameRounding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_FrameBorderSize") :: ImGuiStyleVar_FrameBorderSize
  integer(c_int), bind(C,name="const_ImGuiStyleVar_ItemSpacing") :: ImGuiStyleVar_ItemSpacing
  integer(c_int), bind(C,name="const_ImGuiStyleVar_ItemInnerSpacing") :: ImGuiStyleVar_ItemInnerSpacing
  integer(c_int), bind(C,name="const_ImGuiStyleVar_IndentSpacing") :: ImGuiStyleVar_IndentSpacing
  integer(c_int), bind(C,name="const_ImGuiStyleVar_CellPadding") :: ImGuiStyleVar_CellPadding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_ScrollbarSize") :: ImGuiStyleVar_ScrollbarSize
  integer(c_int), bind(C,name="const_ImGuiStyleVar_ScrollbarRounding") :: ImGuiStyleVar_ScrollbarRounding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_GrabMinSize") :: ImGuiStyleVar_GrabMinSize
  integer(c_int), bind(C,name="const_ImGuiStyleVar_GrabRounding") :: ImGuiStyleVar_GrabRounding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_TabRounding") :: ImGuiStyleVar_TabRounding
  integer(c_int), bind(C,name="const_ImGuiStyleVar_ButtonTextAlign") :: ImGuiStyleVar_ButtonTextAlign
  integer(c_int), bind(C,name="const_ImGuiStyleVar_SelectableTextAlign") :: ImGuiStyleVar_SelectableTextAlign
  integer(c_int), bind(C,name="const_ImGuiStyleVar_COUNT") :: ImGuiStyleVar_COUNT
  ! enum ImGuiButtonFlags_
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_None") :: ImGuiButtonFlags_None
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_MouseButtonLeft") :: ImGuiButtonFlags_MouseButtonLeft
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_MouseButtonRight") :: ImGuiButtonFlags_MouseButtonRight
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_MouseButtonMiddle") :: ImGuiButtonFlags_MouseButtonMiddle
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_MouseButtonMask_") :: ImGuiButtonFlags_MouseButtonMask_
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_MouseButtonDefault_") :: ImGuiButtonFlags_MouseButtonDefault_
  ! enum ImGuiColorEditFlags_
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_None") :: ImGuiColorEditFlags_None
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoAlpha") :: ImGuiColorEditFlags_NoAlpha
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoPicker") :: ImGuiColorEditFlags_NoPicker
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoOptions") :: ImGuiColorEditFlags_NoOptions
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoSmallPreview") :: ImGuiColorEditFlags_NoSmallPreview
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoInputs") :: ImGuiColorEditFlags_NoInputs
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoTooltip") :: ImGuiColorEditFlags_NoTooltip
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoLabel") :: ImGuiColorEditFlags_NoLabel
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoSidePreview") :: ImGuiColorEditFlags_NoSidePreview
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoDragDrop") :: ImGuiColorEditFlags_NoDragDrop
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_NoBorder") :: ImGuiColorEditFlags_NoBorder
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_AlphaBar") :: ImGuiColorEditFlags_AlphaBar
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_AlphaPreview") :: ImGuiColorEditFlags_AlphaPreview
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_AlphaPreviewHalf") :: ImGuiColorEditFlags_AlphaPreviewHalf
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_HDR") :: ImGuiColorEditFlags_HDR
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_DisplayRGB") :: ImGuiColorEditFlags_DisplayRGB
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_DisplayHSV") :: ImGuiColorEditFlags_DisplayHSV
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_DisplayHex") :: ImGuiColorEditFlags_DisplayHex
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_Uint8") :: ImGuiColorEditFlags_Uint8
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_Float") :: ImGuiColorEditFlags_Float
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_PickerHueBar") :: ImGuiColorEditFlags_PickerHueBar
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_PickerHueWheel") :: ImGuiColorEditFlags_PickerHueWheel
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_InputRGB") :: ImGuiColorEditFlags_InputRGB
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_InputHSV") :: ImGuiColorEditFlags_InputHSV
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_DefaultOptions_") :: ImGuiColorEditFlags_DefaultOptions_
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_DisplayMask_") :: ImGuiColorEditFlags_DisplayMask_
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_DataTypeMask_") :: ImGuiColorEditFlags_DataTypeMask_
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_PickerMask_") :: ImGuiColorEditFlags_PickerMask_
  integer(c_int), bind(C,name="const_ImGuiColorEditFlags_InputMask_") :: ImGuiColorEditFlags_InputMask_
  ! enum ImGuiSliderFlags_
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_None") :: ImGuiSliderFlags_None
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_AlwaysClamp") :: ImGuiSliderFlags_AlwaysClamp
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_Logarithmic") :: ImGuiSliderFlags_Logarithmic
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_NoRoundToFormat") :: ImGuiSliderFlags_NoRoundToFormat
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_NoInput") :: ImGuiSliderFlags_NoInput
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_InvalidMask_") :: ImGuiSliderFlags_InvalidMask_
  ! enum ImGuiMouseButton_
  integer(c_int), bind(C,name="const_ImGuiMouseButton_Left") :: ImGuiMouseButton_Left
  integer(c_int), bind(C,name="const_ImGuiMouseButton_Right") :: ImGuiMouseButton_Right
  integer(c_int), bind(C,name="const_ImGuiMouseButton_Middle") :: ImGuiMouseButton_Middle
  integer(c_int), bind(C,name="const_ImGuiMouseButton_COUNT") :: ImGuiMouseButton_COUNT
  ! enum ImGuiMouseCursor_
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_None") :: ImGuiMouseCursor_None
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_Arrow") :: ImGuiMouseCursor_Arrow
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_TextInput") :: ImGuiMouseCursor_TextInput
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_ResizeAll") :: ImGuiMouseCursor_ResizeAll
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_ResizeNS") :: ImGuiMouseCursor_ResizeNS
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_ResizeEW") :: ImGuiMouseCursor_ResizeEW
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_ResizeNESW") :: ImGuiMouseCursor_ResizeNESW
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_ResizeNWSE") :: ImGuiMouseCursor_ResizeNWSE
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_Hand") :: ImGuiMouseCursor_Hand
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_NotAllowed") :: ImGuiMouseCursor_NotAllowed
  integer(c_int), bind(C,name="const_ImGuiMouseCursor_COUNT") :: ImGuiMouseCursor_COUNT
  ! enum ImGuiCond_
  integer(c_int), bind(C,name="const_ImGuiCond_None") :: ImGuiCond_None
  integer(c_int), bind(C,name="const_ImGuiCond_Always") :: ImGuiCond_Always
  integer(c_int), bind(C,name="const_ImGuiCond_Once") :: ImGuiCond_Once
  integer(c_int), bind(C,name="const_ImGuiCond_FirstUseEver") :: ImGuiCond_FirstUseEver
  integer(c_int), bind(C,name="const_ImGuiCond_Appearing") :: ImGuiCond_Appearing
  ! enum ImDrawFlags_
  integer(c_int), bind(C,name="const_ImDrawFlags_None") :: ImDrawFlags_None
  integer(c_int), bind(C,name="const_ImDrawFlags_Closed") :: ImDrawFlags_Closed
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersTopLeft") :: ImDrawFlags_RoundCornersTopLeft
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersTopRight") :: ImDrawFlags_RoundCornersTopRight
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersBottomLeft") :: ImDrawFlags_RoundCornersBottomLeft
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersBottomRight") :: ImDrawFlags_RoundCornersBottomRight
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersNone") :: ImDrawFlags_RoundCornersNone
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersTop") :: ImDrawFlags_RoundCornersTop
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersBottom") :: ImDrawFlags_RoundCornersBottom
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersLeft") :: ImDrawFlags_RoundCornersLeft
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersRight") :: ImDrawFlags_RoundCornersRight
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersAll") :: ImDrawFlags_RoundCornersAll
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersDefault_") :: ImDrawFlags_RoundCornersDefault_
  integer(c_int), bind(C,name="const_ImDrawFlags_RoundCornersMask_") :: ImDrawFlags_RoundCornersMask_
  ! enum ImDrawListFlags_
  integer(c_int), bind(C,name="const_ImDrawListFlags_None") :: ImDrawListFlags_None
  integer(c_int), bind(C,name="const_ImDrawListFlags_AntiAliasedLines") :: ImDrawListFlags_AntiAliasedLines
  integer(c_int), bind(C,name="const_ImDrawListFlags_AntiAliasedLinesUseTex") :: ImDrawListFlags_AntiAliasedLinesUseTex
  integer(c_int), bind(C,name="const_ImDrawListFlags_AntiAliasedFill") :: ImDrawListFlags_AntiAliasedFill
  integer(c_int), bind(C,name="const_ImDrawListFlags_AllowVtxOffset") :: ImDrawListFlags_AllowVtxOffset
  ! enum ImFontAtlasFlags_
  integer(c_int), bind(C,name="const_ImFontAtlasFlags_None") :: ImFontAtlasFlags_None
  integer(c_int), bind(C,name="const_ImFontAtlasFlags_NoPowerOfTwoHeight") :: ImFontAtlasFlags_NoPowerOfTwoHeight
  integer(c_int), bind(C,name="const_ImFontAtlasFlags_NoMouseCursors") :: ImFontAtlasFlags_NoMouseCursors
  integer(c_int), bind(C,name="const_ImFontAtlasFlags_NoBakedLines") :: ImFontAtlasFlags_NoBakedLines
  ! enum ImGuiViewportFlags_
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_None") :: ImGuiViewportFlags_None
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_IsPlatformWindow") :: ImGuiViewportFlags_IsPlatformWindow
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_IsPlatformMonitor") :: ImGuiViewportFlags_IsPlatformMonitor
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_OwnedByApp") :: ImGuiViewportFlags_OwnedByApp
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_NoDecoration") :: ImGuiViewportFlags_NoDecoration
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_NoTaskBarIcon") :: ImGuiViewportFlags_NoTaskBarIcon
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_NoFocusOnAppearing") :: ImGuiViewportFlags_NoFocusOnAppearing
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_NoFocusOnClick") :: ImGuiViewportFlags_NoFocusOnClick
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_NoInputs") :: ImGuiViewportFlags_NoInputs
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_NoRendererClear") :: ImGuiViewportFlags_NoRendererClear
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_TopMost") :: ImGuiViewportFlags_TopMost
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_Minimized") :: ImGuiViewportFlags_Minimized
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_NoAutoMerge") :: ImGuiViewportFlags_NoAutoMerge
  integer(c_int), bind(C,name="const_ImGuiViewportFlags_CanHostOtherWindows") :: ImGuiViewportFlags_CanHostOtherWindows
  ! enum ImGuiItemFlags_
  integer(c_int), bind(C,name="const_ImGuiItemFlags_None") :: ImGuiItemFlags_None
  integer(c_int), bind(C,name="const_ImGuiItemFlags_NoTabStop") :: ImGuiItemFlags_NoTabStop
  integer(c_int), bind(C,name="const_ImGuiItemFlags_ButtonRepeat") :: ImGuiItemFlags_ButtonRepeat
  integer(c_int), bind(C,name="const_ImGuiItemFlags_Disabled") :: ImGuiItemFlags_Disabled
  integer(c_int), bind(C,name="const_ImGuiItemFlags_NoNav") :: ImGuiItemFlags_NoNav
  integer(c_int), bind(C,name="const_ImGuiItemFlags_NoNavDefaultFocus") :: ImGuiItemFlags_NoNavDefaultFocus
  integer(c_int), bind(C,name="const_ImGuiItemFlags_SelectableDontClosePopup") :: ImGuiItemFlags_SelectableDontClosePopup
  integer(c_int), bind(C,name="const_ImGuiItemFlags_MixedValue") :: ImGuiItemFlags_MixedValue
  integer(c_int), bind(C,name="const_ImGuiItemFlags_ReadOnly") :: ImGuiItemFlags_ReadOnly
  integer(c_int), bind(C,name="const_ImGuiItemFlags_Inputable") :: ImGuiItemFlags_Inputable
  ! enum ImGuiItemStatusFlags_
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_None") :: ImGuiItemStatusFlags_None
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_HoveredRect") :: ImGuiItemStatusFlags_HoveredRect
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_HasDisplayRect") :: ImGuiItemStatusFlags_HasDisplayRect
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_Edited") :: ImGuiItemStatusFlags_Edited
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_ToggledSelection") :: ImGuiItemStatusFlags_ToggledSelection
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_ToggledOpen") :: ImGuiItemStatusFlags_ToggledOpen
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_HasDeactivated") :: ImGuiItemStatusFlags_HasDeactivated
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_Deactivated") :: ImGuiItemStatusFlags_Deactivated
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_HoveredWindow") :: ImGuiItemStatusFlags_HoveredWindow
  integer(c_int), bind(C,name="const_ImGuiItemStatusFlags_FocusedByTabbing") :: ImGuiItemStatusFlags_FocusedByTabbing
  ! enum ImGuiInputTextFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_Multiline") :: ImGuiInputTextFlags_Multiline
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_NoMarkEdited") :: ImGuiInputTextFlags_NoMarkEdited
  integer(c_int), bind(C,name="const_ImGuiInputTextFlags_MergedItem") :: ImGuiInputTextFlags_MergedItem
  ! enum ImGuiButtonFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnClick") :: ImGuiButtonFlags_PressedOnClick
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnClickRelease") :: ImGuiButtonFlags_PressedOnClickRelease
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnClickReleaseAnywhere") :: &
     ImGuiButtonFlags_PressedOnClickReleaseAnywhere
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnRelease") :: ImGuiButtonFlags_PressedOnRelease
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnDoubleClick") :: ImGuiButtonFlags_PressedOnDoubleClick
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnDragDropHold") :: ImGuiButtonFlags_PressedOnDragDropHold
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_Repeat") :: ImGuiButtonFlags_Repeat
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_FlattenChildren") :: ImGuiButtonFlags_FlattenChildren
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_AllowItemOverlap") :: ImGuiButtonFlags_AllowItemOverlap
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_DontClosePopups") :: ImGuiButtonFlags_DontClosePopups
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_AlignTextBaseLine") :: ImGuiButtonFlags_AlignTextBaseLine
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_NoKeyModifiers") :: ImGuiButtonFlags_NoKeyModifiers
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_NoHoldingActiveId") :: ImGuiButtonFlags_NoHoldingActiveId
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_NoNavFocus") :: ImGuiButtonFlags_NoNavFocus
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_NoHoveredOnFocus") :: ImGuiButtonFlags_NoHoveredOnFocus
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnMask_") :: ImGuiButtonFlags_PressedOnMask_
  integer(c_int), bind(C,name="const_ImGuiButtonFlags_PressedOnDefault_") :: ImGuiButtonFlags_PressedOnDefault_
  ! enum ImGuiComboFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiComboFlags_CustomPreview") :: ImGuiComboFlags_CustomPreview
  ! enum ImGuiSliderFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_Vertical") :: ImGuiSliderFlags_Vertical
  integer(c_int), bind(C,name="const_ImGuiSliderFlags_ReadOnly") :: ImGuiSliderFlags_ReadOnly
  ! enum ImGuiSelectableFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_NoHoldingActiveID") :: ImGuiSelectableFlags_NoHoldingActiveID
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_SelectOnNav") :: ImGuiSelectableFlags_SelectOnNav
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_SelectOnClick") :: ImGuiSelectableFlags_SelectOnClick
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_SelectOnRelease") :: ImGuiSelectableFlags_SelectOnRelease
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_SpanAvailWidth") :: ImGuiSelectableFlags_SpanAvailWidth
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_DrawHoveredWhenHeld") :: ImGuiSelectableFlags_DrawHoveredWhenHeld
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_SetNavIdOnHover") :: ImGuiSelectableFlags_SetNavIdOnHover
  integer(c_int), bind(C,name="const_ImGuiSelectableFlags_NoPadWithHalfSpacing") :: ImGuiSelectableFlags_NoPadWithHalfSpacing
  ! enum ImGuiTreeNodeFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiTreeNodeFlags_ClipLabelForTrailingButton") :: &
     ImGuiTreeNodeFlags_ClipLabelForTrailingButton
  ! enum ImGuiSeparatorFlags_
  integer(c_int), bind(C,name="const_ImGuiSeparatorFlags_None") :: ImGuiSeparatorFlags_None
  integer(c_int), bind(C,name="const_ImGuiSeparatorFlags_Horizontal") :: ImGuiSeparatorFlags_Horizontal
  integer(c_int), bind(C,name="const_ImGuiSeparatorFlags_Vertical") :: ImGuiSeparatorFlags_Vertical
  integer(c_int), bind(C,name="const_ImGuiSeparatorFlags_SpanAllColumns") :: ImGuiSeparatorFlags_SpanAllColumns
  ! enum ImGuiTextFlags_
  integer(c_int), bind(C,name="const_ImGuiTextFlags_None") :: ImGuiTextFlags_None
  integer(c_int), bind(C,name="const_ImGuiTextFlags_NoWidthForLargeClippedText") :: ImGuiTextFlags_NoWidthForLargeClippedText
  ! enum ImGuiTooltipFlags_
  integer(c_int), bind(C,name="const_ImGuiTooltipFlags_None") :: ImGuiTooltipFlags_None
  integer(c_int), bind(C,name="const_ImGuiTooltipFlags_OverridePreviousTooltip") :: ImGuiTooltipFlags_OverridePreviousTooltip
  ! enum ImGuiLayoutType_
  integer(c_int), bind(C,name="const_ImGuiLayoutType_Horizontal") :: ImGuiLayoutType_Horizontal
  integer(c_int), bind(C,name="const_ImGuiLayoutType_Vertical") :: ImGuiLayoutType_Vertical
  ! enum ImGuiLogType
  integer(c_int), bind(C,name="const_ImGuiLogType_None") :: ImGuiLogType_None
  integer(c_int), bind(C,name="const_ImGuiLogType_TTY") :: ImGuiLogType_TTY
  integer(c_int), bind(C,name="const_ImGuiLogType_File") :: ImGuiLogType_File
  integer(c_int), bind(C,name="const_ImGuiLogType_Buffer") :: ImGuiLogType_Buffer
  integer(c_int), bind(C,name="const_ImGuiLogType_Clipboard") :: ImGuiLogType_Clipboard
  ! enum ImGuiAxis
  integer(c_int), bind(C,name="const_ImGuiAxis_None") :: ImGuiAxis_None
  integer(c_int), bind(C,name="const_ImGuiAxis_X") :: ImGuiAxis_X
  integer(c_int), bind(C,name="const_ImGuiAxis_Y") :: ImGuiAxis_Y
  ! enum ImGuiPlotType
  integer(c_int), bind(C,name="const_ImGuiPlotType_Lines") :: ImGuiPlotType_Lines
  integer(c_int), bind(C,name="const_ImGuiPlotType_Histogram") :: ImGuiPlotType_Histogram
  ! enum ImGuiPopupPositionPolicy
  integer(c_int), bind(C,name="const_ImGuiPopupPositionPolicy_Default") :: ImGuiPopupPositionPolicy_Default
  integer(c_int), bind(C,name="const_ImGuiPopupPositionPolicy_ComboBox") :: ImGuiPopupPositionPolicy_ComboBox
  integer(c_int), bind(C,name="const_ImGuiPopupPositionPolicy_Tooltip") :: ImGuiPopupPositionPolicy_Tooltip
  ! enum ImGuiDataTypePrivate_
  integer(c_int), bind(C,name="const_ImGuiDataType_String") :: ImGuiDataType_String
  integer(c_int), bind(C,name="const_ImGuiDataType_Pointer") :: ImGuiDataType_Pointer
  integer(c_int), bind(C,name="const_ImGuiDataType_ID") :: ImGuiDataType_ID
  ! enum ImGuiNextWindowDataFlags_
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_None") :: ImGuiNextWindowDataFlags_None
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasPos") :: ImGuiNextWindowDataFlags_HasPos
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasSize") :: ImGuiNextWindowDataFlags_HasSize
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasContentSize") :: ImGuiNextWindowDataFlags_HasContentSize
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasCollapsed") :: ImGuiNextWindowDataFlags_HasCollapsed
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasSizeConstraint") :: ImGuiNextWindowDataFlags_HasSizeConstraint
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasFocus") :: ImGuiNextWindowDataFlags_HasFocus
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasBgAlpha") :: ImGuiNextWindowDataFlags_HasBgAlpha
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasScroll") :: ImGuiNextWindowDataFlags_HasScroll
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasViewport") :: ImGuiNextWindowDataFlags_HasViewport
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasDock") :: ImGuiNextWindowDataFlags_HasDock
  integer(c_int), bind(C,name="const_ImGuiNextWindowDataFlags_HasWindowClass") :: ImGuiNextWindowDataFlags_HasWindowClass
  ! enum ImGuiNextItemDataFlags_
  integer(c_int), bind(C,name="const_ImGuiNextItemDataFlags_None") :: ImGuiNextItemDataFlags_None
  integer(c_int), bind(C,name="const_ImGuiNextItemDataFlags_HasWidth") :: ImGuiNextItemDataFlags_HasWidth
  integer(c_int), bind(C,name="const_ImGuiNextItemDataFlags_HasOpen") :: ImGuiNextItemDataFlags_HasOpen
  ! enum ImGuiKeyPrivate_
  integer(c_int), bind(C,name="const_ImGuiKey_LegacyNativeKey_BEGIN") :: ImGuiKey_LegacyNativeKey_BEGIN
  integer(c_int), bind(C,name="const_ImGuiKey_LegacyNativeKey_END") :: ImGuiKey_LegacyNativeKey_END
  integer(c_int), bind(C,name="const_ImGuiKey_Gamepad_BEGIN") :: ImGuiKey_Gamepad_BEGIN
  integer(c_int), bind(C,name="const_ImGuiKey_Gamepad_END") :: ImGuiKey_Gamepad_END
  ! enum ImGuiInputEventType
  integer(c_int), bind(C,name="const_ImGuiInputEventType_None") :: ImGuiInputEventType_None
  integer(c_int), bind(C,name="const_ImGuiInputEventType_MousePos") :: ImGuiInputEventType_MousePos
  integer(c_int), bind(C,name="const_ImGuiInputEventType_MouseWheel") :: ImGuiInputEventType_MouseWheel
  integer(c_int), bind(C,name="const_ImGuiInputEventType_MouseButton") :: ImGuiInputEventType_MouseButton
  integer(c_int), bind(C,name="const_ImGuiInputEventType_MouseViewport") :: ImGuiInputEventType_MouseViewport
  integer(c_int), bind(C,name="const_ImGuiInputEventType_Key") :: ImGuiInputEventType_Key
  integer(c_int), bind(C,name="const_ImGuiInputEventType_Text") :: ImGuiInputEventType_Text
  integer(c_int), bind(C,name="const_ImGuiInputEventType_Focus") :: ImGuiInputEventType_Focus
  integer(c_int), bind(C,name="const_ImGuiInputEventType_COUNT") :: ImGuiInputEventType_COUNT
  ! enum ImGuiInputSource
  integer(c_int), bind(C,name="const_ImGuiInputSource_None") :: ImGuiInputSource_None
  integer(c_int), bind(C,name="const_ImGuiInputSource_Mouse") :: ImGuiInputSource_Mouse
  integer(c_int), bind(C,name="const_ImGuiInputSource_Keyboard") :: ImGuiInputSource_Keyboard
  integer(c_int), bind(C,name="const_ImGuiInputSource_Gamepad") :: ImGuiInputSource_Gamepad
  integer(c_int), bind(C,name="const_ImGuiInputSource_Clipboard") :: ImGuiInputSource_Clipboard
  integer(c_int), bind(C,name="const_ImGuiInputSource_Nav") :: ImGuiInputSource_Nav
  integer(c_int), bind(C,name="const_ImGuiInputSource_COUNT") :: ImGuiInputSource_COUNT
  ! enum ImGuiNavReadMode
  integer(c_int), bind(C,name="const_ImGuiNavReadMode_Down") :: ImGuiNavReadMode_Down
  integer(c_int), bind(C,name="const_ImGuiNavReadMode_Pressed") :: ImGuiNavReadMode_Pressed
  integer(c_int), bind(C,name="const_ImGuiNavReadMode_Released") :: ImGuiNavReadMode_Released
  integer(c_int), bind(C,name="const_ImGuiNavReadMode_Repeat") :: ImGuiNavReadMode_Repeat
  integer(c_int), bind(C,name="const_ImGuiNavReadMode_RepeatSlow") :: ImGuiNavReadMode_RepeatSlow
  integer(c_int), bind(C,name="const_ImGuiNavReadMode_RepeatFast") :: ImGuiNavReadMode_RepeatFast
  ! enum ImGuiActivateFlags_
  integer(c_int), bind(C,name="const_ImGuiActivateFlags_None") :: ImGuiActivateFlags_None
  integer(c_int), bind(C,name="const_ImGuiActivateFlags_PreferInput") :: ImGuiActivateFlags_PreferInput
  integer(c_int), bind(C,name="const_ImGuiActivateFlags_PreferTweak") :: ImGuiActivateFlags_PreferTweak
  integer(c_int), bind(C,name="const_ImGuiActivateFlags_TryToPreserveState") :: ImGuiActivateFlags_TryToPreserveState
  ! enum ImGuiScrollFlags_
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_None") :: ImGuiScrollFlags_None
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_KeepVisibleEdgeX") :: ImGuiScrollFlags_KeepVisibleEdgeX
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_KeepVisibleEdgeY") :: ImGuiScrollFlags_KeepVisibleEdgeY
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_KeepVisibleCenterX") :: ImGuiScrollFlags_KeepVisibleCenterX
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_KeepVisibleCenterY") :: ImGuiScrollFlags_KeepVisibleCenterY
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_AlwaysCenterX") :: ImGuiScrollFlags_AlwaysCenterX
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_AlwaysCenterY") :: ImGuiScrollFlags_AlwaysCenterY
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_NoScrollParent") :: ImGuiScrollFlags_NoScrollParent
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_MaskX_") :: ImGuiScrollFlags_MaskX_
  integer(c_int), bind(C,name="const_ImGuiScrollFlags_MaskY_") :: ImGuiScrollFlags_MaskY_
  ! enum ImGuiNavHighlightFlags_
  integer(c_int), bind(C,name="const_ImGuiNavHighlightFlags_None") :: ImGuiNavHighlightFlags_None
  integer(c_int), bind(C,name="const_ImGuiNavHighlightFlags_TypeDefault") :: ImGuiNavHighlightFlags_TypeDefault
  integer(c_int), bind(C,name="const_ImGuiNavHighlightFlags_TypeThin") :: ImGuiNavHighlightFlags_TypeThin
  integer(c_int), bind(C,name="const_ImGuiNavHighlightFlags_AlwaysDraw") :: ImGuiNavHighlightFlags_AlwaysDraw
  integer(c_int), bind(C,name="const_ImGuiNavHighlightFlags_NoRounding") :: ImGuiNavHighlightFlags_NoRounding
  ! enum ImGuiNavDirSourceFlags_
  integer(c_int), bind(C,name="const_ImGuiNavDirSourceFlags_None") :: ImGuiNavDirSourceFlags_None
  integer(c_int), bind(C,name="const_ImGuiNavDirSourceFlags_RawKeyboard") :: ImGuiNavDirSourceFlags_RawKeyboard
  integer(c_int), bind(C,name="const_ImGuiNavDirSourceFlags_Keyboard") :: ImGuiNavDirSourceFlags_Keyboard
  integer(c_int), bind(C,name="const_ImGuiNavDirSourceFlags_PadDPad") :: ImGuiNavDirSourceFlags_PadDPad
  integer(c_int), bind(C,name="const_ImGuiNavDirSourceFlags_PadLStick") :: ImGuiNavDirSourceFlags_PadLStick
  ! enum ImGuiNavMoveFlags_
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_None") :: ImGuiNavMoveFlags_None
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_LoopX") :: ImGuiNavMoveFlags_LoopX
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_LoopY") :: ImGuiNavMoveFlags_LoopY
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_WrapX") :: ImGuiNavMoveFlags_WrapX
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_WrapY") :: ImGuiNavMoveFlags_WrapY
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_AllowCurrentNavId") :: ImGuiNavMoveFlags_AllowCurrentNavId
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_AlsoScoreVisibleSet") :: ImGuiNavMoveFlags_AlsoScoreVisibleSet
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_ScrollToEdgeY") :: ImGuiNavMoveFlags_ScrollToEdgeY
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_Forwarded") :: ImGuiNavMoveFlags_Forwarded
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_DebugNoResult") :: ImGuiNavMoveFlags_DebugNoResult
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_FocusApi") :: ImGuiNavMoveFlags_FocusApi
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_Tabbing") :: ImGuiNavMoveFlags_Tabbing
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_Activate") :: ImGuiNavMoveFlags_Activate
  integer(c_int), bind(C,name="const_ImGuiNavMoveFlags_DontSetNavHighlight") :: ImGuiNavMoveFlags_DontSetNavHighlight
  ! enum ImGuiNavLayer
  integer(c_int), bind(C,name="const_ImGuiNavLayer_Main") :: ImGuiNavLayer_Main
  integer(c_int), bind(C,name="const_ImGuiNavLayer_Menu") :: ImGuiNavLayer_Menu
  integer(c_int), bind(C,name="const_ImGuiNavLayer_COUNT") :: ImGuiNavLayer_COUNT
  ! enum ImGuiOldColumnFlags_
  integer(c_int), bind(C,name="const_ImGuiOldColumnFlags_None") :: ImGuiOldColumnFlags_None
  integer(c_int), bind(C,name="const_ImGuiOldColumnFlags_NoBorder") :: ImGuiOldColumnFlags_NoBorder
  integer(c_int), bind(C,name="const_ImGuiOldColumnFlags_NoResize") :: ImGuiOldColumnFlags_NoResize
  integer(c_int), bind(C,name="const_ImGuiOldColumnFlags_NoPreserveWidths") :: ImGuiOldColumnFlags_NoPreserveWidths
  integer(c_int), bind(C,name="const_ImGuiOldColumnFlags_NoForceWithinWindow") :: ImGuiOldColumnFlags_NoForceWithinWindow
  integer(c_int), bind(C,name="const_ImGuiOldColumnFlags_GrowParentContentsSize") :: ImGuiOldColumnFlags_GrowParentContentsSize
  ! enum ImGuiDockNodeFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_DockSpace") :: ImGuiDockNodeFlags_DockSpace
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_CentralNode") :: ImGuiDockNodeFlags_CentralNode
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoTabBar") :: ImGuiDockNodeFlags_NoTabBar
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_HiddenTabBar") :: ImGuiDockNodeFlags_HiddenTabBar
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoWindowMenuButton") :: ImGuiDockNodeFlags_NoWindowMenuButton
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoCloseButton") :: ImGuiDockNodeFlags_NoCloseButton
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoDocking") :: ImGuiDockNodeFlags_NoDocking
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoDockingSplitMe") :: ImGuiDockNodeFlags_NoDockingSplitMe
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoDockingSplitOther") :: ImGuiDockNodeFlags_NoDockingSplitOther
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoDockingOverMe") :: ImGuiDockNodeFlags_NoDockingOverMe
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoDockingOverOther") :: ImGuiDockNodeFlags_NoDockingOverOther
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoDockingOverEmpty") :: ImGuiDockNodeFlags_NoDockingOverEmpty
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoResizeX") :: ImGuiDockNodeFlags_NoResizeX
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoResizeY") :: ImGuiDockNodeFlags_NoResizeY
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_SharedFlagsInheritMask_") :: ImGuiDockNodeFlags_SharedFlagsInheritMask_
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_NoResizeFlagsMask_") :: ImGuiDockNodeFlags_NoResizeFlagsMask_
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_LocalFlagsMask_") :: ImGuiDockNodeFlags_LocalFlagsMask_
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_LocalFlagsTransferMask_") :: ImGuiDockNodeFlags_LocalFlagsTransferMask_
  integer(c_int), bind(C,name="const_ImGuiDockNodeFlags_SavedFlagsMask_") :: ImGuiDockNodeFlags_SavedFlagsMask_
  ! enum ImGuiDataAuthority_
  integer(c_int), bind(C,name="const_ImGuiDataAuthority_Auto") :: ImGuiDataAuthority_Auto
  integer(c_int), bind(C,name="const_ImGuiDataAuthority_DockNode") :: ImGuiDataAuthority_DockNode
  integer(c_int), bind(C,name="const_ImGuiDataAuthority_Window") :: ImGuiDataAuthority_Window
  ! enum ImGuiDockNodeState
  integer(c_int), bind(C,name="const_ImGuiDockNodeState_Unknown") :: ImGuiDockNodeState_Unknown
  integer(c_int), bind(C,name="const_ImGuiDockNodeState_HostWindowHiddenBecauseSingleWindow") :: &
     ImGuiDockNodeState_HostWindowHiddenBecauseSingleWindow
  integer(c_int), bind(C,name="const_ImGuiDockNodeState_HostWindowHiddenBecauseWindowsAreResizing") :: &
     ImGuiDockNodeState_HostWindowHiddenBecauseWindowsAreResizing
  integer(c_int), bind(C,name="const_ImGuiDockNodeState_HostWindowVisible") :: ImGuiDockNodeState_HostWindowVisible
  ! enum ImGuiWindowDockStyleCol
  integer(c_int), bind(C,name="const_ImGuiWindowDockStyleCol_Text") :: ImGuiWindowDockStyleCol_Text
  integer(c_int), bind(C,name="const_ImGuiWindowDockStyleCol_Tab") :: ImGuiWindowDockStyleCol_Tab
  integer(c_int), bind(C,name="const_ImGuiWindowDockStyleCol_TabHovered") :: ImGuiWindowDockStyleCol_TabHovered
  integer(c_int), bind(C,name="const_ImGuiWindowDockStyleCol_TabActive") :: ImGuiWindowDockStyleCol_TabActive
  integer(c_int), bind(C,name="const_ImGuiWindowDockStyleCol_TabUnfocused") :: ImGuiWindowDockStyleCol_TabUnfocused
  integer(c_int), bind(C,name="const_ImGuiWindowDockStyleCol_TabUnfocusedActive") :: ImGuiWindowDockStyleCol_TabUnfocusedActive
  integer(c_int), bind(C,name="const_ImGuiWindowDockStyleCol_COUNT") :: ImGuiWindowDockStyleCol_COUNT
  ! enum ImGuiDebugLogFlags_
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_None") :: ImGuiDebugLogFlags_None
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventActiveId") :: ImGuiDebugLogFlags_EventActiveId
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventFocus") :: ImGuiDebugLogFlags_EventFocus
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventPopup") :: ImGuiDebugLogFlags_EventPopup
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventNav") :: ImGuiDebugLogFlags_EventNav
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventIO") :: ImGuiDebugLogFlags_EventIO
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventDocking") :: ImGuiDebugLogFlags_EventDocking
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventViewport") :: ImGuiDebugLogFlags_EventViewport
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_EventMask_") :: ImGuiDebugLogFlags_EventMask_
  integer(c_int), bind(C,name="const_ImGuiDebugLogFlags_OutputToTTY") :: ImGuiDebugLogFlags_OutputToTTY
  ! enum ImGuiContextHookType
  integer(c_int), bind(C,name="const_ImGuiContextHookType_NewFramePre") :: ImGuiContextHookType_NewFramePre
  integer(c_int), bind(C,name="const_ImGuiContextHookType_NewFramePost") :: ImGuiContextHookType_NewFramePost
  integer(c_int), bind(C,name="const_ImGuiContextHookType_EndFramePre") :: ImGuiContextHookType_EndFramePre
  integer(c_int), bind(C,name="const_ImGuiContextHookType_EndFramePost") :: ImGuiContextHookType_EndFramePost
  integer(c_int), bind(C,name="const_ImGuiContextHookType_RenderPre") :: ImGuiContextHookType_RenderPre
  integer(c_int), bind(C,name="const_ImGuiContextHookType_RenderPost") :: ImGuiContextHookType_RenderPost
  integer(c_int), bind(C,name="const_ImGuiContextHookType_Shutdown") :: ImGuiContextHookType_Shutdown
  integer(c_int), bind(C,name="const_ImGuiContextHookType_PendingRemoval_") :: ImGuiContextHookType_PendingRemoval_
  ! enum ImGuiTabBarFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_DockNode") :: ImGuiTabBarFlags_DockNode
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_IsFocused") :: ImGuiTabBarFlags_IsFocused
  integer(c_int), bind(C,name="const_ImGuiTabBarFlags_SaveSettings") :: ImGuiTabBarFlags_SaveSettings
  ! enum ImGuiTabItemFlagsPrivate_
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_SectionMask_") :: ImGuiTabItemFlags_SectionMask_
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_NoCloseButton") :: ImGuiTabItemFlags_NoCloseButton
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_Button") :: ImGuiTabItemFlags_Button
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_Unsorted") :: ImGuiTabItemFlags_Unsorted
  integer(c_int), bind(C,name="const_ImGuiTabItemFlags_Preview") :: ImGuiTabItemFlags_Preview
  ! enum ImGuiFreeTypeBuilderFlags
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_NoHinting") :: ImGuiFreeTypeBuilderFlags_NoHinting
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_NoAutoHint") :: ImGuiFreeTypeBuilderFlags_NoAutoHint
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_ForceAutoHint") :: ImGuiFreeTypeBuilderFlags_ForceAutoHint
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_LightHinting") :: ImGuiFreeTypeBuilderFlags_LightHinting
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_MonoHinting") :: ImGuiFreeTypeBuilderFlags_MonoHinting
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_Bold") :: ImGuiFreeTypeBuilderFlags_Bold
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_Oblique") :: ImGuiFreeTypeBuilderFlags_Oblique
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_Monochrome") :: ImGuiFreeTypeBuilderFlags_Monochrome
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_LoadColor") :: ImGuiFreeTypeBuilderFlags_LoadColor
  integer(c_int), bind(C,name="const_ImGuiFreeTypeBuilderFlags_Bitmap") :: ImGuiFreeTypeBuilderFlags_Bitmap

  !! filedialog
  ! enum ImGuiFileDialogFlags
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_None") :: ImGuiFileDialogFlags_None
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_ConfirmOverwrite") :: &
     ImGuiFileDialogFlags_ConfirmOverwrite
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_DontShowHiddenFiles") :: &
     ImGuiFileDialogFlags_DontShowHiddenFiles
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_DisableCreateDirectoryButton") :: &
     ImGuiFileDialogFlags_DisableCreateDirectoryButton
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_HideColumnType") :: ImGuiFileDialogFlags_HideColumnType
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_HideColumnSize") :: ImGuiFileDialogFlags_HideColumnSize
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_HideColumnDate") :: ImGuiFileDialogFlags_HideColumnDate
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_NoDialog") :: ImGuiFileDialogFlags_NoDialog
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_ReadOnlyFileNameField") :: &
     ImGuiFileDialogFlags_ReadOnlyFileNameField
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_CaseInsensitiveExtention") :: &
     ImGuiFileDialogFlags_CaseInsensitiveExtention
  integer(c_int), bind(C,name="const_ImGuiFileDialogFlags_Modal") :: ImGuiFileDialogFlags_Modal
  ! enum IGFD_FileStyle
  integer(c_int), bind(C,name="const_IGFD_FileStyle_None") :: IGFD_FileStyle_None
  integer(c_int), bind(C,name="const_IGFD_FileStyleByTypeFile") :: IGFD_FileStyleByTypeFile
  integer(c_int), bind(C,name="const_IGFD_FileStyleByTypeDir") :: IGFD_FileStyleByTypeDir
  integer(c_int), bind(C,name="const_IGFD_FileStyleByTypeLink") :: IGFD_FileStyleByTypeLink
  integer(c_int), bind(C,name="const_IGFD_FileStyleByExtention") :: IGFD_FileStyleByExtention
  integer(c_int), bind(C,name="const_IGFD_FileStyleByFullName") :: IGFD_FileStyleByFullName
  integer(c_int), bind(C,name="const_IGFD_FileStyleByContainedInFullName") :: IGFD_FileStyleByContainedInFullName

  ! default font for the critic2 GUI (POINTER)
  type(c_ptr), bind(C,name="const_myfont_ttf_compressed_data_base85_ptr") :: myfont_ttf_compressed_data_base85_ptr
