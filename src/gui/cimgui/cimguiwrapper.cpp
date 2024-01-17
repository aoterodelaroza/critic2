// Copyright (c) 2015 Alberto Otero de la Roza
// <aoterodelaroza@gmail.com>,
// Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
// <victor@fluor.quimica.uniovi.es>.
//
// critic2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at
// your option) any later version.
//
// critic2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

#include "cimguiwrapper.h"

#include "../imgui/imgui_internal.h"
#include "../imgui/imgui.h"
#include "cimgui.h"
#include "../imgui/imgui_freetype.h"

// enum ImGuiWindowFlags_
const int const_ImGuiWindowFlags_None = ImGuiWindowFlags_None;
const int const_ImGuiWindowFlags_NoTitleBar = ImGuiWindowFlags_NoTitleBar;
const int const_ImGuiWindowFlags_NoResize = ImGuiWindowFlags_NoResize;
const int const_ImGuiWindowFlags_NoMove = ImGuiWindowFlags_NoMove;
const int const_ImGuiWindowFlags_NoScrollbar = ImGuiWindowFlags_NoScrollbar;
const int const_ImGuiWindowFlags_NoScrollWithMouse = ImGuiWindowFlags_NoScrollWithMouse;
const int const_ImGuiWindowFlags_NoCollapse = ImGuiWindowFlags_NoCollapse;
const int const_ImGuiWindowFlags_AlwaysAutoResize = ImGuiWindowFlags_AlwaysAutoResize;
const int const_ImGuiWindowFlags_NoBackground = ImGuiWindowFlags_NoBackground;
const int const_ImGuiWindowFlags_NoSavedSettings = ImGuiWindowFlags_NoSavedSettings;
const int const_ImGuiWindowFlags_NoMouseInputs = ImGuiWindowFlags_NoMouseInputs;
const int const_ImGuiWindowFlags_MenuBar = ImGuiWindowFlags_MenuBar;
const int const_ImGuiWindowFlags_HorizontalScrollbar = ImGuiWindowFlags_HorizontalScrollbar;
const int const_ImGuiWindowFlags_NoFocusOnAppearing = ImGuiWindowFlags_NoFocusOnAppearing;
const int const_ImGuiWindowFlags_NoBringToFrontOnFocus = ImGuiWindowFlags_NoBringToFrontOnFocus;
const int const_ImGuiWindowFlags_AlwaysVerticalScrollbar = ImGuiWindowFlags_AlwaysVerticalScrollbar;
const int const_ImGuiWindowFlags_AlwaysHorizontalScrollbar = ImGuiWindowFlags_AlwaysHorizontalScrollbar;
const int const_ImGuiWindowFlags_AlwaysUseWindowPadding = ImGuiWindowFlags_AlwaysUseWindowPadding;
const int const_ImGuiWindowFlags_NoNavInputs = ImGuiWindowFlags_NoNavInputs;
const int const_ImGuiWindowFlags_NoNavFocus = ImGuiWindowFlags_NoNavFocus;
const int const_ImGuiWindowFlags_UnsavedDocument = ImGuiWindowFlags_UnsavedDocument;
const int const_ImGuiWindowFlags_NoDocking = ImGuiWindowFlags_NoDocking;
const int const_ImGuiWindowFlags_NoNav = ImGuiWindowFlags_NoNav;
const int const_ImGuiWindowFlags_NoDecoration = ImGuiWindowFlags_NoDecoration;
const int const_ImGuiWindowFlags_NoInputs = ImGuiWindowFlags_NoInputs;
const int const_ImGuiWindowFlags_NavFlattened = ImGuiWindowFlags_NavFlattened;
const int const_ImGuiWindowFlags_ChildWindow = ImGuiWindowFlags_ChildWindow;
const int const_ImGuiWindowFlags_Tooltip = ImGuiWindowFlags_Tooltip;
const int const_ImGuiWindowFlags_Popup = ImGuiWindowFlags_Popup;
const int const_ImGuiWindowFlags_Modal = ImGuiWindowFlags_Modal;
const int const_ImGuiWindowFlags_ChildMenu = ImGuiWindowFlags_ChildMenu;
const int const_ImGuiWindowFlags_DockNodeHost = ImGuiWindowFlags_DockNodeHost;
// enum ImGuiInputTextFlags_
const int const_ImGuiInputTextFlags_None = ImGuiInputTextFlags_None;
const int const_ImGuiInputTextFlags_CharsDecimal = ImGuiInputTextFlags_CharsDecimal;
const int const_ImGuiInputTextFlags_CharsHexadecimal = ImGuiInputTextFlags_CharsHexadecimal;
const int const_ImGuiInputTextFlags_CharsUppercase = ImGuiInputTextFlags_CharsUppercase;
const int const_ImGuiInputTextFlags_CharsNoBlank = ImGuiInputTextFlags_CharsNoBlank;
const int const_ImGuiInputTextFlags_AutoSelectAll = ImGuiInputTextFlags_AutoSelectAll;
const int const_ImGuiInputTextFlags_EnterReturnsTrue = ImGuiInputTextFlags_EnterReturnsTrue;
const int const_ImGuiInputTextFlags_CallbackCompletion = ImGuiInputTextFlags_CallbackCompletion;
const int const_ImGuiInputTextFlags_CallbackHistory = ImGuiInputTextFlags_CallbackHistory;
const int const_ImGuiInputTextFlags_CallbackAlways = ImGuiInputTextFlags_CallbackAlways;
const int const_ImGuiInputTextFlags_CallbackCharFilter = ImGuiInputTextFlags_CallbackCharFilter;
const int const_ImGuiInputTextFlags_AllowTabInput = ImGuiInputTextFlags_AllowTabInput;
const int const_ImGuiInputTextFlags_CtrlEnterForNewLine = ImGuiInputTextFlags_CtrlEnterForNewLine;
const int const_ImGuiInputTextFlags_NoHorizontalScroll = ImGuiInputTextFlags_NoHorizontalScroll;
const int const_ImGuiInputTextFlags_AlwaysOverwrite = ImGuiInputTextFlags_AlwaysOverwrite;
const int const_ImGuiInputTextFlags_ReadOnly = ImGuiInputTextFlags_ReadOnly;
const int const_ImGuiInputTextFlags_Password = ImGuiInputTextFlags_Password;
const int const_ImGuiInputTextFlags_NoUndoRedo = ImGuiInputTextFlags_NoUndoRedo;
const int const_ImGuiInputTextFlags_CharsScientific = ImGuiInputTextFlags_CharsScientific;
const int const_ImGuiInputTextFlags_CallbackResize = ImGuiInputTextFlags_CallbackResize;
const int const_ImGuiInputTextFlags_CallbackEdit = ImGuiInputTextFlags_CallbackEdit;
// enum ImGuiTreeNodeFlags_;
const int const_ImGuiTreeNodeFlags_None = ImGuiTreeNodeFlags_None;
const int const_ImGuiTreeNodeFlags_Selected = ImGuiTreeNodeFlags_Selected;
const int const_ImGuiTreeNodeFlags_Framed = ImGuiTreeNodeFlags_Framed;
const int const_ImGuiTreeNodeFlags_AllowItemOverlap = ImGuiTreeNodeFlags_AllowItemOverlap;
const int const_ImGuiTreeNodeFlags_NoTreePushOnOpen = ImGuiTreeNodeFlags_NoTreePushOnOpen;
const int const_ImGuiTreeNodeFlags_NoAutoOpenOnLog = ImGuiTreeNodeFlags_NoAutoOpenOnLog;
const int const_ImGuiTreeNodeFlags_DefaultOpen = ImGuiTreeNodeFlags_DefaultOpen;
const int const_ImGuiTreeNodeFlags_OpenOnDoubleClick = ImGuiTreeNodeFlags_OpenOnDoubleClick;
const int const_ImGuiTreeNodeFlags_OpenOnArrow = ImGuiTreeNodeFlags_OpenOnArrow;
const int const_ImGuiTreeNodeFlags_Leaf = ImGuiTreeNodeFlags_Leaf;
const int const_ImGuiTreeNodeFlags_Bullet = ImGuiTreeNodeFlags_Bullet;
const int const_ImGuiTreeNodeFlags_FramePadding = ImGuiTreeNodeFlags_FramePadding;
const int const_ImGuiTreeNodeFlags_SpanAvailWidth = ImGuiTreeNodeFlags_SpanAvailWidth;
const int const_ImGuiTreeNodeFlags_SpanFullWidth = ImGuiTreeNodeFlags_SpanFullWidth;
const int const_ImGuiTreeNodeFlags_NavLeftJumpsBackHere = ImGuiTreeNodeFlags_NavLeftJumpsBackHere;
const int const_ImGuiTreeNodeFlags_CollapsingHeader = ImGuiTreeNodeFlags_CollapsingHeader;
// enum ImGuiPopupFlags_
const int const_ImGuiPopupFlags_None = ImGuiPopupFlags_None;
const int const_ImGuiPopupFlags_MouseButtonLeft = ImGuiPopupFlags_MouseButtonLeft;
const int const_ImGuiPopupFlags_MouseButtonRight = ImGuiPopupFlags_MouseButtonRight;
const int const_ImGuiPopupFlags_MouseButtonMiddle = ImGuiPopupFlags_MouseButtonMiddle;
const int const_ImGuiPopupFlags_MouseButtonMask_ = ImGuiPopupFlags_MouseButtonMask_;
const int const_ImGuiPopupFlags_MouseButtonDefault_ = ImGuiPopupFlags_MouseButtonDefault_;
const int const_ImGuiPopupFlags_NoOpenOverExistingPopup = ImGuiPopupFlags_NoOpenOverExistingPopup;
const int const_ImGuiPopupFlags_NoOpenOverItems = ImGuiPopupFlags_NoOpenOverItems;
const int const_ImGuiPopupFlags_AnyPopupId = ImGuiPopupFlags_AnyPopupId;
const int const_ImGuiPopupFlags_AnyPopupLevel = ImGuiPopupFlags_AnyPopupLevel;
const int const_ImGuiPopupFlags_AnyPopup = ImGuiPopupFlags_AnyPopup;
// enum ImGuiSelectableFlags_
const int const_ImGuiSelectableFlags_None = ImGuiSelectableFlags_None;
const int const_ImGuiSelectableFlags_DontClosePopups = ImGuiSelectableFlags_DontClosePopups;
const int const_ImGuiSelectableFlags_SpanAllColumns = ImGuiSelectableFlags_SpanAllColumns;
const int const_ImGuiSelectableFlags_AllowDoubleClick = ImGuiSelectableFlags_AllowDoubleClick;
const int const_ImGuiSelectableFlags_Disabled = ImGuiSelectableFlags_Disabled;
const int const_ImGuiSelectableFlags_AllowItemOverlap = ImGuiSelectableFlags_AllowItemOverlap;
// enum ImGuiComboFlags_
const int const_ImGuiComboFlags_None = ImGuiComboFlags_None;
const int const_ImGuiComboFlags_PopupAlignLeft = ImGuiComboFlags_PopupAlignLeft;
const int const_ImGuiComboFlags_HeightSmall = ImGuiComboFlags_HeightSmall;
const int const_ImGuiComboFlags_HeightRegular = ImGuiComboFlags_HeightRegular;
const int const_ImGuiComboFlags_HeightLarge = ImGuiComboFlags_HeightLarge;
const int const_ImGuiComboFlags_HeightLargest = ImGuiComboFlags_HeightLargest;
const int const_ImGuiComboFlags_NoArrowButton = ImGuiComboFlags_NoArrowButton;
const int const_ImGuiComboFlags_NoPreview = ImGuiComboFlags_NoPreview;
const int const_ImGuiComboFlags_HeightMask_ = ImGuiComboFlags_HeightMask_;
// enum ImGuiTabBarFlags_
const int const_ImGuiTabBarFlags_None = ImGuiTabBarFlags_None;
const int const_ImGuiTabBarFlags_Reorderable = ImGuiTabBarFlags_Reorderable;
const int const_ImGuiTabBarFlags_AutoSelectNewTabs = ImGuiTabBarFlags_AutoSelectNewTabs;
const int const_ImGuiTabBarFlags_TabListPopupButton = ImGuiTabBarFlags_TabListPopupButton;
const int const_ImGuiTabBarFlags_NoCloseWithMiddleMouseButton = ImGuiTabBarFlags_NoCloseWithMiddleMouseButton;
const int const_ImGuiTabBarFlags_NoTabListScrollingButtons = ImGuiTabBarFlags_NoTabListScrollingButtons;
const int const_ImGuiTabBarFlags_NoTooltip = ImGuiTabBarFlags_NoTooltip;
const int const_ImGuiTabBarFlags_FittingPolicyResizeDown = ImGuiTabBarFlags_FittingPolicyResizeDown;
const int const_ImGuiTabBarFlags_FittingPolicyScroll = ImGuiTabBarFlags_FittingPolicyScroll;
const int const_ImGuiTabBarFlags_FittingPolicyMask_ = ImGuiTabBarFlags_FittingPolicyMask_;
const int const_ImGuiTabBarFlags_FittingPolicyDefault_ = ImGuiTabBarFlags_FittingPolicyDefault_;
// enum ImGuiTabItemFlags_
const int const_ImGuiTabItemFlags_None = ImGuiTabItemFlags_None;
const int const_ImGuiTabItemFlags_UnsavedDocument = ImGuiTabItemFlags_UnsavedDocument;
const int const_ImGuiTabItemFlags_SetSelected = ImGuiTabItemFlags_SetSelected;
const int const_ImGuiTabItemFlags_NoCloseWithMiddleMouseButton = ImGuiTabItemFlags_NoCloseWithMiddleMouseButton;
const int const_ImGuiTabItemFlags_NoPushId = ImGuiTabItemFlags_NoPushId;
const int const_ImGuiTabItemFlags_NoTooltip = ImGuiTabItemFlags_NoTooltip;
const int const_ImGuiTabItemFlags_NoReorder = ImGuiTabItemFlags_NoReorder;
const int const_ImGuiTabItemFlags_Leading = ImGuiTabItemFlags_Leading;
const int const_ImGuiTabItemFlags_Trailing = ImGuiTabItemFlags_Trailing;
// enum ImGuiTableFlags_
const int const_ImGuiTableFlags_None = ImGuiTableFlags_None;
const int const_ImGuiTableFlags_Resizable = ImGuiTableFlags_Resizable;
const int const_ImGuiTableFlags_Reorderable = ImGuiTableFlags_Reorderable;
const int const_ImGuiTableFlags_Hideable = ImGuiTableFlags_Hideable;
const int const_ImGuiTableFlags_Sortable = ImGuiTableFlags_Sortable;
const int const_ImGuiTableFlags_NoSavedSettings = ImGuiTableFlags_NoSavedSettings;
const int const_ImGuiTableFlags_ContextMenuInBody = ImGuiTableFlags_ContextMenuInBody;
const int const_ImGuiTableFlags_RowBg = ImGuiTableFlags_RowBg;
const int const_ImGuiTableFlags_BordersInnerH = ImGuiTableFlags_BordersInnerH;
const int const_ImGuiTableFlags_BordersOuterH = ImGuiTableFlags_BordersOuterH;
const int const_ImGuiTableFlags_BordersInnerV = ImGuiTableFlags_BordersInnerV;
const int const_ImGuiTableFlags_BordersOuterV = ImGuiTableFlags_BordersOuterV;
const int const_ImGuiTableFlags_BordersH = ImGuiTableFlags_BordersH;
const int const_ImGuiTableFlags_BordersV = ImGuiTableFlags_BordersV;
const int const_ImGuiTableFlags_BordersInner = ImGuiTableFlags_BordersInner;
const int const_ImGuiTableFlags_BordersOuter = ImGuiTableFlags_BordersOuter;
const int const_ImGuiTableFlags_Borders = ImGuiTableFlags_Borders;
const int const_ImGuiTableFlags_NoBordersInBody = ImGuiTableFlags_NoBordersInBody;
const int const_ImGuiTableFlags_NoBordersInBodyUntilResize = ImGuiTableFlags_NoBordersInBodyUntilResize;
const int const_ImGuiTableFlags_SizingFixedFit = ImGuiTableFlags_SizingFixedFit;
const int const_ImGuiTableFlags_SizingFixedSame = ImGuiTableFlags_SizingFixedSame;
const int const_ImGuiTableFlags_SizingStretchProp = ImGuiTableFlags_SizingStretchProp;
const int const_ImGuiTableFlags_SizingStretchSame = ImGuiTableFlags_SizingStretchSame;
const int const_ImGuiTableFlags_NoHostExtendX = ImGuiTableFlags_NoHostExtendX;
const int const_ImGuiTableFlags_NoHostExtendY = ImGuiTableFlags_NoHostExtendY;
const int const_ImGuiTableFlags_NoKeepColumnsVisible = ImGuiTableFlags_NoKeepColumnsVisible;
const int const_ImGuiTableFlags_PreciseWidths = ImGuiTableFlags_PreciseWidths;
const int const_ImGuiTableFlags_NoClip = ImGuiTableFlags_NoClip;
const int const_ImGuiTableFlags_PadOuterX = ImGuiTableFlags_PadOuterX;
const int const_ImGuiTableFlags_NoPadOuterX = ImGuiTableFlags_NoPadOuterX;
const int const_ImGuiTableFlags_NoPadInnerX = ImGuiTableFlags_NoPadInnerX;
const int const_ImGuiTableFlags_ScrollX = ImGuiTableFlags_ScrollX;
const int const_ImGuiTableFlags_ScrollY = ImGuiTableFlags_ScrollY;
const int const_ImGuiTableFlags_SortMulti = ImGuiTableFlags_SortMulti;
const int const_ImGuiTableFlags_SortTristate = ImGuiTableFlags_SortTristate;
const int const_ImGuiTableFlags_SizingMask_ = ImGuiTableFlags_SizingMask_;
// enum ImGuiTableColumnFlags_
const int const_ImGuiTableColumnFlags_None = ImGuiTableColumnFlags_None;
const int const_ImGuiTableColumnFlags_Disabled = ImGuiTableColumnFlags_Disabled;
const int const_ImGuiTableColumnFlags_DefaultHide = ImGuiTableColumnFlags_DefaultHide;
const int const_ImGuiTableColumnFlags_DefaultSort = ImGuiTableColumnFlags_DefaultSort;
const int const_ImGuiTableColumnFlags_WidthStretch = ImGuiTableColumnFlags_WidthStretch;
const int const_ImGuiTableColumnFlags_WidthFixed = ImGuiTableColumnFlags_WidthFixed;
const int const_ImGuiTableColumnFlags_NoResize = ImGuiTableColumnFlags_NoResize;
const int const_ImGuiTableColumnFlags_NoReorder = ImGuiTableColumnFlags_NoReorder;
const int const_ImGuiTableColumnFlags_NoHide = ImGuiTableColumnFlags_NoHide;
const int const_ImGuiTableColumnFlags_NoClip = ImGuiTableColumnFlags_NoClip;
const int const_ImGuiTableColumnFlags_NoSort = ImGuiTableColumnFlags_NoSort;
const int const_ImGuiTableColumnFlags_NoSortAscending = ImGuiTableColumnFlags_NoSortAscending;
const int const_ImGuiTableColumnFlags_NoSortDescending = ImGuiTableColumnFlags_NoSortDescending;
const int const_ImGuiTableColumnFlags_NoHeaderLabel = ImGuiTableColumnFlags_NoHeaderLabel;
const int const_ImGuiTableColumnFlags_NoHeaderWidth = ImGuiTableColumnFlags_NoHeaderWidth;
const int const_ImGuiTableColumnFlags_PreferSortAscending = ImGuiTableColumnFlags_PreferSortAscending;
const int const_ImGuiTableColumnFlags_PreferSortDescending = ImGuiTableColumnFlags_PreferSortDescending;
const int const_ImGuiTableColumnFlags_IndentEnable = ImGuiTableColumnFlags_IndentEnable;
const int const_ImGuiTableColumnFlags_IndentDisable = ImGuiTableColumnFlags_IndentDisable;
const int const_ImGuiTableColumnFlags_IsEnabled = ImGuiTableColumnFlags_IsEnabled;
const int const_ImGuiTableColumnFlags_IsVisible = ImGuiTableColumnFlags_IsVisible;
const int const_ImGuiTableColumnFlags_IsSorted = ImGuiTableColumnFlags_IsSorted;
const int const_ImGuiTableColumnFlags_IsHovered = ImGuiTableColumnFlags_IsHovered;
const int const_ImGuiTableColumnFlags_WidthMask_ = ImGuiTableColumnFlags_WidthMask_;
const int const_ImGuiTableColumnFlags_IndentMask_ = ImGuiTableColumnFlags_IndentMask_;
const int const_ImGuiTableColumnFlags_StatusMask_ = ImGuiTableColumnFlags_StatusMask_;
const int const_ImGuiTableColumnFlags_NoDirectResize_ = ImGuiTableColumnFlags_NoDirectResize_;
// enum ImGuiTableRowFlags_
const int const_ImGuiTableRowFlags_None = ImGuiTableRowFlags_None;
const int const_ImGuiTableRowFlags_Headers = ImGuiTableRowFlags_Headers;
// enum ImGuiTableBgTarget_
const int const_ImGuiTableBgTarget_None = ImGuiTableBgTarget_None;
const int const_ImGuiTableBgTarget_RowBg0 = ImGuiTableBgTarget_RowBg0;
const int const_ImGuiTableBgTarget_RowBg1 = ImGuiTableBgTarget_RowBg1;
const int const_ImGuiTableBgTarget_CellBg = ImGuiTableBgTarget_CellBg;
// enum ImGuiFocusedFlags_
const int const_ImGuiFocusedFlags_None = ImGuiFocusedFlags_None;
const int const_ImGuiFocusedFlags_ChildWindows = ImGuiFocusedFlags_ChildWindows;
const int const_ImGuiFocusedFlags_RootWindow = ImGuiFocusedFlags_RootWindow;
const int const_ImGuiFocusedFlags_AnyWindow = ImGuiFocusedFlags_AnyWindow;
const int const_ImGuiFocusedFlags_NoPopupHierarchy = ImGuiFocusedFlags_NoPopupHierarchy;
const int const_ImGuiFocusedFlags_DockHierarchy = ImGuiFocusedFlags_DockHierarchy;
const int const_ImGuiFocusedFlags_RootAndChildWindows = ImGuiFocusedFlags_RootAndChildWindows;
// enum ImGuiHoveredFlags_
const int const_ImGuiHoveredFlags_None = ImGuiHoveredFlags_None;
const int const_ImGuiHoveredFlags_ChildWindows = ImGuiHoveredFlags_ChildWindows;
const int const_ImGuiHoveredFlags_RootWindow = ImGuiHoveredFlags_RootWindow;
const int const_ImGuiHoveredFlags_AnyWindow = ImGuiHoveredFlags_AnyWindow;
const int const_ImGuiHoveredFlags_NoPopupHierarchy = ImGuiHoveredFlags_NoPopupHierarchy;
const int const_ImGuiHoveredFlags_DockHierarchy = ImGuiHoveredFlags_DockHierarchy;
const int const_ImGuiHoveredFlags_AllowWhenBlockedByPopup = ImGuiHoveredFlags_AllowWhenBlockedByPopup;
const int const_ImGuiHoveredFlags_AllowWhenBlockedByActiveItem = ImGuiHoveredFlags_AllowWhenBlockedByActiveItem;
const int const_ImGuiHoveredFlags_AllowWhenOverlapped = ImGuiHoveredFlags_AllowWhenOverlapped;
const int const_ImGuiHoveredFlags_AllowWhenDisabled = ImGuiHoveredFlags_AllowWhenDisabled;
const int const_ImGuiHoveredFlags_NoNavOverride = ImGuiHoveredFlags_NoNavOverride;
const int const_ImGuiHoveredFlags_RectOnly = ImGuiHoveredFlags_RectOnly;
const int const_ImGuiHoveredFlags_RootAndChildWindows = ImGuiHoveredFlags_RootAndChildWindows;
// enum ImGuiDockNodeFlags_
const int const_ImGuiDockNodeFlags_None = ImGuiDockNodeFlags_None;
const int const_ImGuiDockNodeFlags_KeepAliveOnly = ImGuiDockNodeFlags_KeepAliveOnly;
const int const_ImGuiDockNodeFlags_NoDockingInCentralNode = ImGuiDockNodeFlags_NoDockingInCentralNode;
const int const_ImGuiDockNodeFlags_PassthruCentralNode = ImGuiDockNodeFlags_PassthruCentralNode;
const int const_ImGuiDockNodeFlags_NoSplit = ImGuiDockNodeFlags_NoSplit;
const int const_ImGuiDockNodeFlags_NoResize = ImGuiDockNodeFlags_NoResize;
const int const_ImGuiDockNodeFlags_AutoHideTabBar = ImGuiDockNodeFlags_AutoHideTabBar;
// enum ImGuiDragDropFlags_
const int const_ImGuiDragDropFlags_None = ImGuiDragDropFlags_None;
const int const_ImGuiDragDropFlags_SourceNoPreviewTooltip = ImGuiDragDropFlags_SourceNoPreviewTooltip;
const int const_ImGuiDragDropFlags_SourceNoDisableHover = ImGuiDragDropFlags_SourceNoDisableHover;
const int const_ImGuiDragDropFlags_SourceNoHoldToOpenOthers = ImGuiDragDropFlags_SourceNoHoldToOpenOthers;
const int const_ImGuiDragDropFlags_SourceAllowNullID = ImGuiDragDropFlags_SourceAllowNullID;
const int const_ImGuiDragDropFlags_SourceExtern = ImGuiDragDropFlags_SourceExtern;
const int const_ImGuiDragDropFlags_SourceAutoExpirePayload = ImGuiDragDropFlags_SourceAutoExpirePayload;
const int const_ImGuiDragDropFlags_AcceptBeforeDelivery = ImGuiDragDropFlags_AcceptBeforeDelivery;
const int const_ImGuiDragDropFlags_AcceptNoDrawDefaultRect = ImGuiDragDropFlags_AcceptNoDrawDefaultRect;
const int const_ImGuiDragDropFlags_AcceptNoPreviewTooltip = ImGuiDragDropFlags_AcceptNoPreviewTooltip;
const int const_ImGuiDragDropFlags_AcceptPeekOnly = ImGuiDragDropFlags_AcceptPeekOnly;
// enum ImGuiDataType_
const int const_ImGuiDataType_S8 = ImGuiDataType_S8;
const int const_ImGuiDataType_U8 = ImGuiDataType_U8;
const int const_ImGuiDataType_S16 = ImGuiDataType_S16;
const int const_ImGuiDataType_U16 = ImGuiDataType_U16;
const int const_ImGuiDataType_S32 = ImGuiDataType_S32;
const int const_ImGuiDataType_U32 = ImGuiDataType_U32;
const int const_ImGuiDataType_S64 = ImGuiDataType_S64;
const int const_ImGuiDataType_U64 = ImGuiDataType_U64;
const int const_ImGuiDataType_Float = ImGuiDataType_Float;
const int const_ImGuiDataType_Double = ImGuiDataType_Double;
const int const_ImGuiDataType_COUNT = ImGuiDataType_COUNT;
// enum ImGuiDir_
const int const_ImGuiDir_None = ImGuiDir_None;
const int const_ImGuiDir_Left = ImGuiDir_Left;
const int const_ImGuiDir_Right = ImGuiDir_Right;
const int const_ImGuiDir_Up = ImGuiDir_Up;
const int const_ImGuiDir_Down = ImGuiDir_Down;
const int const_ImGuiDir_COUNT = ImGuiDir_COUNT;
// enum ImGuiSortDirection_
const int const_ImGuiSortDirection_None = ImGuiSortDirection_None;
const int const_ImGuiSortDirection_Ascending = ImGuiSortDirection_Ascending;
const int const_ImGuiSortDirection_Descending = ImGuiSortDirection_Descending;
// enum ImGuiKey_
const int const_ImGuiKey_None = ImGuiKey_None;
const int const_ImGuiKey_Tab = ImGuiKey_Tab;
const int const_ImGuiKey_LeftArrow = ImGuiKey_LeftArrow;
const int const_ImGuiKey_RightArrow = ImGuiKey_RightArrow;
const int const_ImGuiKey_UpArrow = ImGuiKey_UpArrow;
const int const_ImGuiKey_DownArrow = ImGuiKey_DownArrow;
const int const_ImGuiKey_PageUp = ImGuiKey_PageUp;
const int const_ImGuiKey_PageDown = ImGuiKey_PageDown;
const int const_ImGuiKey_Home = ImGuiKey_Home;
const int const_ImGuiKey_End = ImGuiKey_End;
const int const_ImGuiKey_Insert = ImGuiKey_Insert;
const int const_ImGuiKey_Delete = ImGuiKey_Delete;
const int const_ImGuiKey_Backspace = ImGuiKey_Backspace;
const int const_ImGuiKey_Space = ImGuiKey_Space;
const int const_ImGuiKey_Enter = ImGuiKey_Enter;
const int const_ImGuiKey_Escape = ImGuiKey_Escape;
const int const_ImGuiKey_LeftCtrl = ImGuiKey_LeftCtrl;
const int const_ImGuiKey_LeftShift = ImGuiKey_LeftShift;
const int const_ImGuiKey_LeftAlt = ImGuiKey_LeftAlt;
const int const_ImGuiKey_LeftSuper = ImGuiKey_LeftSuper;
const int const_ImGuiKey_RightCtrl = ImGuiKey_RightCtrl;
const int const_ImGuiKey_RightShift = ImGuiKey_RightShift;
const int const_ImGuiKey_RightAlt = ImGuiKey_RightAlt;
const int const_ImGuiKey_RightSuper = ImGuiKey_RightSuper;
const int const_ImGuiKey_Menu = ImGuiKey_Menu;
const int const_ImGuiKey_0 = ImGuiKey_0;
const int const_ImGuiKey_1 = ImGuiKey_1;
const int const_ImGuiKey_2 = ImGuiKey_2;
const int const_ImGuiKey_3 = ImGuiKey_3;
const int const_ImGuiKey_4 = ImGuiKey_4;
const int const_ImGuiKey_5 = ImGuiKey_5;
const int const_ImGuiKey_6 = ImGuiKey_6;
const int const_ImGuiKey_7 = ImGuiKey_7;
const int const_ImGuiKey_8 = ImGuiKey_8;
const int const_ImGuiKey_9 = ImGuiKey_9;
const int const_ImGuiKey_A = ImGuiKey_A;
const int const_ImGuiKey_B = ImGuiKey_B;
const int const_ImGuiKey_C = ImGuiKey_C;
const int const_ImGuiKey_D = ImGuiKey_D;
const int const_ImGuiKey_E = ImGuiKey_E;
const int const_ImGuiKey_F = ImGuiKey_F;
const int const_ImGuiKey_G = ImGuiKey_G;
const int const_ImGuiKey_H = ImGuiKey_H;
const int const_ImGuiKey_I = ImGuiKey_I;
const int const_ImGuiKey_J = ImGuiKey_J;
const int const_ImGuiKey_K = ImGuiKey_K;
const int const_ImGuiKey_L = ImGuiKey_L;
const int const_ImGuiKey_M = ImGuiKey_M;
const int const_ImGuiKey_N = ImGuiKey_N;
const int const_ImGuiKey_O = ImGuiKey_O;
const int const_ImGuiKey_P = ImGuiKey_P;
const int const_ImGuiKey_Q = ImGuiKey_Q;
const int const_ImGuiKey_R = ImGuiKey_R;
const int const_ImGuiKey_S = ImGuiKey_S;
const int const_ImGuiKey_T = ImGuiKey_T;
const int const_ImGuiKey_U = ImGuiKey_U;
const int const_ImGuiKey_V = ImGuiKey_V;
const int const_ImGuiKey_W = ImGuiKey_W;
const int const_ImGuiKey_X = ImGuiKey_X;
const int const_ImGuiKey_Y = ImGuiKey_Y;
const int const_ImGuiKey_Z = ImGuiKey_Z;
const int const_ImGuiKey_F1 = ImGuiKey_F1;
const int const_ImGuiKey_F2 = ImGuiKey_F2;
const int const_ImGuiKey_F3 = ImGuiKey_F3;
const int const_ImGuiKey_F4 = ImGuiKey_F4;
const int const_ImGuiKey_F5 = ImGuiKey_F5;
const int const_ImGuiKey_F6 = ImGuiKey_F6;
const int const_ImGuiKey_F7 = ImGuiKey_F7;
const int const_ImGuiKey_F8 = ImGuiKey_F8;
const int const_ImGuiKey_F9 = ImGuiKey_F9;
const int const_ImGuiKey_F10 = ImGuiKey_F10;
const int const_ImGuiKey_F11 = ImGuiKey_F11;
const int const_ImGuiKey_F12 = ImGuiKey_F12;
const int const_ImGuiKey_Apostrophe = ImGuiKey_Apostrophe;
const int const_ImGuiKey_Comma = ImGuiKey_Comma;
const int const_ImGuiKey_Minus = ImGuiKey_Minus;
const int const_ImGuiKey_Period = ImGuiKey_Period;
const int const_ImGuiKey_Slash = ImGuiKey_Slash;
const int const_ImGuiKey_Semicolon = ImGuiKey_Semicolon;
const int const_ImGuiKey_Equal = ImGuiKey_Equal;
const int const_ImGuiKey_LeftBracket = ImGuiKey_LeftBracket;
const int const_ImGuiKey_Backslash = ImGuiKey_Backslash;
const int const_ImGuiKey_RightBracket = ImGuiKey_RightBracket;
const int const_ImGuiKey_GraveAccent = ImGuiKey_GraveAccent;
const int const_ImGuiKey_CapsLock = ImGuiKey_CapsLock;
const int const_ImGuiKey_ScrollLock = ImGuiKey_ScrollLock;
const int const_ImGuiKey_NumLock = ImGuiKey_NumLock;
const int const_ImGuiKey_PrintScreen = ImGuiKey_PrintScreen;
const int const_ImGuiKey_Pause = ImGuiKey_Pause;
const int const_ImGuiKey_Keypad0 = ImGuiKey_Keypad0;
const int const_ImGuiKey_Keypad1 = ImGuiKey_Keypad1;
const int const_ImGuiKey_Keypad2 = ImGuiKey_Keypad2;
const int const_ImGuiKey_Keypad3 = ImGuiKey_Keypad3;
const int const_ImGuiKey_Keypad4 = ImGuiKey_Keypad4;
const int const_ImGuiKey_Keypad5 = ImGuiKey_Keypad5;
const int const_ImGuiKey_Keypad6 = ImGuiKey_Keypad6;
const int const_ImGuiKey_Keypad7 = ImGuiKey_Keypad7;
const int const_ImGuiKey_Keypad8 = ImGuiKey_Keypad8;
const int const_ImGuiKey_Keypad9 = ImGuiKey_Keypad9;
const int const_ImGuiKey_KeypadDecimal = ImGuiKey_KeypadDecimal;
const int const_ImGuiKey_KeypadDivide = ImGuiKey_KeypadDivide;
const int const_ImGuiKey_KeypadMultiply = ImGuiKey_KeypadMultiply;
const int const_ImGuiKey_KeypadSubtract = ImGuiKey_KeypadSubtract;
const int const_ImGuiKey_KeypadAdd = ImGuiKey_KeypadAdd;
const int const_ImGuiKey_KeypadEnter = ImGuiKey_KeypadEnter;
const int const_ImGuiKey_KeypadEqual = ImGuiKey_KeypadEqual;
const int const_ImGuiKey_GamepadStart = ImGuiKey_GamepadStart;
const int const_ImGuiKey_GamepadBack = ImGuiKey_GamepadBack;
const int const_ImGuiKey_GamepadFaceUp = ImGuiKey_GamepadFaceUp;
const int const_ImGuiKey_GamepadFaceDown = ImGuiKey_GamepadFaceDown;
const int const_ImGuiKey_GamepadFaceLeft = ImGuiKey_GamepadFaceLeft;
const int const_ImGuiKey_GamepadFaceRight = ImGuiKey_GamepadFaceRight;
const int const_ImGuiKey_GamepadDpadUp = ImGuiKey_GamepadDpadUp;
const int const_ImGuiKey_GamepadDpadDown = ImGuiKey_GamepadDpadDown;
const int const_ImGuiKey_GamepadDpadLeft = ImGuiKey_GamepadDpadLeft;
const int const_ImGuiKey_GamepadDpadRight = ImGuiKey_GamepadDpadRight;
const int const_ImGuiKey_GamepadL1 = ImGuiKey_GamepadL1;
const int const_ImGuiKey_GamepadR1 = ImGuiKey_GamepadR1;
const int const_ImGuiKey_GamepadL2 = ImGuiKey_GamepadL2;
const int const_ImGuiKey_GamepadR2 = ImGuiKey_GamepadR2;
const int const_ImGuiKey_GamepadL3 = ImGuiKey_GamepadL3;
const int const_ImGuiKey_GamepadR3 = ImGuiKey_GamepadR3;
const int const_ImGuiKey_GamepadLStickUp = ImGuiKey_GamepadLStickUp;
const int const_ImGuiKey_GamepadLStickDown = ImGuiKey_GamepadLStickDown;
const int const_ImGuiKey_GamepadLStickLeft = ImGuiKey_GamepadLStickLeft;
const int const_ImGuiKey_GamepadLStickRight = ImGuiKey_GamepadLStickRight;
const int const_ImGuiKey_GamepadRStickUp = ImGuiKey_GamepadRStickUp;
const int const_ImGuiKey_GamepadRStickDown = ImGuiKey_GamepadRStickDown;
const int const_ImGuiKey_GamepadRStickLeft = ImGuiKey_GamepadRStickLeft;
const int const_ImGuiKey_GamepadRStickRight = ImGuiKey_GamepadRStickRight;
const int const_ImGuiKey_ModCtrl = ImGuiKey_ModCtrl;
const int const_ImGuiKey_ModShift = ImGuiKey_ModShift;
const int const_ImGuiKey_ModAlt = ImGuiKey_ModAlt;
const int const_ImGuiKey_ModSuper = ImGuiKey_ModSuper;
// ImGuiKey_COUNT
const int const_ImGuiKey_NamedKey_BEGIN = ImGuiKey_NamedKey_BEGIN;
const int const_ImGuiKey_NamedKey_END = ImGuiKey_NamedKey_END;
// const int const_ImGuiKey_NamedKey_COUNT = ImGuiKey_NamedKey_COUNT;
// ImGuiKey_KeysData_SIZE
const int const_ImGuiKey_KeysData_OFFSET = ImGuiKey_KeysData_OFFSET;
// enum ImGuiModFlags_
const int const_ImGuiModFlags_None = ImGuiModFlags_None;
const int const_ImGuiModFlags_Ctrl = ImGuiModFlags_Ctrl;
const int const_ImGuiModFlags_Shift = ImGuiModFlags_Shift;
const int const_ImGuiModFlags_Alt = ImGuiModFlags_Alt;
const int const_ImGuiModFlags_Super = ImGuiModFlags_Super;
// enum ImGuiNavInput_
const int const_ImGuiNavInput_Activate = ImGuiNavInput_Activate;
const int const_ImGuiNavInput_Cancel = ImGuiNavInput_Cancel;
const int const_ImGuiNavInput_Input = ImGuiNavInput_Input;
const int const_ImGuiNavInput_Menu = ImGuiNavInput_Menu;
const int const_ImGuiNavInput_DpadLeft = ImGuiNavInput_DpadLeft;
const int const_ImGuiNavInput_DpadRight = ImGuiNavInput_DpadRight;
const int const_ImGuiNavInput_DpadUp = ImGuiNavInput_DpadUp;
const int const_ImGuiNavInput_DpadDown = ImGuiNavInput_DpadDown;
const int const_ImGuiNavInput_LStickLeft = ImGuiNavInput_LStickLeft;
const int const_ImGuiNavInput_LStickRight = ImGuiNavInput_LStickRight;
const int const_ImGuiNavInput_LStickUp = ImGuiNavInput_LStickUp;
const int const_ImGuiNavInput_LStickDown = ImGuiNavInput_LStickDown;
const int const_ImGuiNavInput_FocusPrev = ImGuiNavInput_FocusPrev;
const int const_ImGuiNavInput_FocusNext = ImGuiNavInput_FocusNext;
const int const_ImGuiNavInput_TweakSlow = ImGuiNavInput_TweakSlow;
const int const_ImGuiNavInput_TweakFast = ImGuiNavInput_TweakFast;
const int const_ImGuiNavInput_KeyLeft_ = ImGuiNavInput_KeyLeft_;
const int const_ImGuiNavInput_KeyRight_ = ImGuiNavInput_KeyRight_;
const int const_ImGuiNavInput_KeyUp_ = ImGuiNavInput_KeyUp_;
const int const_ImGuiNavInput_KeyDown_ = ImGuiNavInput_KeyDown_;
// ImGuiNavInput_COUNT
// enum ImGuiConfigFlags_
const int const_ImGuiConfigFlags_None = ImGuiConfigFlags_None;
const int const_ImGuiConfigFlags_NavEnableKeyboard = ImGuiConfigFlags_NavEnableKeyboard;
const int const_ImGuiConfigFlags_NavEnableGamepad = ImGuiConfigFlags_NavEnableGamepad;
const int const_ImGuiConfigFlags_NavEnableSetMousePos = ImGuiConfigFlags_NavEnableSetMousePos;
const int const_ImGuiConfigFlags_NavNoCaptureKeyboard = ImGuiConfigFlags_NavNoCaptureKeyboard;
const int const_ImGuiConfigFlags_NoMouse = ImGuiConfigFlags_NoMouse;
const int const_ImGuiConfigFlags_NoMouseCursorChange = ImGuiConfigFlags_NoMouseCursorChange;
const int const_ImGuiConfigFlags_DockingEnable = ImGuiConfigFlags_DockingEnable;
const int const_ImGuiConfigFlags_ViewportsEnable = ImGuiConfigFlags_ViewportsEnable;
const int const_ImGuiConfigFlags_DpiEnableScaleViewports = ImGuiConfigFlags_DpiEnableScaleViewports;
const int const_ImGuiConfigFlags_DpiEnableScaleFonts = ImGuiConfigFlags_DpiEnableScaleFonts;
const int const_ImGuiConfigFlags_IsSRGB = ImGuiConfigFlags_IsSRGB;
const int const_ImGuiConfigFlags_IsTouchScreen = ImGuiConfigFlags_IsTouchScreen;
// enum ImGuiBackendFlags_
const int const_ImGuiBackendFlags_None = ImGuiBackendFlags_None;
const int const_ImGuiBackendFlags_HasGamepad = ImGuiBackendFlags_HasGamepad;
const int const_ImGuiBackendFlags_HasMouseCursors = ImGuiBackendFlags_HasMouseCursors;
const int const_ImGuiBackendFlags_HasSetMousePos = ImGuiBackendFlags_HasSetMousePos;
const int const_ImGuiBackendFlags_RendererHasVtxOffset = ImGuiBackendFlags_RendererHasVtxOffset;
const int const_ImGuiBackendFlags_PlatformHasViewports = ImGuiBackendFlags_PlatformHasViewports;
const int const_ImGuiBackendFlags_HasMouseHoveredViewport = ImGuiBackendFlags_HasMouseHoveredViewport;
const int const_ImGuiBackendFlags_RendererHasViewports = ImGuiBackendFlags_RendererHasViewports;
// enum ImGuiCol_
const int const_ImGuiCol_Text = ImGuiCol_Text;
const int const_ImGuiCol_TextDisabled = ImGuiCol_TextDisabled;
const int const_ImGuiCol_WindowBg = ImGuiCol_WindowBg;
const int const_ImGuiCol_ChildBg = ImGuiCol_ChildBg;
const int const_ImGuiCol_PopupBg = ImGuiCol_PopupBg;
const int const_ImGuiCol_Border = ImGuiCol_Border;
const int const_ImGuiCol_BorderShadow = ImGuiCol_BorderShadow;
const int const_ImGuiCol_FrameBg = ImGuiCol_FrameBg;
const int const_ImGuiCol_FrameBgHovered = ImGuiCol_FrameBgHovered;
const int const_ImGuiCol_FrameBgActive = ImGuiCol_FrameBgActive;
const int const_ImGuiCol_TitleBg = ImGuiCol_TitleBg;
const int const_ImGuiCol_TitleBgActive = ImGuiCol_TitleBgActive;
const int const_ImGuiCol_TitleBgCollapsed = ImGuiCol_TitleBgCollapsed;
const int const_ImGuiCol_MenuBarBg = ImGuiCol_MenuBarBg;
const int const_ImGuiCol_ScrollbarBg = ImGuiCol_ScrollbarBg;
const int const_ImGuiCol_ScrollbarGrab = ImGuiCol_ScrollbarGrab;
const int const_ImGuiCol_ScrollbarGrabHovered = ImGuiCol_ScrollbarGrabHovered;
const int const_ImGuiCol_ScrollbarGrabActive = ImGuiCol_ScrollbarGrabActive;
const int const_ImGuiCol_CheckMark = ImGuiCol_CheckMark;
const int const_ImGuiCol_SliderGrab = ImGuiCol_SliderGrab;
const int const_ImGuiCol_SliderGrabActive = ImGuiCol_SliderGrabActive;
const int const_ImGuiCol_Button = ImGuiCol_Button;
const int const_ImGuiCol_ButtonHovered = ImGuiCol_ButtonHovered;
const int const_ImGuiCol_ButtonActive = ImGuiCol_ButtonActive;
const int const_ImGuiCol_Header = ImGuiCol_Header;
const int const_ImGuiCol_HeaderHovered = ImGuiCol_HeaderHovered;
const int const_ImGuiCol_HeaderActive = ImGuiCol_HeaderActive;
const int const_ImGuiCol_Separator = ImGuiCol_Separator;
const int const_ImGuiCol_SeparatorHovered = ImGuiCol_SeparatorHovered;
const int const_ImGuiCol_SeparatorActive = ImGuiCol_SeparatorActive;
const int const_ImGuiCol_ResizeGrip = ImGuiCol_ResizeGrip;
const int const_ImGuiCol_ResizeGripHovered = ImGuiCol_ResizeGripHovered;
const int const_ImGuiCol_ResizeGripActive = ImGuiCol_ResizeGripActive;
const int const_ImGuiCol_Tab = ImGuiCol_Tab;
const int const_ImGuiCol_TabHovered = ImGuiCol_TabHovered;
const int const_ImGuiCol_TabActive = ImGuiCol_TabActive;
const int const_ImGuiCol_TabUnfocused = ImGuiCol_TabUnfocused;
const int const_ImGuiCol_TabUnfocusedActive = ImGuiCol_TabUnfocusedActive;
const int const_ImGuiCol_DockingPreview = ImGuiCol_DockingPreview;
const int const_ImGuiCol_DockingEmptyBg = ImGuiCol_DockingEmptyBg;
const int const_ImGuiCol_PlotLines = ImGuiCol_PlotLines;
const int const_ImGuiCol_PlotLinesHovered = ImGuiCol_PlotLinesHovered;
const int const_ImGuiCol_PlotHistogram = ImGuiCol_PlotHistogram;
const int const_ImGuiCol_PlotHistogramHovered = ImGuiCol_PlotHistogramHovered;
const int const_ImGuiCol_TableHeaderBg = ImGuiCol_TableHeaderBg;
const int const_ImGuiCol_TableBorderStrong = ImGuiCol_TableBorderStrong;
const int const_ImGuiCol_TableBorderLight = ImGuiCol_TableBorderLight;
const int const_ImGuiCol_TableRowBg = ImGuiCol_TableRowBg;
const int const_ImGuiCol_TableRowBgAlt = ImGuiCol_TableRowBgAlt;
const int const_ImGuiCol_TextSelectedBg = ImGuiCol_TextSelectedBg;
const int const_ImGuiCol_DragDropTarget = ImGuiCol_DragDropTarget;
const int const_ImGuiCol_NavHighlight = ImGuiCol_NavHighlight;
const int const_ImGuiCol_NavWindowingHighlight = ImGuiCol_NavWindowingHighlight;
const int const_ImGuiCol_NavWindowingDimBg = ImGuiCol_NavWindowingDimBg;
const int const_ImGuiCol_ModalWindowDimBg = ImGuiCol_ModalWindowDimBg;
// const int const_ImGuiCol_COUNT = ImGuiCol_COUNT;
// enum ImGuiStyleVar_
const int const_ImGuiStyleVar_Alpha = ImGuiStyleVar_Alpha;
const int const_ImGuiStyleVar_DisabledAlpha = ImGuiStyleVar_DisabledAlpha;
const int const_ImGuiStyleVar_WindowPadding = ImGuiStyleVar_WindowPadding;
const int const_ImGuiStyleVar_WindowRounding = ImGuiStyleVar_WindowRounding;
const int const_ImGuiStyleVar_WindowBorderSize = ImGuiStyleVar_WindowBorderSize;
const int const_ImGuiStyleVar_WindowMinSize = ImGuiStyleVar_WindowMinSize;
const int const_ImGuiStyleVar_WindowTitleAlign = ImGuiStyleVar_WindowTitleAlign;
const int const_ImGuiStyleVar_ChildRounding = ImGuiStyleVar_ChildRounding;
const int const_ImGuiStyleVar_ChildBorderSize = ImGuiStyleVar_ChildBorderSize;
const int const_ImGuiStyleVar_PopupRounding = ImGuiStyleVar_PopupRounding;
const int const_ImGuiStyleVar_PopupBorderSize = ImGuiStyleVar_PopupBorderSize;
const int const_ImGuiStyleVar_FramePadding = ImGuiStyleVar_FramePadding;
const int const_ImGuiStyleVar_FrameRounding = ImGuiStyleVar_FrameRounding;
const int const_ImGuiStyleVar_FrameBorderSize = ImGuiStyleVar_FrameBorderSize;
const int const_ImGuiStyleVar_ItemSpacing = ImGuiStyleVar_ItemSpacing;
const int const_ImGuiStyleVar_ItemInnerSpacing = ImGuiStyleVar_ItemInnerSpacing;
const int const_ImGuiStyleVar_IndentSpacing = ImGuiStyleVar_IndentSpacing;
const int const_ImGuiStyleVar_CellPadding = ImGuiStyleVar_CellPadding;
const int const_ImGuiStyleVar_ScrollbarSize = ImGuiStyleVar_ScrollbarSize;
const int const_ImGuiStyleVar_ScrollbarRounding = ImGuiStyleVar_ScrollbarRounding;
const int const_ImGuiStyleVar_GrabMinSize = ImGuiStyleVar_GrabMinSize;
const int const_ImGuiStyleVar_GrabRounding = ImGuiStyleVar_GrabRounding;
const int const_ImGuiStyleVar_TabRounding = ImGuiStyleVar_TabRounding;
const int const_ImGuiStyleVar_ButtonTextAlign = ImGuiStyleVar_ButtonTextAlign;
const int const_ImGuiStyleVar_SelectableTextAlign = ImGuiStyleVar_SelectableTextAlign;
const int const_ImGuiStyleVar_COUNT = ImGuiStyleVar_COUNT;
// enum ImGuiButtonFlags_
const int const_ImGuiButtonFlags_None = ImGuiButtonFlags_None;
const int const_ImGuiButtonFlags_MouseButtonLeft = ImGuiButtonFlags_MouseButtonLeft;
const int const_ImGuiButtonFlags_MouseButtonRight = ImGuiButtonFlags_MouseButtonRight;
const int const_ImGuiButtonFlags_MouseButtonMiddle = ImGuiButtonFlags_MouseButtonMiddle;
const int const_ImGuiButtonFlags_MouseButtonMask_ = ImGuiButtonFlags_MouseButtonMask_;
const int const_ImGuiButtonFlags_MouseButtonDefault_ = ImGuiButtonFlags_MouseButtonDefault_;
// enum ImGuiColorEditFlags_
const int const_ImGuiColorEditFlags_None = ImGuiColorEditFlags_None;
const int const_ImGuiColorEditFlags_NoAlpha = ImGuiColorEditFlags_NoAlpha;
const int const_ImGuiColorEditFlags_NoPicker = ImGuiColorEditFlags_NoPicker;
const int const_ImGuiColorEditFlags_NoOptions = ImGuiColorEditFlags_NoOptions;
const int const_ImGuiColorEditFlags_NoSmallPreview = ImGuiColorEditFlags_NoSmallPreview;
const int const_ImGuiColorEditFlags_NoInputs = ImGuiColorEditFlags_NoInputs;
const int const_ImGuiColorEditFlags_NoTooltip = ImGuiColorEditFlags_NoTooltip;
const int const_ImGuiColorEditFlags_NoLabel = ImGuiColorEditFlags_NoLabel;
const int const_ImGuiColorEditFlags_NoSidePreview = ImGuiColorEditFlags_NoSidePreview;
const int const_ImGuiColorEditFlags_NoDragDrop = ImGuiColorEditFlags_NoDragDrop;
const int const_ImGuiColorEditFlags_NoBorder = ImGuiColorEditFlags_NoBorder;
const int const_ImGuiColorEditFlags_AlphaBar = ImGuiColorEditFlags_AlphaBar;
const int const_ImGuiColorEditFlags_AlphaPreview = ImGuiColorEditFlags_AlphaPreview;
const int const_ImGuiColorEditFlags_AlphaPreviewHalf = ImGuiColorEditFlags_AlphaPreviewHalf;
const int const_ImGuiColorEditFlags_HDR = ImGuiColorEditFlags_HDR;
const int const_ImGuiColorEditFlags_DisplayRGB = ImGuiColorEditFlags_DisplayRGB;
const int const_ImGuiColorEditFlags_DisplayHSV = ImGuiColorEditFlags_DisplayHSV;
const int const_ImGuiColorEditFlags_DisplayHex = ImGuiColorEditFlags_DisplayHex;
const int const_ImGuiColorEditFlags_Uint8 = ImGuiColorEditFlags_Uint8;
const int const_ImGuiColorEditFlags_Float = ImGuiColorEditFlags_Float;
const int const_ImGuiColorEditFlags_PickerHueBar = ImGuiColorEditFlags_PickerHueBar;
const int const_ImGuiColorEditFlags_PickerHueWheel = ImGuiColorEditFlags_PickerHueWheel;
const int const_ImGuiColorEditFlags_InputRGB = ImGuiColorEditFlags_InputRGB;
const int const_ImGuiColorEditFlags_InputHSV = ImGuiColorEditFlags_InputHSV;
const int const_ImGuiColorEditFlags_DefaultOptions_ = ImGuiColorEditFlags_DefaultOptions_;
const int const_ImGuiColorEditFlags_DisplayMask_ = ImGuiColorEditFlags_DisplayMask_;
const int const_ImGuiColorEditFlags_DataTypeMask_ = ImGuiColorEditFlags_DataTypeMask_;
const int const_ImGuiColorEditFlags_PickerMask_ = ImGuiColorEditFlags_PickerMask_;
const int const_ImGuiColorEditFlags_InputMask_ = ImGuiColorEditFlags_InputMask_;
// enum ImGuiSliderFlags_
const int const_ImGuiSliderFlags_None = ImGuiSliderFlags_None;
const int const_ImGuiSliderFlags_AlwaysClamp = ImGuiSliderFlags_AlwaysClamp;
const int const_ImGuiSliderFlags_Logarithmic = ImGuiSliderFlags_Logarithmic;
const int const_ImGuiSliderFlags_NoRoundToFormat = ImGuiSliderFlags_NoRoundToFormat;
const int const_ImGuiSliderFlags_NoInput = ImGuiSliderFlags_NoInput;
const int const_ImGuiSliderFlags_InvalidMask_ = ImGuiSliderFlags_InvalidMask_;
// enum ImGuiMouseButton_
const int const_ImGuiMouseButton_Left = ImGuiMouseButton_Left;
const int const_ImGuiMouseButton_Right = ImGuiMouseButton_Right;
const int const_ImGuiMouseButton_Middle = ImGuiMouseButton_Middle;
const int const_ImGuiMouseButton_COUNT = ImGuiMouseButton_COUNT;
// enum ImGuiMouseCursor_
const int const_ImGuiMouseCursor_None = ImGuiMouseCursor_None;
const int const_ImGuiMouseCursor_Arrow = ImGuiMouseCursor_Arrow;
const int const_ImGuiMouseCursor_TextInput = ImGuiMouseCursor_TextInput;
const int const_ImGuiMouseCursor_ResizeAll = ImGuiMouseCursor_ResizeAll;
const int const_ImGuiMouseCursor_ResizeNS = ImGuiMouseCursor_ResizeNS;
const int const_ImGuiMouseCursor_ResizeEW = ImGuiMouseCursor_ResizeEW;
const int const_ImGuiMouseCursor_ResizeNESW = ImGuiMouseCursor_ResizeNESW;
const int const_ImGuiMouseCursor_ResizeNWSE = ImGuiMouseCursor_ResizeNWSE;
const int const_ImGuiMouseCursor_Hand = ImGuiMouseCursor_Hand;
const int const_ImGuiMouseCursor_NotAllowed = ImGuiMouseCursor_NotAllowed;
const int const_ImGuiMouseCursor_COUNT = ImGuiMouseCursor_COUNT;
// enum ImGuiCond_
const int const_ImGuiCond_None = ImGuiCond_None;
const int const_ImGuiCond_Always = ImGuiCond_Always;
const int const_ImGuiCond_Once = ImGuiCond_Once;
const int const_ImGuiCond_FirstUseEver = ImGuiCond_FirstUseEver;
const int const_ImGuiCond_Appearing = ImGuiCond_Appearing;
// enum ImDrawFlags_
const int const_ImDrawFlags_None = ImDrawFlags_None;
const int const_ImDrawFlags_Closed = ImDrawFlags_Closed;
const int const_ImDrawFlags_RoundCornersTopLeft = ImDrawFlags_RoundCornersTopLeft;
const int const_ImDrawFlags_RoundCornersTopRight = ImDrawFlags_RoundCornersTopRight;
const int const_ImDrawFlags_RoundCornersBottomLeft = ImDrawFlags_RoundCornersBottomLeft;
const int const_ImDrawFlags_RoundCornersBottomRight = ImDrawFlags_RoundCornersBottomRight;
const int const_ImDrawFlags_RoundCornersNone = ImDrawFlags_RoundCornersNone;
const int const_ImDrawFlags_RoundCornersTop = ImDrawFlags_RoundCornersTop;
const int const_ImDrawFlags_RoundCornersBottom = ImDrawFlags_RoundCornersBottom;
const int const_ImDrawFlags_RoundCornersLeft = ImDrawFlags_RoundCornersLeft;
const int const_ImDrawFlags_RoundCornersRight = ImDrawFlags_RoundCornersRight;
const int const_ImDrawFlags_RoundCornersAll = ImDrawFlags_RoundCornersAll;
const int const_ImDrawFlags_RoundCornersDefault_ = ImDrawFlags_RoundCornersDefault_;
const int const_ImDrawFlags_RoundCornersMask_ = ImDrawFlags_RoundCornersMask_;
// enum ImDrawListFlags_
const int const_ImDrawListFlags_None = ImDrawListFlags_None;
const int const_ImDrawListFlags_AntiAliasedLines = ImDrawListFlags_AntiAliasedLines;
const int const_ImDrawListFlags_AntiAliasedLinesUseTex = ImDrawListFlags_AntiAliasedLinesUseTex;
const int const_ImDrawListFlags_AntiAliasedFill = ImDrawListFlags_AntiAliasedFill;
const int const_ImDrawListFlags_AllowVtxOffset = ImDrawListFlags_AllowVtxOffset;
// enum ImFontAtlasFlags_
const int const_ImFontAtlasFlags_None = ImFontAtlasFlags_None;
const int const_ImFontAtlasFlags_NoPowerOfTwoHeight = ImFontAtlasFlags_NoPowerOfTwoHeight;
const int const_ImFontAtlasFlags_NoMouseCursors = ImFontAtlasFlags_NoMouseCursors;
const int const_ImFontAtlasFlags_NoBakedLines = ImFontAtlasFlags_NoBakedLines;
// enum ImGuiViewportFlags_
const int const_ImGuiViewportFlags_None = ImGuiViewportFlags_None;
const int const_ImGuiViewportFlags_IsPlatformWindow = ImGuiViewportFlags_IsPlatformWindow;
const int const_ImGuiViewportFlags_IsPlatformMonitor = ImGuiViewportFlags_IsPlatformMonitor;
const int const_ImGuiViewportFlags_OwnedByApp = ImGuiViewportFlags_OwnedByApp;
const int const_ImGuiViewportFlags_NoDecoration = ImGuiViewportFlags_NoDecoration;
const int const_ImGuiViewportFlags_NoTaskBarIcon = ImGuiViewportFlags_NoTaskBarIcon;
const int const_ImGuiViewportFlags_NoFocusOnAppearing = ImGuiViewportFlags_NoFocusOnAppearing;
const int const_ImGuiViewportFlags_NoFocusOnClick = ImGuiViewportFlags_NoFocusOnClick;
const int const_ImGuiViewportFlags_NoInputs = ImGuiViewportFlags_NoInputs;
const int const_ImGuiViewportFlags_NoRendererClear = ImGuiViewportFlags_NoRendererClear;
const int const_ImGuiViewportFlags_TopMost = ImGuiViewportFlags_TopMost;
const int const_ImGuiViewportFlags_Minimized = ImGuiViewportFlags_Minimized;
const int const_ImGuiViewportFlags_NoAutoMerge = ImGuiViewportFlags_NoAutoMerge;
const int const_ImGuiViewportFlags_CanHostOtherWindows = ImGuiViewportFlags_CanHostOtherWindows;
// enum ImGuiItemFlags_
const int const_ImGuiItemFlags_None = ImGuiItemFlags_None;
const int const_ImGuiItemFlags_NoTabStop = ImGuiItemFlags_NoTabStop;
const int const_ImGuiItemFlags_ButtonRepeat = ImGuiItemFlags_ButtonRepeat;
const int const_ImGuiItemFlags_Disabled = ImGuiItemFlags_Disabled;
const int const_ImGuiItemFlags_NoNav = ImGuiItemFlags_NoNav;
const int const_ImGuiItemFlags_NoNavDefaultFocus = ImGuiItemFlags_NoNavDefaultFocus;
const int const_ImGuiItemFlags_SelectableDontClosePopup = ImGuiItemFlags_SelectableDontClosePopup;
const int const_ImGuiItemFlags_MixedValue = ImGuiItemFlags_MixedValue;
const int const_ImGuiItemFlags_ReadOnly = ImGuiItemFlags_ReadOnly;
const int const_ImGuiItemFlags_Inputable = ImGuiItemFlags_Inputable;
// enum ImGuiItemStatusFlags_
const int const_ImGuiItemStatusFlags_None = ImGuiItemStatusFlags_None;
const int const_ImGuiItemStatusFlags_HoveredRect = ImGuiItemStatusFlags_HoveredRect;
const int const_ImGuiItemStatusFlags_HasDisplayRect = ImGuiItemStatusFlags_HasDisplayRect;
const int const_ImGuiItemStatusFlags_Edited = ImGuiItemStatusFlags_Edited;
const int const_ImGuiItemStatusFlags_ToggledSelection = ImGuiItemStatusFlags_ToggledSelection;
const int const_ImGuiItemStatusFlags_ToggledOpen = ImGuiItemStatusFlags_ToggledOpen;
const int const_ImGuiItemStatusFlags_HasDeactivated = ImGuiItemStatusFlags_HasDeactivated;
const int const_ImGuiItemStatusFlags_Deactivated = ImGuiItemStatusFlags_Deactivated;
const int const_ImGuiItemStatusFlags_HoveredWindow = ImGuiItemStatusFlags_HoveredWindow;
const int const_ImGuiItemStatusFlags_FocusedByTabbing = ImGuiItemStatusFlags_FocusedByTabbing;
// enum ImGuiInputTextFlagsPrivate_
const int const_ImGuiInputTextFlags_Multiline = ImGuiInputTextFlags_Multiline;
const int const_ImGuiInputTextFlags_NoMarkEdited = ImGuiInputTextFlags_NoMarkEdited;
const int const_ImGuiInputTextFlags_MergedItem = ImGuiInputTextFlags_MergedItem;
// enum ImGuiButtonFlagsPrivate_
const int const_ImGuiButtonFlags_PressedOnClick = ImGuiButtonFlags_PressedOnClick;
const int const_ImGuiButtonFlags_PressedOnClickRelease = ImGuiButtonFlags_PressedOnClickRelease;
const int const_ImGuiButtonFlags_PressedOnClickReleaseAnywhere = ImGuiButtonFlags_PressedOnClickReleaseAnywhere;
const int const_ImGuiButtonFlags_PressedOnRelease = ImGuiButtonFlags_PressedOnRelease;
const int const_ImGuiButtonFlags_PressedOnDoubleClick = ImGuiButtonFlags_PressedOnDoubleClick;
const int const_ImGuiButtonFlags_PressedOnDragDropHold = ImGuiButtonFlags_PressedOnDragDropHold;
const int const_ImGuiButtonFlags_Repeat = ImGuiButtonFlags_Repeat;
const int const_ImGuiButtonFlags_FlattenChildren = ImGuiButtonFlags_FlattenChildren;
const int const_ImGuiButtonFlags_AllowItemOverlap = ImGuiButtonFlags_AllowItemOverlap;
const int const_ImGuiButtonFlags_DontClosePopups = ImGuiButtonFlags_DontClosePopups;
const int const_ImGuiButtonFlags_AlignTextBaseLine = ImGuiButtonFlags_AlignTextBaseLine;
const int const_ImGuiButtonFlags_NoKeyModifiers = ImGuiButtonFlags_NoKeyModifiers;
const int const_ImGuiButtonFlags_NoHoldingActiveId = ImGuiButtonFlags_NoHoldingActiveId;
const int const_ImGuiButtonFlags_NoNavFocus = ImGuiButtonFlags_NoNavFocus;
const int const_ImGuiButtonFlags_NoHoveredOnFocus = ImGuiButtonFlags_NoHoveredOnFocus;
const int const_ImGuiButtonFlags_PressedOnMask_ = ImGuiButtonFlags_PressedOnMask_;
const int const_ImGuiButtonFlags_PressedOnDefault_ = ImGuiButtonFlags_PressedOnDefault_;
// enum ImGuiComboFlagsPrivate_
const int const_ImGuiComboFlags_CustomPreview = ImGuiComboFlags_CustomPreview;
// enum ImGuiSliderFlagsPrivate_
const int const_ImGuiSliderFlags_Vertical = ImGuiSliderFlags_Vertical;
const int const_ImGuiSliderFlags_ReadOnly = ImGuiSliderFlags_ReadOnly;
// enum ImGuiSelectableFlagsPrivate_
const int const_ImGuiSelectableFlags_NoHoldingActiveID = ImGuiSelectableFlags_NoHoldingActiveID;
const int const_ImGuiSelectableFlags_SelectOnNav = ImGuiSelectableFlags_SelectOnNav;
const int const_ImGuiSelectableFlags_SelectOnClick = ImGuiSelectableFlags_SelectOnClick;
const int const_ImGuiSelectableFlags_SelectOnRelease = ImGuiSelectableFlags_SelectOnRelease;
const int const_ImGuiSelectableFlags_SpanAvailWidth = ImGuiSelectableFlags_SpanAvailWidth;
const int const_ImGuiSelectableFlags_DrawHoveredWhenHeld = ImGuiSelectableFlags_DrawHoveredWhenHeld;
const int const_ImGuiSelectableFlags_SetNavIdOnHover = ImGuiSelectableFlags_SetNavIdOnHover;
const int const_ImGuiSelectableFlags_NoPadWithHalfSpacing = ImGuiSelectableFlags_NoPadWithHalfSpacing;
// enum ImGuiTreeNodeFlagsPrivate_
const int const_ImGuiTreeNodeFlags_ClipLabelForTrailingButton = ImGuiTreeNodeFlags_ClipLabelForTrailingButton;
// enum ImGuiSeparatorFlags_
const int const_ImGuiSeparatorFlags_None = ImGuiSeparatorFlags_None;
const int const_ImGuiSeparatorFlags_Horizontal = ImGuiSeparatorFlags_Horizontal;
const int const_ImGuiSeparatorFlags_Vertical = ImGuiSeparatorFlags_Vertical;
const int const_ImGuiSeparatorFlags_SpanAllColumns = ImGuiSeparatorFlags_SpanAllColumns;
// enum ImGuiTextFlags_
const int const_ImGuiTextFlags_None = ImGuiTextFlags_None;
const int const_ImGuiTextFlags_NoWidthForLargeClippedText = ImGuiTextFlags_NoWidthForLargeClippedText;
// enum ImGuiTooltipFlags_
const int const_ImGuiTooltipFlags_None = ImGuiTooltipFlags_None;
const int const_ImGuiTooltipFlags_OverridePreviousTooltip = ImGuiTooltipFlags_OverridePreviousTooltip;
// enum ImGuiLayoutType_
const int const_ImGuiLayoutType_Horizontal = ImGuiLayoutType_Horizontal;
const int const_ImGuiLayoutType_Vertical = ImGuiLayoutType_Vertical;
// enum ImGuiLogType
const int const_ImGuiLogType_None = ImGuiLogType_None;
const int const_ImGuiLogType_TTY = ImGuiLogType_TTY;
const int const_ImGuiLogType_File = ImGuiLogType_File;
const int const_ImGuiLogType_Buffer = ImGuiLogType_Buffer;
const int const_ImGuiLogType_Clipboard = ImGuiLogType_Clipboard;
// enum ImGuiAxis
const int const_ImGuiAxis_None = ImGuiAxis_None;
const int const_ImGuiAxis_X = ImGuiAxis_X;
const int const_ImGuiAxis_Y = ImGuiAxis_Y;
// enum ImGuiPlotType
const int const_ImGuiPlotType_Lines = ImGuiPlotType_Lines;
const int const_ImGuiPlotType_Histogram = ImGuiPlotType_Histogram;
// enum ImGuiPopupPositionPolicy
const int const_ImGuiPopupPositionPolicy_Default = ImGuiPopupPositionPolicy_Default;
const int const_ImGuiPopupPositionPolicy_ComboBox = ImGuiPopupPositionPolicy_ComboBox;
const int const_ImGuiPopupPositionPolicy_Tooltip = ImGuiPopupPositionPolicy_Tooltip;
// enum ImGuiDataTypePrivate_
const int const_ImGuiDataType_String = ImGuiDataType_String;
const int const_ImGuiDataType_Pointer = ImGuiDataType_Pointer;
const int const_ImGuiDataType_ID = ImGuiDataType_ID;
// enum ImGuiNextWindowDataFlags_
const int const_ImGuiNextWindowDataFlags_None = ImGuiNextWindowDataFlags_None;
const int const_ImGuiNextWindowDataFlags_HasPos = ImGuiNextWindowDataFlags_HasPos;
const int const_ImGuiNextWindowDataFlags_HasSize = ImGuiNextWindowDataFlags_HasSize;
const int const_ImGuiNextWindowDataFlags_HasContentSize = ImGuiNextWindowDataFlags_HasContentSize;
const int const_ImGuiNextWindowDataFlags_HasCollapsed = ImGuiNextWindowDataFlags_HasCollapsed;
const int const_ImGuiNextWindowDataFlags_HasSizeConstraint = ImGuiNextWindowDataFlags_HasSizeConstraint;
const int const_ImGuiNextWindowDataFlags_HasFocus = ImGuiNextWindowDataFlags_HasFocus;
const int const_ImGuiNextWindowDataFlags_HasBgAlpha = ImGuiNextWindowDataFlags_HasBgAlpha;
const int const_ImGuiNextWindowDataFlags_HasScroll = ImGuiNextWindowDataFlags_HasScroll;
const int const_ImGuiNextWindowDataFlags_HasViewport = ImGuiNextWindowDataFlags_HasViewport;
const int const_ImGuiNextWindowDataFlags_HasDock = ImGuiNextWindowDataFlags_HasDock;
const int const_ImGuiNextWindowDataFlags_HasWindowClass = ImGuiNextWindowDataFlags_HasWindowClass;
// enum ImGuiNextItemDataFlags_
const int const_ImGuiNextItemDataFlags_None = ImGuiNextItemDataFlags_None;
const int const_ImGuiNextItemDataFlags_HasWidth = ImGuiNextItemDataFlags_HasWidth;
const int const_ImGuiNextItemDataFlags_HasOpen = ImGuiNextItemDataFlags_HasOpen;
// enum ImGuiKeyPrivate_
const int const_ImGuiKey_LegacyNativeKey_BEGIN = ImGuiKey_LegacyNativeKey_BEGIN;
const int const_ImGuiKey_LegacyNativeKey_END = ImGuiKey_LegacyNativeKey_END;
const int const_ImGuiKey_Gamepad_BEGIN = ImGuiKey_Gamepad_BEGIN;
const int const_ImGuiKey_Gamepad_END = ImGuiKey_Gamepad_END;
// enum ImGuiInputEventType
const int const_ImGuiInputEventType_None = ImGuiInputEventType_None;
const int const_ImGuiInputEventType_MousePos = ImGuiInputEventType_MousePos;
const int const_ImGuiInputEventType_MouseWheel = ImGuiInputEventType_MouseWheel;
const int const_ImGuiInputEventType_MouseButton = ImGuiInputEventType_MouseButton;
const int const_ImGuiInputEventType_MouseViewport = ImGuiInputEventType_MouseViewport;
const int const_ImGuiInputEventType_Key = ImGuiInputEventType_Key;
const int const_ImGuiInputEventType_Text = ImGuiInputEventType_Text;
const int const_ImGuiInputEventType_Focus = ImGuiInputEventType_Focus;
const int const_ImGuiInputEventType_COUNT = ImGuiInputEventType_COUNT;
// enum ImGuiInputSource
const int const_ImGuiInputSource_None = ImGuiInputSource_None;
const int const_ImGuiInputSource_Mouse = ImGuiInputSource_Mouse;
const int const_ImGuiInputSource_Keyboard = ImGuiInputSource_Keyboard;
const int const_ImGuiInputSource_Gamepad = ImGuiInputSource_Gamepad;
const int const_ImGuiInputSource_Clipboard = ImGuiInputSource_Clipboard;
const int const_ImGuiInputSource_Nav = ImGuiInputSource_Nav;
const int const_ImGuiInputSource_COUNT = ImGuiInputSource_COUNT;
// enum ImGuiNavReadMode
const int const_ImGuiNavReadMode_Down = ImGuiNavReadMode_Down;
const int const_ImGuiNavReadMode_Pressed = ImGuiNavReadMode_Pressed;
const int const_ImGuiNavReadMode_Released = ImGuiNavReadMode_Released;
const int const_ImGuiNavReadMode_Repeat = ImGuiNavReadMode_Repeat;
const int const_ImGuiNavReadMode_RepeatSlow = ImGuiNavReadMode_RepeatSlow;
const int const_ImGuiNavReadMode_RepeatFast = ImGuiNavReadMode_RepeatFast;
// enum ImGuiActivateFlags_
const int const_ImGuiActivateFlags_None = ImGuiActivateFlags_None;
const int const_ImGuiActivateFlags_PreferInput = ImGuiActivateFlags_PreferInput;
const int const_ImGuiActivateFlags_PreferTweak = ImGuiActivateFlags_PreferTweak;
const int const_ImGuiActivateFlags_TryToPreserveState = ImGuiActivateFlags_TryToPreserveState;
// enum ImGuiScrollFlags_
const int const_ImGuiScrollFlags_None = ImGuiScrollFlags_None;
const int const_ImGuiScrollFlags_KeepVisibleEdgeX = ImGuiScrollFlags_KeepVisibleEdgeX;
const int const_ImGuiScrollFlags_KeepVisibleEdgeY = ImGuiScrollFlags_KeepVisibleEdgeY;
const int const_ImGuiScrollFlags_KeepVisibleCenterX = ImGuiScrollFlags_KeepVisibleCenterX;
const int const_ImGuiScrollFlags_KeepVisibleCenterY = ImGuiScrollFlags_KeepVisibleCenterY;
const int const_ImGuiScrollFlags_AlwaysCenterX = ImGuiScrollFlags_AlwaysCenterX;
const int const_ImGuiScrollFlags_AlwaysCenterY = ImGuiScrollFlags_AlwaysCenterY;
const int const_ImGuiScrollFlags_NoScrollParent = ImGuiScrollFlags_NoScrollParent;
const int const_ImGuiScrollFlags_MaskX_ = ImGuiScrollFlags_MaskX_;
const int const_ImGuiScrollFlags_MaskY_ = ImGuiScrollFlags_MaskY_;
// enum ImGuiNavHighlightFlags_
const int const_ImGuiNavHighlightFlags_None = ImGuiNavHighlightFlags_None;
const int const_ImGuiNavHighlightFlags_TypeDefault = ImGuiNavHighlightFlags_TypeDefault;
const int const_ImGuiNavHighlightFlags_TypeThin = ImGuiNavHighlightFlags_TypeThin;
const int const_ImGuiNavHighlightFlags_AlwaysDraw = ImGuiNavHighlightFlags_AlwaysDraw;
const int const_ImGuiNavHighlightFlags_NoRounding = ImGuiNavHighlightFlags_NoRounding;
// enum ImGuiNavDirSourceFlags_
const int const_ImGuiNavDirSourceFlags_None = ImGuiNavDirSourceFlags_None;
const int const_ImGuiNavDirSourceFlags_RawKeyboard = ImGuiNavDirSourceFlags_RawKeyboard;
const int const_ImGuiNavDirSourceFlags_Keyboard = ImGuiNavDirSourceFlags_Keyboard;
const int const_ImGuiNavDirSourceFlags_PadDPad = ImGuiNavDirSourceFlags_PadDPad;
const int const_ImGuiNavDirSourceFlags_PadLStick = ImGuiNavDirSourceFlags_PadLStick;
// enum ImGuiNavMoveFlags_
const int const_ImGuiNavMoveFlags_None = ImGuiNavMoveFlags_None;
const int const_ImGuiNavMoveFlags_LoopX = ImGuiNavMoveFlags_LoopX;
const int const_ImGuiNavMoveFlags_LoopY = ImGuiNavMoveFlags_LoopY;
const int const_ImGuiNavMoveFlags_WrapX = ImGuiNavMoveFlags_WrapX;
const int const_ImGuiNavMoveFlags_WrapY = ImGuiNavMoveFlags_WrapY;
const int const_ImGuiNavMoveFlags_AllowCurrentNavId = ImGuiNavMoveFlags_AllowCurrentNavId;
const int const_ImGuiNavMoveFlags_AlsoScoreVisibleSet = ImGuiNavMoveFlags_AlsoScoreVisibleSet;
const int const_ImGuiNavMoveFlags_ScrollToEdgeY = ImGuiNavMoveFlags_ScrollToEdgeY;
const int const_ImGuiNavMoveFlags_Forwarded = ImGuiNavMoveFlags_Forwarded;
const int const_ImGuiNavMoveFlags_DebugNoResult = ImGuiNavMoveFlags_DebugNoResult;
const int const_ImGuiNavMoveFlags_FocusApi = ImGuiNavMoveFlags_FocusApi;
const int const_ImGuiNavMoveFlags_Tabbing = ImGuiNavMoveFlags_Tabbing;
const int const_ImGuiNavMoveFlags_Activate = ImGuiNavMoveFlags_Activate;
const int const_ImGuiNavMoveFlags_DontSetNavHighlight = ImGuiNavMoveFlags_DontSetNavHighlight;
// enum ImGuiNavLayer
const int const_ImGuiNavLayer_Main = ImGuiNavLayer_Main;
const int const_ImGuiNavLayer_Menu = ImGuiNavLayer_Menu;
// const int const_ImGuiNavLayer_COUNT = ImGuiNavLayer_COUNT;
// enum ImGuiOldColumnFlags_
const int const_ImGuiOldColumnFlags_None = ImGuiOldColumnFlags_None;
const int const_ImGuiOldColumnFlags_NoBorder = ImGuiOldColumnFlags_NoBorder;
const int const_ImGuiOldColumnFlags_NoResize = ImGuiOldColumnFlags_NoResize;
const int const_ImGuiOldColumnFlags_NoPreserveWidths = ImGuiOldColumnFlags_NoPreserveWidths;
const int const_ImGuiOldColumnFlags_NoForceWithinWindow = ImGuiOldColumnFlags_NoForceWithinWindow;
const int const_ImGuiOldColumnFlags_GrowParentContentsSize = ImGuiOldColumnFlags_GrowParentContentsSize;
// enum ImGuiDockNodeFlagsPrivate_
const int const_ImGuiDockNodeFlags_DockSpace = ImGuiDockNodeFlags_DockSpace;
const int const_ImGuiDockNodeFlags_CentralNode = ImGuiDockNodeFlags_CentralNode;
const int const_ImGuiDockNodeFlags_NoTabBar = ImGuiDockNodeFlags_NoTabBar;
const int const_ImGuiDockNodeFlags_HiddenTabBar = ImGuiDockNodeFlags_HiddenTabBar;
const int const_ImGuiDockNodeFlags_NoWindowMenuButton = ImGuiDockNodeFlags_NoWindowMenuButton;
const int const_ImGuiDockNodeFlags_NoCloseButton = ImGuiDockNodeFlags_NoCloseButton;
const int const_ImGuiDockNodeFlags_NoDocking = ImGuiDockNodeFlags_NoDocking;
const int const_ImGuiDockNodeFlags_NoDockingSplitMe = ImGuiDockNodeFlags_NoDockingSplitMe;
const int const_ImGuiDockNodeFlags_NoDockingSplitOther = ImGuiDockNodeFlags_NoDockingSplitOther;
const int const_ImGuiDockNodeFlags_NoDockingOverMe = ImGuiDockNodeFlags_NoDockingOverMe;
const int const_ImGuiDockNodeFlags_NoDockingOverOther = ImGuiDockNodeFlags_NoDockingOverOther;
const int const_ImGuiDockNodeFlags_NoDockingOverEmpty = ImGuiDockNodeFlags_NoDockingOverEmpty;
const int const_ImGuiDockNodeFlags_NoResizeX = ImGuiDockNodeFlags_NoResizeX;
const int const_ImGuiDockNodeFlags_NoResizeY = ImGuiDockNodeFlags_NoResizeY;
const int const_ImGuiDockNodeFlags_SharedFlagsInheritMask_ = ImGuiDockNodeFlags_SharedFlagsInheritMask_;
const int const_ImGuiDockNodeFlags_NoResizeFlagsMask_ = ImGuiDockNodeFlags_NoResizeFlagsMask_;
const int const_ImGuiDockNodeFlags_LocalFlagsMask_ = ImGuiDockNodeFlags_LocalFlagsMask_;
const int const_ImGuiDockNodeFlags_LocalFlagsTransferMask_ = ImGuiDockNodeFlags_LocalFlagsTransferMask_;
const int const_ImGuiDockNodeFlags_SavedFlagsMask_ = ImGuiDockNodeFlags_SavedFlagsMask_;
// enum ImGuiDataAuthority_
const int const_ImGuiDataAuthority_Auto = ImGuiDataAuthority_Auto;
const int const_ImGuiDataAuthority_DockNode = ImGuiDataAuthority_DockNode;
const int const_ImGuiDataAuthority_Window = ImGuiDataAuthority_Window;
// enum ImGuiDockNodeState
const int const_ImGuiDockNodeState_Unknown = ImGuiDockNodeState_Unknown;
const int const_ImGuiDockNodeState_HostWindowHiddenBecauseSingleWindow = ImGuiDockNodeState_HostWindowHiddenBecauseSingleWindow;
const int const_ImGuiDockNodeState_HostWindowHiddenBecauseWindowsAreResizing = ImGuiDockNodeState_HostWindowHiddenBecauseWindowsAreResizing;
const int const_ImGuiDockNodeState_HostWindowVisible = ImGuiDockNodeState_HostWindowVisible;
// enum ImGuiWindowDockStyleCol
const int const_ImGuiWindowDockStyleCol_Text = ImGuiWindowDockStyleCol_Text;
const int const_ImGuiWindowDockStyleCol_Tab = ImGuiWindowDockStyleCol_Tab;
const int const_ImGuiWindowDockStyleCol_TabHovered = ImGuiWindowDockStyleCol_TabHovered;
const int const_ImGuiWindowDockStyleCol_TabActive = ImGuiWindowDockStyleCol_TabActive;
const int const_ImGuiWindowDockStyleCol_TabUnfocused = ImGuiWindowDockStyleCol_TabUnfocused;
const int const_ImGuiWindowDockStyleCol_TabUnfocusedActive = ImGuiWindowDockStyleCol_TabUnfocusedActive;
// const int const_ImGuiWindowDockStyleCol_COUNT = ImGuiWindowDockStyleCol_COUNT;
// enum ImGuiDebugLogFlags_
const int const_ImGuiDebugLogFlags_None = ImGuiDebugLogFlags_None;
const int const_ImGuiDebugLogFlags_EventActiveId = ImGuiDebugLogFlags_EventActiveId;
const int const_ImGuiDebugLogFlags_EventFocus = ImGuiDebugLogFlags_EventFocus;
const int const_ImGuiDebugLogFlags_EventPopup = ImGuiDebugLogFlags_EventPopup;
const int const_ImGuiDebugLogFlags_EventNav = ImGuiDebugLogFlags_EventNav;
const int const_ImGuiDebugLogFlags_EventIO = ImGuiDebugLogFlags_EventIO;
const int const_ImGuiDebugLogFlags_EventDocking = ImGuiDebugLogFlags_EventDocking;
const int const_ImGuiDebugLogFlags_EventViewport = ImGuiDebugLogFlags_EventViewport;
const int const_ImGuiDebugLogFlags_EventMask_ = ImGuiDebugLogFlags_EventMask_;
const int const_ImGuiDebugLogFlags_OutputToTTY = ImGuiDebugLogFlags_OutputToTTY;
// enum ImGuiContextHookType
const int const_ImGuiContextHookType_NewFramePre = ImGuiContextHookType_NewFramePre;
const int const_ImGuiContextHookType_NewFramePost = ImGuiContextHookType_NewFramePost;
const int const_ImGuiContextHookType_EndFramePre = ImGuiContextHookType_EndFramePre;
const int const_ImGuiContextHookType_EndFramePost = ImGuiContextHookType_EndFramePost;
const int const_ImGuiContextHookType_RenderPre = ImGuiContextHookType_RenderPre;
const int const_ImGuiContextHookType_RenderPost = ImGuiContextHookType_RenderPost;
const int const_ImGuiContextHookType_Shutdown = ImGuiContextHookType_Shutdown;
const int const_ImGuiContextHookType_PendingRemoval_ = ImGuiContextHookType_PendingRemoval_;
// enum ImGuiTabBarFlagsPrivate_
const int const_ImGuiTabBarFlags_DockNode = ImGuiTabBarFlags_DockNode;
const int const_ImGuiTabBarFlags_IsFocused = ImGuiTabBarFlags_IsFocused;
const int const_ImGuiTabBarFlags_SaveSettings = ImGuiTabBarFlags_SaveSettings;
// enum ImGuiTabItemFlagsPrivate_
const int const_ImGuiTabItemFlags_SectionMask_ = ImGuiTabItemFlags_SectionMask_;
const int const_ImGuiTabItemFlags_NoCloseButton = ImGuiTabItemFlags_NoCloseButton;
const int const_ImGuiTabItemFlags_Button = ImGuiTabItemFlags_Button;
const int const_ImGuiTabItemFlags_Unsorted = ImGuiTabItemFlags_Unsorted;
const int const_ImGuiTabItemFlags_Preview = ImGuiTabItemFlags_Preview;
#ifdef IMGUI_ENABLE_FREETYPE
// enum ImGuiFreeTypeBuilderFlags
const int const_ImGuiFreeTypeBuilderFlags_NoHinting = ImGuiFreeTypeBuilderFlags_NoHinting;
const int const_ImGuiFreeTypeBuilderFlags_NoAutoHint = ImGuiFreeTypeBuilderFlags_NoAutoHint;
const int const_ImGuiFreeTypeBuilderFlags_ForceAutoHint = ImGuiFreeTypeBuilderFlags_ForceAutoHint;
const int const_ImGuiFreeTypeBuilderFlags_LightHinting = ImGuiFreeTypeBuilderFlags_LightHinting;
const int const_ImGuiFreeTypeBuilderFlags_MonoHinting = ImGuiFreeTypeBuilderFlags_MonoHinting;
const int const_ImGuiFreeTypeBuilderFlags_Bold = ImGuiFreeTypeBuilderFlags_Bold;
const int const_ImGuiFreeTypeBuilderFlags_Oblique = ImGuiFreeTypeBuilderFlags_Oblique;
const int const_ImGuiFreeTypeBuilderFlags_Monochrome = ImGuiFreeTypeBuilderFlags_Monochrome;
const int const_ImGuiFreeTypeBuilderFlags_LoadColor = ImGuiFreeTypeBuilderFlags_LoadColor;
const int const_ImGuiFreeTypeBuilderFlags_Bitmap = ImGuiFreeTypeBuilderFlags_Bitmap;
#endif
