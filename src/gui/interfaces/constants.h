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

// Preprocessor-constants for Fortran.

#ifndef CONSTANTS_H
#define CONSTANTS_H

// threads
/* Function return values */
extern "C" const int const_thrd_error;
extern "C" const int const_thrd_success;
extern "C" const int const_thrd_timedout;
extern "C" const int const_thrd_busy;
extern "C" const int const_thrd_nomem;
/* Mutex types */
extern "C" const int const_mtx_plain;
extern "C" const int const_mtx_timed;
extern "C" const int const_mtx_recursive;

// OpenGL
extern "C" const int const_GL_COLOR_BUFFER_BIT;

extern "C" // GLFW;
extern "C" const int const_GLFW_TRUE;
extern "C" const int const_GLFW_FALSE;
extern "C" const int const_GLFW_SAMPLES;
extern "C" const int const_GLFW_CONTEXT_VERSION_MAJOR;
extern "C" const int const_GLFW_CONTEXT_VERSION_MINOR;
extern "C" const int const_GLFW_OPENGL_PROFILE;
extern "C" const int const_GLFW_OPENGL_FORWARD_COMPAT;
extern "C" const int const_GLFW_OPENGL_CORE_PROFILE;
extern "C" const int const_GLFW_STICKY_KEYS;

// cimgui
// enum ImGuiWindowFlags_
extern "C" const int const_ImGuiWindowFlags_None;
extern "C" const int const_ImGuiWindowFlags_NoTitleBar;
extern "C" const int const_ImGuiWindowFlags_NoResize;
extern "C" const int const_ImGuiWindowFlags_NoMove;
extern "C" const int const_ImGuiWindowFlags_NoScrollbar;
extern "C" const int const_ImGuiWindowFlags_NoScrollWithMouse;
extern "C" const int const_ImGuiWindowFlags_NoCollapse;
extern "C" const int const_ImGuiWindowFlags_AlwaysAutoResize;
extern "C" const int const_ImGuiWindowFlags_NoBackground;
extern "C" const int const_ImGuiWindowFlags_NoSavedSettings;
extern "C" const int const_ImGuiWindowFlags_NoMouseInputs;
extern "C" const int const_ImGuiWindowFlags_MenuBar;
extern "C" const int const_ImGuiWindowFlags_HorizontalScrollbar;
extern "C" const int const_ImGuiWindowFlags_NoFocusOnAppearing;
extern "C" const int const_ImGuiWindowFlags_NoBringToFrontOnFocus;
extern "C" const int const_ImGuiWindowFlags_AlwaysVerticalScrollbar;
extern "C" const int const_ImGuiWindowFlags_AlwaysHorizontalScrollbar;
extern "C" const int const_ImGuiWindowFlags_AlwaysUseWindowPadding;
extern "C" const int const_ImGuiWindowFlags_NoNavInputs;
extern "C" const int const_ImGuiWindowFlags_NoNavFocus;
extern "C" const int const_ImGuiWindowFlags_UnsavedDocument;
extern "C" const int const_ImGuiWindowFlags_NoDocking;
extern "C" const int const_ImGuiWindowFlags_NoNav;
extern "C" const int const_ImGuiWindowFlags_NoDecoration;
extern "C" const int const_ImGuiWindowFlags_NoInputs;
extern "C" const int const_ImGuiWindowFlags_NavFlattened;
extern "C" const int const_ImGuiWindowFlags_ChildWindow;
extern "C" const int const_ImGuiWindowFlags_Tooltip;
extern "C" const int const_ImGuiWindowFlags_Popup;
extern "C" const int const_ImGuiWindowFlags_Modal;
extern "C" const int const_ImGuiWindowFlags_ChildMenu;
extern "C" const int const_ImGuiWindowFlags_DockNodeHost;
// enum ImGuiInputTextFlags_
extern "C" const int const_ImGuiInputTextFlags_None;
extern "C" const int const_ImGuiInputTextFlags_CharsDecimal;
extern "C" const int const_ImGuiInputTextFlags_CharsHexadecimal;
extern "C" const int const_ImGuiInputTextFlags_CharsUppercase;
extern "C" const int const_ImGuiInputTextFlags_CharsNoBlank;
extern "C" const int const_ImGuiInputTextFlags_AutoSelectAll;
extern "C" const int const_ImGuiInputTextFlags_EnterReturnsTrue;
extern "C" const int const_ImGuiInputTextFlags_CallbackCompletion;
extern "C" const int const_ImGuiInputTextFlags_CallbackHistory;
extern "C" const int const_ImGuiInputTextFlags_CallbackAlways;
extern "C" const int const_ImGuiInputTextFlags_CallbackCharFilter;
extern "C" const int const_ImGuiInputTextFlags_AllowTabInput;
extern "C" const int const_ImGuiInputTextFlags_CtrlEnterForNewLine;
extern "C" const int const_ImGuiInputTextFlags_NoHorizontalScroll;
extern "C" const int const_ImGuiInputTextFlags_AlwaysOverwrite;
extern "C" const int const_ImGuiInputTextFlags_ReadOnly;
extern "C" const int const_ImGuiInputTextFlags_Password;
extern "C" const int const_ImGuiInputTextFlags_NoUndoRedo;
extern "C" const int const_ImGuiInputTextFlags_CharsScientific;
extern "C" const int const_ImGuiInputTextFlags_CallbackResize;
extern "C" const int const_ImGuiInputTextFlags_CallbackEdit;
// enum ImGuiTreeNodeFlags_;
extern "C" const int const_ImGuiTreeNodeFlags_None;
extern "C" const int const_ImGuiTreeNodeFlags_Selected;
extern "C" const int const_ImGuiTreeNodeFlags_Framed;
extern "C" const int const_ImGuiTreeNodeFlags_AllowItemOverlap;
extern "C" const int const_ImGuiTreeNodeFlags_NoTreePushOnOpen;
extern "C" const int const_ImGuiTreeNodeFlags_NoAutoOpenOnLog;
extern "C" const int const_ImGuiTreeNodeFlags_DefaultOpen;
extern "C" const int const_ImGuiTreeNodeFlags_OpenOnDoubleClick;
extern "C" const int const_ImGuiTreeNodeFlags_OpenOnArrow;
extern "C" const int const_ImGuiTreeNodeFlags_Leaf;
extern "C" const int const_ImGuiTreeNodeFlags_Bullet;
extern "C" const int const_ImGuiTreeNodeFlags_FramePadding;
extern "C" const int const_ImGuiTreeNodeFlags_SpanAvailWidth;
extern "C" const int const_ImGuiTreeNodeFlags_SpanFullWidth;
extern "C" const int const_ImGuiTreeNodeFlags_NavLeftJumpsBackHere;
extern "C" const int const_ImGuiTreeNodeFlags_CollapsingHeader;
// enum ImGuiPopupFlags_
extern "C" const int const_ImGuiPopupFlags_None;
extern "C" const int const_ImGuiPopupFlags_MouseButtonLeft;
extern "C" const int const_ImGuiPopupFlags_MouseButtonRight;
extern "C" const int const_ImGuiPopupFlags_MouseButtonMiddle;
extern "C" const int const_ImGuiPopupFlags_MouseButtonMask_;
extern "C" const int const_ImGuiPopupFlags_MouseButtonDefault_;
extern "C" const int const_ImGuiPopupFlags_NoOpenOverExistingPopup;
extern "C" const int const_ImGuiPopupFlags_NoOpenOverItems;
extern "C" const int const_ImGuiPopupFlags_AnyPopupId;
extern "C" const int const_ImGuiPopupFlags_AnyPopupLevel;
extern "C" const int const_ImGuiPopupFlags_AnyPopup;
// enum ImGuiSelectableFlags_
extern "C" const int const_ImGuiSelectableFlags_None;
extern "C" const int const_ImGuiSelectableFlags_DontClosePopups;
extern "C" const int const_ImGuiSelectableFlags_SpanAllColumns;
extern "C" const int const_ImGuiSelectableFlags_AllowDoubleClick;
extern "C" const int const_ImGuiSelectableFlags_Disabled;
extern "C" const int const_ImGuiSelectableFlags_AllowItemOverlap;
// enum ImGuiComboFlags_
extern "C" const int const_ImGuiComboFlags_None;
extern "C" const int const_ImGuiComboFlags_PopupAlignLeft;
extern "C" const int const_ImGuiComboFlags_HeightSmall;
extern "C" const int const_ImGuiComboFlags_HeightRegular;
extern "C" const int const_ImGuiComboFlags_HeightLarge;
extern "C" const int const_ImGuiComboFlags_HeightLargest;
extern "C" const int const_ImGuiComboFlags_NoArrowButton;
extern "C" const int const_ImGuiComboFlags_NoPreview;
extern "C" const int const_ImGuiComboFlags_HeightMask_;
// enum ImGuiTabBarFlags_
extern "C" const int const_ImGuiTabBarFlags_None;
extern "C" const int const_ImGuiTabBarFlags_Reorderable;
extern "C" const int const_ImGuiTabBarFlags_AutoSelectNewTabs;
extern "C" const int const_ImGuiTabBarFlags_TabListPopupButton;
extern "C" const int const_ImGuiTabBarFlags_NoCloseWithMiddleMouseButton;
extern "C" const int const_ImGuiTabBarFlags_NoTabListScrollingButtons;
extern "C" const int const_ImGuiTabBarFlags_NoTooltip;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyResizeDown;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyScroll;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyMask_;
extern "C" const int const_ImGuiTabBarFlags_FittingPolicyDefault_;
// enum ImGuiTabItemFlags_
extern "C" const int const_ImGuiTabItemFlags_None;
extern "C" const int const_ImGuiTabItemFlags_UnsavedDocument;
extern "C" const int const_ImGuiTabItemFlags_SetSelected;
extern "C" const int const_ImGuiTabItemFlags_NoCloseWithMiddleMouseButton;
extern "C" const int const_ImGuiTabItemFlags_NoPushId;
extern "C" const int const_ImGuiTabItemFlags_NoTooltip;
extern "C" const int const_ImGuiTabItemFlags_NoReorder;
extern "C" const int const_ImGuiTabItemFlags_Leading;
extern "C" const int const_ImGuiTabItemFlags_Trailing;
// enum ImGuiTableFlags_
extern "C" const int const_ImGuiTableFlags_None;
extern "C" const int const_ImGuiTableFlags_Resizable;
extern "C" const int const_ImGuiTableFlags_Reorderable;
extern "C" const int const_ImGuiTableFlags_Hideable;
extern "C" const int const_ImGuiTableFlags_Sortable;
extern "C" const int const_ImGuiTableFlags_NoSavedSettings;
extern "C" const int const_ImGuiTableFlags_ContextMenuInBody;
extern "C" const int const_ImGuiTableFlags_RowBg;
extern "C" const int const_ImGuiTableFlags_BordersInnerH;
extern "C" const int const_ImGuiTableFlags_BordersOuterH;
extern "C" const int const_ImGuiTableFlags_BordersInnerV;
extern "C" const int const_ImGuiTableFlags_BordersOuterV;
extern "C" const int const_ImGuiTableFlags_BordersH;
extern "C" const int const_ImGuiTableFlags_BordersV;
extern "C" const int const_ImGuiTableFlags_BordersInner;
extern "C" const int const_ImGuiTableFlags_BordersOuter;
extern "C" const int const_ImGuiTableFlags_Borders;
extern "C" const int const_ImGuiTableFlags_NoBordersInBody;
extern "C" const int const_ImGuiTableFlags_NoBordersInBodyUntilResize;
extern "C" const int const_ImGuiTableFlags_SizingFixedFit;
extern "C" const int const_ImGuiTableFlags_SizingFixedSame;
extern "C" const int const_ImGuiTableFlags_SizingStretchProp;
extern "C" const int const_ImGuiTableFlags_SizingStretchSame;
extern "C" const int const_ImGuiTableFlags_NoHostExtendX;
extern "C" const int const_ImGuiTableFlags_NoHostExtendY;
extern "C" const int const_ImGuiTableFlags_NoKeepColumnsVisible;
extern "C" const int const_ImGuiTableFlags_PreciseWidths;
extern "C" const int const_ImGuiTableFlags_NoClip;
extern "C" const int const_ImGuiTableFlags_PadOuterX;
extern "C" const int const_ImGuiTableFlags_NoPadOuterX;
extern "C" const int const_ImGuiTableFlags_NoPadInnerX;
extern "C" const int const_ImGuiTableFlags_ScrollX;
extern "C" const int const_ImGuiTableFlags_ScrollY;
extern "C" const int const_ImGuiTableFlags_SortMulti;
extern "C" const int const_ImGuiTableFlags_SortTristate;
extern "C" const int const_ImGuiTableFlags_SizingMask_;
// enum ImGuiTableColumnFlags_
extern "C" const int const_ImGuiTableColumnFlags_None;
extern "C" const int const_ImGuiTableColumnFlags_Disabled;
extern "C" const int const_ImGuiTableColumnFlags_DefaultHide;
extern "C" const int const_ImGuiTableColumnFlags_DefaultSort;
extern "C" const int const_ImGuiTableColumnFlags_WidthStretch;
extern "C" const int const_ImGuiTableColumnFlags_WidthFixed;
extern "C" const int const_ImGuiTableColumnFlags_NoResize;
extern "C" const int const_ImGuiTableColumnFlags_NoReorder;
extern "C" const int const_ImGuiTableColumnFlags_NoHide;
extern "C" const int const_ImGuiTableColumnFlags_NoClip;
extern "C" const int const_ImGuiTableColumnFlags_NoSort;
extern "C" const int const_ImGuiTableColumnFlags_NoSortAscending;
extern "C" const int const_ImGuiTableColumnFlags_NoSortDescending;
extern "C" const int const_ImGuiTableColumnFlags_NoHeaderLabel;
extern "C" const int const_ImGuiTableColumnFlags_NoHeaderWidth;
extern "C" const int const_ImGuiTableColumnFlags_PreferSortAscending;
extern "C" const int const_ImGuiTableColumnFlags_PreferSortDescending;
extern "C" const int const_ImGuiTableColumnFlags_IndentEnable;
extern "C" const int const_ImGuiTableColumnFlags_IndentDisable;
extern "C" const int const_ImGuiTableColumnFlags_IsEnabled;
extern "C" const int const_ImGuiTableColumnFlags_IsVisible;
extern "C" const int const_ImGuiTableColumnFlags_IsSorted;
extern "C" const int const_ImGuiTableColumnFlags_IsHovered;
extern "C" const int const_ImGuiTableColumnFlags_WidthMask_;
extern "C" const int const_ImGuiTableColumnFlags_IndentMask_;
extern "C" const int const_ImGuiTableColumnFlags_StatusMask_;
extern "C" const int const_ImGuiTableColumnFlags_NoDirectResize_;
// enum ImGuiTableRowFlags_
extern "C" const int const_ImGuiTableRowFlags_None;
extern "C" const int const_ImGuiTableRowFlags_Headers;
// enum ImGuiTableBgTarget_
extern "C" const int const_ImGuiTableBgTarget_None;
extern "C" const int const_ImGuiTableBgTarget_RowBg0;
extern "C" const int const_ImGuiTableBgTarget_RowBg1;
extern "C" const int const_ImGuiTableBgTarget_CellBg;
// enum ImGuiFocusedFlags_
extern "C" const int const_ImGuiFocusedFlags_None;
extern "C" const int const_ImGuiFocusedFlags_ChildWindows;
extern "C" const int const_ImGuiFocusedFlags_RootWindow;
extern "C" const int const_ImGuiFocusedFlags_AnyWindow;
extern "C" const int const_ImGuiFocusedFlags_NoPopupHierarchy;
extern "C" const int const_ImGuiFocusedFlags_DockHierarchy;
extern "C" const int const_ImGuiFocusedFlags_RootAndChildWindows;
// enum ImGuiHoveredFlags_
extern "C" const int const_ImGuiHoveredFlags_None;
extern "C" const int const_ImGuiHoveredFlags_ChildWindows;
extern "C" const int const_ImGuiHoveredFlags_RootWindow;
extern "C" const int const_ImGuiHoveredFlags_AnyWindow;
extern "C" const int const_ImGuiHoveredFlags_NoPopupHierarchy;
extern "C" const int const_ImGuiHoveredFlags_DockHierarchy;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenBlockedByPopup;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenBlockedByActiveItem;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenOverlapped;
extern "C" const int const_ImGuiHoveredFlags_AllowWhenDisabled;
extern "C" const int const_ImGuiHoveredFlags_NoNavOverride;
extern "C" const int const_ImGuiHoveredFlags_RectOnly;
extern "C" const int const_ImGuiHoveredFlags_RootAndChildWindows;
// enum ImGuiDockNodeFlags_
extern "C" const int const_ImGuiDockNodeFlags_None;
extern "C" const int const_ImGuiDockNodeFlags_KeepAliveOnly;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingInCentralNode;
extern "C" const int const_ImGuiDockNodeFlags_PassthruCentralNode;
extern "C" const int const_ImGuiDockNodeFlags_NoSplit;
extern "C" const int const_ImGuiDockNodeFlags_NoResize;
extern "C" const int const_ImGuiDockNodeFlags_AutoHideTabBar;
// enum ImGuiDragDropFlags_
extern "C" const int const_ImGuiDragDropFlags_None;
extern "C" const int const_ImGuiDragDropFlags_SourceNoPreviewTooltip;
extern "C" const int const_ImGuiDragDropFlags_SourceNoDisableHover;
extern "C" const int const_ImGuiDragDropFlags_SourceNoHoldToOpenOthers;
extern "C" const int const_ImGuiDragDropFlags_SourceAllowNullID;
extern "C" const int const_ImGuiDragDropFlags_SourceExtern;
extern "C" const int const_ImGuiDragDropFlags_SourceAutoExpirePayload;
extern "C" const int const_ImGuiDragDropFlags_AcceptBeforeDelivery;
extern "C" const int const_ImGuiDragDropFlags_AcceptNoDrawDefaultRect;
extern "C" const int const_ImGuiDragDropFlags_AcceptNoPreviewTooltip;
extern "C" const int const_ImGuiDragDropFlags_AcceptPeekOnly;
// enum ImGuiDataType_
extern "C" const int const_ImGuiDataType_S8;
extern "C" const int const_ImGuiDataType_U8;
extern "C" const int const_ImGuiDataType_S16;
extern "C" const int const_ImGuiDataType_U16;
extern "C" const int const_ImGuiDataType_S32;
extern "C" const int const_ImGuiDataType_U32;
extern "C" const int const_ImGuiDataType_S64;
extern "C" const int const_ImGuiDataType_U64;
extern "C" const int const_ImGuiDataType_Float;
extern "C" const int const_ImGuiDataType_Double;
extern "C" const int const_ImGuiDataType_COUNT;
// enum ImGuiDir_
extern "C" const int const_ImGuiDir_None;
extern "C" const int const_ImGuiDir_Left;
extern "C" const int const_ImGuiDir_Right;
extern "C" const int const_ImGuiDir_Up;
extern "C" const int const_ImGuiDir_Down;
extern "C" const int const_ImGuiDir_COUNT;
// enum ImGuiSortDirection_
extern "C" const int const_ImGuiSortDirection_None;
extern "C" const int const_ImGuiSortDirection_Ascending;
extern "C" const int const_ImGuiSortDirection_Descending;
// enum ImGuiKey_
extern "C" const int const_ImGuiKey_None;
extern "C" const int const_ImGuiKey_Tab;
extern "C" const int const_ImGuiKey_LeftArrow;
extern "C" const int const_ImGuiKey_RightArrow;
extern "C" const int const_ImGuiKey_UpArrow;
extern "C" const int const_ImGuiKey_DownArrow;
extern "C" const int const_ImGuiKey_PageUp;
extern "C" const int const_ImGuiKey_PageDown;
extern "C" const int const_ImGuiKey_Home;
extern "C" const int const_ImGuiKey_End;
extern "C" const int const_ImGuiKey_Insert;
extern "C" const int const_ImGuiKey_Delete;
extern "C" const int const_ImGuiKey_Backspace;
extern "C" const int const_ImGuiKey_Space;
extern "C" const int const_ImGuiKey_Enter;
extern "C" const int const_ImGuiKey_Escape;
extern "C" const int const_ImGuiKey_LeftCtrl;
extern "C" const int const_ImGuiKey_LeftShift;
extern "C" const int const_ImGuiKey_LeftAlt;
extern "C" const int const_ImGuiKey_LeftSuper;
extern "C" const int const_ImGuiKey_RightCtrl;
extern "C" const int const_ImGuiKey_RightShift;
extern "C" const int const_ImGuiKey_RightAlt;
extern "C" const int const_ImGuiKey_RightSuper;
extern "C" const int const_ImGuiKey_Menu;
extern "C" const int const_ImGuiKey_0;
extern "C" const int const_ImGuiKey_1;
extern "C" const int const_ImGuiKey_2;
extern "C" const int const_ImGuiKey_3;
extern "C" const int const_ImGuiKey_4;
extern "C" const int const_ImGuiKey_5;
extern "C" const int const_ImGuiKey_6;
extern "C" const int const_ImGuiKey_7;
extern "C" const int const_ImGuiKey_8;
extern "C" const int const_ImGuiKey_9;
extern "C" const int const_ImGuiKey_A;
extern "C" const int const_ImGuiKey_B;
extern "C" const int const_ImGuiKey_C;
extern "C" const int const_ImGuiKey_D;
extern "C" const int const_ImGuiKey_E;
extern "C" const int const_ImGuiKey_F;
extern "C" const int const_ImGuiKey_G;
extern "C" const int const_ImGuiKey_H;
extern "C" const int const_ImGuiKey_I;
extern "C" const int const_ImGuiKey_J;
extern "C" const int const_ImGuiKey_K;
extern "C" const int const_ImGuiKey_L;
extern "C" const int const_ImGuiKey_M;
extern "C" const int const_ImGuiKey_N;
extern "C" const int const_ImGuiKey_O;
extern "C" const int const_ImGuiKey_P;
extern "C" const int const_ImGuiKey_Q;
extern "C" const int const_ImGuiKey_R;
extern "C" const int const_ImGuiKey_S;
extern "C" const int const_ImGuiKey_T;
extern "C" const int const_ImGuiKey_U;
extern "C" const int const_ImGuiKey_V;
extern "C" const int const_ImGuiKey_W;
extern "C" const int const_ImGuiKey_X;
extern "C" const int const_ImGuiKey_Y;
extern "C" const int const_ImGuiKey_Z;
extern "C" const int const_ImGuiKey_F1;
extern "C" const int const_ImGuiKey_F2;
extern "C" const int const_ImGuiKey_F3;
extern "C" const int const_ImGuiKey_F4;
extern "C" const int const_ImGuiKey_F5;
extern "C" const int const_ImGuiKey_F6;
extern "C" const int const_ImGuiKey_F7;
extern "C" const int const_ImGuiKey_F8;
extern "C" const int const_ImGuiKey_F9;
extern "C" const int const_ImGuiKey_F10;
extern "C" const int const_ImGuiKey_F11;
extern "C" const int const_ImGuiKey_F12;
extern "C" const int const_ImGuiKey_Apostrophe;
extern "C" const int const_ImGuiKey_Comma;
extern "C" const int const_ImGuiKey_Minus;
extern "C" const int const_ImGuiKey_Period;
extern "C" const int const_ImGuiKey_Slash;
extern "C" const int const_ImGuiKey_Semicolon;
extern "C" const int const_ImGuiKey_Equal;
extern "C" const int const_ImGuiKey_LeftBracket;
extern "C" const int const_ImGuiKey_Backslash;
extern "C" const int const_ImGuiKey_RightBracket;
extern "C" const int const_ImGuiKey_GraveAccent;
extern "C" const int const_ImGuiKey_CapsLock;
extern "C" const int const_ImGuiKey_ScrollLock;
extern "C" const int const_ImGuiKey_NumLock;
extern "C" const int const_ImGuiKey_PrintScreen;
extern "C" const int const_ImGuiKey_Pause;
extern "C" const int const_ImGuiKey_Keypad0;
extern "C" const int const_ImGuiKey_Keypad1;
extern "C" const int const_ImGuiKey_Keypad2;
extern "C" const int const_ImGuiKey_Keypad3;
extern "C" const int const_ImGuiKey_Keypad4;
extern "C" const int const_ImGuiKey_Keypad5;
extern "C" const int const_ImGuiKey_Keypad6;
extern "C" const int const_ImGuiKey_Keypad7;
extern "C" const int const_ImGuiKey_Keypad8;
extern "C" const int const_ImGuiKey_Keypad9;
extern "C" const int const_ImGuiKey_KeypadDecimal;
extern "C" const int const_ImGuiKey_KeypadDivide;
extern "C" const int const_ImGuiKey_KeypadMultiply;
extern "C" const int const_ImGuiKey_KeypadSubtract;
extern "C" const int const_ImGuiKey_KeypadAdd;
extern "C" const int const_ImGuiKey_KeypadEnter;
extern "C" const int const_ImGuiKey_KeypadEqual;
extern "C" const int const_ImGuiKey_GamepadStart;
extern "C" const int const_ImGuiKey_GamepadBack;
extern "C" const int const_ImGuiKey_GamepadFaceUp;
extern "C" const int const_ImGuiKey_GamepadFaceDown;
extern "C" const int const_ImGuiKey_GamepadFaceLeft;
extern "C" const int const_ImGuiKey_GamepadFaceRight;
extern "C" const int const_ImGuiKey_GamepadDpadUp;
extern "C" const int const_ImGuiKey_GamepadDpadDown;
extern "C" const int const_ImGuiKey_GamepadDpadLeft;
extern "C" const int const_ImGuiKey_GamepadDpadRight;
extern "C" const int const_ImGuiKey_GamepadL1;
extern "C" const int const_ImGuiKey_GamepadR1;
extern "C" const int const_ImGuiKey_GamepadL2;
extern "C" const int const_ImGuiKey_GamepadR2;
extern "C" const int const_ImGuiKey_GamepadL3;
extern "C" const int const_ImGuiKey_GamepadR3;
extern "C" const int const_ImGuiKey_GamepadLStickUp;
extern "C" const int const_ImGuiKey_GamepadLStickDown;
extern "C" const int const_ImGuiKey_GamepadLStickLeft;
extern "C" const int const_ImGuiKey_GamepadLStickRight;
extern "C" const int const_ImGuiKey_GamepadRStickUp;
extern "C" const int const_ImGuiKey_GamepadRStickDown;
extern "C" const int const_ImGuiKey_GamepadRStickLeft;
extern "C" const int const_ImGuiKey_GamepadRStickRight;
extern "C" const int const_ImGuiKey_ModCtrl;
extern "C" const int const_ImGuiKey_ModShift;
extern "C" const int const_ImGuiKey_ModAlt;
extern "C" const int const_ImGuiKey_ModSuper;
// ImGuiKey_COUNT
extern "C" const int const_ImGuiKey_NamedKey_BEGIN;
extern "C" const int const_ImGuiKey_NamedKey_END;
// extern "C" const int const_ImGuiKey_NamedKey_COUNT;
// ImGuiKey_KeysData_SIZE
extern "C" const int const_ImGuiKey_KeysData_OFFSET;
// enum ImGuiModFlags_
extern "C" const int const_ImGuiModFlags_None;
extern "C" const int const_ImGuiModFlags_Ctrl;
extern "C" const int const_ImGuiModFlags_Shift;
extern "C" const int const_ImGuiModFlags_Alt;
extern "C" const int const_ImGuiModFlags_Super;
// enum ImGuiNavInput_
extern "C" const int const_ImGuiNavInput_Activate;
extern "C" const int const_ImGuiNavInput_Cancel;
extern "C" const int const_ImGuiNavInput_Input;
extern "C" const int const_ImGuiNavInput_Menu;
extern "C" const int const_ImGuiNavInput_DpadLeft;
extern "C" const int const_ImGuiNavInput_DpadRight;
extern "C" const int const_ImGuiNavInput_DpadUp;
extern "C" const int const_ImGuiNavInput_DpadDown;
extern "C" const int const_ImGuiNavInput_LStickLeft;
extern "C" const int const_ImGuiNavInput_LStickRight;
extern "C" const int const_ImGuiNavInput_LStickUp;
extern "C" const int const_ImGuiNavInput_LStickDown;
extern "C" const int const_ImGuiNavInput_FocusPrev;
extern "C" const int const_ImGuiNavInput_FocusNext;
extern "C" const int const_ImGuiNavInput_TweakSlow;
extern "C" const int const_ImGuiNavInput_TweakFast;
extern "C" const int const_ImGuiNavInput_KeyLeft_;
extern "C" const int const_ImGuiNavInput_KeyRight_;
extern "C" const int const_ImGuiNavInput_KeyUp_;
extern "C" const int const_ImGuiNavInput_KeyDown_;
// ImGuiNavInput_COUNT
// enum ImGuiConfigFlags_
extern "C" const int const_ImGuiConfigFlags_None;
extern "C" const int const_ImGuiConfigFlags_NavEnableKeyboard;
extern "C" const int const_ImGuiConfigFlags_NavEnableGamepad;
extern "C" const int const_ImGuiConfigFlags_NavEnableSetMousePos;
extern "C" const int const_ImGuiConfigFlags_NavNoCaptureKeyboard;
extern "C" const int const_ImGuiConfigFlags_NoMouse;
extern "C" const int const_ImGuiConfigFlags_NoMouseCursorChange;
extern "C" const int const_ImGuiConfigFlags_DockingEnable;
extern "C" const int const_ImGuiConfigFlags_ViewportsEnable;
extern "C" const int const_ImGuiConfigFlags_DpiEnableScaleViewports;
extern "C" const int const_ImGuiConfigFlags_DpiEnableScaleFonts;
extern "C" const int const_ImGuiConfigFlags_IsSRGB;
extern "C" const int const_ImGuiConfigFlags_IsTouchScreen;
// enum ImGuiBackendFlags_
extern "C" const int const_ImGuiBackendFlags_None;
extern "C" const int const_ImGuiBackendFlags_HasGamepad;
extern "C" const int const_ImGuiBackendFlags_HasMouseCursors;
extern "C" const int const_ImGuiBackendFlags_HasSetMousePos;
extern "C" const int const_ImGuiBackendFlags_RendererHasVtxOffset;
extern "C" const int const_ImGuiBackendFlags_PlatformHasViewports;
extern "C" const int const_ImGuiBackendFlags_HasMouseHoveredViewport;
extern "C" const int const_ImGuiBackendFlags_RendererHasViewports;
// enum ImGuiCol_
extern "C" const int const_ImGuiCol_Text;
extern "C" const int const_ImGuiCol_TextDisabled;
extern "C" const int const_ImGuiCol_WindowBg;
extern "C" const int const_ImGuiCol_ChildBg;
extern "C" const int const_ImGuiCol_PopupBg;
extern "C" const int const_ImGuiCol_Border;
extern "C" const int const_ImGuiCol_BorderShadow;
extern "C" const int const_ImGuiCol_FrameBg;
extern "C" const int const_ImGuiCol_FrameBgHovered;
extern "C" const int const_ImGuiCol_FrameBgActive;
extern "C" const int const_ImGuiCol_TitleBg;
extern "C" const int const_ImGuiCol_TitleBgActive;
extern "C" const int const_ImGuiCol_TitleBgCollapsed;
extern "C" const int const_ImGuiCol_MenuBarBg;
extern "C" const int const_ImGuiCol_ScrollbarBg;
extern "C" const int const_ImGuiCol_ScrollbarGrab;
extern "C" const int const_ImGuiCol_ScrollbarGrabHovered;
extern "C" const int const_ImGuiCol_ScrollbarGrabActive;
extern "C" const int const_ImGuiCol_CheckMark;
extern "C" const int const_ImGuiCol_SliderGrab;
extern "C" const int const_ImGuiCol_SliderGrabActive;
extern "C" const int const_ImGuiCol_Button;
extern "C" const int const_ImGuiCol_ButtonHovered;
extern "C" const int const_ImGuiCol_ButtonActive;
extern "C" const int const_ImGuiCol_Header;
extern "C" const int const_ImGuiCol_HeaderHovered;
extern "C" const int const_ImGuiCol_HeaderActive;
extern "C" const int const_ImGuiCol_Separator;
extern "C" const int const_ImGuiCol_SeparatorHovered;
extern "C" const int const_ImGuiCol_SeparatorActive;
extern "C" const int const_ImGuiCol_ResizeGrip;
extern "C" const int const_ImGuiCol_ResizeGripHovered;
extern "C" const int const_ImGuiCol_ResizeGripActive;
extern "C" const int const_ImGuiCol_Tab;
extern "C" const int const_ImGuiCol_TabHovered;
extern "C" const int const_ImGuiCol_TabActive;
extern "C" const int const_ImGuiCol_TabUnfocused;
extern "C" const int const_ImGuiCol_TabUnfocusedActive;
extern "C" const int const_ImGuiCol_DockingPreview;
extern "C" const int const_ImGuiCol_DockingEmptyBg;
extern "C" const int const_ImGuiCol_PlotLines;
extern "C" const int const_ImGuiCol_PlotLinesHovered;
extern "C" const int const_ImGuiCol_PlotHistogram;
extern "C" const int const_ImGuiCol_PlotHistogramHovered;
extern "C" const int const_ImGuiCol_TableHeaderBg;
extern "C" const int const_ImGuiCol_TableBorderStrong;
extern "C" const int const_ImGuiCol_TableBorderLight;
extern "C" const int const_ImGuiCol_TableRowBg;
extern "C" const int const_ImGuiCol_TableRowBgAlt;
extern "C" const int const_ImGuiCol_TextSelectedBg;
extern "C" const int const_ImGuiCol_DragDropTarget;
extern "C" const int const_ImGuiCol_NavHighlight;
extern "C" const int const_ImGuiCol_NavWindowingHighlight;
extern "C" const int const_ImGuiCol_NavWindowingDimBg;
extern "C" const int const_ImGuiCol_ModalWindowDimBg;
// extern "C" const int const_ImGuiCol_COUNT;
// enum ImGuiStyleVar_
extern "C" const int const_ImGuiStyleVar_Alpha;
extern "C" const int const_ImGuiStyleVar_DisabledAlpha;
extern "C" const int const_ImGuiStyleVar_WindowPadding;
extern "C" const int const_ImGuiStyleVar_WindowRounding;
extern "C" const int const_ImGuiStyleVar_WindowBorderSize;
extern "C" const int const_ImGuiStyleVar_WindowMinSize;
extern "C" const int const_ImGuiStyleVar_WindowTitleAlign;
extern "C" const int const_ImGuiStyleVar_ChildRounding;
extern "C" const int const_ImGuiStyleVar_ChildBorderSize;
extern "C" const int const_ImGuiStyleVar_PopupRounding;
extern "C" const int const_ImGuiStyleVar_PopupBorderSize;
extern "C" const int const_ImGuiStyleVar_FramePadding;
extern "C" const int const_ImGuiStyleVar_FrameRounding;
extern "C" const int const_ImGuiStyleVar_FrameBorderSize;
extern "C" const int const_ImGuiStyleVar_ItemSpacing;
extern "C" const int const_ImGuiStyleVar_ItemInnerSpacing;
extern "C" const int const_ImGuiStyleVar_IndentSpacing;
extern "C" const int const_ImGuiStyleVar_CellPadding;
extern "C" const int const_ImGuiStyleVar_ScrollbarSize;
extern "C" const int const_ImGuiStyleVar_ScrollbarRounding;
extern "C" const int const_ImGuiStyleVar_GrabMinSize;
extern "C" const int const_ImGuiStyleVar_GrabRounding;
extern "C" const int const_ImGuiStyleVar_TabRounding;
extern "C" const int const_ImGuiStyleVar_ButtonTextAlign;
extern "C" const int const_ImGuiStyleVar_SelectableTextAlign;
extern "C" const int const_ImGuiStyleVar_COUNT;
// enum ImGuiButtonFlags_
extern "C" const int const_ImGuiButtonFlags_None;
extern "C" const int const_ImGuiButtonFlags_MouseButtonLeft;
extern "C" const int const_ImGuiButtonFlags_MouseButtonRight;
extern "C" const int const_ImGuiButtonFlags_MouseButtonMiddle;
extern "C" const int const_ImGuiButtonFlags_MouseButtonMask_;
extern "C" const int const_ImGuiButtonFlags_MouseButtonDefault_;
// enum ImGuiColorEditFlags_
extern "C" const int const_ImGuiColorEditFlags_None;
extern "C" const int const_ImGuiColorEditFlags_NoAlpha;
extern "C" const int const_ImGuiColorEditFlags_NoPicker;
extern "C" const int const_ImGuiColorEditFlags_NoOptions;
extern "C" const int const_ImGuiColorEditFlags_NoSmallPreview;
extern "C" const int const_ImGuiColorEditFlags_NoInputs;
extern "C" const int const_ImGuiColorEditFlags_NoTooltip;
extern "C" const int const_ImGuiColorEditFlags_NoLabel;
extern "C" const int const_ImGuiColorEditFlags_NoSidePreview;
extern "C" const int const_ImGuiColorEditFlags_NoDragDrop;
extern "C" const int const_ImGuiColorEditFlags_NoBorder;
extern "C" const int const_ImGuiColorEditFlags_AlphaBar;
extern "C" const int const_ImGuiColorEditFlags_AlphaPreview;
extern "C" const int const_ImGuiColorEditFlags_AlphaPreviewHalf;
extern "C" const int const_ImGuiColorEditFlags_HDR;
extern "C" const int const_ImGuiColorEditFlags_DisplayRGB;
extern "C" const int const_ImGuiColorEditFlags_DisplayHSV;
extern "C" const int const_ImGuiColorEditFlags_DisplayHex;
extern "C" const int const_ImGuiColorEditFlags_Uint8;
extern "C" const int const_ImGuiColorEditFlags_Float;
extern "C" const int const_ImGuiColorEditFlags_PickerHueBar;
extern "C" const int const_ImGuiColorEditFlags_PickerHueWheel;
extern "C" const int const_ImGuiColorEditFlags_InputRGB;
extern "C" const int const_ImGuiColorEditFlags_InputHSV;
extern "C" const int const_ImGuiColorEditFlags_DefaultOptions_;
extern "C" const int const_ImGuiColorEditFlags_DisplayMask_;
extern "C" const int const_ImGuiColorEditFlags_DataTypeMask_;
extern "C" const int const_ImGuiColorEditFlags_PickerMask_;
extern "C" const int const_ImGuiColorEditFlags_InputMask_;
// enum ImGuiSliderFlags_
extern "C" const int const_ImGuiSliderFlags_None;
extern "C" const int const_ImGuiSliderFlags_AlwaysClamp;
extern "C" const int const_ImGuiSliderFlags_Logarithmic;
extern "C" const int const_ImGuiSliderFlags_NoRoundToFormat;
extern "C" const int const_ImGuiSliderFlags_NoInput;
extern "C" const int const_ImGuiSliderFlags_InvalidMask_;
// enum ImGuiMouseButton_
extern "C" const int const_ImGuiMouseButton_Left;
extern "C" const int const_ImGuiMouseButton_Right;
extern "C" const int const_ImGuiMouseButton_Middle;
extern "C" const int const_ImGuiMouseButton_COUNT;
// enum ImGuiMouseCursor_
extern "C" const int const_ImGuiMouseCursor_None;
extern "C" const int const_ImGuiMouseCursor_Arrow;
extern "C" const int const_ImGuiMouseCursor_TextInput;
extern "C" const int const_ImGuiMouseCursor_ResizeAll;
extern "C" const int const_ImGuiMouseCursor_ResizeNS;
extern "C" const int const_ImGuiMouseCursor_ResizeEW;
extern "C" const int const_ImGuiMouseCursor_ResizeNESW;
extern "C" const int const_ImGuiMouseCursor_ResizeNWSE;
extern "C" const int const_ImGuiMouseCursor_Hand;
extern "C" const int const_ImGuiMouseCursor_NotAllowed;
extern "C" const int const_ImGuiMouseCursor_COUNT;
// enum ImGuiCond_
extern "C" const int const_ImGuiCond_None;
extern "C" const int const_ImGuiCond_Always;
extern "C" const int const_ImGuiCond_Once;
extern "C" const int const_ImGuiCond_FirstUseEver;
extern "C" const int const_ImGuiCond_Appearing;
// enum ImDrawFlags_
extern "C" const int const_ImDrawFlags_None;
extern "C" const int const_ImDrawFlags_Closed;
extern "C" const int const_ImDrawFlags_RoundCornersTopLeft;
extern "C" const int const_ImDrawFlags_RoundCornersTopRight;
extern "C" const int const_ImDrawFlags_RoundCornersBottomLeft;
extern "C" const int const_ImDrawFlags_RoundCornersBottomRight;
extern "C" const int const_ImDrawFlags_RoundCornersNone;
extern "C" const int const_ImDrawFlags_RoundCornersTop;
extern "C" const int const_ImDrawFlags_RoundCornersBottom;
extern "C" const int const_ImDrawFlags_RoundCornersLeft;
extern "C" const int const_ImDrawFlags_RoundCornersRight;
extern "C" const int const_ImDrawFlags_RoundCornersAll;
extern "C" const int const_ImDrawFlags_RoundCornersDefault_;
extern "C" const int const_ImDrawFlags_RoundCornersMask_;
// enum ImDrawListFlags_
extern "C" const int const_ImDrawListFlags_None;
extern "C" const int const_ImDrawListFlags_AntiAliasedLines;
extern "C" const int const_ImDrawListFlags_AntiAliasedLinesUseTex;
extern "C" const int const_ImDrawListFlags_AntiAliasedFill;
extern "C" const int const_ImDrawListFlags_AllowVtxOffset;
// enum ImFontAtlasFlags_
extern "C" const int const_ImFontAtlasFlags_None;
extern "C" const int const_ImFontAtlasFlags_NoPowerOfTwoHeight;
extern "C" const int const_ImFontAtlasFlags_NoMouseCursors;
extern "C" const int const_ImFontAtlasFlags_NoBakedLines;
// enum ImGuiViewportFlags_
extern "C" const int const_ImGuiViewportFlags_None;
extern "C" const int const_ImGuiViewportFlags_IsPlatformWindow;
extern "C" const int const_ImGuiViewportFlags_IsPlatformMonitor;
extern "C" const int const_ImGuiViewportFlags_OwnedByApp;
extern "C" const int const_ImGuiViewportFlags_NoDecoration;
extern "C" const int const_ImGuiViewportFlags_NoTaskBarIcon;
extern "C" const int const_ImGuiViewportFlags_NoFocusOnAppearing;
extern "C" const int const_ImGuiViewportFlags_NoFocusOnClick;
extern "C" const int const_ImGuiViewportFlags_NoInputs;
extern "C" const int const_ImGuiViewportFlags_NoRendererClear;
extern "C" const int const_ImGuiViewportFlags_TopMost;
extern "C" const int const_ImGuiViewportFlags_Minimized;
extern "C" const int const_ImGuiViewportFlags_NoAutoMerge;
extern "C" const int const_ImGuiViewportFlags_CanHostOtherWindows;
// enum ImGuiItemFlags_
extern "C" const int const_ImGuiItemFlags_None;
extern "C" const int const_ImGuiItemFlags_NoTabStop;
extern "C" const int const_ImGuiItemFlags_ButtonRepeat;
extern "C" const int const_ImGuiItemFlags_Disabled;
extern "C" const int const_ImGuiItemFlags_NoNav;
extern "C" const int const_ImGuiItemFlags_NoNavDefaultFocus;
extern "C" const int const_ImGuiItemFlags_SelectableDontClosePopup;
extern "C" const int const_ImGuiItemFlags_MixedValue;
extern "C" const int const_ImGuiItemFlags_ReadOnly;
extern "C" const int const_ImGuiItemFlags_Inputable;
// enum ImGuiItemStatusFlags_
extern "C" const int const_ImGuiItemStatusFlags_None;
extern "C" const int const_ImGuiItemStatusFlags_HoveredRect;
extern "C" const int const_ImGuiItemStatusFlags_HasDisplayRect;
extern "C" const int const_ImGuiItemStatusFlags_Edited;
extern "C" const int const_ImGuiItemStatusFlags_ToggledSelection;
extern "C" const int const_ImGuiItemStatusFlags_ToggledOpen;
extern "C" const int const_ImGuiItemStatusFlags_HasDeactivated;
extern "C" const int const_ImGuiItemStatusFlags_Deactivated;
extern "C" const int const_ImGuiItemStatusFlags_HoveredWindow;
extern "C" const int const_ImGuiItemStatusFlags_FocusedByTabbing;
// enum ImGuiInputTextFlagsPrivate_
extern "C" const int const_ImGuiInputTextFlags_Multiline;
extern "C" const int const_ImGuiInputTextFlags_NoMarkEdited;
extern "C" const int const_ImGuiInputTextFlags_MergedItem;
// enum ImGuiButtonFlagsPrivate_
extern "C" const int const_ImGuiButtonFlags_PressedOnClick;
extern "C" const int const_ImGuiButtonFlags_PressedOnClickRelease;
extern "C" const int const_ImGuiButtonFlags_PressedOnClickReleaseAnywhere;
extern "C" const int const_ImGuiButtonFlags_PressedOnRelease;
extern "C" const int const_ImGuiButtonFlags_PressedOnDoubleClick;
extern "C" const int const_ImGuiButtonFlags_PressedOnDragDropHold;
extern "C" const int const_ImGuiButtonFlags_Repeat;
extern "C" const int const_ImGuiButtonFlags_FlattenChildren;
extern "C" const int const_ImGuiButtonFlags_AllowItemOverlap;
extern "C" const int const_ImGuiButtonFlags_DontClosePopups;
extern "C" const int const_ImGuiButtonFlags_AlignTextBaseLine;
extern "C" const int const_ImGuiButtonFlags_NoKeyModifiers;
extern "C" const int const_ImGuiButtonFlags_NoHoldingActiveId;
extern "C" const int const_ImGuiButtonFlags_NoNavFocus;
extern "C" const int const_ImGuiButtonFlags_NoHoveredOnFocus;
extern "C" const int const_ImGuiButtonFlags_PressedOnMask_;
extern "C" const int const_ImGuiButtonFlags_PressedOnDefault_;
// enum ImGuiComboFlagsPrivate_
extern "C" const int const_ImGuiComboFlags_CustomPreview;
// enum ImGuiSliderFlagsPrivate_
extern "C" const int const_ImGuiSliderFlags_Vertical;
extern "C" const int const_ImGuiSliderFlags_ReadOnly;
// enum ImGuiSelectableFlagsPrivate_
extern "C" const int const_ImGuiSelectableFlags_NoHoldingActiveID;
extern "C" const int const_ImGuiSelectableFlags_SelectOnNav;
extern "C" const int const_ImGuiSelectableFlags_SelectOnClick;
extern "C" const int const_ImGuiSelectableFlags_SelectOnRelease;
extern "C" const int const_ImGuiSelectableFlags_SpanAvailWidth;
extern "C" const int const_ImGuiSelectableFlags_DrawHoveredWhenHeld;
extern "C" const int const_ImGuiSelectableFlags_SetNavIdOnHover;
extern "C" const int const_ImGuiSelectableFlags_NoPadWithHalfSpacing;
// enum ImGuiTreeNodeFlagsPrivate_
extern "C" const int const_ImGuiTreeNodeFlags_ClipLabelForTrailingButton;
// enum ImGuiSeparatorFlags_
extern "C" const int const_ImGuiSeparatorFlags_None;
extern "C" const int const_ImGuiSeparatorFlags_Horizontal;
extern "C" const int const_ImGuiSeparatorFlags_Vertical;
extern "C" const int const_ImGuiSeparatorFlags_SpanAllColumns;
// enum ImGuiTextFlags_
extern "C" const int const_ImGuiTextFlags_None;
extern "C" const int const_ImGuiTextFlags_NoWidthForLargeClippedText;
// enum ImGuiTooltipFlags_
extern "C" const int const_ImGuiTooltipFlags_None;
extern "C" const int const_ImGuiTooltipFlags_OverridePreviousTooltip;
// enum ImGuiLayoutType_
extern "C" const int const_ImGuiLayoutType_Horizontal;
extern "C" const int const_ImGuiLayoutType_Vertical;
// enum ImGuiLogType
extern "C" const int const_ImGuiLogType_None;
extern "C" const int const_ImGuiLogType_TTY;
extern "C" const int const_ImGuiLogType_File;
extern "C" const int const_ImGuiLogType_Buffer;
extern "C" const int const_ImGuiLogType_Clipboard;
// enum ImGuiAxis
extern "C" const int const_ImGuiAxis_None;
extern "C" const int const_ImGuiAxis_X;
extern "C" const int const_ImGuiAxis_Y;
// enum ImGuiPlotType
extern "C" const int const_ImGuiPlotType_Lines;
extern "C" const int const_ImGuiPlotType_Histogram;
// enum ImGuiPopupPositionPolicy
extern "C" const int const_ImGuiPopupPositionPolicy_Default;
extern "C" const int const_ImGuiPopupPositionPolicy_ComboBox;
extern "C" const int const_ImGuiPopupPositionPolicy_Tooltip;
// enum ImGuiDataTypePrivate_
extern "C" const int const_ImGuiDataType_String;
extern "C" const int const_ImGuiDataType_Pointer;
extern "C" const int const_ImGuiDataType_ID;
// enum ImGuiNextWindowDataFlags_
extern "C" const int const_ImGuiNextWindowDataFlags_None;
extern "C" const int const_ImGuiNextWindowDataFlags_HasPos;
extern "C" const int const_ImGuiNextWindowDataFlags_HasSize;
extern "C" const int const_ImGuiNextWindowDataFlags_HasContentSize;
extern "C" const int const_ImGuiNextWindowDataFlags_HasCollapsed;
extern "C" const int const_ImGuiNextWindowDataFlags_HasSizeConstraint;
extern "C" const int const_ImGuiNextWindowDataFlags_HasFocus;
extern "C" const int const_ImGuiNextWindowDataFlags_HasBgAlpha;
extern "C" const int const_ImGuiNextWindowDataFlags_HasScroll;
extern "C" const int const_ImGuiNextWindowDataFlags_HasViewport;
extern "C" const int const_ImGuiNextWindowDataFlags_HasDock;
extern "C" const int const_ImGuiNextWindowDataFlags_HasWindowClass;
// enum ImGuiNextItemDataFlags_
extern "C" const int const_ImGuiNextItemDataFlags_None;
extern "C" const int const_ImGuiNextItemDataFlags_HasWidth;
extern "C" const int const_ImGuiNextItemDataFlags_HasOpen;
// enum ImGuiKeyPrivate_
extern "C" const int const_ImGuiKey_LegacyNativeKey_BEGIN;
extern "C" const int const_ImGuiKey_LegacyNativeKey_END;
extern "C" const int const_ImGuiKey_Gamepad_BEGIN;
extern "C" const int const_ImGuiKey_Gamepad_END;
// enum ImGuiInputEventType
extern "C" const int const_ImGuiInputEventType_None;
extern "C" const int const_ImGuiInputEventType_MousePos;
extern "C" const int const_ImGuiInputEventType_MouseWheel;
extern "C" const int const_ImGuiInputEventType_MouseButton;
extern "C" const int const_ImGuiInputEventType_MouseViewport;
extern "C" const int const_ImGuiInputEventType_Key;
extern "C" const int const_ImGuiInputEventType_Text;
extern "C" const int const_ImGuiInputEventType_Focus;
extern "C" const int const_ImGuiInputEventType_COUNT;
// enum ImGuiInputSource
extern "C" const int const_ImGuiInputSource_None;
extern "C" const int const_ImGuiInputSource_Mouse;
extern "C" const int const_ImGuiInputSource_Keyboard;
extern "C" const int const_ImGuiInputSource_Gamepad;
extern "C" const int const_ImGuiInputSource_Clipboard;
extern "C" const int const_ImGuiInputSource_Nav;
extern "C" const int const_ImGuiInputSource_COUNT;
// enum ImGuiNavReadMode
extern "C" const int const_ImGuiNavReadMode_Down;
extern "C" const int const_ImGuiNavReadMode_Pressed;
extern "C" const int const_ImGuiNavReadMode_Released;
extern "C" const int const_ImGuiNavReadMode_Repeat;
extern "C" const int const_ImGuiNavReadMode_RepeatSlow;
extern "C" const int const_ImGuiNavReadMode_RepeatFast;
// enum ImGuiActivateFlags_
extern "C" const int const_ImGuiActivateFlags_None;
extern "C" const int const_ImGuiActivateFlags_PreferInput;
extern "C" const int const_ImGuiActivateFlags_PreferTweak;
extern "C" const int const_ImGuiActivateFlags_TryToPreserveState;
// enum ImGuiScrollFlags_
extern "C" const int const_ImGuiScrollFlags_None;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleEdgeX;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleEdgeY;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleCenterX;
extern "C" const int const_ImGuiScrollFlags_KeepVisibleCenterY;
extern "C" const int const_ImGuiScrollFlags_AlwaysCenterX;
extern "C" const int const_ImGuiScrollFlags_AlwaysCenterY;
extern "C" const int const_ImGuiScrollFlags_NoScrollParent;
extern "C" const int const_ImGuiScrollFlags_MaskX_;
extern "C" const int const_ImGuiScrollFlags_MaskY_;
// enum ImGuiNavHighlightFlags_
extern "C" const int const_ImGuiNavHighlightFlags_None;
extern "C" const int const_ImGuiNavHighlightFlags_TypeDefault;
extern "C" const int const_ImGuiNavHighlightFlags_TypeThin;
extern "C" const int const_ImGuiNavHighlightFlags_AlwaysDraw;
extern "C" const int const_ImGuiNavHighlightFlags_NoRounding;
// enum ImGuiNavDirSourceFlags_
extern "C" const int const_ImGuiNavDirSourceFlags_None;
extern "C" const int const_ImGuiNavDirSourceFlags_RawKeyboard;
extern "C" const int const_ImGuiNavDirSourceFlags_Keyboard;
extern "C" const int const_ImGuiNavDirSourceFlags_PadDPad;
extern "C" const int const_ImGuiNavDirSourceFlags_PadLStick;
// enum ImGuiNavMoveFlags_
extern "C" const int const_ImGuiNavMoveFlags_None;
extern "C" const int const_ImGuiNavMoveFlags_LoopX;
extern "C" const int const_ImGuiNavMoveFlags_LoopY;
extern "C" const int const_ImGuiNavMoveFlags_WrapX;
extern "C" const int const_ImGuiNavMoveFlags_WrapY;
extern "C" const int const_ImGuiNavMoveFlags_AllowCurrentNavId;
extern "C" const int const_ImGuiNavMoveFlags_AlsoScoreVisibleSet;
extern "C" const int const_ImGuiNavMoveFlags_ScrollToEdgeY;
extern "C" const int const_ImGuiNavMoveFlags_Forwarded;
extern "C" const int const_ImGuiNavMoveFlags_DebugNoResult;
extern "C" const int const_ImGuiNavMoveFlags_FocusApi;
extern "C" const int const_ImGuiNavMoveFlags_Tabbing;
extern "C" const int const_ImGuiNavMoveFlags_Activate;
extern "C" const int const_ImGuiNavMoveFlags_DontSetNavHighlight;
// enum ImGuiNavLayer
extern "C" const int const_ImGuiNavLayer_Main;
extern "C" const int const_ImGuiNavLayer_Menu;
extern "C" const int const_ImGuiNavLayer_COUNT;
// enum ImGuiOldColumnFlags_
extern "C" const int const_ImGuiOldColumnFlags_None;
extern "C" const int const_ImGuiOldColumnFlags_NoBorder;
extern "C" const int const_ImGuiOldColumnFlags_NoResize;
extern "C" const int const_ImGuiOldColumnFlags_NoPreserveWidths;
extern "C" const int const_ImGuiOldColumnFlags_NoForceWithinWindow;
extern "C" const int const_ImGuiOldColumnFlags_GrowParentContentsSize;
// enum ImGuiDockNodeFlagsPrivate_
extern "C" const int const_ImGuiDockNodeFlags_DockSpace;
extern "C" const int const_ImGuiDockNodeFlags_CentralNode;
extern "C" const int const_ImGuiDockNodeFlags_NoTabBar;
extern "C" const int const_ImGuiDockNodeFlags_HiddenTabBar;
extern "C" const int const_ImGuiDockNodeFlags_NoWindowMenuButton;
extern "C" const int const_ImGuiDockNodeFlags_NoCloseButton;
extern "C" const int const_ImGuiDockNodeFlags_NoDocking;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingSplitMe;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingSplitOther;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingOverMe;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingOverOther;
extern "C" const int const_ImGuiDockNodeFlags_NoDockingOverEmpty;
extern "C" const int const_ImGuiDockNodeFlags_NoResizeX;
extern "C" const int const_ImGuiDockNodeFlags_NoResizeY;
extern "C" const int const_ImGuiDockNodeFlags_SharedFlagsInheritMask_;
extern "C" const int const_ImGuiDockNodeFlags_NoResizeFlagsMask_;
extern "C" const int const_ImGuiDockNodeFlags_LocalFlagsMask_;
extern "C" const int const_ImGuiDockNodeFlags_LocalFlagsTransferMask_;
extern "C" const int const_ImGuiDockNodeFlags_SavedFlagsMask_;
// enum ImGuiDataAuthority_
extern "C" const int const_ImGuiDataAuthority_Auto;
extern "C" const int const_ImGuiDataAuthority_DockNode;
extern "C" const int const_ImGuiDataAuthority_Window;
// enum ImGuiDockNodeState
extern "C" const int const_ImGuiDockNodeState_Unknown;
extern "C" const int const_ImGuiDockNodeState_HostWindowHiddenBecauseSingleWindow;
extern "C" const int const_ImGuiDockNodeState_HostWindowHiddenBecauseWindowsAreResizing;
extern "C" const int const_ImGuiDockNodeState_HostWindowVisible;
// enum ImGuiWindowDockStyleCol
extern "C" const int const_ImGuiWindowDockStyleCol_Text;
extern "C" const int const_ImGuiWindowDockStyleCol_Tab;
extern "C" const int const_ImGuiWindowDockStyleCol_TabHovered;
extern "C" const int const_ImGuiWindowDockStyleCol_TabActive;
extern "C" const int const_ImGuiWindowDockStyleCol_TabUnfocused;
extern "C" const int const_ImGuiWindowDockStyleCol_TabUnfocusedActive;
extern "C" const int const_ImGuiWindowDockStyleCol_COUNT;
// enum ImGuiDebugLogFlags_
extern "C" const int const_ImGuiDebugLogFlags_None;
extern "C" const int const_ImGuiDebugLogFlags_EventActiveId;
extern "C" const int const_ImGuiDebugLogFlags_EventFocus;
extern "C" const int const_ImGuiDebugLogFlags_EventPopup;
extern "C" const int const_ImGuiDebugLogFlags_EventNav;
extern "C" const int const_ImGuiDebugLogFlags_EventIO;
extern "C" const int const_ImGuiDebugLogFlags_EventDocking;
extern "C" const int const_ImGuiDebugLogFlags_EventViewport;
extern "C" const int const_ImGuiDebugLogFlags_EventMask_;
extern "C" const int const_ImGuiDebugLogFlags_OutputToTTY;
// enum ImGuiContextHookType
extern "C" const int const_ImGuiContextHookType_NewFramePre;
extern "C" const int const_ImGuiContextHookType_NewFramePost;
extern "C" const int const_ImGuiContextHookType_EndFramePre;
extern "C" const int const_ImGuiContextHookType_EndFramePost;
extern "C" const int const_ImGuiContextHookType_RenderPre;
extern "C" const int const_ImGuiContextHookType_RenderPost;
extern "C" const int const_ImGuiContextHookType_Shutdown;
extern "C" const int const_ImGuiContextHookType_PendingRemoval_;
// enum ImGuiTabBarFlagsPrivate_
extern "C" const int const_ImGuiTabBarFlags_DockNode;
extern "C" const int const_ImGuiTabBarFlags_IsFocused;
extern "C" const int const_ImGuiTabBarFlags_SaveSettings;
// enum ImGuiTabItemFlagsPrivate_
extern "C" const int const_ImGuiTabItemFlags_SectionMask_;
extern "C" const int const_ImGuiTabItemFlags_NoCloseButton;
extern "C" const int const_ImGuiTabItemFlags_Button;
extern "C" const int const_ImGuiTabItemFlags_Unsorted;
extern "C" const int const_ImGuiTabItemFlags_Preview;
// enum ImGuiFreeTypeBuilderFlags
extern "C" const int const_ImGuiFreeTypeBuilderFlags_NoHinting;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_NoAutoHint;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_ForceAutoHint;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_LightHinting;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_MonoHinting;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Bold;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Oblique;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Monochrome;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_LoadColor;
extern "C" const int const_ImGuiFreeTypeBuilderFlags_Bitmap;

//// filedialog
// enum ImGuiFileDialogFlags;
extern "C" const int const_ImGuiFileDialogFlags_None;
extern "C" const int const_ImGuiFileDialogFlags_ConfirmOverwrite;
extern "C" const int const_ImGuiFileDialogFlags_DontShowHiddenFiles;
extern "C" const int const_ImGuiFileDialogFlags_DisableCreateDirectoryButton;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnType;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnSize;
extern "C" const int const_ImGuiFileDialogFlags_HideColumnDate;
extern "C" const int const_ImGuiFileDialogFlags_NoDialog;
extern "C" const int const_ImGuiFileDialogFlags_ReadOnlyFileNameField;
extern "C" const int const_ImGuiFileDialogFlags_CaseInsensitiveExtention;
extern "C" const int const_ImGuiFileDialogFlags_Modal;
// enum IGFD_FileStyle
extern "C" const int const_IGFD_FileStyle_None;
extern "C" const int const_IGFD_FileStyleByTypeFile;
extern "C" const int const_IGFD_FileStyleByTypeDir;
extern "C" const int const_IGFD_FileStyleByTypeLink;
extern "C" const int const_IGFD_FileStyleByExtention;
extern "C" const int const_IGFD_FileStyleByFullName;
extern "C" const int const_IGFD_FileStyleByContainedInFullName;

// The default font for the critic2 GUI
extern "C" const char *const_myfont_ttf_compressed_data_base85_ptr;

#endif CONSTANTS_H
