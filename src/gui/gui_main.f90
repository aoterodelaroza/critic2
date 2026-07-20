! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! Structure class and routines for basic crystallography computations
module gui_main
  use iso_c_binding, only: c_ptr, c_float, c_int, c_null_ptr, c_bool
  use interfaces_cimgui, only: ImGuiIO, ImGuiContext, ImVec4, ImVec2, ImGuiViewport,&
     ImFontAtlas
  use param, only: maxzat0
  implicit none

  private

  ! DATA TYPES used in the critic2 GUI
  ! + Real variables declared as c_float:
  !   - Parameters associated with ImGui (scrolling step, tooltip delay, font sizes,...)
  !   - Mouse positions and screen positions
  !   - Colors
  !   - Scene graphics parameters (atom scale and radii, border thickness, offsets,...)
  !   - Transformation matrices and related variables (fov, reset distance,...)
  !   - Draw list variables (dl_sphere, dl_cylinder,...)
  ! + Real variables declared as *8:
  !   - Time flags for events
  !   - System coordinates and distances (bond distance, bond factor,...)
  ! + For now, c_int and integer are not distinguished, but should the need arise
  ! this need to be looked into. At some point, the whole program should be converted
  ! to variable type usage conforming to the Fortran standard.

  ! variables to GUI's structures & data
  type(ImGuiIO), pointer, public :: io ! pointer to ImGui's IO object
  type(ImFontAtlas), pointer, public :: fonts ! pointer to IO%Fonts
  type(ImGuiContext), pointer, public :: g ! pointer to ImGui's context
  type(ImGuiViewport), pointer, public :: mainvwp ! pointer to main viewport
  type(c_ptr), public :: rootwin ! the root window pointer (GLFWwindow*)
  type(ImVec2), public :: fontsize ! font size (sensitive to scaling)
  real(c_float), parameter, public :: fontbakesize = 16._c_float ! normal bake size (for GUI)
  real(c_float), parameter, public :: fontbakesize_large = 128._c_float ! large bake size (for rendering)
  type(c_ptr), public :: font_normal, font_large ! GUI and rendering font pointers

  ! GUI control parameters
  ! integer(c_int), parameter, public :: ms_samples = 1 ! number of samples in multisamples
  logical, public :: tooltip_enabled = .true. ! whether tooltips are enabled
  real(c_float), public :: tooltip_delay = 0.5_c_float ! tooltip delay, in seconds
  real(c_float), public :: tooltip_wrap_factor = 40._c_float ! tooltip wrap factor (fontsize)
  logical, parameter, public :: reuse_mid_empty_systems = .false. ! whether to reuse the empty systems in the middle
  logical, public :: tree_select_updates_inpcon = .true. ! selecting in tree chooses system in input console
  logical, public :: tree_select_updates_view = .true. ! selecting in tree chooses system in view
  integer, public :: lockbehavior = 1 ! 0=no-lock, 1=only SCF, 2=all systems

  ! GUI colors, defaults
  real(c_float), parameter, public :: ColorTableCellBg_def(4,0:8) = reshape((/& ! tree table name cell
     0.65_c_float,0.05_c_float,0.05_c_float,0.4_c_float,&  ! 3d periodic (dark red)
     0.90_c_float,0.35_c_float,0.10_c_float,0.4_c_float,&  ! 3d periodic, layered (red-orange)
     1.00_c_float,0.60_c_float,0.10_c_float,0.4_c_float,&  ! 3d periodic, chain (orange)
     1.00_c_float,0.80_c_float,0.25_c_float,0.4_c_float,&  ! 3d periodic, molecular (yellow-orange)
     0.00_c_float,0.75_c_float,0.65_c_float,0.4_c_float,&  ! slab (teal)
     0.25_c_float,0.80_c_float,0.30_c_float,0.4_c_float,&  ! chain (green)
     0.45_c_float,0.55_c_float,0.70_c_float,0.4_c_float,&  ! molecule in a box (blue-grey)
     0.25_c_float,0.45_c_float,0.95_c_float,0.4_c_float,&  ! single molecule (blue)
     0.65_c_float,0.35_c_float,0.95_c_float,0.4_c_float/),&! molecular cluster (purple)
     shape(ColorTableCellBg_def))
  real(c_float), parameter, public :: ColorHighlightScene_def(4) = (/1.0_c_float,1.0_c_float,0.5_c_float,0.7_c_float/) ! hover highlight of atoms in view
  real(c_float), parameter, public :: ColorHighlightSelectScene_def(4) = (/0.2_c_float,0.64_c_float,0.9_c_float,0.7_c_float/) ! edit geometry selection color
  real(c_float), parameter, public :: ColorHighlightBondScene(4) = (/0.3_c_float,0.85_c_float,0.3_c_float,0.7_c_float/) ! bonds tab: bonded-neighbor highlight
  real(c_float), parameter, public :: ColorHighlightBondScene2(4) = (/0.85_c_float,0.3_c_float,0.3_c_float,0.7_c_float/) ! bonds tab: selected bonded-neighbor highlight
  real(c_float), parameter, public :: ColorMeasureSelect_def(4,4) = reshape((/&
     1._c_float,  0.4_c_float, 0.4_c_float, 0.5_c_float,&
     0.4_c_float, 1._c_float,  0.4_c_float, 0.5_c_float,&
     0.4_c_float, 0.4_c_float, 1._c_float, 0.5_c_float,&
     0.9_c_float, 0.7_c_float, 0.4_c_float, 0.5_c_float/),(/4,4/))

  ! GUI colors, actual values
  real(c_float), public :: ColorTableCellBg(4,0:8)
  type(ImVec4), parameter, public :: ColorDialogDir = ImVec4(0.9, 0.9, 0.5, 1.0) ! directories in the dialog
  type(ImVec4), parameter, public :: ColorDialogFile = ImVec4(1.0, 1.0, 1.0, 1.0) ! files in the dialog
  type(ImVec4), parameter, public :: ColorHighlightText = ImVec4(0.2, 0.64, 0.9, 1.0) ! highlighted text
  type(ImVec4), parameter, public :: ColorDangerButton = ImVec4(0.50, 0.08, 0.08, 1.0) ! important button
  type(ImVec4), parameter, public :: ColorDangerText = ImVec4(0.80, 0.08, 0.08, 1.0) ! important text
  type(ImVec4), parameter, public :: ColorWaitBg = ImVec4(0.80, 0.80, 0.80, 0.6) ! dim the background while waiting
  type(ImVec4), parameter, public :: ColorFrameBgAlt = ImVec4(0.29,0.16,0.48,0.54) ! alternate framebg
  type(ImVec4), parameter, public :: ColorFieldSelected = ImVec4(0.91,1.00,0.00,0.31) ! selected field
  type(ImVec4), parameter, public :: ColorTableHighlightRow = ImVec4(1._c_float,0.8_c_float,0.1_c_float,0.5_c_float) ! selectable highlight color
  real(c_float), public :: ColorHighlightScene(4)
  real(c_float), public :: ColorHighlightSelectScene(4)
  real(c_float), public :: ColorMeasureSelect(4,4)

  ! generic black/white constants
  type(ImVec4), parameter, public :: ColorBlack = ImVec4(0.0, 0.0, 0.0, 1.0) ! opaque black
  type(ImVec4), parameter, public :: ColorWhite = ImVec4(1.0, 1.0, 1.0, 1.0) ! opaque white

  ! representation default colors
  real(c_float), parameter, public :: ColorAtomBorder_def(3) = (/0._c_float,0._c_float,0._c_float/) ! atom border (black)
  real(c_float), parameter, public :: ColorOccEmpty_def(3) = (/0.80_c_float,0.80_c_float,0.80_c_float/) ! empty occupancy sector (light gray)
  real(c_float), parameter, public :: ColorBond_def(3) = (/0._c_float,0._c_float,0._c_float/) ! bond (black)
  real(c_float), parameter, public :: ColorBondBorder_def(3) = (/0._c_float,0._c_float,0._c_float/) ! bond border (black)
  real(c_float), parameter, public :: ColorVdwContacts_def(3) = (/0.51_c_float,0.83_c_float,0.11_c_float/) ! vdw contacts bond (green)
  real(c_float), parameter, public :: ColorHbonds_def(3) = (/0.11_c_float,0.44_c_float,0.83_c_float/) ! hydrogen bond (blue)
  real(c_float), parameter, public :: ColorHbondStrong_def(3) = (/0.00_c_float,0.00_c_float,1.00_c_float/) ! strong H-bond (blue)
  real(c_float), parameter, public :: ColorHbondModerate_def(3) = (/0.00_c_float,0.70_c_float,0.00_c_float/) ! moderate H-bond (green)
  real(c_float), parameter, public :: ColorHbondWeak_def(3) = (/1.00_c_float,0.00_c_float,0.00_c_float/) ! weak H-bond (red)
  real(c_float), parameter, public :: ColorLabel_def(3) = (/0._c_float,0._c_float,0._c_float/) ! atom label (black)
  real(c_float), parameter, public :: ColorRotaxis_def(3) = (/0._c_float,0._c_float,0._c_float/) ! rotation axis (black)
  real(c_float), parameter, public :: ColorAxes_def(3,3) = reshape((/& ! canonical axis colors
     1._c_float,0._c_float,0._c_float,&  ! x = red
     0._c_float,1._c_float,0._c_float,&  ! y = green
     0._c_float,0._c_float,1._c_float/),&! z = blue
     shape(ColorAxes_def))

  ! scene default colors
  real(c_float), parameter, public :: ColorSceneBg_def(3) = (/1._c_float,1._c_float,1._c_float/) ! scene background (white)

  ! non-palette color-related numerics
  real(c_float), parameter, public :: ColorRootBg(4) = (/0.45_c_float,0.55_c_float,0.60_c_float,1.00_c_float/) ! root window clear color
  real(c_float), parameter, public :: ColorClearTransparent(4) = (/0._c_float,0._c_float,0._c_float,0._c_float/) ! FBO transparent clear
  real(c_float), parameter, public :: ColorButtonHoverFactor = 1.2_c_float ! hovered button brighten factor
  real(c_float), parameter, public :: ColorButtonActiveFactor = 0.8_c_float ! active button darken factor
  real(c_float), parameter, public :: lumweights(3) = (/0.299_c_float,0.587_c_float,0.114_c_float/) ! perceived-luminance weights

  ! element colors
  real(c_float), public :: ColorElement(3,0:maxzat0)

  ! flags to control main's behavior
  integer, public :: force_run_commands = 0 ! execute commands from the input console (0=no,1=only selected,2=all)
  logical, public, volatile :: force_quit_threads = .false. ! set to true to force all threads to quit as soon as possible (volatile: written by the main thread, polled by the worker)

  ! public procedures
  public :: gui_start
  public :: set_default_ui_settings
  public :: show_tools_menu

  interface
     module subroutine gui_start()
     end subroutine gui_start
     module subroutine set_default_ui_settings()
     end subroutine set_default_ui_settings
     module subroutine show_tools_menu(isys,idparent,ttshown,launchgeometry)
       integer, intent(in) :: isys
       integer, intent(in) :: idparent
       logical, intent(inout) :: ttshown
       logical, intent(inout), optional :: launchgeometry
     end subroutine show_tools_menu
  end interface

end module gui_main
