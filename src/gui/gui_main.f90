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
     0.80_c_float,0.00_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic
     1.00_c_float,0.43_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic, layered
     1.00_c_float,0.79_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic, chain
     0.58_c_float,1.00_c_float,0.00_c_float,0.4_c_float,&  ! 3d periodic, molecular
     0.00_c_float,1.00_c_float,0.96_c_float,0.4_c_float,&  ! slab
     0.00_c_float,0.28_c_float,1.00_c_float,0.4_c_float,&  ! chain
     0.51_c_float,0.51_c_float,0.51_c_float,0.4_c_float,&  ! molecule in a box
     0.56_c_float,0.00_c_float,1.00_c_float,0.4_c_float,&  ! single molecule
     1.00_c_float,0.00_c_float,0.66_c_float,0.4_c_float/),&! molecular cluster
     shape(ColorTableCellBg_def))
  real(c_float), parameter, public :: ColorHighlightScene_def(4) = (/1.0_c_float,1.0_c_float,0.5_c_float,0.7_c_float/) ! hover highlight of atoms in view
  real(c_float), parameter, public :: ColorHighlightSelectScene_def(4) = (/0.2_c_float,0.64_c_float,0.9_c_float,0.7_c_float/) ! edit geometry selection color
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

  ! element colors
  real(c_float), public :: ColorElement(3,0:maxzat0)

  ! flags to control main's behavior
  integer, public :: force_run_commands = 0 ! execute commands from the input console (0=no,1=only selected,2=all)
  logical, public :: force_quit_threads = .false. ! set to true to force all threads to quit as soon as possible

  ! public procedures
  public :: gui_start
  public :: set_default_ui_settings

  interface
     module subroutine gui_start()
     end subroutine gui_start
     module subroutine set_default_ui_settings()
     end subroutine set_default_ui_settings
  end interface

end module gui_main
