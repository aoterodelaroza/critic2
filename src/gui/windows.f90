! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! The class to handle ImGui windows.
module windows
  use iso_c_binding
  use representations, only: representation
  use scenes, only: scene
  use interfaces_cimgui, only: ImVec2
  use global, only: rborder_def
  use param, only: isformat_r_unknown, eye, mlen
  implicit none

  private

  ! the buffer for the output console
  character(kind=c_char,len=:), allocatable, target :: outputb
  integer(c_size_t) :: lob = 0
  integer(c_size_t), parameter :: maxlob = 2000000
  real*8 :: timelast_output_written = 0d0 ! time the output was last written

  ! the buffer for the input console
  character(kind=c_char,len=:), allocatable, target :: inputb
  integer(c_size_t), parameter :: maxlib = 40000

  ! view modes (positive = normal, user-selectable; negative = forced).
  integer, parameter, public :: vm_builder_bondorder = -11 ! forced by builder: cycle the bond order (persistent)
  integer, parameter, public :: vm_builder_bondremove = -10 ! forced by builder: remove bonds (persistent)
  integer, parameter, public :: vm_builder_bondh = -9 ! forced by builder: create bonds, dropping a hydrogen (persistent)
  integer, parameter, public :: vm_builder_bond = -8 ! forced by builder: create bonds (persistent)
  integer, parameter, public :: vm_builder_trim = -7 ! forced by builder: trim branches (persistent)
  integer, parameter, public :: vm_builder_addfragment = -6 ! forced by builder: add fragments (persistent)
  integer, parameter, public :: vm_builder_addatom = -5 ! forced by builder: add atoms (persistent)
  integer, parameter, public :: vm_builder_remove = -4 ! forced by builder: remove atoms (persistent)
  integer, parameter, public :: vm_builder_valence = -3 ! forced by builder: change valence (persistent)
  integer, parameter, public :: vm_mdinteract = -2 ! forced during interactive dynamics
  integer, parameter, public :: vm_pick_atom = -1 ! forced by a window awaiting an atom pick
  integer, parameter, public :: vm_navigate  = 0
  integer, parameter, public :: vm_select    = 1
  integer, parameter, public :: vm_movemol   = 2
  integer, parameter, public :: vm_moveatom  = 3
  integer, parameter, public :: vm_NUM = 3 ! highest user-selectable mode (combo)
  integer, parameter, public :: vm_builder_lo = vm_builder_bondorder ! lower bound of the builder-mode range
  integer, parameter, public :: vm_builder_hi = vm_builder_valence ! upper bound of the builder-mode range

  character(len=17), parameter, public :: vmnames(vm_builder_lo:vm_NUM) = (/&
     "Bond Order       ",& ! vm_builder_bondorder
     "Remove Bonds     ",& ! vm_builder_bondremove
     "Create Bonds (-H)",& ! vm_builder_bondh
     "Create Bonds     ",& ! vm_builder_bond
     "Trim Branches    ",& ! vm_builder_trim
     "Add Fragments    ",& ! vm_builder_addfragment
     "Add Atoms        ",& ! vm_builder_addatom
     "Remove Atoms     ",& ! vm_builder_remove
     "Change Valence   ",& ! vm_builder_valence
     "Interact (MD)    ",& ! vm_mdinteract
     "Pick Atoms       ",& ! vm_pick_atom
     "Navigate         ",& ! vm_navigate
     "Select           ",& ! vm_select
     "Move Molecules   ",& ! vm_movemol
     "Move Atoms       "&  ! vm_moveatom
     /)

  ! A staged atom pick for the operations that need two atoms chosen
  ! one after the other (create bond, geometry add-bond, editrep bond
  ! anchor). Holds the first pick and the time it was staged (or the
  ! pick session was armed), so a geometry change in between
  ! invalidates the stored cell-atom index.
  type pairpick
     integer :: idx(4) = 0 ! staged atom (cell index + lattice vector; 0 = none)
     real*8 :: time = 0d0 ! time the atom was staged or the pick was armed
   contains
     procedure :: arm => pairpick_arm ! stamp the time without staging an atom
     procedure :: stage => pairpick_stage ! store the atom and stamp the time
     procedure :: clear => pairpick_clear ! forget the staged atom and the time
     procedure :: is_staged => pairpick_is_staged ! an atom is staged
     procedure :: is_stale => pairpick_is_stale ! the geometry changed since the stamp
     procedure :: same => pairpick_same ! the staged atom equals this pick
  end type pairpick

  ! view mode data structure for window_forced modes
  type viewmode_data
     character(len=:), allocatable :: msg ! caller-supplied prompt shown in the view bar (pick-atom mode only)
     integer :: tooltip_iz = 0 ! atomic number of the element shown at the cursor (add-atoms mode; <= 0 = none)
     character(len=:), allocatable :: tooltip_frag ! fragment name shown at the cursor (add-fragments mode)
     logical :: frag_isligand = .false. ! the fragment bonds to the clicked atom instead of replacing it (add-fragments mode)
     integer(c_int) :: idx(4) = 0 ! atom identifier under mouse position (cell index + lattice vector)
     integer(c_int) :: bidx(5) = 0 ! bond identifier under mouse position (two cell indices + lattice vector)
     integer :: flag = 0 ! pick bind that fired (1 = main, 2 = alternate)
     real(c_float) :: xpos(2) = 0._c_float ! texture position of the click (add-atoms mode)
     integer :: owner = 0 ! owner window ID
  end type viewmode_data

  ! user data for the file open dialog
  type :: dialog_userdata
     type(c_ptr) :: dptr = c_null_ptr ! the pointer for the file dialog
     integer(c_int) :: mol = -1 ! -1 = auto, 0 = crystal, 1 = molecule
     logical :: showhidden = .false._c_bool ! show hidden files
     integer(c_int) :: isformat = isformat_r_unknown ! force structure format
     logical :: readlastonly = .false._c_bool ! read only the last structure
     integer(c_int) :: purpose ! the purpose of the dialog
     logical :: molcubic = .false. ! whether to read the cell as cubic in a molecule
     real(c_float) :: rborder = real(rborder_def,c_float) ! border for the cell in a molecule
  end type dialog_userdata

  ! Wrapper class to handle ImGui windows
  type window
     ! global window parameters
     logical :: isinit = .false. ! whether this window has been initialized
     logical :: permanent = .false. ! permanent = do not steal window's ID if it is closed
     logical :: firstpass = .true. ! true if this is the first draw() call
     logical(c_bool) :: isopen ! whether the window is open
     integer :: type ! the window type
     integer(c_int) :: id ! internal ID for this window (index from win(:))
     integer :: idparent = 0 ! internal ID (from win(:)) for the caller window
     integer :: itoken = 0 ! parent token for some dialogs functions
     integer :: purpose ! purpose of the window at creation (dialogs: wpurp_dialog_*; views: wpurp_view_*)
     integer(c_int) :: flags ! window flags
     character(kind=c_char,len=:), allocatable :: name ! name of the window
     character(kind=c_char,len=:), allocatable :: errmsg ! error message in the window
     type(c_ptr) :: ptr ! ImGuiWindow* pointer (use only after Begin())
     type(c_ptr) :: dptr ! ImGuiFileDialog* pointer for dialogs
     integer :: isys = 1 ! the system on which the window operates
     integer :: irep = 0 ! the representation on which the window operates
     logical :: tied_to_tree = .false. ! whether the system in this window is tied to the tree system
     real(c_float) :: pos(2) = (/0._c_float,0._c_float/) ! the position of the window's top left corner
     logical :: isdocked = .false. ! whether the window is docked
     integer, allocatable :: iord(:) ! table order (multiple windows)
     integer :: lastselected = 0 ! selectable, last element selected (multiple windows)
     character(len=:), allocatable :: tabselected ! which tab is selected (multiple windows)
     integer(c_int) :: sortcid = 0 ! sort the window's table by this column id (tree, geometry)
     integer(c_int) :: sortdir = 1 ! sort the window's table with this direction (tree, geometry)
     real*8 :: timelast_sort = 0d0 ! time the window's table was last sorted (tree, geometry)
     real*8 :: timelast_assign = 0d0 ! time this window was last assigned a system (tree, view, inpcon)
     real*8 :: timelast_focused = 0d0 ! time this window was last focused (close-dialog bind)
     ! tree table parameters
     integer :: forceselect = 0 ! make the tree select this system in the next pass
     real*8 :: timelast_tree_update = 0d0 ! time the tree was last updated
     real*8 :: timelast_tree_resize = 0d0 ! time the tree columnes were last resized
     ! view parameters
     logical :: ismain ! whether this is the main view or an alternate
     type(scene), pointer :: sc ! pointer to the view scene
     integer(c_int) :: FBO ! draw framebuffer
     integer(c_int) :: FBOtex ! draw framebuffer -> texture
     integer(c_int) :: FBOdepth ! draw framebuffer -> depth buffer
     integer(c_int) :: FBOpick ! object pick framebuffer
     integer(c_int) :: FBOdepthp ! object pick framebuffer -> depth buffer
     integer(c_int) :: FBOrgba ! object pick framebuffer -> rgba buffer
     integer(c_int) :: FBOside ! side of the render texture (pixels)
     type(ImVec2) :: v_rmin, v_rmax ! view image rectangle
     logical :: forcerender = .true. ! force render of the scene
     logical :: lowresrender = .false. ! last render was at reduced (interactive) resolution
     logical :: pickbonds_last = .false. ! whether the last pick render included the bonds
     integer :: viewmode = vm_navigate ! view mode (see vm_* above)
     logical :: viewmode_transient = .false. ! true if view mode is transient (resets every frame)
     type(viewmode_data) :: vmdata ! data associated with window_forced view modes
     type(ImVec2) :: mousepos_lastpick ! mouse position at the last atom pick
     integer(c_int) :: mousepos_idx(5) ! identifier for the atom under mouse position
     integer(c_int) :: mousepos_bidx(5) = 0 ! identifier for the bond under mouse position (bond pick modes only)
     type(ImVec2) :: mposlast ! mouse parameters ----v
     real(c_float) :: mpos0_r(3), mpos0_l(3), mpos0_m(3), cpos0_l(3), cpos0_m(3)
     real(c_float) :: oldview(4,4)
     real(c_float) :: mpos0_s
     integer :: ilock = 0 ! mouse parameters -^
     integer :: moveobj_icel = 0 ! cell atom in a move (molecules/atoms) drag
     integer :: moveobj_imol = 0 ! molecule in a move drag
     logical :: moveobj_isdiscrete = .false. ! whether the move fragment is discrete
     logical :: selrect_active = .false. ! rubber-band selection drag in progress (vm_select)
     type(ImVec2) :: press_p0 ! press position (mouse/screen coords): click-vs-drag test and rubber-band anchor
     integer :: measure_pend = 0 ! pending press capture (0=none, 1=measure add, 2=measure delete, 3=forced-mode pick, 4=alternate pick)
     integer :: measure_pend_idx(5) = 0 ! atom captured at the press for the pending measurement/pick
     integer :: measure_pend_bidx(5) = 0 ! bond captured at the press for the pending pick
     real*8 :: timelast_view_getpixel = 0d0 ! time the pick buffer was last queried for atom ID
     ! dialog parameters
     type(dialog_userdata) :: dialog_data ! for the side pane callback
     ! input console parameters
     ! output console parameters
     ! new structure from library & save image
     character(kind=c_char,len=:), allocatable :: okfile ! ok file
     character(kind=c_char,len=:), allocatable :: okfilter ! ok filter
     logical :: okfile_set = .false. ! whether the library file has been set by the user
     logical :: okfile_read = .false. ! whether the structure list should be re-read from the lib
     integer(c_int) :: okfile_format = 0 ! the file format
     ! load field parameters
     ! scf plot and tree plot parameters
     real(c_double) :: ymin, ymax ! y-end of the plot
     integer :: plotn ! number of plot data
     real(c_double), allocatable :: plotx(:), ploty(:) ! plot data
     ! edit representation parameters
     type(representation), pointer :: rep => NULL() ! the representation on which the e.r. window operates
     real*8 :: timelast_plot_update = 0d0 ! time the plot was last updaed
     integer :: editrep_pick_item = 0 ! text/measurement item waiting for an atom pick (0 = idle)
     integer :: editrep_pick_slot = 0 ! anchor or measurement atom the pick will fill
     type(pairpick) :: editrep_pick ! pick session stamp + staged first bond atom (committed when the pair completes)
     ! export image parameters
     integer(c_int) :: nsample ! number of samples for anti-aliasing
     integer(c_int) :: jpgquality ! jpg quality
     logical :: exportview ! export viewport or whole texture
     integer(c_int) :: npixel ! number of pixels in the export buffer
     logical :: transparentbg ! transparent background
     ! vibrations parameters
     integer(c_int) :: ifrequnit = 0 ! frequency unit (0 = cm-1, 1 = THz)
     integer(c_int) :: iqptunit = 0 ! qpt unit (0 = fract, 1 = Cartesian (1/bohr), 2 = Cartesian (1/ang))
     ! geometry parameters
     integer(c_int) :: geometry_atomtype = 1 ! coord type for atoms [0=spc, 1=nneq, 2=ncel(frac), 3=ncel(cart)]
     integer(c_int) :: geometry_moltype = 3 ! coord type for molecules [0=spc, 1=nneq, 2=ncel(frac), 3=ncel(cart)]
     logical :: geometry_forcewyc = .true. ! force respecting symmetry when displacing nneq
     logical :: geometry_keepbonding = .true. ! preserve bonding through atom-preserving transformations
     integer :: geometry_euler_drag_mol = 0 ! molecule whose Euler angles are being dragged (0 = none)
     real*8 :: geometry_euler_drag_val(3) = 0d0 ! unwrapped Euler angles (degrees) during an active drag
     real(c_float) :: geometry_select_rgba(4) ! highlight color
     real*8 :: geometry_input_coord(3) = 0d0 ! coordinates for the new atom in add button
     integer :: geometry_input_species = 1 ! species for the new atom in add button
     type(pairpick) :: geometry_addbond ! staged atom waiting for an add-bond pick
     integer :: geometry_addbond_iview = 0 ! view window commanded for the add-bond pick
     ! builder parameters
     integer :: builder_vm = 0 ! forced mode commanded to the parent view (0 = idle, else one of the vm_builder_* modes)
     integer :: builder_isys = 0 ! system latched for the builder picks (0 = no mode active)
     real*8 :: builder_time = 0d0 ! time of the last click-free poll (stale-click guard)
     type(pairpick) :: builder_bond ! create bonds: staged first atom
     integer :: builder_addatom_z = 6 ! add atoms: selected element (Z)
     integer :: builder_addatom_ig = 6 ! add atoms: local geometry (0-based combo index, 6 = tetrahedral)
     character(len=:), allocatable :: builder_frag_name ! add fragments: name of the selected fragment
     integer :: builder_frag_nat = 0 ! add fragments: number of atoms in the fragment
     integer, allocatable :: builder_frag_z(:) ! add fragments: atomic numbers
     real*8, allocatable :: builder_frag_x(:,:) ! add fragments: Cartesian coordinates (bohr)
     integer :: builder_frag_ianchor = 0 ! add fragments: anchor atom (attaches to the structure)
     integer :: builder_frag_iattach = 0 ! add fragments: placeholder atom marking the attachment direction
     real*8 :: builder_frag_radius = 0d0 ! add fragments: fictitious covalent radius, bohr (ligand; <= 0 = substituent)
     logical :: edit_pending = .false. ! keybinding request to toggle an edit session
     integer :: edit_kind = 0 ! active edit session and its number of atoms: 0 = none, 2 = distance, 3 = angle, 4 = dihedral
     integer :: edit_isys = 0 ! system latched for the edit session
     integer :: edit_idx(4,4) = 0 ! latched atoms (cell atom + lattice vector); for angles, column 2 is the vertex; for dihedrals, columns 2-3 are the axis
     integer :: edit_imove(4) = 0 ! per-atom move mode (0 = fixed, 1 = atom, 2 = fragment/group, 3 = dihedral half, moving the 2-3 substituents)
     logical :: edit_fragok(4) = .false. ! whether the fragment/group move option is available per atom
     integer :: edit_nfrag(4) = 0 ! number of atoms in the masked fragment of each atom
     integer, allocatable :: edit_frag(:,:) ! masked-fragment cell atoms of the latched atoms; for dihedrals, columns 2-3 hold the two halves of the severed 2-3 bond
     logical :: edit_dirty = .false. ! in-place moves applied but not yet committed (rebuild pending)
     real*8 :: edit_time = 0d0 ! time of the last edit-free poll (external-change guard)
     character(len=:), allocatable :: geometry_expression ! expression for column in atoms table
     logical :: geometry_expression_ok = .false. ! is the expression ok?
     character(len=:), allocatable :: geometry_expr_error ! error for expression
     real*8 :: timelast_geometry_clearhighlights = 0d0 ! time the highlights were last cleared
     logical :: geometry_cell_simple = .true. ! simple (integer multiples of a/b/c) vs. full matrix transformation
     integer(c_int) :: geometry_cell_nrep(3) = (/1_c_int,1_c_int,1_c_int/) ! a/b/c multiples for the simple transformation
     integer(c_int) :: geometry_cell_intmat(3,3) = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/)) ! integer part of the full transformation matrix
     integer :: geometry_cell_cen(3) = 1 ! per-row centering index added to the full transformation matrix (1 = zero centering)
     real*8 :: geometry_cell_origin(3) = 0d0 ! origin shift for the matrix transformation
     integer(c_int) :: geometry_cell_inice = 10 ! maximum supercell size for the nice cell search
     real*8, allocatable :: geometry_cell_nice_rmax(:) ! inscribed-sphere radii from the nice search
     real*8, allocatable :: geometry_cell_nice_mmax(:,:,:) ! transformation matrices from the nice search
     character(len=mlen), allocatable :: geometry_sym_ops(:) ! cached symmetry operations (crystallographic notation)
     character(len=mlen), allocatable :: geometry_sym_hm(:) ! cached Hermann-Mauguin symbols of the operations
     real*8, allocatable :: geometry_sym_axes(:,:) ! cached rotation axes (crystallographic coordinates)
     logical, allocatable :: geometry_sym_sel(:) ! per-operation selection flags (symmetry tab)
     integer :: geometry_sym_selgen = 0 ! generation counter, bumped when the symmetry selection changes
     real*8, allocatable :: geometry_sym_analyze_eps(:) ! tolerances from the symmetry-vs-epsilon analysis
     character(len=11), allocatable :: geometry_sym_analyze_sym(:) ! space group symbols from the analysis
     integer, allocatable :: geometry_sym_analyze_num(:) ! space group numbers from the analysis
     ! preferences parameters
     logical :: color_preferences_reset_reps = .true. ! whether changing the element colors resets current representations
     ! water cluster demonstration parameters
     integer(c_int) :: wc_nwat = 12 ! number of water molecules to generate
     integer(c_int) :: wc_placement = 3 ! initial placement of the monomers (0 = random, 1 = row, 2 = ring, 3 = flat ring)
     integer(c_int) :: wc_mode = 1 ! run mode (0 = continuous, 1 = timed)
     integer(c_int) :: wc_time = 60 ! time limit in the timed mode (s)
     logical :: wc_started = .false. ! whether the relaxation has been auto-started for the cluster (isys)
     logical :: wc_clockrun = .false. ! timed mode: whether the clock is ticking (settling included)
     logical :: wc_timeup = .false. ! timed mode: whether the player's time is over
     real*8 :: wc_clock0 = 0d0 ! timed mode: time at which the clock started (s)
     real*8 :: wc_elapsed = 0d0 ! timed mode: elapsed time shown on the clock (s)
     character(len=:), allocatable :: wc_name ! participant name, shown on the scene near the score
   contains
     procedure :: init => window_init ! initialize the window
     procedure :: end => window_end ! finalize the window
     procedure :: focused => window_focused ! whether the root window is focused (even if not current)
     procedure :: draw => window_draw ! draw the window, calls one of the draw commands below
     ! tree procedures
     procedure :: draw_tree
     procedure :: remap_tree ! update the system information shown by the tree
     procedure :: sort_tree ! sort the systems in the tree
     procedure :: reassign_tree ! reassign the currently selected system in the tree
     procedure :: select_system_tree ! select a system in the tree
     ! treeplot procedures
     procedure :: draw_treeplot
     ! view procedures
     procedure :: draw_view
     procedure :: create_texture_view ! create the texture for the view
     procedure :: delete_texture_view ! delete the texture for the view
     procedure :: viewmode_set_mode ! set the viewmode based on user keypresses
     procedure :: viewmode_set_forced ! enter a window-forced pick view mode (pick an atom, builder)
     procedure :: viewmode_release_forced ! cancel the forced view mode if owned by the caller
     procedure :: viewmode_exit_forced ! exit any forced view mode, back to navigation
     procedure :: viewmode_bar_display ! bar display for the current view mode
     procedure :: viewmode_activate_picking ! activate picking by view mode
     procedure :: viewmode_process_events ! process mouse events according to view mode
     procedure :: select_view ! select the system to show in a view
     procedure :: draw_cursor_overlay ! draw the overlay at the mouse cursor (mode icon + measurement)
     procedure :: mousepos_to_texpos ! mouse position to texture position
     procedure :: texpos_to_mousepos ! texture position to mouse position
     procedure :: getpixel ! get depth or color from texture position
     procedure :: view_to_texpos ! view coordinates to texture position
     procedure :: texpos_to_view ! texture position to view coordinates
     procedure :: world_to_texpos ! world coordinates to texture position
     procedure :: texpos_to_world ! texture position to world coordinates
     procedure :: export_to_image ! export current scene to an image file
     procedure :: export_to_png_simple ! export to a png file with the same name as system
     ! dialog procedures
     procedure :: draw_dialog
     ! input console procedures
     procedure :: draw_ci
     procedure :: run_commands_ci ! run the commands in the input buffer
     procedure :: block_gui_ci ! block the GUI while the input buffer commands are run
     procedure :: get_input_details_ci ! get the system & field strings for current input
     ! output console procedures
     procedure :: update_co
     procedure :: draw_co
     ! about window
     procedure :: draw_about
     ! new structure procedures
     procedure :: draw_new_struct
     ! new structure from library procedures
     procedure :: draw_new_struct_library
     ! load field
     procedure :: draw_load_field
     ! scf plot
     procedure :: update_scfplot
     procedure :: draw_scfplot
     ! edit representation
     procedure :: update_editrep
     procedure :: draw_editrep
     procedure :: draw_editrep_atoms
     procedure :: draw_editrep_unitcell
     procedure :: draw_editrep_axes
     procedure :: draw_editrep_symelem
     procedure :: draw_editrep_text
     procedure :: draw_editrep_measure
     ! export image
     procedure :: draw_exportimage
     ! vibrations
     procedure :: draw_vibrations
     ! dynamics
     procedure :: draw_dynamics
     ! water cluster demonstration
     procedure :: draw_water_cluster
     ! geometry
     procedure :: update_geometry
     procedure :: draw_geometry
     ! preferences
     procedure :: draw_preferences
     ! builder window
     procedure :: draw_builder
     procedure :: edit_stop
  end type window
  public :: window

  ! the window stack and named windows
  integer, public :: nwin
  type(window), allocatable, target, public :: win(:)
  integer, public :: iwin_tree
  integer, public :: iwin_console_input
  integer, public :: iwin_console_output
  integer, public :: iwin_view
  integer, public :: iwin_about
  public :: windows_init
  public :: vm_is_forcedpick
  public :: vm_exits_on_empty
  public :: draw_ff_backend_combo
  public :: paste_clipboard_fragment

  ! window types
  integer, parameter, public :: wintype_tree = 1
  integer, parameter, public :: wintype_view = 2
  integer, parameter, public :: wintype_console_input = 3
  integer, parameter, public :: wintype_console_output = 4
  integer, parameter, public :: wintype_dialog = 5
  integer, parameter, public :: wintype_new_struct = 6
  integer, parameter, public :: wintype_new_struct_library = 7
  integer, parameter, public :: wintype_load_field = 8
  integer, parameter, public :: wintype_scfplot = 9
  integer, parameter, public :: wintype_editrep = 10
  integer, parameter, public :: wintype_exportimage = 11
  integer, parameter, public :: wintype_vibrations = 12
  integer, parameter, public :: wintype_about = 13
  integer, parameter, public :: wintype_preferences = 14
  integer, parameter, public :: wintype_treeplot = 15
  integer, parameter, public :: wintype_geometry = 16
  integer, parameter, public :: wintype_builder = 17
  integer, parameter, public :: wintype_dynamics = 18
  integer, parameter, public :: wintype_water_cluster = 19

  ! window purposes
  integer, parameter, public :: wpurp_unknown = 0
  integer, parameter, public :: wpurp_dialog_openfiles = 1
  integer, parameter, public :: wpurp_dialog_savelogfile = 2
  integer, parameter, public :: wpurp_dialog_openlibraryfile = 3
  integer, parameter, public :: wpurp_dialog_openfieldfile = 4
  integer, parameter, public :: wpurp_dialog_openonefilemodal = 5
  integer, parameter, public :: wpurp_dialog_saveimagefile = 6
  integer, parameter, public :: wpurp_dialog_openvibfile = 7
  integer, parameter, public :: wpurp_view_main = 8
  integer, parameter, public :: wpurp_view_alternate = 9

  ! column ids for the table in the tree widget
  integer(c_int), parameter, public :: ic_tree_closebutton = 0
  integer(c_int), parameter, public :: ic_tree_expandbutton = 1
  integer(c_int), parameter, public :: ic_tree_id = 2
  integer(c_int), parameter, public :: ic_tree_format = 3
  integer(c_int), parameter, public :: ic_tree_props = 4
  integer(c_int), parameter, public :: ic_tree_name = 5
  integer(c_int), parameter, public :: ic_tree_spg = 6
  integer(c_int), parameter, public :: ic_tree_v = 7
  integer(c_int), parameter, public :: ic_tree_vmol = 8
  integer(c_int), parameter, public :: ic_tree_nneq = 9
  integer(c_int), parameter, public :: ic_tree_ncel = 10
  integer(c_int), parameter, public :: ic_tree_nmol = 11
  integer(c_int), parameter, public :: ic_tree_zprime = 12
  integer(c_int), parameter, public :: ic_tree_a = 13
  integer(c_int), parameter, public :: ic_tree_b = 14
  integer(c_int), parameter, public :: ic_tree_c = 15
  integer(c_int), parameter, public :: ic_tree_alpha = 16
  integer(c_int), parameter, public :: ic_tree_beta = 17
  integer(c_int), parameter, public :: ic_tree_gamma = 18
  integer(c_int), parameter, public :: ic_tree_e = 19
  integer(c_int), parameter, public :: ic_tree_emol = 20
  integer(c_int), parameter, public :: ic_tree_p = 21
  integer(c_int), parameter, public :: ic_tree_NUMCOLUMNS = 22 ! keep up to date

  ! routines to manipulate the window stack
  public :: stack_realloc_maybe
  public :: stack_create_window
  public :: regenerate_window_pointers
  public :: read_output_uout
  public :: fill_input_ci

  !xx! Interfaces
  interface
     !xx! proc submodule !xx!
     module subroutine pairpick_arm(pp)
       class(pairpick), intent(inout) :: pp
     end subroutine pairpick_arm
     module subroutine pairpick_stage(pp,idx)
       class(pairpick), intent(inout) :: pp
       integer, intent(in) :: idx(4)
     end subroutine pairpick_stage
     module subroutine pairpick_clear(pp)
       class(pairpick), intent(inout) :: pp
     end subroutine pairpick_clear
     pure module function pairpick_is_staged(pp)
       class(pairpick), intent(in) :: pp
       logical :: pairpick_is_staged
     end function pairpick_is_staged
     pure module function pairpick_is_stale(pp,timelastchange)
       class(pairpick), intent(in) :: pp
       real*8, intent(in) :: timelastchange
       logical :: pairpick_is_stale
     end function pairpick_is_stale
     pure module function pairpick_same(pp,idx)
       class(pairpick), intent(in) :: pp
       integer, intent(in) :: idx(4)
       logical :: pairpick_same
     end function pairpick_same
     pure module function vm_is_forcedpick(mode)
       integer, intent(in) :: mode
       logical :: vm_is_forcedpick
     end function vm_is_forcedpick
     pure module function vm_exits_on_empty(mode)
       integer, intent(in) :: mode
       logical :: vm_exits_on_empty
     end function vm_exits_on_empty
     module subroutine windows_init()
     end subroutine windows_init
     module subroutine stack_realloc_maybe()
     end subroutine stack_realloc_maybe
     module function stack_create_window(type,isopen,purpose,isys,irep,idparent,itoken,permanent,orraise)
       integer, intent(in) :: type
       logical, intent(in) :: isopen
       integer, intent(in), optional :: purpose
       integer, intent(in), optional :: isys
       integer, intent(in), optional :: irep
       integer, intent(in), optional :: idparent
       integer, intent(in), optional :: itoken
       logical, intent(in), optional :: permanent
       integer, intent(in), optional :: orraise
       integer :: stack_create_window
     end function stack_create_window
     module subroutine regenerate_window_pointers()
     end subroutine regenerate_window_pointers
     module subroutine window_init(w,type,isopen,id,purpose,isys,irep,idparent,itoken)
       class(window), intent(inout), target :: w
       integer, intent(in) :: type
       logical, intent(in) :: isopen
       integer, intent(in) :: id
       integer, intent(in), optional :: purpose
       integer, intent(in), optional :: isys
       integer, intent(in), optional :: irep
       integer, intent(in), optional :: idparent
       integer, intent(in), optional :: itoken
     end subroutine window_init
     module subroutine window_end(w)
       class(window), intent(inout), target :: w
     end subroutine window_end
     module function window_focused(w)
       class(window), intent(inout) :: w
       logical :: window_focused
     end function window_focused
     module subroutine window_draw(w)
       class(window), intent(inout), target :: w
     end subroutine window_draw
     !xx! tree submodule !xx!
     module subroutine draw_tree(w)
       class(window), intent(inout), target :: w
     end subroutine draw_tree
     module subroutine remap_tree(w)
       class(window), intent(inout) :: w
     end subroutine remap_tree
     module subroutine sort_tree(w)
       class(window), intent(inout) :: w
     end subroutine sort_tree
     module subroutine reassign_tree(w,cfilter)
       class(window), intent(inout) :: w
       type(c_ptr), intent(inout) :: cfilter
     end subroutine reassign_tree
     module subroutine select_system_tree(w,idx)
       class(window), intent(inout) :: w
       integer, intent(in) :: idx
     end subroutine select_system_tree
     !xx! treeplot submodule !xx!
     module subroutine draw_treeplot(w)
       class(window), intent(inout), target :: w
     end subroutine draw_treeplot
     !xx! view submodule !xx!
     module subroutine draw_view(w)
       class(window), intent(inout), target :: w
     end subroutine draw_view
     module subroutine create_texture_view(w,atex)
       class(window), intent(inout), target :: w
       integer, intent(in) :: atex
     end subroutine create_texture_view
     module subroutine delete_texture_view(w)
       class(window), intent(inout), target :: w
     end subroutine delete_texture_view
     module subroutine select_view(w,isys)
       class(window), intent(inout), target :: w
       integer, intent(in) :: isys
     end subroutine select_view
     module subroutine viewmode_set_mode(w)
       class(window), intent(inout), target :: w
     end subroutine viewmode_set_mode
     module subroutine viewmode_set_forced(w,mode,message,idcaller)
       class(window), intent(inout), target :: w
       integer, intent(in) :: mode
       character(len=*), intent(in), optional :: message
       integer, intent(in) :: idcaller
     end subroutine viewmode_set_forced
     module subroutine viewmode_exit_forced(w)
       class(window), intent(inout), target :: w
     end subroutine viewmode_exit_forced
     module subroutine viewmode_release_forced(w,idcaller,mode)
       class(window), intent(inout), target :: w
       integer, intent(in) :: idcaller
       integer, intent(in), optional :: mode
     end subroutine viewmode_release_forced
     module subroutine viewmode_bar_display(w)
       class(window), intent(inout), target :: w
     end subroutine viewmode_bar_display
     module function viewmode_activate_picking(w,hover)
       class(window), intent(inout), target :: w
       logical, intent(in) :: hover
       logical :: viewmode_activate_picking
     end function viewmode_activate_picking
     module subroutine viewmode_process_events(w,hover)
       class(window), intent(inout), target :: w
       logical, intent(in) :: hover
     end subroutine viewmode_process_events
     module subroutine draw_cursor_overlay(w,idx)
       class(window), intent(inout), target :: w
       integer(c_int), intent(in) :: idx(5)
     end subroutine draw_cursor_overlay
     module subroutine mousepos_to_texpos(w,pos)
       class(window), intent(inout), target :: w
       type(ImVec2), intent(inout) :: pos
     end subroutine mousepos_to_texpos
     module subroutine texpos_to_mousepos(w,pos)
       class(window), intent(inout), target :: w
       type(ImVec2), intent(inout) :: pos
     end subroutine texpos_to_mousepos
     module subroutine getpixel(w,pos,depth,rgba)
       class(window), intent(inout), target :: w
       type(ImVec2), intent(inout) :: pos
       real(c_float), intent(out), optional :: depth
       real(c_float), intent(out), optional :: rgba(4)
     end subroutine getpixel
     module subroutine view_to_texpos(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine view_to_texpos
     module subroutine texpos_to_view(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine texpos_to_view
     module subroutine world_to_texpos(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine world_to_texpos
     module subroutine texpos_to_world(w,pos)
       class(window), intent(inout), target :: w
       real(c_float), intent(inout) :: pos(3)
     end subroutine texpos_to_world
     module subroutine export_to_image(w,file,fformat,nsample,npixel,transparentbg,&
        exportview,jpgquality,errmsg)
       class(window), intent(inout), target :: w
       character(kind=c_char,len=*), intent(in) :: file
       character(kind=c_char,len=*), intent(in) :: fformat
       integer(c_int), intent(in) :: nsample, jpgquality, npixel
       logical, intent(in) :: exportview, transparentbg
       character(kind=c_char,len=:), allocatable, intent(inout) :: errmsg
     end subroutine export_to_image
     module subroutine export_to_png_simple(w)
       class(window), intent(inout), target :: w
     end subroutine export_to_png_simple
     !xx! dialog submodule !xx!
     module subroutine draw_dialog(w)
       class(window), intent(inout), target :: w
     end subroutine draw_dialog
     !xx! ci submodule !xx!
     module subroutine draw_ci(w)
       class(window), intent(inout), target :: w
     end subroutine draw_ci
     module subroutine run_commands_ci(w)
       class(window), intent(inout), target :: w
     end subroutine run_commands_ci
     module subroutine block_gui_ci(w,allsys)
       class(window), intent(inout), target :: w
       logical, intent(in) :: allsys
     end subroutine block_gui_ci
     module subroutine get_input_details_ci(w,csystem,cfield)
       class(window), intent(inout), target :: w
       character(len=:), allocatable, intent(inout) :: csystem, cfield
     end subroutine get_input_details_ci
     module subroutine fill_input_ci(str)
       character(len=*), intent(in) :: str
     end subroutine fill_input_ci
     !xx! co submodule !xx!
     module subroutine update_co(w)
       class(window), intent(inout), target :: w
     end subroutine update_co
     module subroutine draw_co(w)
       class(window), intent(inout), target :: w
     end subroutine draw_co
     module function read_output_uout(iscom,cominfo)
       logical, intent(in) :: iscom
       character*(*), intent(in), optional :: cominfo
       logical :: read_output_uout
     end function read_output_uout
     !xx! about submodule !xx!
     module subroutine draw_about(w)
       class(window), intent(inout), target :: w
     end subroutine draw_about
     !xx! preferences submodule !xx!
     module subroutine draw_preferences(w)
       class(window), intent(inout), target :: w
     end subroutine draw_preferences
     !xx! new_struct submodule !xx!
     module subroutine draw_new_struct(w)
       class(window), intent(inout), target :: w
     end subroutine draw_new_struct
     !xx! new_struct_library submodule !xx!
     module subroutine draw_new_struct_library(w)
       class(window), intent(inout), target :: w
     end subroutine draw_new_struct_library
     !xx! load_field submodule !xx!
     module subroutine draw_load_field(w)
       class(window), intent(inout), target :: w
     end subroutine draw_load_field
     !xx! scfplot submodule !xx!
     module subroutine update_scfplot(w)
       class(window), intent(inout), target :: w
     end subroutine update_scfplot
     module subroutine draw_scfplot(w)
       class(window), intent(inout), target :: w
     end subroutine draw_scfplot
     !xx! editrep submodule !xx!
     module subroutine update_editrep(w)
       class(window), intent(inout), target :: w
     end subroutine update_editrep
     module subroutine draw_editrep(w)
       class(window), intent(inout), target :: w
     end subroutine draw_editrep
     module function draw_editrep_atoms(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical :: changed
     end function draw_editrep_atoms
     module function draw_editrep_unitcell(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical :: changed
     end function draw_editrep_unitcell
     module function draw_editrep_axes(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical :: changed
     end function draw_editrep_axes
     module function draw_editrep_symelem(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical :: changed
     end function draw_editrep_symelem
     module function draw_editrep_text(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical :: changed
     end function draw_editrep_text
     module function draw_editrep_measure(w,ttshown) result(changed)
       class(window), intent(inout), target :: w
       logical, intent(inout) :: ttshown
       logical :: changed
     end function draw_editrep_measure
     !xx! exportimage submodule !xx!
     module subroutine draw_exportimage(w)
       class(window), intent(inout), target :: w
     end subroutine draw_exportimage
     !xx! vibrations submodule !xx!
     module subroutine draw_vibrations(w)
       class(window), intent(inout), target :: w
     end subroutine draw_vibrations
     !xx! dynamics submodule !xx!
     module subroutine draw_dynamics(w)
       class(window), intent(inout), target :: w
     end subroutine draw_dynamics
     !xx! water_cluster submodule !xx!
     module subroutine draw_water_cluster(w)
       class(window), intent(inout), target :: w
     end subroutine draw_water_cluster
     !xx! geometry submodule !xx!
     module subroutine update_geometry(w)
       class(window), intent(inout), target :: w
     end subroutine update_geometry
     module subroutine draw_geometry(w)
       class(window), intent(inout), target :: w
     end subroutine draw_geometry
     !xx! builder submodule !xx!
     module subroutine draw_builder(w)
       class(window), intent(inout), target :: w
     end subroutine draw_builder
     module subroutine paste_clipboard_fragment(isys,iview,xpos,icel,atcenter)
       integer, intent(in) :: isys, iview, icel
       real(c_float), intent(in) :: xpos(2)
       logical, intent(in), optional :: atcenter
     end subroutine paste_clipboard_fragment
     module subroutine edit_stop(w)
       class(window), intent(inout) :: w
     end subroutine edit_stop
     module subroutine draw_ff_backend_combo(isys,strid,nchars)
       integer, intent(in) :: isys, nchars
       character(len=*), intent(in) :: strid
     end subroutine draw_ff_backend_combo
     module subroutine addatom_geom_paint(ig,p0,side,dl,iz,bgrgb)
       integer, intent(in) :: ig
       type(ImVec2), intent(in) :: p0
       real(c_float), intent(in) :: side
       type(c_ptr), intent(in) :: dl
       integer, intent(in), optional :: iz
       real(c_float), intent(in), optional :: bgrgb(3)
     end subroutine addatom_geom_paint
  end interface

end module windows
