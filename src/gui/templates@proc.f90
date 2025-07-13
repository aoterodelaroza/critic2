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

! Some utilities for building the GUI (e.g. wrappers around ImGui routines).
submodule (templates) proc
  use iso_c_binding
  use hashmod, only: hash
  use param, only: newline
  implicit none

  ! the keywords
  !- structural tools
  integer, parameter :: ikeyw_none = 0
  integer, parameter :: ikeyw_amd = 1              ! AMD
  integer, parameter :: ikeyw_atomlabel = 2        ! ATOMLABEL
  integer, parameter :: ikeyw_bz = 3               ! BZ
  integer, parameter :: ikeyw_compare = 4          ! COMPARE
  integer, parameter :: ikeyw_comparevc = 5        ! COMPAREVC
  integer, parameter :: ikeyw_coord = 6            ! COORD
  integer, parameter :: ikeyw_econ = 7             ! ECON
  integer, parameter :: ikeyw_environ = 8          ! ENVIRON
  integer, parameter :: ikeyw_ewald = 9            ! EWALD
  integer, parameter :: ikeyw_identify = 10        ! IDENTIFY
  integer, parameter :: ikeyw_kpoints = 11         ! KPOINTS
  integer, parameter :: ikeyw_molmove = 12         ! MOLMOVE
  integer, parameter :: ikeyw_molreorder = 13      ! MOLREORDER
  integer, parameter :: ikeyw_newcell = 14         ! NEWCELL
  integer, parameter :: ikeyw_packing = 15         ! PACKING
  integer, parameter :: ikeyw_polyhedra = 16       ! POLYHEDRA
  integer, parameter :: ikeyw_powder = 17          ! POWDER
  integer, parameter :: ikeyw_rdf = 18             ! RDF
  integer, parameter :: ikeyw_spg = 19             ! SPG
  integer, parameter :: ikeyw_sym = 20             ! SYM/SYMM/NOSYM/NOSYMM
  integer, parameter :: ikeyw_vdw = 21             ! VDW
  !- fields
  integer, parameter :: ikeyw_load = 22            ! LOAD
  integer, parameter :: ikeyw_reference = 23       ! REFERENCE
  integer, parameter :: ikeyw_setfield = 24        ! SETFIELD
  integer, parameter :: ikeyw_unload = 25          ! UNLOAD
  !- read and write files
  integer, parameter :: ikeyw_makemolsnc = 26      ! MAKEMOLSNC
  integer, parameter :: ikeyw_write = 27           ! WRITE
  !- misc
  integer, parameter :: ikeyw_benchmark = 28       ! BENCHMARK
  integer, parameter :: ikeyw_libxc = 29           ! LIBXC
  integer, parameter :: ikeyw_molcell = 30         ! MOLCELL
  integer, parameter :: ikeyw_root = 31            ! ROOT
  integer, parameter :: ikeyw_units = 32           ! UNITS
  integer, parameter :: ikeyw_zpsp = 33            ! ZPSP/Q/QAT/NOCORE
  integer, parameter :: ikeyw_addopts = 34         ! additional options
  !- variables
  integer, parameter :: ikeyw_clear = 35           ! CLEAR
  integer, parameter :: ikeyw_list = 36            ! LIST
  !- field evaluation
  integer, parameter :: ikeyw_point = 37           ! POINT
  integer, parameter :: ikeyw_line = 38            ! LINE
  integer, parameter :: ikeyw_plane = 39           ! PLANE
  integer, parameter :: ikeyw_cube = 40            ! CUBE
  !- critical points
  integer, parameter :: ikeyw_auto = 41            ! AUTO
  integer, parameter :: ikeyw_cpreport = 42        ! CPREPORT
  integer, parameter :: ikeyw_pointprop = 43       ! POINTPROP
  !- qtaim plots
  integer, parameter :: ikeyw_grdvec = 44          ! GRDVEC
  integer, parameter :: ikeyw_fluxprint = 45       ! FLUXPRINT
  integer, parameter :: ikeyw_basinplot = 46       ! BASINPLOT
  integer, parameter :: ikeyw_bundleplot = 47      ! BUNDLEPLOT
  !- integration
  integer, parameter :: ikeyw_integrable = 48      ! INTEGRABLE
  integer, parameter :: ikeyw_integrals = 49       ! INTEGRALS
  integer, parameter :: ikeyw_sphereintegrals = 50 ! SPHEREINTEGRALS
  integer, parameter :: ikeyw_qtree = 51           ! QTREE
  integer, parameter :: ikeyw_yt = 52              ! YT
  integer, parameter :: ikeyw_bader = 53           ! BADER
  integer, parameter :: ikeyw_isosurface = 54      ! ISOSURFACE
  integer, parameter :: ikeyw_hirshfeld = 55       ! HIRSHFELD
  integer, parameter :: ikeyw_voronoi = 56         ! VORONOI
  integer, parameter :: ikeyw_molcalc = 57         ! MOLCALC
  !- nciplot
  integer, parameter :: ikeyw_nciplot = 58         ! NCIPLOT
  !- stm
  integer, parameter :: ikeyw_stm = 59             ! STM
  !- sigma hole
  integer, parameter :: ikeyw_sigmahole = 60       ! SIGMAHOLE
  !- xdm
  integer, parameter :: ikeyw_xdm = 61             ! XDM
  integer, parameter :: ikeyw_NUM = 61

  ! keyword sections (need to be sequential)
  integer, parameter :: isection_none = 0
  integer, parameter :: isection_structural_tools = 1 ! structural tools
  integer, parameter :: isection_fields = 2           ! fields
  integer, parameter :: isection_readwrite = 3        ! read & write files
  integer, parameter :: isection_misc = 4             ! miscellaneous
  integer, parameter :: isection_variables = 5        ! variables
  integer, parameter :: isection_field_evaluation = 6 ! field evaluation
  integer, parameter :: isection_cps = 7              ! critical points
  integer, parameter :: isection_aimplot = 8          ! aim plots
  integer, parameter :: isection_integrals = 9        ! integrals
  integer, parameter :: isection_nciplot = 10         ! nciplot
  integer, parameter :: isection_stm = 11             ! stm
  integer, parameter :: isection_sigmahole = 12       ! sigmahole
  integer, parameter :: isection_xdm = 13             ! xdm
  integer, parameter :: isection_NUM = 13
  integer, parameter :: ikeyw_section(ikeyw_NUM) = (/&
     isection_structural_tools,& ! AMD
     isection_structural_tools,& ! ATOMLABEL
     isection_structural_tools,& ! BZ
     isection_structural_tools,& ! COMPARE
     isection_structural_tools,& ! COMPAREVC
     isection_structural_tools,& ! COORD
     isection_structural_tools,& ! ECON
     isection_structural_tools,& ! ENVIRON
     isection_structural_tools,& ! EWALD
     isection_structural_tools,& ! IDENTIFY
     isection_structural_tools,& ! KPOINTS
     isection_structural_tools,& ! MOLMOVE
     isection_structural_tools,& ! MOLREORDER
     isection_structural_tools,& ! NEWCELL
     isection_structural_tools,& ! PACKING
     isection_structural_tools,& ! POLYHEDRA
     isection_structural_tools,& ! POWDER
     isection_structural_tools,& ! RDF
     isection_structural_tools,& ! SPG
     isection_structural_tools,& ! SYM/SYMM/NOSYM/NOSYMM
     isection_structural_tools,& ! VDW
     isection_fields,&           ! LOAD
     isection_fields,&           ! SETFIELD
     isection_fields,&           ! UNLOAD
     isection_fields,&           ! REFERENCE
     isection_readwrite,&        ! MAKEMOLSNC
     isection_readwrite,&        ! WRITE
     isection_misc,&             ! BENCHMARK
     isection_misc,&             ! LIBXC
     isection_misc,&             ! MOLCELL
     isection_misc,&             ! ROOT
     isection_misc,&             ! UNITS
     isection_misc,&             ! ZPSP/Q/QAT/NOCORE
     isection_misc,&             ! additional options
     isection_variables,&        ! CLEAR
     isection_variables,&        ! LIST
     isection_field_evaluation,& ! POINT
     isection_field_evaluation,& ! LINE
     isection_field_evaluation,& ! PLANE
     isection_field_evaluation,& ! CUBE
     isection_cps,&              ! AUTO
     isection_cps,&              ! CPREPORT
     isection_cps,&              ! POINTPROP
     isection_aimplot,&          ! GRDVEC
     isection_aimplot,&          ! FLUXPRINT
     isection_aimplot,&          ! BASINPLOT
     isection_aimplot,&          ! BUNDLEPLOT
     isection_integrals,&        ! INTEGRABLE
     isection_integrals,&        ! INTEGRALS
     isection_integrals,&        ! SPHEREINTEGRALS
     isection_integrals,&        ! QTREE
     isection_integrals,&        ! YT
     isection_integrals,&        ! BADER
     isection_integrals,&        ! ISOSURFACE
     isection_integrals,&        ! HIRSHFELD
     isection_integrals,&        ! VORONOI
     isection_integrals,&        ! MOLCALC
     isection_nciplot,&          ! NCIPLOT
     isection_stm,&              ! STM
     isection_sigmahole,&        ! SIGMAHOLE
     isection_xdm&               ! XDM
     /)

  ! keyword titles
  character(len=*,kind=c_char), parameter :: keyword_titles(ikeyw_NUM) = (/&
     "AMD (average minimum distances)             ",& ! AMD
     "ATOMLABEL (relabel atoms)                   ",& ! ATOMLABEL
     "BZ (print Brillouin zone)                   ",& ! BZ
     "COMPARE (compare structures)                ",& ! COMPARE
     "COMPAREVC (variable-cell compare structures)",& ! COMPAREVC
     "COORD (2- and 3-coordination numbers)       ",& ! COORD
     "ECON (effective coordination numbers)       ",& ! ECON
     "ENVIRON (calculate atomic environments)     ",& ! ENVIRON
     "EWALD (Ewald potential and energy)          ",& ! EWALD
     "IDENTIFY (identify atomic positions)        ",& ! IDENTIFY
     "KPOINTS (calculate k-point grid sizes)      ",& ! KPOINTS
     "MOLMOVE (move atoms in a molecular crystal) ",& ! MOLMOVE
     "MOLREORDER (reorder atoms in molecules)     ",& ! MOLREORDER
     "NEWCELL (change unit cell)                  ",& ! NEWCELL
     "PACKING (packing ratio)                     ",& ! PACKING
     "POLYHEDRA (coordination polyhedra)          ",& ! POLYHEDRA
     "POWDER (powder diffraction patterns)        ",& ! POWDER
     "RDF (radial distribution functions)         ",& ! RDF
     "SPG (list space group types)                ",& ! SPG
     "SYM (symmetry analysis & refinement)        ",& ! SYM
     "VDW (van der Waals volume)                  ",& ! VDW
     "LOAD (load a field)                         ",& ! LOAD
     "REFERENCE (set the reference field)         ",& ! REFERENCE
     "SETFIELD (set field options)                ",& ! SETFIELD
     "UNLOAD (unload a field)                     ",& ! UNLOAD
     "MAKEMOLSNC (write DMACRYS mols file)        ",& ! MAKEMOLSNC
     "WRITE (write structure to a file)           ",& ! WRITE
     "BENCHMARK (test field evaluation)           ",& ! BENCHMARK
     "LIBXC (list xc functionals)                 ",& ! LIBXC
     "MOLCELL (change molecular cell)             ",& ! MOLCELL
     "ROOT (change new files prefix)              ",& ! ROOT
     "UNITS (change distance units in output)     ",& ! UNITS
     "ZPSP/Q (atomic & pseudo charges)            ",& ! ZPSP/Q/QAT/NOCORE
     "Additional critic2 options                  ",& ! Additional options
     "CLEAR (clear variables)                     ",& ! CLEAR
     "LIST (list variables)                       ",& ! LIST
     "POINT (evaluate field at a point)           ",& ! POINT
     "LINE (evaluate field on a line)             ",& ! LINE
     "PLANE (evaluate field on a plane)           ",& ! PLANE
     "CUBE (evaluate field on a 3D grid)          ",& ! CUBE
     "AUTO (locate critical points)               ",& ! AUTO
     "CPREPORT (report critical points)           ",& ! CPREPORT
     "POINTPROP (properties calculcated at CPs)   ",& ! POINTPROP
     "GRDVEC (plot gradient paths in a plane)     ",& ! GRDVEC
     "FLUXPRINT (plot gradient paths in 3D)       ",& ! FLUXPRINT
     "BASINPLOT (plot atomic basins)              ",& ! BASINPLOT
     "BUNDLEPLOT (plot gradient path bundles)     ",& ! BUNDLEPLOT
     "INTEGRABLE (change integrated properties)   ",& ! INTEGRABLE
     "INTEGRALS (atomic integrals, bisection)     ",& ! INTEGRALS
     "SPHEREINTEGRALS (integrate atomic spheres)  ",& ! SPHEREINTEGRALS
     "QTREE (atomic integrals, quad-tree)         ",& ! QTREE
     "YT (atomic integrals, Yu-Trinkle, grids)    ",& ! YT
     "BADER (atomic integrals, Henkelman, grids)  ",& ! BADER
     "ISOSURFACE (isosurface integrals, grids)    ",& ! ISOSURFACE
     "HIRSHFELD (Hirshfeld atomic integrals)      ",& ! HIRSHFELD
     "VORONOI (Voronoi atomic integrals, grids)   ",& ! VORONOI
     "MOLCALC (mesh integration)                  ",& ! MOLCALC
     "NCIPLOT (non-covalent interaction plots)    ",& ! NCIPLOT
     "STM (scanning tunneling microscopy plots)   ",& ! STM
     "SIGMAHOLE (characterize sigma holes)        ",& ! SIGMAHOLE
     "XDM (exchange-hole dipole moment dispersion)"&  ! XDM
     /)

  ! section titles
  character(len=*,kind=c_char), parameter :: section_titles(isection_NUM) = (/&
     "Structural Tools  ",& ! structural tools
     "Load/Unload Fields",& ! fields
     "Read/Write Files  ",& ! read/write files
     "Miscellaneous     ",& ! miscellaneous
     "Variables         ",& ! variables
     "Field Evaluation  ",& ! field evaluation
     "Critical Points   ",& ! critical points
     "QTAIM Plots       ",& ! aim plots
     "Integrals         ",& ! integrals
     "NCI plots         ",& ! nciplots
     "STM plots         ",& ! stm plots
     "Sigma holes       ",& ! sigma hole
     "XDM dispersion    "&  ! xdm
     /)

  ! section ranges
  integer, parameter :: section_ranges(2,isection_NUM) = reshape((/&
     1, 21,& ! structural tools
     22,25,& ! fields
     26,27,& ! read/write files
     28,34,& ! miscellaneous
     35,36,& ! variables
     37,40,& ! field evaluation
     41,43,& ! critical points
     44,47,& ! qtaim plots
     48,57,& ! integration
     58,58,& ! nciplot
     59,59,& ! stm
     60,60,& ! sigmahole
     61,61&  ! xdm
     /),shape(section_ranges))

  ! template (keyw) files
  character*(*), parameter :: template_file(ikeyw_NUM) = (/&
     "amd            ",& ! AMD
     "atomlabel      ",& ! ATOMLABEL
     "bz             ",& ! BZ
     "compare        ",& ! COMPARE
     "comparevc      ",& ! COMPAREVC
     "coord          ",& ! COORD
     "econ           ",& ! ECON
     "environ        ",& ! ENVIRON
     "ewald          ",& ! EWALD
     "identify       ",& ! IDENTIFY
     "kpoints        ",& ! KPOINTS
     "molmove        ",& ! MOLMOVE
     "molreorder     ",& ! MOLREORDER
     "newcell        ",& ! NEWCELL
     "packing        ",& ! PACKING
     "polyhedra      ",& ! POLYHEDRA
     "powder         ",& ! POWDER
     "rdf            ",& ! RDF
     "spg            ",& ! SPG
     "sym            ",& ! SYM
     "vdw            ",& ! VDW
     "load           ",& ! LOAD
     "reference      ",& ! REFERENCE
     "setfield       ",& ! SETFIELD
     "unload         ",& ! UNLOAD
     "makemolsnc     ",& ! MAKEMOLSNC
     "write          ",& ! WRITE
     "benchmark      ",& ! BENCHMARK
     "libxc          ",& ! LIBXC
     "molcell        ",& ! MOLCELL
     "root           ",& ! ROOT
     "units          ",& ! UNITS
     "zpsp           ",& ! ZPSP/Q/QAT/NOCORE
     "addopts        ",& ! additional options
     "clear          ",& ! CLEAR
     "list           ",& ! LIST
     "point          ",& ! POINT
     "line           ",& ! LINE
     "plane          ",& ! PLANE
     "cube           ",& ! CUBE
     "auto           ",& ! AUTO
     "cpreport       ",& ! CPREPORT
     "pointprop      ",& ! POINTPROP
     "grdvec         ",& ! GRDVEC
     "fluxprint      ",& ! FLUXPRINT
     "basinplot      ",& ! BASINPLOT
     "bundleplot     ",& ! BUNDLEPLOT
     "integrable     ",& ! INTEGRABLE
     "integrals      ",& ! INTEGRALS
     "sphereintegrals",& ! SPHEREINTEGRALS
     "qtree          ",& ! QTREE
     "yt             ",& ! YT
     "bader          ",& ! BADER
     "isosurface     ",& ! ISOSURFACE
     "hirshfeld      ",& ! HIRSHFELD
     "voronoi        ",& ! VORONOI
     "molcalc        ",& ! MOLCALC
     "nciplot        ",& ! NCIPLOT
     "stm            ",& ! STM
     "sigmahole      ",& ! SIGMAHOLE
     "xdm            "&  ! XDM
     /)

  ! documentation (md) files
  character*(*), parameter :: doclink(ikeyw_NUM) = (/&
     "structure/#c2-amd         ",& ! AMD
     "structure/#c2-atomlabel   ",& ! ATOMLABEL
     "structure/#c2-bz          ",& ! BZ
     "structure/#c2-compare     ",& ! COMPARE
     "structure/#c2-comparevc   ",& ! COMPAREVC
     "structure/#c2-coord       ",& ! COORD
     "structure/#c2-econ        ",& ! ECON
     "structure/#c2-environ     ",& ! ENVIRON
     "structure/#c2-ewald       ",& ! EWALD
     "structure/#c2-identify    ",& ! IDENTIFY
     "structure/#c2-kpoints     ",& ! KPOINTS
     "structure/#c2-molmove     ",& ! MOLMOVE
     "structure/#c2-molreorder  ",& ! MOLREORDER
     "structure/#c2-newcell     ",& ! NEWCELL
     "structure/#c2-packing     ",& ! PACKING
     "structure/#c2-polyhedra   ",& ! POLYHEDRA
     "structure/#c2-powder      ",& ! POWDER
     "structure/#c2-rdf         ",& ! RDF
     "crystal/#c2-spg           ",& ! SPG
     "crystal/#c2-symm          ",& ! SYM
     "structure/#c2-vdw         ",& ! VDW
     "fields/#c2-load           ",& ! LOAD
     "fields/#c2-reference      ",& ! REFERENCE
     "fields/#c2-setfield       ",& ! SETFIELD
     "fields/#c2-unload         ",& ! UNLOAD
     "write/#c2-makemolsnc      ",& ! MAKEMOLSNC
     "write/#c2-write           ",& ! WRITE
     "misc/#c2-benchmark        ",& ! BENCHMARK
     "arithmetics/#libxc        ",& ! LIBXC
     "molecule/#c2-molcell      ",& ! MOLCELL
     "misc/#c2-root             ",& ! ROOT
     "inputoutput/#c2-units     ",& ! UNITS
     "crystal/#c2-charge        ",& ! ZPSP/Q/QAT/NOCORE
     "misc/#c2-control          ",& ! additional options
     "arithmetics/#c2-clear     ",& ! CLEAR
     "arithmetics/#c2-list      ",& ! LIST
     "graphics/#c2-point        ",& ! POINT
     "graphics/#c2-line         ",& ! LINE
     "graphics/#c2-plane        ",& ! PLANE
     "graphics/#c2-cube         ",& ! CUBE
     "cpsearch/#c2-auto         ",& ! AUTO
     "cpsearch/#c2-cpreport     ",& ! CPREPORT
     "cpsearch/#c2-pointprop    ",& ! POINTPROP
     "gradientpath/#c2-grdvec   ",& ! GRDVEC
     "gradientpath/#c2-fluxprint",& ! FLUXPRINT
     "basinplot/#c2-basinplot   ",& ! BASINPLOT
     "basinplot/#c2-bundleplot  ",& ! BUNDLEPLOT
     "integrate/#c2-integrable  ",& ! INTEGRABLE
     "integrate/#c2-integrals   ",& ! INTEGRALS
     "integrate/#c2-integrals   ",& ! SPHEREINTEGRALS
     "integrate/#c2-qtree       ",& ! QTREE
     "integrate/#c2-yt          ",& ! YT
     "integrate/#c2-bader       ",& ! BADER
     "integrate/#c2-isosurface  ",& ! ISOSURFACE
     "integrate/#c2-hirshfeld   ",& ! HIRSHFELD
     "integrate/#c2-voronoi     ",& ! VORONOI
     "misc/#c2-molcalc          ",& ! MOLCALC
     "nciplot/                  ",& ! NCIPLOT
     "stm/                      ",& ! STM
     "misc/#c2-sigmahole        ",& ! SIGMAHOLE
     "misc/#c2-xdm              "&  ! XDM
     /)

  ! template hash
  type(hash) :: thash

  !xx! private procedures
  ! subroutine launch_keyword_action(textinsert,ikeyw)
  ! function get_keyword_template(ikeyw)

contains

  !> Draw the keyword context menu in the templates and help buttons
  !> of the consolie input window. If textinsert is true, insert the
  !> tempalte, otherwise bring up the help.
  module subroutine draw_keyword_context_menu(textinsert)
    use interfaces_cimgui
    use utils, only: iw_menuitem
    logical, intent(in) :: textinsert

    character(kind=c_char,len=:), allocatable, target :: str1, str2

    integer :: i, j

    do i = 1, isection_NUM
       str1 = trim(section_titles(i)) // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          do j = section_ranges(1,i), section_ranges(2,i)
             ! str2 = trim(keyword_titles(j)) // c_null_char
             if (iw_menuitem(trim(keyword_titles(j)))) &
                call launch_keyword_action(textinsert,j)
          end do

          call igEndMenu()
       end if
    end do
    call igEndPopup()

  end subroutine draw_keyword_context_menu

  !xx! private procedures

  !> Launch the keyword action for the given keyword ikeyw. If
  !> textinsert, insert the template into the console
  !> input. Otherwise, open the documentation.
  subroutine launch_keyword_action(textinsert,ikeyw)
    use windows, only: win, iwin_console_input
    use interfaces_cimgui
    logical, intent(in) :: textinsert
    integer, intent(in) :: ikeyw

    character(kind=c_char,len=:), allocatable, target :: str

    if (textinsert) then
       str = get_keyword_template(ikeyw)
       call win(iwin_console_input)%fill_input_ci(str)
    else
       str = "https://aoterodelaroza.github.io/critic2/manual/" // trim(doclink(ikeyw)) //&
          c_null_char
       call openLink(c_loc(str))
    end if

  end subroutine launch_keyword_action

  !> Returns the template for the given keyword.
  function get_keyword_template(ikeyw)
    use global, only: critic_home
    use tools_io, only: fopen_read, fclose, getline_raw
    use param, only: dirsep, newline
    integer, intent(in) :: ikeyw
    character(kind=c_char,len=:), allocatable :: get_keyword_template

    character(kind=c_char,len=:), allocatable :: file, keyw, line
    logical :: exist
    integer :: lu

    keyw = trim(template_file(ikeyw))
    if (.not.thash%iskey(keyw)) then
       ! check the file exists
       file = trim(critic_home) // dirsep // "helpdoc" // dirsep // keyw // ".keyw"
       inquire(file=file,exist=exist)
       if (.not.exist) then
          get_keyword_template = &
             "## Template file not found. Please make sure the CRITIC_HOME is set or" // newline //&
             "## that critic2 is installed propertly."
          return
       end if

       ! read the file
       get_keyword_template = ""
       lu = fopen_read(file)
       do while (getline_raw(lu,line,.false.))
          get_keyword_template = get_keyword_template // line // newline
       end do
       call fclose(lu)

       ! save it to the hash
       call thash%put(keyw,get_keyword_template)
    else
       ! get it from the hash
       get_keyword_template = thash%get(keyw,'a')
    end if

  end function get_keyword_template

end submodule proc
