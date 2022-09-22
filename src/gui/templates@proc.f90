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
  integer, parameter :: ikeyw_coord = 5            ! COORD
  integer, parameter :: ikeyw_econ = 6             ! ECON
  integer, parameter :: ikeyw_environ = 7          ! ENVIRON
  integer, parameter :: ikeyw_ewald = 8            ! EWALD
  integer, parameter :: ikeyw_identify = 9         ! IDENTIFY
  integer, parameter :: ikeyw_kpoints = 10          ! KPOINTS
  integer, parameter :: ikeyw_molmove = 11         ! MOLMOVE
  integer, parameter :: ikeyw_molreorder = 12      ! MOLREORDER
  integer, parameter :: ikeyw_newcell = 13         ! NEWCELL
  integer, parameter :: ikeyw_packing = 14         ! PACKING
  integer, parameter :: ikeyw_polyhedra = 15       ! POLYHEDRA
  integer, parameter :: ikeyw_powder = 16          ! POWDER
  integer, parameter :: ikeyw_rdf = 17             ! RDF
  integer, parameter :: ikeyw_spg = 18             ! SPG
  integer, parameter :: ikeyw_sym = 19             ! SYM/SYMM/NOSYM/NOSYMM
  integer, parameter :: ikeyw_vdw = 20             ! VDW
  !- fields
  integer, parameter :: ikeyw_load = 21            ! LOAD
  integer, parameter :: ikeyw_reference = 22       ! REFERENCE
  integer, parameter :: ikeyw_setfield = 23        ! SETFIELD
  integer, parameter :: ikeyw_unload = 24          ! UNLOAD
  !- read and write files
  integer, parameter :: ikeyw_makemolsnc = 25      ! MAKEMOLSNC
  integer, parameter :: ikeyw_write = 26           ! WRITE
  !- misc
  integer, parameter :: ikeyw_benchmark = 27       ! BENCHMARK
  integer, parameter :: ikeyw_libxc = 28           ! LIBXC
  integer, parameter :: ikeyw_molcell = 29         ! MOLCELL
  integer, parameter :: ikeyw_root = 30            ! ROOT
  integer, parameter :: ikeyw_units = 31           ! UNITS
  integer, parameter :: ikeyw_zpsp = 32            ! ZPSP/Q/QAT/NOCORE
  integer, parameter :: ikeyw_addopts = 33         ! additional options
  !- variables
  integer, parameter :: ikeyw_clear = 34           ! CLEAR
  integer, parameter :: ikeyw_list = 35            ! LIST
  !- field evaluation
  integer, parameter :: ikeyw_point = 36           ! POINT
  integer, parameter :: ikeyw_line = 37            ! LINE
  integer, parameter :: ikeyw_plane = 38           ! PLANE
  integer, parameter :: ikeyw_cube = 39            ! CUBE
  !- critical points
  integer, parameter :: ikeyw_auto = 40            ! AUTO
  integer, parameter :: ikeyw_cpreport = 41        ! CPREPORT
  integer, parameter :: ikeyw_pointprop = 42       ! POINTPROP
  !- qtaim plots
  integer, parameter :: ikeyw_grdvec = 43          ! GRDVEC
  integer, parameter :: ikeyw_fluxprint = 44       ! FLUXPRINT
  integer, parameter :: ikeyw_basinplot = 45       ! BASINPLOT
  integer, parameter :: ikeyw_bundleplot = 46      ! BUNDLEPLOT
  !- integration
  integer, parameter :: ikeyw_integrable = 47      ! INTEGRABLE
  integer, parameter :: ikeyw_integrals = 48       ! INTEGRALS
  integer, parameter :: ikeyw_sphereintegrals = 49 ! SPHEREINTEGRALS
  integer, parameter :: ikeyw_qtree = 50           ! QTREE
  integer, parameter :: ikeyw_yt = 51              ! YT
  integer, parameter :: ikeyw_bader = 52           ! BADER
  integer, parameter :: ikeyw_isosurface = 53      ! ISOSURFACE
  integer, parameter :: ikeyw_hirshfeld = 54       ! HIRSHFELD
  integer, parameter :: ikeyw_voronoi = 55         ! VORONOI
  integer, parameter :: ikeyw_molcalc = 56         ! MOLCALC
  !- nciplot
  integer, parameter :: ikeyw_nciplot = 57         ! NCIPLOT
  !- stm
  integer, parameter :: ikeyw_stm = 58             ! STM
  !- sigma hole
  integer, parameter :: ikeyw_sigmahole = 59       ! SIGMAHOLE
  !- xdm
  integer, parameter :: ikeyw_xdm = 60             ! XDM
  integer, parameter :: ikeyw_NUM = 60

  ! keyword sections (need to be sequential)
  integer, parameter :: isection_none = 0
  integer, parameter :: isection_structural_tools = 1 ! structural tools
  integer, parameter :: isection_readwrite = 2        ! read & write files
  integer, parameter :: isection_fields = 3           ! fields
  integer, parameter :: isection_field_evaluation = 4 ! field evaluation
  integer, parameter :: isection_cps = 5              ! critical points
  integer, parameter :: isection_aimplot = 6          ! aim plots
  integer, parameter :: isection_integrals = 7        ! integrals
  integer, parameter :: isection_stm = 8              ! stm
  integer, parameter :: isection_nciplot = 9          ! nciplot
  integer, parameter :: isection_xdm = 10             ! xdm
  integer, parameter :: isection_sigmahole = 11       ! sigmahole
  integer, parameter :: isection_variables = 12       ! variables
  integer, parameter :: isection_misc = 13            ! miscellaneous
  integer, parameter :: isection_NUM = 13
  integer, parameter :: ikeyw_section(ikeyw_NUM) = (/&
     isection_structural_tools,& ! AMD
     isection_structural_tools,& ! ATOMLABEL
     isection_structural_tools,& ! BZ
     isection_structural_tools,& ! COMPARE
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
     "Read/Write Files  ",& ! read/write files
     "Load/Unload Fields",& ! fields
     "Field Evaluation  ",& ! field evaluation
     "Critical Points   ",& ! critical points
     "QTAIM Plots       ",& ! aim plots
     "Integrals         ",& ! integrals
     "STM plots         ",& ! stm plots
     "NCI plots         ",& ! nciplots
     "XDM dispersion    ",& ! xdm
     "Sigma holes       ",& ! sigma hole
     "Variables         ",& ! variables
     "Miscellaneous     "& ! miscellaneous
     /)

  ! section ranges
  integer, parameter :: section_ranges(2,isection_NUM) = reshape((/&
     1,19,&  ! structural tools
     24,25,& ! read/write files
     20,23,& ! fields
     35,38,& ! field evaluation
     39,41,& ! critical points
     42,45,& ! aim plots
     46,55,& ! integrals
     57,57,& ! stm plots
     56,56,& ! nciplots
     59,59,& ! xdm
     58,58,& ! sigma hole
     33,34,& ! variables
     26,32& ! miscellaneous
     /),shape(section_ranges))

  ! template (keyw) files
  character*(*), parameter :: template_file(ikeyw_NUM) = (/&
     "amd            ",& ! AMD
     "atomlabel      ",& ! ATOMLABEL
     "bz             ",& ! BZ
     "compare        ",& ! COMPARE
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
     "makemolsnc     ",& ! WRITE
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
     "fluxprint      ",& ! GRDVEC
     "basinplot      ",& ! GRDVEC
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

contains

  !> Draw the keyword context menu in the templates and help buttons
  !> of the consolie input window. If textinsert is true, insert the
  !> tempalte, otherwise bring up the help.
  module subroutine draw_keyword_context_menu(textinsert)
    use interfaces_cimgui
    logical, intent(in) :: textinsert

    character(kind=c_char,len=:), allocatable, target :: str1, str2

    integer :: i, j

    do i = 1, isection_NUM
       str1 = trim(section_titles(i)) // c_null_char
       if (igBeginMenu(c_loc(str1),.true._c_bool)) then
          do j = section_ranges(1,i), section_ranges(2,i)
             str2 = trim(keyword_titles(j)) // c_null_char
             if (igMenuItem_Bool(c_loc(str2),c_null_ptr,.false._c_bool,.true._c_bool)) then
                call launch_keyword_action(textinsert,j)
             end if
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
    use windows, only: win, iwin_console_input, stack_create_window
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
