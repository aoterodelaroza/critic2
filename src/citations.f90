! Copyright (c) 2007-2026 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Registry of literature citations. Every "Please cite" block the
! program writes comes from this table, in one format, and the GUI
! help window reads the same entries (by key) to show them with their
! DOI links. Add a new paper here (index parameter + table entry) and
! refer to it by index from the code or by key from the help files.
module citations
  implicit none

  private

  public :: cite_print
  public :: cite_find
  public :: cite_line
  public :: cite_doi_url
  public :: cite_bibtex

  !> A literature citation: a short key (used by the GUI help files),
  !> the BibTeX entry kind (article or misc), the fields, and the bare
  !> DOI (blank if the work has none). The author list is in BibTeX
  !> form (A and B and C). All fields are plain ASCII: they go to the
  !> CLI output and the start-up banner, which Windows consoles do not
  !> read as UTF-8.
  type citation
     character(len=24) :: key
     character(len=8) :: kind
     character(len=256) :: author
     character(len=256) :: title
     character(len=64) :: journal ! for misc: how it was published (name and URL)
     character(len=8) :: volume
     character(len=4) :: year
     character(len=16) :: pages
     character(len=80) :: doi
  end type citation
  public :: citation

  ! citation indices
  integer, parameter, public :: cite_critic2_2014 = 1
  integer, parameter, public :: cite_critic2_2009 = 2
  integer, parameter, public :: cite_xdm_2007 = 3
  integer, parameter, public :: cite_xdm_2012 = 4
  integer, parameter, public :: cite_xdm_2013 = 5
  integer, parameter, public :: cite_nci_2010 = 6
  integer, parameter, public :: cite_nci_2012 = 7
  integer, parameter, public :: cite_bader_2006 = 8
  integer, parameter, public :: cite_bader_2007 = 9
  integer, parameter, public :: cite_bader_2009 = 10
  integer, parameter, public :: cite_yt_2011 = 11
  integer, parameter, public :: cite_hirshfeld_i_2007 = 12
  integer, parameter, public :: cite_becke_1988 = 13
  integer, parameter, public :: cite_lebedev_1999 = 14
  integer, parameter, public :: cite_horton_2017 = 15
  integer, parameter, public :: cite_smoothrho_2022 = 16
  integer, parameter, public :: cite_peach_2008 = 17
  integer, parameter, public :: cite_sigmahole_2016 = 18
  integer, parameter, public :: cite_keast_1986 = 19
  integer, parameter, public :: cite_cubpack_2003 = 20
  integer, parameter, public :: cite_rodriguez_2009 = 21
  integer, parameter, public :: cite_blatov_2004 = 22
  integer, parameter, public :: cite_blatov_2014 = 23
  integer, parameter, public :: cite_amd_2022 = 24
  integer, parameter, public :: cite_gpwdf_2024 = 25
  integer, parameter, public :: cite_degelder_2001 = 26
  integer, parameter, public :: cite_emd_2000 = 27
  integer, parameter, public :: cite_vcpwdf_2022 = 28
  integer, parameter, public :: cite_econ_2016 = 29
  integer, parameter, public :: cite_econ_1979 = 30
  integer, parameter, public :: cite_hinuma_2017 = 31
  integer, parameter, public :: cite_xrpd_bg_2001 = 32
  integer, parameter :: ncite = 32

  !> The citations of the program itself (start-up banner, About window)
  integer, parameter, public :: cite_program(2) = (/cite_critic2_2014,cite_critic2_2009/)

  !> The citation table, indexed by the cite_* parameters.
  type(citation), parameter, public :: citelist(ncite) = (/&
     citation("critic2_2014","article",&
        "A. Otero-de-la-Roza and E. R. Johnson and V. Luana",&
        "Critic2: A program for real-space analysis of quantum chemical interactions in solids",&
        "Comput. Phys. Commun.","185","2014","1007-1018",&
        "10.1016/j.cpc.2013.10.026"),&
     citation("critic2_2009","article",&
        "A. Otero-de-la-Roza and M. A. Blanco and A. Martin Pendas and V. Luana",&
        "Critic: a new program for the topological analysis of solid-state electron densities",&
        "Comput. Phys. Commun.","180","2009","157-166",&
        "10.1016/j.cpc.2008.07.018"),&
     citation("xdm_2007","article",&
        "A. D. Becke and E. R. Johnson",&
        "Exchange-hole dipole moment and the dispersion interaction revisited",&
        "J. Chem. Phys.","127","2007","154108",&
        "10.1063/1.2795701"),&
     citation("xdm_2012","article",&
        "A. Otero-de-la-Roza and E. R. Johnson",&
        "Van der Waals interactions in solids using the exchange-hole dipole moment model",&
        "J. Chem. Phys.","136","2012","174109",&
        "10.1063/1.4705760"),&
     citation("xdm_2013","article",&
        "A. Otero-de-la-Roza and E. R. Johnson",&
        "Non-covalent interactions and thermochemistry using XDM-corrected hybrid and range-separated hybrid density functionals",&
        "J. Chem. Phys.","138","2013","204109",&
        "10.1063/1.4807330"),&
     citation("nci_2010","article",&
        "E. R. Johnson and S. Keinan and P. Mori-Sanchez and J. Contreras-Garcia and A. J. Cohen and W. Yang",&
        "Revealing Noncovalent Interactions",&
        "J. Am. Chem. Soc.","132","2010","6498-6506",&
        "10.1021/ja100936w"),&
     citation("nci_2012","article",&
        "A. Otero-de-la-Roza and E. R. Johnson and J. Contreras-Garcia",&
        "Revealing non-covalent interactions in solids: NCI plots revisited",&
        "Phys. Chem. Chem. Phys.","14","2012","12165-12172",&
        "10.1039/C2CP41395G"),&
     citation("bader_2006","article",&
        "G. Henkelman and A. Arnaldsson and H. Jonsson",&
        "A fast and robust algorithm for Bader decomposition of charge density",&
        "Comput. Mater. Sci.","36","2006","354-360",&
        "10.1016/j.commatsci.2005.04.010"),&
     citation("bader_2007","article",&
        "E. Sanville and S. D. Kenny and R. Smith and G. Henkelman",&
        "Improved grid-based algorithm for Bader charge allocation",&
        "J. Comput. Chem.","28","2007","899-908",&
        "10.1002/jcc.20575"),&
     citation("bader_2009","article",&
        "W. Tang and E. Sanville and G. Henkelman",&
        "A grid-based Bader analysis algorithm without lattice bias",&
        "J. Phys.: Condens. Matter","21","2009","084204",&
        "10.1088/0953-8984/21/8/084204"),&
     citation("yt_2011","article",&
        "M. Yu and D. R. Trinkle",&
        "Accurate and efficient algorithm for Bader charge integration",&
        "J. Chem. Phys.","134","2011","064111",&
        "10.1063/1.3553716"),&
     citation("hirshfeld_i_2007","article",&
        "P. Bultinck and C. Van Alsenoy and P. W. Ayers and R. Carbo-Dorca",&
        "Critical analysis and extension of the Hirshfeld atoms in molecules",&
        "J. Chem. Phys.","126","2007","144111",&
        "10.1063/1.2715563"),&
     citation("becke_1988","article",&
        "A. D. Becke",&
        "A multicenter numerical integration scheme for polyatomic molecules",&
        "J. Chem. Phys.","88","1988","2547-2553",&
        "10.1063/1.454033"),&
     citation("lebedev_1999","article",&
        "V. I. Lebedev and D. N. Laikov",&
        "A quadrature formula for the sphere of the 131st algebraic order of accuracy",&
        "Dokl. Math.","59","1999","477-481",&
        ""),&
     citation("horton_2017","misc",&
        "T. Verstraelen and P. Tecmer and F. Heidar-Zadeh and C. E. Gonzalez-Espinoza and M. Chan and T. D. Kim and K. Boguslawski and S. Fias and S. Vandenbrande and D. Berrocal and P. W. Ayers",&
        "HORTON 2.1.0",&
        "HORTON 2.1.0, http://theochem.github.com/horton/","","2017","",&
        ""),&
     citation("smoothrho_2022","article",&
        "A. Otero-de-la-Roza",&
        "Finding critical points and reconstruction of electron densities on grids",&
        "J. Chem. Phys.","156","2022","224116",&
        "10.1063/5.0090232"),&
     citation("peach_2008","article",&
        "M. J. G. Peach and P. Benfield and T. Helgaker and D. J. Tozer",&
        "Excitation energies in density functional theory: An evaluation and a diagnostic test",&
        "J. Chem. Phys.","128","2008","044118",&
        "10.1063/1.2831900"),&
     citation("sigmahole_2016","article",&
        "M. H. Kolar and P. Hobza",&
        "Computer Modeling of Halogen Bonds and Other sigma-Hole Interactions",&
        "Chem. Rev.","116","2016","5155-5187",&
        "10.1021/acs.chemrev.5b00560"),&
     citation("keast_1986","article",&
        "P. Keast",&
        "Moderate-degree tetrahedral quadrature formulas",&
        "Comput. Methods Appl. Mech. Eng.","55","1986","339-348",&
        "10.1016/0045-7825(86)90059-9"),&
     citation("cubpack_2003","article",&
        "R. Cools and A. Haegemans",&
        "Algorithm 824: CUBPACK: a package for automatic cubature; framework description",&
        "ACM Trans. Math. Softw.","29","2003","287-296",&
        "10.1145/838250.838253"),&
     citation("rodriguez_2009","article",&
        "J. I. Rodriguez and A. M. Koster and P. W. Ayers and A. Santos-Valle and A. Vela and G. Merino",&
        "An efficient grid-based scheme to compute QTAIM atomic properties without explicit calculation of zero-flux surfaces",&
        "J. Comput. Chem.","30","2009","1082-1092",&
        "10.1002/jcc.21134"),&
     citation("blatov_2004","article",&
        "V. A. Blatov",&
        "Voronoi-dirichlet polyhedra in crystal chemistry: theory and applications",&
        "Crystallogr. Rev.","10","2004","249-318",&
        "10.1080/08893110412331323170"),&
     citation("blatov_2014","article",&
        "V. A. Blatov and A. P. Shevchenko and D. M. Proserpio",&
        "Applied Topological Analysis of Crystal Structures with the Program Package ToposPro",&
        "Cryst. Growth Des.","14","2014","3576-3586",&
        "10.1021/cg500498k"),&
     citation("amd_2022","article",&
        "D. Widdowson and M. M. Mosca and A. Pulido and A. I. Cooper and V. Kurlin",&
        "Average minimum distances of periodic point sets - foundational invariants for mapping periodic crystals",&
        "MATCH Commun. Math. Comput. Chem.","87","2022","529-559",&
        "10.46793/match.87-3.529W"),&
     citation("gpwdf_2024","article",&
        "A. Otero-de-la-Roza",&
        "Powder-diffraction-based structural comparison for crystal structure prediction without prior indexing",&
        "J. Appl. Cryst.","57","2024","1401-1414",&
        "10.1107/S1600576724007489"),&
     citation("degelder_2001","article",&
        "R. de Gelder and R. Wehrens and J. A. Hageman",&
        "A generalized expression for the similarity of spectra: application to powder diffraction pattern classification",&
        "J. Comput. Chem.","22","2001","273-289",&
        "10.1002/1096-987X(200102)22:3<273::AID-JCC1001>3.0.CO;2-0"),&
     citation("emd_2000","article",&
        "Y. Rubner and C. Tomasi and L. J. Guibas",&
        "The Earth Mover's Distance as a Metric for Image Retrieval",&
        "Int. J. Comput. Vis.","40","2000","99-121",&
        "10.1023/A:1026543900054"),&
     citation("vcpwdf_2022","article",&
        "R. A. Mayo and A. Otero-de-la-Roza and E. R. Johnson",&
        "Development and assessment of an improved powder-diffraction-based method for molecular crystal structure similarity",&
        "CrystEngComm","24","2022","8326-8338",&
        "10.1039/D2CE01080A"),&
     citation("econ_2016","article",&
        "M. Nespolo",&
        "Charge distribution as a tool to investigate structural details. IV. A new route to heteroligand polyhedra",&
        "Acta Cryst. B","72","2016","51-66",&
        "10.1107/S2052520615019472"),&
     citation("econ_1979","article",&
        "R. Hoppe",&
        "Effective coordination numbers (ECoN) and mean fictive ionic radii (MEFIR)",&
        "Z. Kristallogr.","150","1979","23-52",&
        "10.1524/zkri.1979.150.14.23"),&
     citation("hinuma_2017","article",&
        "Y. Hinuma and G. Pizzi and Y. Kumagai and F. Oba and I. Tanaka",&
        "Band structure diagram paths based on crystallography",&
        "Comput. Mater. Sci.","128","2017","140-184",&
        "10.1016/j.commatsci.2016.10.015"),&
     citation("xrpd_bg_2001","article",&
        "W. I. F. David and D. S. Sivia",&
        "Background estimation using a robust Bayesian analysis",&
        "J. Appl. Cryst.","34","2001","318-324",&
        "10.1107/S0021889801004332")&
     /)

  interface
     module subroutine cite_print(ids,header,indent,comment)
       integer, intent(in) :: ids(:)
       character*(*), intent(in), optional :: header
       integer, intent(in), optional :: indent
       logical, intent(in), optional :: comment
     end subroutine cite_print
     module function cite_find(key) result(id)
       character*(*), intent(in) :: key
       integer :: id
     end function cite_find
     module function cite_line(id) result(str)
       integer, intent(in) :: id
       character(len=:), allocatable :: str
     end function cite_line
     module function cite_doi_url(id) result(str)
       integer, intent(in) :: id
       character(len=:), allocatable :: str
     end function cite_doi_url
     module function cite_bibtex(id) result(str)
       integer, intent(in) :: id
       character(len=:), allocatable :: str
     end function cite_bibtex
  end interface

end module citations
