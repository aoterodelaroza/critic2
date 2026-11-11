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

  !> A literature citation: a short key (used by the GUI help files),
  !> the reference text (authors, journal, volume (year) pages) and the
  !> bare DOI (blank if the work has none). The text is plain ASCII:
  !> it goes to the CLI output and the start-up banner, which Windows
  !> consoles do not read as UTF-8.
  type citation
     character(len=24) :: key
     character(len=256) :: text
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
     citation("critic2_2014","A. Otero-de-la-Roza, E. R. Johnson, and V. Luana, Comput. Phys. Commun. 185 (2014) 1007-1018.",&
     "10.1016/j.cpc.2013.10.026"),&
     citation("critic2_2009","A. Otero-de-la-Roza, M. A. Blanco, A. Martin Pendas, and V. Luana, Comput. Phys. Commun. 180 (2009) 157-166.",&
     "10.1016/j.cpc.2008.07.018"),&
     citation("xdm_2007","A. D. Becke and E. R. Johnson, J. Chem. Phys. 127 (2007) 154108.",&
     "10.1063/1.2795701"),&
     citation("xdm_2012","A. Otero-de-la-Roza and E. R. Johnson, J. Chem. Phys. 136 (2012) 174109.",&
     "10.1063/1.4705760"),&
     citation("xdm_2013","A. Otero-de-la-Roza and E. R. Johnson, J. Chem. Phys. 138 (2013) 204109.",&
     "10.1063/1.4807330"),&
     citation("nci_2010","E. R. Johnson, S. Keinan, P. Mori-Sanchez, J. Contreras-Garcia, A. J. Cohen, and W. Yang, J. Am. Chem. Soc. 132 (2010) 6498-6506.",&
     "10.1021/ja100936w"),&
     citation("nci_2012","A. Otero-de-la-Roza, E. R. Johnson, and J. Contreras-Garcia, Phys. Chem. Chem. Phys. 14 (2012) 12165-12172.",&
     "10.1039/C2CP41395G"),&
     citation("bader_2006","G. Henkelman, A. Arnaldsson, and H. Jonsson, Comput. Mater. Sci. 36 (2006) 354-360.",&
     "10.1016/j.commatsci.2005.04.010"),&
     citation("bader_2007","E. Sanville, S. D. Kenny, R. Smith, and G. Henkelman, J. Comput. Chem. 28 (2007) 899-908.",&
     "10.1002/jcc.20575"),&
     citation("bader_2009","W. Tang, E. Sanville, and G. Henkelman, J. Phys.: Condens. Matter 21 (2009) 084204.",&
     "10.1088/0953-8984/21/8/084204"),&
     citation("yt_2011","M. Yu and D. R. Trinkle, J. Chem. Phys. 134 (2011) 064111.",&
     "10.1063/1.3553716"),&
     citation("hirshfeld_i_2007","P. Bultinck, C. Van Alsenoy, P. W. Ayers, and R. Carbo-Dorca, J. Chem. Phys. 126 (2007) 144111.",&
     "10.1063/1.2715563"),&
     citation("becke_1988","A. D. Becke, J. Chem. Phys. 88 (1988) 2547-2553.",&
     "10.1063/1.454033"),&
     citation("lebedev_1999","V. I. Lebedev and D. N. Laikov, Dokl. Math. 59 (1999) 477-481.",&
     ""),&
     citation("horton_2017","T. Verstraelen, P. Tecmer, F. Heidar-Zadeh, C. E. Gonzalez-Espinoza, M. Chan, T. D. Kim, K. Boguslawski, S. Fias, S. Vandenbrande, D. Berrocal, and P. W. Ayers, HORTON 2.1.0, http://theochem.github.com/horton/ (2017).",&
     ""),&
     citation("smoothrho_2022","A. Otero-de-la-Roza, J. Chem. Phys. 156 (2022) 224116.",&
     "10.1063/5.0090232"),&
     citation("peach_2008","M. J. G. Peach, P. Benfield, T. Helgaker, and D. J. Tozer, J. Chem. Phys. 128 (2008) 044118.",&
     "10.1063/1.2831900"),&
     citation("sigmahole_2016","M. H. Kolar and P. Hobza, Chem. Rev. 116 (2016) 5155-5187.",&
     "10.1021/acs.chemrev.5b00560"),&
     citation("keast_1986","P. Keast, Comput. Methods Appl. Mech. Eng. 55 (1986) 339-348.",&
     "10.1016/0045-7825(86)90059-9"),&
     citation("cubpack_2003","R. Cools and A. Haegemans, ACM Trans. Math. Softw. 29 (2003) 287-296.",&
     "10.1145/838250.838253"),&
     citation("rodriguez_2009","J. I. Rodriguez, A. M. Koster, P. W. Ayers, A. Santos-Valle, A. Vela, and G. Merino, J. Comput. Chem. 30 (2009) 1082-1092.",&
     "10.1002/jcc.21134"),&
     citation("blatov_2004","V. A. Blatov, Crystallogr. Rev. 10 (2004) 249-318.",&
     "10.1080/08893110412331323170"),&
     citation("blatov_2014","V. A. Blatov, A. P. Shevchenko, and D. M. Proserpio, Cryst. Growth Des. 14 (2014) 3576-3586.",&
     "10.1021/cg500498k"),&
     citation("amd_2022","D. Widdowson, M. M. Mosca, A. Pulido, A. I. Cooper, and V. Kurlin, MATCH Commun. Math. Comput. Chem. 87 (2022) 529-559.",&
     "10.46793/match.87-3.529W"),&
     citation("gpwdf_2024","A. Otero-de-la-Roza, J. Appl. Cryst. 57 (2024) 1401-1414.",&
     "10.1107/S1600576724007489"),&
     citation("degelder_2001","R. de Gelder, R. Wehrens, and J. A. Hageman, J. Comput. Chem. 22 (2001) 273-289.",&
     "10.1002/1096-987X(200102)22:3<273::AID-JCC1001>3.0.CO;2-0"),&
     citation("emd_2000","Y. Rubner, C. Tomasi, and L. J. Guibas, Int. J. Comput. Vis. 40 (2000) 99-121.",&
     "10.1023/A:1026543900054"),&
     citation("vcpwdf_2022","R. A. Mayo, A. Otero-de-la-Roza, and E. R. Johnson, CrystEngComm 24 (2022) 8326-8338.",&
     "10.1039/D2CE01080A"),&
     citation("econ_2016","M. Nespolo, Acta Cryst. B 72 (2016) 51-66.",&
     "10.1107/S2052520615019472"),&
     citation("econ_1979","R. Hoppe, Z. Kristallogr. 150 (1979) 23-52.",&
     "10.1524/zkri.1979.150.14.23"),&
     citation("hinuma_2017","Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, and I. Tanaka, Comput. Mater. Sci. 128 (2017) 140-184.",&
     "10.1016/j.commatsci.2016.10.015"),&
     citation("xrpd_bg_2001","W. I. F. David and D. S. Sivia, J. Appl. Cryst. 34 (2001) 318-324.",&
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
  end interface

end module citations
