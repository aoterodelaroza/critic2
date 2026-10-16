! Copyright (c) 2015-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Evaluation of arithmetic expressions. This entire module is
! THREAD-SAFE (or should be, at least).
module arithmetic
  use iso_c_binding, only: c_ptr
  implicit none

  ! token type
  type token
     integer :: type = 0 ! type of token (see list in implementation module)
     real*8 :: fval = 0d0 ! associated real value
     integer :: ival = 0 ! associated integer value
     character(len=:), allocatable :: sval ! associated string
     character*10 :: fder = "" ! associated modifier
  end type token
  public :: token

  public :: eval
  public :: eval_grid
  public :: fields_in_eval
  public :: setvariable
  public :: isvariable
  public :: clearvariable
  public :: clearallvariables
  public :: listvariables
  public :: listlibxc
  public :: pretokenize

  private

  ! enum for operators and functions
  integer, parameter :: fun_unknown     = 0  !< unknown operator
  integer, parameter :: fun_openpar     = 1  !< open parenthesis
  integer, parameter :: fun_closepar    = 2  !< close parenthesis
  integer, parameter :: fun_uplus       = 3  !< unary +
  integer, parameter :: fun_uminus      = 4  !< unary -
  integer, parameter :: fun_abs         = 5  !< abs(.)
  integer, parameter :: fun_exp         = 6  !< exp(.)
  integer, parameter :: fun_sqrt        = 7  !< sqrt(.)
  integer, parameter :: fun_floor       = 8  !< floor(.)
  integer, parameter :: fun_ceiling     = 9  !< ceiling(.)
  integer, parameter :: fun_round       = 10 !< round(.)
  integer, parameter :: fun_log         = 11 !< log(.)
  integer, parameter :: fun_log10       = 12 !< log10(.)
  integer, parameter :: fun_sin         = 13 !< sin(.)
  integer, parameter :: fun_asin        = 14 !< asin(.)
  integer, parameter :: fun_cos         = 15 !< cos(.)
  integer, parameter :: fun_acos        = 16 !< acos(.)
  integer, parameter :: fun_tan         = 17 !< tan(.)
  integer, parameter :: fun_atan        = 18 !< atan(.)
  integer, parameter :: fun_atan2       = 19 !< atan2(.)
  integer, parameter :: fun_sinh        = 20 !< sinh(.)
  integer, parameter :: fun_cosh        = 21 !< cosh(.)
  integer, parameter :: fun_erf         = 22 !< erf(.)
  integer, parameter :: fun_erfc        = 23 !< erfc(.)
  integer, parameter :: fun_min         = 24 !< min(.)
  integer, parameter :: fun_max         = 25 !< max(.)
  integer, parameter :: fun_power       = 26 !< **
  integer, parameter :: fun_leq         = 27 !< <=
  integer, parameter :: fun_geq         = 28 !< >=
  integer, parameter :: fun_equal       = 29 !< ==
  integer, parameter :: fun_neq         = 30 !< !=
  integer, parameter :: fun_and         = 31 !< &&
  integer, parameter :: fun_or          = 32 !< ||
  integer, parameter :: fun_plus        = 33 !< binary +
  integer, parameter :: fun_minus       = 34 !< binary -
  integer, parameter :: fun_prod        = 35 !< *
  integer, parameter :: fun_div         = 36 !< /
  integer, parameter :: fun_modulo      = 37 !< %
  integer, parameter :: fun_great       = 38 !< >
  integer, parameter :: fun_less        = 39 !< <
  integer, parameter :: fun_xc          = 40 !< xc(.,...)  (X)
  integer, parameter :: fun_gtf         = 41 !< Thomas-Fermi kinetic energy density
  integer, parameter :: fun_vtf         = 42 !< Potential energy density calcd using gtf and the local virial theorem
  integer, parameter :: fun_htf         = 43 !< Total energy density calcd using gtf and the local virial theorem
  integer, parameter :: fun_gtf_kir     = 44 !< Thomas-Fermi ked with Kirzhnits gradient correction
  integer, parameter :: fun_vtf_kir     = 45 !< Potential energy density calcd using gtf_kir and the local virial theorem
  integer, parameter :: fun_htf_kir     = 46 !< Total energy density calcd using gtf_kir and the local virial theorem
  integer, parameter :: fun_gkin        = 47 !< Kinetic energy density, G-version (grho * grho)
  integer, parameter :: fun_kkin        = 48 !< Kinetic energy density, K-version (rho * laprho)
  integer, parameter :: fun_l           = 49 !< Lagrangian density (-1/4 * laprho)
  integer, parameter :: fun_elf         = 50 !< Electron localization function (ELF)
  integer, parameter :: fun_vir         = 51 !< Electronic potential energy density, virial field
  integer, parameter :: fun_he          = 52 !< Energy density, G + V
  integer, parameter :: fun_lol         = 53 !< Localized-orbital locator (LOL)
  integer, parameter :: fun_lol_kir     = 54 !< Localized-orbital locator (LOL), with Kirzhnits G
  integer, parameter :: fun_brhole_a1   = 55 !< BR hole, A prefactor, spin up
  integer, parameter :: fun_brhole_a2   = 56 !< BR hole, A prefactor, spin down
  integer, parameter :: fun_brhole_a    = 57 !< BR hole, A prefactor, spin average
  integer, parameter :: fun_brhole_alf1 = 58 !< BR hole, alpha exponent, spin up
  integer, parameter :: fun_brhole_alf2 = 59 !< BR hole, alpha exponent, spin down
  integer, parameter :: fun_brhole_alf  = 60 !< BR hole, alpha exponent, spin average
  integer, parameter :: fun_brhole_b1   = 61 !< BR hole, hole distance (b), spin up
  integer, parameter :: fun_brhole_b2   = 62 !< BR hole, hole distance (b), spin down
  integer, parameter :: fun_brhole_b    = 63 !< BR hole, hole distance (b), spin average
  integer, parameter :: fun_xhcurv1     = 64 !< exchange-hole curvature, spin up
  integer, parameter :: fun_xhcurv2     = 65 !< exchange-hole curvature, spin down
  integer, parameter :: fun_xhcurv      = 66 !< exchange-hole curvature, spin average
  integer, parameter :: fun_dsigs1      = 67 !< leading coefficient same-spin pair density, spin up
  integer, parameter :: fun_dsigs2      = 68 !< leading coefficient same-spin pair density, spin down
  integer, parameter :: fun_dsigs       = 69 !< leading coefficient same-spin pair density, spin avg.
  integer, parameter :: fun_mep         = 70 !< Molecular electrostatic potential
  integer, parameter :: fun_uslater     = 71 !< Slater potential
  integer, parameter :: fun_nheff       = 72 !< Effective exchange hole normalization
  integer, parameter :: fun_xhole       = 73 !< exchange-hole
  integer, parameter :: fun_rdg         = 74 !< reduced density gradient (RDG)

  ! enum for structural variables
  integer, parameter :: svar_dnuc    = 1  !< Distance to the closest nucleus
  integer, parameter :: svar_xnucx   = 2  !< x of the closest nucleus (crystallographic)
  integer, parameter :: svar_ynucx   = 3  !< y of the closest nucleus (crystallographic)
  integer, parameter :: svar_znucx   = 4  !< z of the closest nucleus (crystallographic)
  integer, parameter :: svar_xnucc   = 5  !< x of the closest nucleus (Cartesian)
  integer, parameter :: svar_ynucc   = 6  !< y of the closest nucleus (Cartesian)
  integer, parameter :: svar_znucc   = 7  !< z of the closest nucleus (Cartesian)
  integer, parameter :: svar_x       = 8  !< x of the evaluation point (default units)
  integer, parameter :: svar_y       = 9  !< y of the evaluation point (default units)
  integer, parameter :: svar_z       = 10 !< z of the evaluation point (default units)
  integer, parameter :: svar_xx      = 11 !< x of the evaluation point (crystallographic)
  integer, parameter :: svar_yx      = 12 !< y of the evaluation point (crystallographic)
  integer, parameter :: svar_zx      = 13 !< z of the evaluation point (crystallographic)
  integer, parameter :: svar_xc      = 14 !< x of the evaluation point (Cartesian)
  integer, parameter :: svar_yc      = 15 !< y of the evaluation point (Cartesian)
  integer, parameter :: svar_zc      = 16 !< z of the evaluation point (Cartesian)
  integer, parameter :: svar_xm      = 17 !< x of the evaluation point (Cartesian molecular)
  integer, parameter :: svar_ym      = 18 !< y of the evaluation point (Cartesian molecular)
  integer, parameter :: svar_zm      = 19 !< z of the evaluation point (Cartesian molecular)
  integer, parameter :: svar_xxr     = 20 !< x of the evaluation point (reduced cryst.)
  integer, parameter :: svar_yxr     = 21 !< y of the evaluation point (reduced cryst.)
  integer, parameter :: svar_zxr     = 22 !< z of the evaluation point (reduced cryst.)
  integer, parameter :: svar_idnuc   = 23 !< complete-list id of the closest nucleus
  integer, parameter :: svar_nidnuc  = 24 !< non-equivalent-list id of the closest nucleus
  integer, parameter :: svar_rho0nuc = 25 !< atomic density contribution from the closest nucleus
  integer, parameter :: svar_spcnuc  = 26 !< species id of the closest nucleus
  integer, parameter :: svar_zatnuc  = 27 !< atomic number of the closest nucleus

  ! The functions and structural variables known to the parser: the
  ! parser reads the names from this table, and the GUI uses it to
  ! list them with their descriptions.
  integer, parameter, public :: arithkind_math = 1 !< function of numbers
  integer, parameter, public :: arithkind_chem = 2 !< chemical function, takes a field id
  integer, parameter, public :: arithkind_svar = 3 !< structural variable (@name)
  integer, parameter, public :: arithreq_none = 0 !< any field
  integer, parameter, public :: arithreq_gkin = 1 !< fields with the kinetic energy density (wavefunctions)
  integer, parameter, public :: arithreq_stress = 2 !< fields with the stress tensor (molecular wavefunctions)
  integer, parameter, public :: arithreq_basis = 3 !< wavefunctions with basis set, and libcint
  integer, parameter, public :: arithreq_rhf = 4 !< restricted molecular wavefunctions
  integer, parameter, public :: arithreq_libxc = 5 !< critic2 compiled with libxc
  type arith_name
     character(len=12) :: name !< name in the expression
     integer :: code !< fun_* or svar_* code
     integer :: kind !< arithkind_*
     integer :: req !< arithreq_*, what the field (or critic2) must provide
     character(len=16) :: args !< the arguments: (x), (id),...; ":n" for svars that take an atom
     character(len=80) :: desc !< description
  end type arith_name
  public :: arith_name
  type(arith_name), parameter, public :: arith_names(*) = (/&
     arith_name("abs",fun_abs,arithkind_math,arithreq_none,"(x)","Absolute value"),&
     arith_name("exp",fun_exp,arithkind_math,arithreq_none,"(x)","Exponential"),&
     arith_name("sqrt",fun_sqrt,arithkind_math,arithreq_none,"(x)","Square root"),&
     arith_name("floor",fun_floor,arithkind_math,arithreq_none,"(x)","Largest integer not above x"),&
     arith_name("ceiling",fun_ceiling,arithkind_math,arithreq_none,"(x)","Smallest integer not below x"),&
     arith_name("ceil",fun_ceiling,arithkind_math,arithreq_none,"(x)","Same as ceiling"),&
     arith_name("round",fun_round,arithkind_math,arithreq_none,"(x)","Nearest integer"),&
     arith_name("log",fun_log,arithkind_math,arithreq_none,"(x)","Natural logarithm"),&
     arith_name("log10",fun_log10,arithkind_math,arithreq_none,"(x)","Base-10 logarithm"),&
     arith_name("sin",fun_sin,arithkind_math,arithreq_none,"(x)","Sine"),&
     arith_name("asin",fun_asin,arithkind_math,arithreq_none,"(x)","Arc sine"),&
     arith_name("cos",fun_cos,arithkind_math,arithreq_none,"(x)","Cosine"),&
     arith_name("acos",fun_acos,arithkind_math,arithreq_none,"(x)","Arc cosine"),&
     arith_name("tan",fun_tan,arithkind_math,arithreq_none,"(x)","Tangent"),&
     arith_name("atan",fun_atan,arithkind_math,arithreq_none,"(x)","Arc tangent"),&
     arith_name("atan2",fun_atan2,arithkind_math,arithreq_none,"(y,x)","Arc tangent of y/x, in the quadrant of (x,y)"),&
     arith_name("sinh",fun_sinh,arithkind_math,arithreq_none,"(x)","Hyperbolic sine"),&
     arith_name("cosh",fun_cosh,arithkind_math,arithreq_none,"(x)","Hyperbolic cosine"),&
     arith_name("erf",fun_erf,arithkind_math,arithreq_none,"(x)","Error function"),&
     arith_name("erfc",fun_erfc,arithkind_math,arithreq_none,"(x)","Complementary error function"),&
     arith_name("min",fun_min,arithkind_math,arithreq_none,"(x,y)","Smaller of x and y"),&
     arith_name("max",fun_max,arithkind_math,arithreq_none,"(x,y)","Larger of x and y"),&
     arith_name("gtf",fun_gtf,arithkind_chem,arithreq_none,"(id)","Thomas-Fermi kinetic energy density"),&
     arith_name("vtf",fun_vtf,arithkind_chem,arithreq_none,"(id)","Potential energy density from gtf and the local virial theorem"),&
     arith_name("htf",fun_htf,arithkind_chem,arithreq_none,"(id)","Total energy density from gtf and the local virial theorem"),&
     arith_name("gtf_kir",fun_gtf_kir,arithkind_chem,arithreq_none,"(id)","Thomas-Fermi kinetic energy density with Kirzhnits correction"),&
     arith_name("vtf_kir",fun_vtf_kir,arithkind_chem,arithreq_none,"(id)","Potential energy density from gtf_kir and the local virial theorem"),&
     arith_name("htf_kir",fun_htf_kir,arithkind_chem,arithreq_none,"(id)","Total energy density from gtf_kir and the local virial theorem"),&
     arith_name("lag",fun_l,arithkind_chem,arithreq_none,"(id)","Lagrangian density (-1/4 of the Laplacian)"),&
     arith_name("lol_kir",fun_lol_kir,arithkind_chem,arithreq_none,"(id)","Localized-orbital locator (LOL) with the Kirzhnits kinetic energy"),&
     arith_name("rdg",fun_rdg,arithkind_chem,arithreq_none,"(id)","Reduced density gradient"),&
     arith_name("gkin",fun_gkin,arithkind_chem,arithreq_gkin,"(id)","Kinetic energy density, G version"),&
     arith_name("kkin",fun_kkin,arithkind_chem,arithreq_gkin,"(id)","Kinetic energy density, K version"),&
     arith_name("vir",fun_vir,arithkind_chem,arithreq_stress,"(id)","Electronic potential energy density (virial field)"),&
     arith_name("he",fun_he,arithkind_chem,arithreq_stress,"(id)","Electronic energy density (vir + gkin)"),&
     arith_name("elf",fun_elf,arithkind_chem,arithreq_gkin,"(id)","Electron localization function (ELF)"),&
     arith_name("lol",fun_lol,arithkind_chem,arithreq_gkin,"(id)","Localized-orbital locator (LOL)"),&
     arith_name("brhole_a1",fun_brhole_a1,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, A prefactor, spin up"),&
     arith_name("brhole_a2",fun_brhole_a2,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, A prefactor, spin down"),&
     arith_name("brhole_a",fun_brhole_a,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, A prefactor, spin average"),&
     arith_name("brhole_alf1",fun_brhole_alf1,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, exponent, spin up"),&
     arith_name("brhole_alf2",fun_brhole_alf2,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, exponent, spin down"),&
     arith_name("brhole_alf",fun_brhole_alf,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, exponent, spin average"),&
     arith_name("brhole_b1",fun_brhole_b1,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, distance b, spin up"),&
     arith_name("brhole_b2",fun_brhole_b2,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, distance b, spin down"),&
     arith_name("brhole_b",fun_brhole_b,arithkind_chem,arithreq_gkin,"(id)","Becke-Roussel hole, distance b, spin average"),&
     arith_name("xhcurv1",fun_xhcurv1,arithkind_chem,arithreq_gkin,"(id)","Exchange-hole curvature, spin up"),&
     arith_name("xhcurv2",fun_xhcurv2,arithkind_chem,arithreq_gkin,"(id)","Exchange-hole curvature, spin down"),&
     arith_name("xhcurv",fun_xhcurv,arithkind_chem,arithreq_gkin,"(id)","Exchange-hole curvature, spin average"),&
     arith_name("dsigs1",fun_dsigs1,arithkind_chem,arithreq_gkin,"(id)","Same-spin pair density, leading coefficient, spin up"),&
     arith_name("dsigs2",fun_dsigs2,arithkind_chem,arithreq_gkin,"(id)","Same-spin pair density, leading coefficient, spin down"),&
     arith_name("dsigs",fun_dsigs,arithkind_chem,arithreq_gkin,"(id)","Same-spin pair density, leading coefficient, spin average"),&
     arith_name("mep",fun_mep,arithkind_chem,arithreq_basis,"(id)","Molecular electrostatic potential"),&
     arith_name("uslater",fun_uslater,arithkind_chem,arithreq_basis,"(id)","Slater exchange potential"),&
     arith_name("nheff",fun_nheff,arithkind_chem,arithreq_basis,"(id)","Reverse Becke-Roussel effective hole normalization"),&
     arith_name("xhole",fun_xhole,arithkind_chem,arithreq_rhf,"(id,x,y,z)","Exchange hole with the reference point at (x,y,z)"),&
     arith_name("xc",fun_xc,arithkind_chem,arithreq_libxc,"(rho,...,idx)","Exchange-correlation energy density of libxc functional idx"),&
     arith_name("dnuc",svar_dnuc,arithkind_svar,arithreq_none,":n","Distance to the nearest atom (to atom n, with :n)"),&
     arith_name("xnucx",svar_xnucx,arithkind_svar,arithreq_none,":n","x of the nearest atom (of atom n, with :n), crystallographic"),&
     arith_name("ynucx",svar_ynucx,arithkind_svar,arithreq_none,":n","y of the nearest atom (of atom n, with :n), crystallographic"),&
     arith_name("znucx",svar_znucx,arithkind_svar,arithreq_none,":n","z of the nearest atom (of atom n, with :n), crystallographic"),&
     arith_name("xnucc",svar_xnucc,arithkind_svar,arithreq_none,":n","x of the nearest atom (of atom n, with :n), Cartesian"),&
     arith_name("ynucc",svar_ynucc,arithkind_svar,arithreq_none,":n","y of the nearest atom (of atom n, with :n), Cartesian"),&
     arith_name("znucc",svar_znucc,arithkind_svar,arithreq_none,":n","z of the nearest atom (of atom n, with :n), Cartesian"),&
     arith_name("x",svar_x,arithkind_svar,arithreq_none,"","x of the point (fractional in crystals, Cartesian in molecules)"),&
     arith_name("y",svar_y,arithkind_svar,arithreq_none,"","y of the point (fractional in crystals, Cartesian in molecules)"),&
     arith_name("z",svar_z,arithkind_svar,arithreq_none,"","z of the point (fractional in crystals, Cartesian in molecules)"),&
     arith_name("xx",svar_xx,arithkind_svar,arithreq_none,"","x of the point, crystallographic"),&
     arith_name("yx",svar_yx,arithkind_svar,arithreq_none,"","y of the point, crystallographic"),&
     arith_name("zx",svar_zx,arithkind_svar,arithreq_none,"","z of the point, crystallographic"),&
     arith_name("xc",svar_xc,arithkind_svar,arithreq_none,"","x of the point, Cartesian (bohr)"),&
     arith_name("yc",svar_yc,arithkind_svar,arithreq_none,"","y of the point, Cartesian (bohr)"),&
     arith_name("zc",svar_zc,arithkind_svar,arithreq_none,"","z of the point, Cartesian (bohr)"),&
     arith_name("xm",svar_xm,arithkind_svar,arithreq_none,"","x of the point, Cartesian (molecular frame and units)"),&
     arith_name("ym",svar_ym,arithkind_svar,arithreq_none,"","y of the point, Cartesian (molecular frame and units)"),&
     arith_name("zm",svar_zm,arithkind_svar,arithreq_none,"","z of the point, Cartesian (molecular frame and units)"),&
     arith_name("xxr",svar_xxr,arithkind_svar,arithreq_none,"","x of the point, crystallographic in the reduced cell"),&
     arith_name("yxr",svar_yxr,arithkind_svar,arithreq_none,"","y of the point, crystallographic in the reduced cell"),&
     arith_name("zxr",svar_zxr,arithkind_svar,arithreq_none,"","z of the point, crystallographic in the reduced cell"),&
     arith_name("idnuc",svar_idnuc,arithkind_svar,arithreq_none,"","Complete-list id of the nearest atom"),&
     arith_name("nidnuc",svar_nidnuc,arithkind_svar,arithreq_none,"","Non-equivalent-list id of the nearest atom"),&
     arith_name("rho0nuc",svar_rho0nuc,arithkind_svar,arithreq_none,":n","Atomic density of the nearest atom (of atom n, with :n)"),&
     arith_name("spcnuc",svar_spcnuc,arithkind_svar,arithreq_none,"","Species id of the nearest atom"),&
     arith_name("zatnuc",svar_zatnuc,arithkind_svar,arithreq_none,"","Atomic number of the nearest atom")/)
#ifdef HAVE_LIBXC
  logical, parameter, public :: arith_have_libxc = .true. !< xc() is available
#else
  logical, parameter, public :: arith_have_libxc = .false. !< xc() is available
#endif

  ! module procedure interfaces
  interface
     recursive module function eval(expr,errmsg,x0,sptr,periodic,toklistin)
       real*8 :: eval
       character(*), intent(in) :: expr
       character(len=:), allocatable, intent(inout) :: errmsg
       real*8, intent(in), optional :: x0(3)
       type(c_ptr), intent(in), optional :: sptr
       logical, intent(in), optional :: periodic
       type(token), intent(in), optional, target :: toklistin(:)
     end function eval
     module subroutine eval_grid(n,expr,sptr,f,iok)
       integer, intent(in) :: n(3)
       character(*), intent(in) :: expr
       type(c_ptr), intent(in) :: sptr
       real*8, intent(out) :: f(n(1),n(2),n(3))
       logical, intent(out) :: iok
     end subroutine eval_grid
     module subroutine fields_in_eval(expr,errmsg,n,idlist,sptr,chemargs)
       use param, only: mlen
       character(*), intent(in) :: expr
       character(len=:), allocatable, intent(inout) :: errmsg
       integer, intent(out) :: n
       character(len=mlen), allocatable, intent(inout) :: idlist(:)
       type(c_ptr), intent(in) :: sptr
       logical, intent(in), optional :: chemargs
     end subroutine fields_in_eval
     module subroutine setvariable(ikey,ival)
       character*(*), intent(in) :: ikey
       real*8, intent(in) :: ival
     end subroutine setvariable
     module function isvariable(ikey,ival)
       character*(*), intent(in) :: ikey
       real*8, intent(out) :: ival
       logical :: isvariable
     end function isvariable
     module subroutine clearvariable(ikey)
       character*(*), intent(in) :: ikey
     end subroutine clearvariable
     module subroutine clearallvariables()
     end subroutine clearallvariables
     module subroutine listvariables()
     end subroutine listvariables
     module subroutine listlibxc(doref,doname,doflags)
       logical, intent(in) :: doref
       logical, intent(in) :: doname
       logical, intent(in) :: doflags
     end subroutine listlibxc
     module subroutine pretokenize(expr,toklist,errmsg,sptr)
       character(*), intent(in) :: expr
       type(token), intent(inout), allocatable :: toklist(:)
       character(len=:), allocatable, intent(inout) :: errmsg
       type(c_ptr), intent(in), optional :: sptr
     end subroutine pretokenize
  end interface

end module arithmetic
