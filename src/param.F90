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

! Mathematical and physical constants, color specifications and factorials.
module param
  use hashmod, only: hash
  implicit none

  public

  ! named constants
  real*8, parameter :: VBIG     = huge(1d0) !< a very big number
  real*8, parameter :: VSMALL   = 1d-80   !< a very small number
  real*8, parameter :: pi       = 3.14159265358979323846d0 !< pi
  real*8, parameter :: log2     = log(2d0) !< ln(2)
  real*8, parameter :: sqpi     = sqrt(pi) !< sqrt(pi)
  real*8, parameter :: sqfp     = sqrt(4d0*pi) !< sqrt(4pi)
  real*8, parameter :: tpi      = 2d0 * pi !< 2pi
  real*8, parameter :: tpi2     = tpi*tpi !< (2pi)^2
  real*8, parameter :: fourpi   = 4d0 * pi !< 4pi
  real*8, parameter :: rad      = pi / 180d0 !< pi / 180
  real*8, parameter :: cte      = 2.71828182845904523536d0 !< e
  real*8, parameter :: ctsq2    = sqrt(2d0) !< sqrt(2)
  real*8, parameter :: ctsq3    = sqrt(3d0) !< sqrt(3)
  real*8, parameter :: ctsq32   = ctsq3 / 2d0 !< sqrt(3) / 2
  real*8, parameter :: cteuler  = 0.57721566490153286061d0 !< gamma
  real*8, parameter :: ctgold   = 1.61803398874989484820d0 !< golden ratio (1+sqrt(5))/2
  real*8, parameter :: bohrtoa  =  0.52917720859d0 !< bohr to angstrom conversion factor
  real*8, parameter :: bohrtocm = 0.52917720859d-8      !< bohr -> cm (nist2006)
  real*8, parameter :: bohrtom = 0.52917720859d-10      !< bohr -> m (nist2006)
  real*8, parameter :: bohrtonm = 0.052917720859d2      !< bohr -> nm (nist2006)
  real*8, parameter :: bohrtopm = 0.52917720859d2       !< bohr -> pm (nist2006)
  real*8, parameter :: hartocm1 =  2.194746313705d5 !< hartree to cm-1 factor
  real*8, parameter :: hartoev  =  27.211386245988d0 !< hartree to ev factor
  real*8, parameter :: autogpa  =  29421.0108037190d0 !< hartree/bohr**3 to gpa
  real*8, parameter :: eva3togpa = 160.2176487028540d0 !< eV/ang**3 to gpa
  real*8, parameter :: pcamu = 1.660538782d-24 !< atomic mass unit [g] (nist2006)
  real*8, parameter :: cm1tothz = 2.99792458d-2 !< cm-1 to THz conversion factor (c*1e-10)
  real*8, parameter :: zero=0d0 !< 0
  real*8, parameter :: one=1d0 !< 1
  real*8, parameter :: two=2d0 !< 2
  real*8, parameter :: three=3d0 !< 3
  real*8, parameter :: four=4d0 !< 4
  real*8, parameter :: five=5d0 !< 5
  real*8, parameter :: half=5d-1 !< 1/2
  real*8, parameter :: third=one/three !< 1/3
  real*8, parameter :: fourth=one/four !< 1/4
  real*8, parameter :: fiveth=one/five !< 1/5
  real*8, parameter :: twothird=two/three !< 2/3
  real*8, parameter :: sixth=one/6d0 !< 1/6
  real*8, parameter :: fivesixth=five/6d0 !< 5/6
  real*8, parameter :: threefourth=three/four !< 3/4
  real*8, parameter :: eye(3,3) = reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),shape(eye))
  real*8, parameter :: eyet(3,4) = reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0/),shape(eyet))
  complex*16, parameter :: img = (0d0,1d0)
  character*1, parameter :: tab = char(9)
  character(*), parameter :: newline = new_line('a')

  ! limits for critic2
  integer, parameter :: maxzat = 118 ! maximum atomic number supported
  integer, parameter :: maxzat0 = 118+5 ! maximum atomic number supported plus critical point types
  integer, parameter :: mlen = 512 ! length for fixed-length strings (filenames, etc.)
  integer, parameter :: mmlen = 2048 ! long length for fixed-length strings

  ! Enumerate for coordinate types
  integer, parameter, public :: icrd_cart = 0 ! Cartesian
  integer, parameter, public :: icrd_crys = 1 ! Crystallographic
  integer, parameter, public :: icrd_rcrys = 2 ! Reduced crystallographic

  ! Enumerate for structure formats
  integer, parameter, public :: isformat_unknown = 0
  integer, parameter, public :: isformat_from_input = 1
  integer, parameter, public :: isformat_from_library = 2
  integer, parameter, public :: isformat_derived = 3
  integer, parameter, public :: isformat_cif = 4
  integer, parameter, public :: isformat_shelx = 5
  integer, parameter, public :: isformat_f21 = 6
  integer, parameter, public :: isformat_cube = 7
  integer, parameter, public :: isformat_bincube = 8
  integer, parameter, public :: isformat_struct = 9
  integer, parameter, public :: isformat_abinit = 10
  integer, parameter, public :: isformat_elk = 11
  integer, parameter, public :: isformat_qein = 12
  integer, parameter, public :: isformat_qeout = 13
  integer, parameter, public :: isformat_crystal = 14
  integer, parameter, public :: isformat_fploout = 15
  integer, parameter, public :: isformat_xyz = 16
  integer, parameter, public :: isformat_wfn = 17
  integer, parameter, public :: isformat_wfx = 18
  integer, parameter, public :: isformat_fchk = 19
  integer, parameter, public :: isformat_molden = 20
  integer, parameter, public :: isformat_gaussian = 21
  integer, parameter, public :: isformat_siesta = 22
  integer, parameter, public :: isformat_xsf = 23
  integer, parameter, public :: isformat_gen = 24
  integer, parameter, public :: isformat_vasp = 25
  integer, parameter, public :: isformat_pwc = 26
  integer, parameter, public :: isformat_axsf = 27
  integer, parameter, public :: isformat_dat = 28
  integer, parameter, public :: isformat_pgout = 29
  integer, parameter, public :: isformat_orca = 30
  integer, parameter, public :: isformat_dmain = 31
  integer, parameter, public :: isformat_aimsin = 32
  integer, parameter, public :: isformat_aimsout = 33
  integer, parameter, public :: isformat_tinkerfrac = 34
  integer, parameter, public :: isformat_gjf = 35
  integer, parameter, public :: isformat_castepcell = 36
  integer, parameter, public :: isformat_castepgeom = 37
  integer, parameter, public :: isformat_mol2 = 38
  integer, parameter, public :: isformat_pdb = 39
  integer, parameter, public :: isformat_zmat = 40

  ! Enumerate for vibration data formats
  integer, parameter, public :: ivformat_unknown = 0
  integer, parameter, public :: ivformat_matdynmodes = 1
  integer, parameter, public :: ivformat_matdyneig = 2
  integer, parameter, public :: ivformat_qedyn = 3
  integer, parameter, public :: ivformat_phonopy_ascii = 4
  integer, parameter, public :: ivformat_phonopy_yaml = 5
  integer, parameter, public :: ivformat_phonopy_hdf5 = 6
  integer, parameter, public :: ivformat_crystal_out = 7
  integer, parameter, public :: ivformat_gaussian_log = 8
  integer, parameter, public :: ivformat_gaussian_fchk = 9

  ! Enumerate for field formats
  integer, parameter, public :: ifformat_unknown = 0
  integer, parameter, public :: ifformat_wien = 1
  integer, parameter, public :: ifformat_elk = 2
  integer, parameter, public :: ifformat_pi = 3
  integer, parameter, public :: ifformat_cube = 4
  integer, parameter, public :: ifformat_bincube = 5
  integer, parameter, public :: ifformat_abinit = 6
  integer, parameter, public :: ifformat_vasp = 7
  integer, parameter, public :: ifformat_vaspnov = 8
  integer, parameter, public :: ifformat_qub = 9
  integer, parameter, public :: ifformat_xsf = 10
  integer, parameter, public :: ifformat_elkgrid = 11
  integer, parameter, public :: ifformat_siestagrid = 12
  integer, parameter, public :: ifformat_fplogrid = 13
  integer, parameter, public :: ifformat_dftb = 14
  integer, parameter, public :: ifformat_pwc = 15
  integer, parameter, public :: ifformat_wfn = 16
  integer, parameter, public :: ifformat_wfx = 17
  integer, parameter, public :: ifformat_fchk = 18
  integer, parameter, public :: ifformat_molden = 19
  integer, parameter, public :: ifformat_fmt = 20
  integer, parameter, public :: ifformat_as = 21
  integer, parameter, public :: ifformat_as_promolecular = 22
  integer, parameter, public :: ifformat_as_core = 23
  integer, parameter, public :: ifformat_as_lap = 24
  integer, parameter, public :: ifformat_as_grad = 25
  integer, parameter, public :: ifformat_as_pot = 26
  integer, parameter, public :: ifformat_as_resample = 27
  integer, parameter, public :: ifformat_as_clm = 28
  integer, parameter, public :: ifformat_as_clm_sub = 29
  integer, parameter, public :: ifformat_as_ghost = 30
  integer, parameter, public :: ifformat_copy = 31
  integer, parameter, public :: ifformat_promolecular = 32
  integer, parameter, public :: ifformat_promolecular_fragment = 33
  integer, parameter, public :: ifformat_as_hxx1 = 34
  integer, parameter, public :: ifformat_as_hxx2 = 35
  integer, parameter, public :: ifformat_as_hxx3 = 36

  ! Enumerate for molecular and crystal properties. These are used
  ! throughout the code as flags for the calculation of scalar fields.
  integer, parameter :: ims = 5    ! number of items in this enumerate
  integer, parameter :: im_null = 0 ! the void property
  integer, parameter :: im_volume = 1 ! volume (1)
  integer, parameter :: im_rho = 2 ! the electron density (rho)
  integer, parameter :: im_gradrho = 3  ! gradient of the electron density
  integer, parameter :: im_gkin = 4  ! kinetic energy density (gradrho * gradrho)
  integer, parameter :: im_b = 5  ! exchange-hole dipole

  ! free atomic polarizabilities from CRC handbook, 88th ed.
  real*8, parameter :: alpha_free(1:maxzat0) = (/  0.6668D0,  0.2051D0, 24.3300D0,  5.6000D0,& ! 1-4
      3.0300D0,  1.7600D0 , 1.1000D0,  0.8020D0,  0.5570D0,  0.3956D0, 24.1100D0, 10.6000D0,& ! 5-12
      6.8000D0,  5.3800D0,  3.6300D0,  2.9000D0,  2.1800D0,  1.6411D0, 43.4000D0, 22.8000D0,& ! 13-20
     17.8000D0, 14.6000D0, 12.4000D0, 11.6000D0,  9.4000D0,  8.4000D0,  7.5000D0,  6.8000D0,& ! 21-28
      6.2000D0,  5.7500D0,  8.1200D0,  6.0700D0,  4.3100D0,  3.7700D0,  3.0500D0,  2.4844D0,& ! 29-36
     47.3000D0, 27.6000D0, 22.7000D0, 17.9000D0, 15.7000D0, 12.8000D0, 11.4000D0,  9.6000D0,& ! 37-44
      8.6000D0,  4.8000D0,  7.2000D0,  7.3600D0, 10.2000D0,  7.7000D0,  6.6000D0,  5.5000D0,& ! 45-52
      5.3500D0,  4.0440D0, 59.4200D0, 39.7000D0, 31.1000D0, 29.6000D0, 28.2000D0, 31.4000D0,& ! 53-60
     30.1000D0, 28.8000D0, 27.7000D0, 23.5000D0, 25.5000D0, 24.5000D0, 23.6000D0, 22.7000D0,& ! 61-68
     21.8000D0, 21.0000D0, 21.9000D0, 16.2000D0, 13.1000D0, 11.1000D0,  9.7000D0,  8.5000D0,& ! 69-76
      7.6000D0,  6.5000D0,  5.8000D0,  5.0200D0,  7.6000D0,  6.8000D0,  7.4000D0,  6.8000D0,& ! 77-84
      6.0000D0,  5.3000D0, 48.6000D0, 38.3000D0, 32.1000D0, 32.1000D0, 25.4000D0, 24.9000D0,& ! 85-92
     24.8000D0, 24.5000D0, 23.3000D0, 23.0000D0, 22.7000D0, 20.5000D0, 19.7000D0, 23.8000D0,& ! 93-100
     18.2000D0, 17.5000D0,    -999d0,    -999d0,    -999d0,    -999d0,    -999d0,    -999d0,& ! 101-108
        -999d0,    -999d0,    -999d0,    -999d0,    -999d0,    -999d0,    -999d0,    -999d0,& ! 109-116
        -999d0,    -999d0,& ! 117-118
        -999d0,    -999d0,    -999d0,    -999d0,    -999d0/) / .52917720859d0**3 ! 119-124

  ! factorials
  integer, parameter :: mfact = 30  !< factorial vector size
  real*8, dimension(0:mfact) :: fact(0:mfact) !< factorial vector, fact(n) = n!

  ! jmol colors (http://jmol.sourceforge.net/jscolors/)
  ! adapted from tessel.
  integer, parameter :: JMLcol(3,0:maxzat0) = reshape((/&
     135,206,250, & ! 000 Bq
     235,235,235, 207,235,235, 204,128,255, & ! 001-003 H,  He, Li
     194,255,000, 255,181,181, 108,108,108, & ! 004-006 Be, B,  C
     048,080,248, 255,013,013, 144,224,080, & ! 007-009 N,  O,  F
     179,227,245, 171,092,242, 138,255,000, & ! 010-012 Ne, Na, Mg
     191,166,166, 240,200,160, 255,128,000, & ! 013-015 Al, Si, P
     255,255,048, 031,240,031, 128,209,227, & ! 016,018 S,  Cl, Ar
     143,064,212, 061,255,000, 230,230,230, & ! 019-021 K,  Ca, Sc
     191,194,199, 166,166,171, 138,153,199, & ! 022-024 Ti, V,  Cr
     156,122,199, 224,102,051, 240,144,160, & ! 025-027 Mn, Fe, Co
     080,208,080, 200,128,051, 125,128,176, & ! 028-030 Ni, Cu, Zn
     194,143,143, 102,143,143, 189,128,227, & ! 031-033 Ga, Ge, As
     255,161,000, 166,041,041, 092,184,209, & ! 034-036 Se, Br, Kr
     112,046,176, 000,255,000, 148,255,255, & ! 037-039 Rb, Sr, Y
     148,224,224, 115,194,201, 084,181,181, & ! 040-042 Zr, Nb, Mo
     059,158,158, 036,143,143, 010,125,140, & ! 043-045 Tc, Ru, Rh
     000,105,133, 192,192,192, 255,217,143, & ! 046-048 Pd, Ag, Cd
     166,117,115, 102,128,128, 158,099,181, & ! 049-051 In, Sn, Sb
     212,122,000, 148,000,148, 066,158,176, & ! 052-054 Te, I,  Xe
     087,023,143, 000,201,000, 112,212,255, & ! 055-057 Cs, Ba, La
     255,255,199, 217,255,199, 199,255,199, & ! 058-060 Ce, Pr, Nd
     163,255,199, 143,255,199, 097,255,199, & ! 061-063 Pm, Sm, Eu
     069,255,199, 048,255,199, 031,255,199, & ! 064-066 Gd, Tb, Dy
     000,255,156, 000,230,117, 000,212,082, & ! 067-069 Ho, Er, Tm
     000,191,056, 000,171,036, 077,194,255, & ! 070-072 Yb, Lu, Hf
     077,166,255, 033,148,214, 038,125,171, & ! 073-075 Ta, W,  Re
     038,102,150, 023,084,135, 208,208,224, & ! 076-078 Os, Ir, Pt
     255,209,035, 184,184,208, 166,084,077, & ! 079-081 Au, Hg, Tl
     087,089,097, 158,079,181, 171,092,000, & ! 082-084 Pb, Bi, Po
     117,079,069, 066,130,150, 066,000,102, & ! 085-087 At, Rn, Fr
     000,125,000, 112,171,250, 000,186,255, & ! 088-090 Ra, Ac, Th
     000,161,255, 000,143,255, 000,128,255, & ! 091-093 Pa, U,  Np
     000,107,255, 084,092,242, 120,092,227, & ! 094-096 Pu, Am, Cm
     138,079,227, 161,054,212, 179,031,212, & ! 097-099 Bk, Cf, Es
     179,031,186, 179,013,166, 189,013,135, & ! 100-102 Fm, Md, No
     199,000,102, 204,000,089, 209,000,079, & ! 103-105 Lr, Rf, Db
     217,000,069, 224,000,056, 230,000,046, & ! 106-108 Sg, Bh, Hs
     235,000,038, 160,000,066, 015,130,015, & ! 109-111 Mt, Ds, Rg
     020,090,255, 200,000,200, 255,180,070, & ! 112-114 Cn, Nh, Fl
     000,220,220, 230,010,010, 140,255,140, & ! 115-117 Mc, Lv, Ts
     112,112,255, &                           ! 118 Og
     072,159,004, 255,217,061, 149,136,255, & ! 119-121 ncp, bcp, rcp
     255,102,087, 044,255,000&                ! 122-123 ccp, xcp
     /),shape(JMLcol)) !< jmol color definitions

  ! jmol colors, slightly darker
  integer, parameter :: JMLcol2(3,0:maxzat0) = reshape((/&
     095,166,210, & ! 000 Bq
     195,195,195, 157,195,195, 144,680,195, & ! 001-003 H,  He, Li
     134,195,000, 195,121,121, 084,084,084, & ! 004-006 Be, B,  C
     000,020,188, 195,000,000, 084,164,020, & ! 007-009 N,  O,  F
     119,167,185, 111,032,182, 078,195,000, & ! 010-012 Ne, Na, Mg
     131,106,106, 180,140,100, 195,068,000, & ! 013-015 Al, Si, P
     195,195,000, 000,180,000, 068,149,167, & ! 016,018 S,  Cl, Ar
     083,004,152, 001,195,000, 170,170,170, & ! 019-021 K,  Ca, Sc
     131,134,139, 106,106,111, 078,093,139, & ! 022-024 Ti, V,  Cr
     096,062,139, 164,042,000, 180,084,100, & ! 025-027 Mn, Fe, Co
     020,148,020, 140,068,000, 065,068,116, & ! 028-030 Ni, Cu, Zn
     134,083,083, 042,083,083, 129,068,167, & ! 031-033 Ga, Ge, As
     195,101,000, 106,000,000, 032,124,149, & ! 034-036 Se, Br, Kr
     052,000,116, 000,195,000, 088,195,195, & ! 037-039 Rb, Sr, Y
     088,164,164, 055,134,141, 024,121,121, & ! 040-042 Zr, Nb, Mo
     000,098,098, 000,083,083, 000,065,080, & ! 043-045 Tc, Ru, Rh
     000,045,073, 132,132,132, 195,157,083, & ! 046-048 Pd, Ag, Cd
     106,057,055, 042,068,068, 098,039,121, & ! 049-051 In, Sn, Sb
     152,062,000, 088,000,088, 006,098,116, & ! 052-054 Te, I,  Xe
     027,000,083, 000,141,000, 052,152,195, & ! 055-057 Cs, Ba, La
     195,195,139, 157,195,139, 139,195,139, & ! 058-060 Ce, Pr, Nd
     103,195,139, 083,195,139, 037,195,139, & ! 061-063 Pm, Sm, Eu
     009,195,139, 000,195,139, 000,195,139, & ! 064-066 Gd, Tb, Dy
     000,195,096, 000,170,057, 000,152,022, & ! 067-069 Ho, Er, Tm
     000,131,000, 000,111,000, 017,134,195, & ! 070-072 Yb, Lu, Hf
     017,106,195, 000,088,154, 000,065,111, & ! 073-075 Ta, W,  Re
     000,042,090, 000,024,075, 148,148,164, & ! 076-078 Os, Ir, Pt
     195,149,000, 124,124,148, 106,024,017, & ! 079-081 Au, Hg, Tl
     027,029,037, 098,019,121, 111,032,000, & ! 082-084 Pb, Bi, Po
     057,019,009, 006,070,090, 006,000,042, & ! 085-087 At, Rn, Fr
     000,065,000, 052,111,190, 000,126,195, & ! 088-090 Ra, Ac, Th
     000,101,195, 000,083,195, 000,068,195, & ! 091-093 Pa, U,  Np
     000,047,195, 024,032,182, 060,032,167, & ! 094-096 Pu, Am, Cm
     078,019,167, 101,000,152, 119,000,152, & ! 097-099 Bk, Cf, Es
     119,000,126, 119,000,106, 129,000,075, & ! 100-102 Fm, Md, No
     139,000,042, 144,000,029, 149,000,019, & ! 103-105 Lr, Rf, Db
     157,000,009, 164,000,000, 170,000,000, & ! 106-108 Sg, Bh, Hs
     175,000,000, 100,000,006, 000,070,000, & ! 109-111 Mt, Ds, Rg
     000,030,195, 140,000,140, 195,120,010, & ! 112-114 Cn, Nh, Fl
     000,160,160, 170,000,000, 080,195,080, & ! 115-117 Mc, Lv, Ts
     112,112,255, &                           ! 118 Uuh
     072,159,004, 255,217,061, 149,136,255, & ! 119-121 ncp, bcp, rcp
     255,102,087, 044,255,000&                ! 122-123 ccp, xcp
     /),shape(JMLcol2)) !< jmol color definitions, slightly darker

  ! Covalent radii in angstrom from Pyykko et al., http://dx.doi.org/10.1021/jp5065819
  ! (single bonds only)
  real*8, parameter :: atmcov0(0:maxzat0) = (/&
     ! 1       2       3       4       5       6       7       8       9       0
     0.60d0,& ! 0
     0.32d0, 0.46d0, 1.33d0, 1.02d0, 0.85d0, 0.75d0, 0.71d0, 0.63d0, 0.64d0, 0.67d0,& ! 1-10
     1.55d0, 1.39d0, 1.26d0, 1.16d0, 1.11d0, 1.03d0, 0.99d0, 0.96d0, 1.96d0, 1.71d0,& ! 11-20
     1.48d0, 1.36d0, 1.34d0, 1.22d0, 1.19d0, 1.16d0, 1.11d0, 1.10d0, 1.12d0, 1.18d0,& ! 21-30
     1.24d0, 1.21d0, 1.21d0, 1.16d0, 1.14d0, 1.17d0, 2.10d0, 1.85d0, 1.63d0, 1.54d0,& ! 31-40
     1.47d0, 1.38d0, 1.28d0, 1.25d0, 1.25d0, 1.20d0, 1.28d0, 1.36d0, 1.42d0, 1.40d0,& ! 41-50
     1.40d0, 1.36d0, 1.33d0, 1.31d0, 2.32d0, 1.96d0, 1.80d0, 1.63d0, 1.76d0, 1.74d0,& ! 51-60
     1.73d0, 1.72d0, 1.68d0, 1.69d0, 1.68d0, 1.67d0, 1.66d0, 1.65d0, 1.64d0, 1.70d0,& ! 61-70
     1.62d0, 1.52d0, 1.46d0, 1.37d0, 1.31d0, 1.29d0, 1.22d0, 1.23d0, 1.24d0, 1.33d0,& ! 71-80
     1.44d0, 1.44d0, 1.51d0, 1.45d0, 1.47d0, 1.42d0, 2.23d0, 2.01d0, 1.86d0, 1.75d0,& ! 81-90
     1.69d0, 1.70d0, 1.71d0, 1.72d0, 1.66d0, 1.66d0, 1.68d0, 1.68d0, 1.65d0, 1.67d0,& ! 91-100
     1.73d0, 1.76d0, 1.61d0, 1.57d0, 1.49d0, 1.43d0, 1.41d0, 1.34d0, 1.29d0, 1.28d0,& ! 101-110
     1.21d0, 1.22d0, 1.36d0, 1.43d0, 1.62d0, 1.75d0, 1.65d0, 1.57d0, 0.00d0, 0.00d0,& ! 111-120
     0.00d0, 0.00d0, 0.00d0/) / bohrtoa ! 121-123
  real*8 :: atmcov(0:maxzat0) = atmcov0

  ! Van der Waals radii follow the CSD:
  ! From Bondi: http://dx.doi.org/10.1021/j100785a001
  ! except H from Rowland and Taylor http://dx.doi.org/10.1021/jp953141
  ! Elements that do not have radius are assigned 2.00.
  real*8, parameter :: atmvdw0(0:maxzat0) = (/&
      0.60d0,& ! 0
      1.09d0, 1.40d0, 1.82d0, 2.00d0, 2.00d0, 1.70d0, 1.55d0, 1.52d0, 1.47d0, 1.54d0,& ! 1-10
      2.27d0, 1.73d0, 2.00d0, 2.10d0, 1.80d0, 1.80d0, 1.75d0, 1.88d0, 2.75d0, 2.00d0,& ! 11-20
      2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 1.63d0, 1.40d0, 1.39d0,& ! 21-30
      1.87d0, 2.00d0, 1.85d0, 1.90d0, 1.85d0, 2.02d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,& ! 31-40
      2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 1.63d0, 1.72d0, 1.58d0, 1.93d0, 2.17d0,& ! 41-50
      2.00d0, 2.06d0, 1.98d0, 2.16d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,& ! 51-60
      2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,& ! 61-70
      2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 1.72d0, 1.66d0, 1.55d0,& ! 71-80
      1.96d0, 2.02d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,& ! 81-90
      2.00d0, 1.86d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,& ! 91-100
      2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0,& ! 101-1100
      2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 2.00d0, 0.00d0, 0.00d0,& ! 111-118
      0.00d0, 0.00d0, 0.00d0/) / bohrtoa                                               ! 119-123
  real*8 :: atmvdw(0:maxzat0) = atmvdw0

  ! standard atomic weights
  real*8, parameter :: atmass(0:maxzat0) = (/&
        0d0,& ! 0
        1.00794d0, 4.002602d0,     6.941d0, 9.012182d0,    10.811d0, 12.0107d0,   14.0067d0, 15.9994d0,& ! 1-8
     18.9984032d0,  20.1797d0, 22.989770d0,  24.3050d0, 26.981538d0, 28.0855d0, 30.973761d0,  32.065d0,& ! 9-16
         35.453d0,   39.948d0,   39.0983d0,   40.078d0, 44.955910d0,  47.867d0,   50.9415d0, 51.9961d0,& ! 17-24
      54.938049d0,   55.845d0, 58.933200d0,  58.6934d0,    63.546d0,  65.409d0,    69.723d0,   72.64d0,& ! 25-32
       74.92160d0,    78.96d0,    79.904d0,   83.798d0,   85.4678d0,   87.62d0,  88.90585d0,  91.224d0,& ! 33-40
       92.90638d0,    95.94d0,  97.90721d0,   101.07d0, 102.90550d0,  106.42d0,  107.8682d0, 112.411d0,& ! 41-48
        114.818d0,  118.710d0,   121.760d0,   127.60d0, 126.90447d0, 131.293d0, 132.90545d0, 137.327d0,& ! 49-56
       138.9055d0,  140.116d0, 140.90765d0,   144.24d0, 144.91276d0,  150.36d0,   151.964d0,  157.25d0,& ! 57-64
      158.92534d0,  162.500d0, 164.93032d0,  167.259d0, 168.93421d0,  173.04d0,   174.967d0,  178.49d0,& ! 65-72
       180.9479d0,   183.84d0,   186.207d0,   190.23d0,   192.217d0, 195.078d0, 196.96655d0,  200.59d0,& ! 73-80
       204.3833d0,    207.2d0, 208.98038d0,      209d0,       210d0,     222d0,       223d0,     226d0,& ! 81-88
            227d0,   232.04d0,    231.04d0,   238.03d0,       237d0,     244d0,       243d0,     247d0,& ! 89-96
            247d0,      251d0,       254d0,      257d0,       258d0,     259d0,       262d0,     267d0,& ! 97-104
            268d0,      271d0,       272d0,      270d0,       276d0,     281d0,       280d0,     285d0,& ! 105-112
            284d0,      289d0,       288d0,      293d0,       292d0,     294d0,&                         ! 113-118
              0d0,        0d0,         0d0,        0d0,         0d0/)                                    ! 119-123

  ! atomic scattering factors
  ! Coefficients for an analytical approximation to the atomic scattering factors.
  ! From ITC, vol. C, 3rd ed. (2004), tables 6.1.1.4 and 6.1.1.5 (p. 578 and ff.).
  ! order is: a1, b1, a2, b2, a3, b3, a4, b4, c
  real*8, parameter :: cscatt(9,94) = reshape((/&
     0.489918d0, 20.6593d0 , 0.262003d0, 7.74039d0 , 0.196767d0, 49.5519d0 , 0.049879d0 , 2.20159d0 , 0.001305d0,& ! H
     0.8734d0,  9.1037d0,  0.6309d0,  3.3568d0,  0.3112d0, 22.9276d0,   0.178d0, 0.9821d0,  0.0064d0,& ! He
     1.1282d0,  3.9546d0,  0.7508d0,  1.0524d0,  0.6175d0, 85.3905d0,  0.4653d0 , 168.261d0,  0.0377d0,& ! Li
     1.5919d0, 43.6427d0,  1.1278d0,  1.8623d0,  0.5391d0, 103.483d0,  0.7029d0,  0.542d0,  0.0385d0,& ! Be
     2.0545d0, 23.2185d0,  1.3326d0,   1.021d0,  1.0979d0, 60.3498d0,  0.7068d0, 0.1403d0, -0.1932d0,& ! B
     2.31d0, 20.8439d0,    1.02d0, 10.2075d0,  1.5886d0,  0.5687d0,   0.865d0 , 51.6512d0,  0.2156d0,& ! C
     12.2126d0,  0.0057d0,  3.1322d0,  9.8933d0,  2.0125d0, 28.9975d0,  1.1663d0, 0.5826d0, -11.529d0,& ! N
     3.0485d0, 13.2771d0,  2.2868d0,  5.7011d0,  1.5463d0,  0.3239d0,   0.867d0 , 32.9089d0,  0.2508d0,& ! O
     3.5392d0, 10.2825d0,  2.6412d0,  4.2944d0,   1.517d0,  0.2615d0,  1.0243d0 , 26.1476d0,  0.2776d0,& ! F
     3.9553d0,  8.4042d0,  3.1125d0,  3.4262d0,  1.4546d0,  0.2306d0,  1.1251d0 , 21.7184d0,  0.3515d0,& ! Ne
     4.7626d0,   3.285d0,  3.1736d0,  8.8422d0,  1.2674d0,  0.3136d0,  1.1128d0 , 129.424d0,   0.676d0,& ! Na
     5.4204d0,  2.8275d0,  2.1735d0, 79.2611d0,  1.2269d0,  0.3808d0,  2.3073d0, 7.1937d0,  0.8584d0,& ! Mg
     6.4202d0,  3.0387d0,  1.9002d0,  0.7426d0,  1.5936d0, 31.5472d0,  1.9646d0 , 85.0886d0,  1.1151d0,& ! Al
     6.2915d0,  2.4386d0,  3.0353d0, 32.3337d0,  1.9891d0,  0.6785d0,   1.541d0 , 81.6937d0,  1.1407d0,& ! Si
     6.4345d0,  1.9067d0,  4.1791d0,  27.157d0,    1.78d0,   0.526d0,  1.4908d0 , 68.1645d0,  1.1149d0,& ! P
     6.9053d0,  1.4679d0,  5.2034d0, 22.2151d0,  1.4379d0,  0.2536d0,  1.5863d0, 56.172d0,  0.8669d0,& ! S
     11.4604d0,  0.0104d0,  7.1964d0,  1.1662d0,  6.2556d0, 18.5194d0,  1.6455d0 , 47.7784d0, -9.5574d0,& ! Cl
     7.4845d0,  0.9072d0,  6.7723d0, 14.8407d0,  0.6539d0, 43.8983d0,  1.6442d0 , 33.3929d0,  1.4445d0,& ! Ar
     8.2186d0, 12.7949d0,  7.4398d0,  0.7748d0,  1.0519d0, 213.187d0,  0.8659d0 , 41.6841d0,  1.4228d0,& ! K
     8.6266d0, 10.4421d0,  7.3873d0,  0.6599d0,  1.5899d0, 85.7484d0,  1.0211d0 , 178.437d0,  1.3751d0,& ! Ca
     9.189d0,  9.0213d0,  7.3679d0,  0.5729d0,  1.6409d0, 136.108d0,   1.468d0 , 51.3531d0,  1.3329d0,& ! Sc
     9.7595d0,  7.8508d0,  7.3558d0,     0.5d0,  1.6991d0, 35.6338d0,  1.9021d0 , 116.105d0,  1.2807d0,& ! Ti
     10.2971d0,  6.8657d0,  7.3511d0,  0.4385d0,  2.0703d0, 26.8938d0,  2.0571d0 , 102.478d0,  1.2199d0,& ! V
     10.6406d0,  6.1038d0,  7.3537d0,   0.392d0,   3.324d0, 20.2626d0,  1.4922d0 , 98.7399d0,  1.1832d0,& ! Cr
     11.2819d0,  5.3409d0,  7.3573d0,  0.3432d0,  3.0193d0, 17.8674d0,  2.2441d0 , 83.7543d0,  1.0896d0,& ! Mn
     11.7695d0,  4.7611d0,  7.3573d0,  0.3072d0,  3.5222d0, 15.3535d0,  2.3045d0 , 76.8805d0,  1.0369d0,& ! Fe
     12.2841d0,  4.2791d0,  7.3409d0,  0.2784d0,  4.0034d0, 13.5359d0,  2.3488d0 , 71.1692d0,  1.0118d0,& ! Co
     12.8376d0,  3.8785d0,   7.292d0,  0.2565d0,  4.4438d0, 12.1763d0,    2.38d0 , 66.3421d0,  1.0341d0,& ! Ni
     13.338d0,  3.5828d0,  7.1676d0,   0.247d0,  5.6158d0, 11.3966d0,  1.6735d0 , 64.8126d0,   1.191d0,& ! Cu
     14.0743d0,  3.2655d0,  7.0318d0,  0.2333d0,  5.1652d0, 10.3163d0,    2.41d0 , 58.7097d0,  1.3041d0,& ! Zn
     15.2354d0,  3.0669d0,  6.7006d0,  0.2412d0,  4.3591d0, 10.7805d0,  2.9623d0 , 61.4135d0,  1.7189d0,& ! Ga
     16.0816d0,  2.8509d0,  6.3747d0,  0.2516d0,  3.7068d0, 11.4468d0,   3.683d0 , 54.7625d0,  2.1313d0,& ! Ge
     16.6723d0,  2.6345d0,  6.0701d0,  0.2647d0,  3.4313d0, 12.9479d0,  4.2779d0 , 47.7972d0,   2.531d0,& ! As
     17.0006d0,  2.4098d0,  5.8196d0,  0.2726d0,  3.9731d0, 15.2372d0,  4.3543d0 , 43.8163d0,  2.8409d0,& ! Se
     17.1789d0,  2.1723d0,  5.2358d0, 16.5796d0,  5.6377d0,  0.2609d0,  3.9851d0 , 41.4328d0,  2.9557d0,& ! Br
     17.3555d0,  1.9384d0,  6.7286d0, 16.5623d0,  5.5493d0,  0.2261d0,  3.5375d0 , 39.3972d0,   2.825d0,& ! Kr
     17.1784d0,  1.7888d0,  9.6435d0, 17.3151d0,  5.1399d0,  0.2748d0,  1.5292d0 , 164.934d0,  3.4873d0,& ! Rb
     17.5663d0,  1.5564d0,  9.8184d0, 14.0988d0,   5.422d0,  0.1664d0,  2.6694d0 , 132.376d0,  2.5064d0,& ! Sr
     17.776d0,  1.4029d0, 10.2946d0, 12.8006d0, 5.72629d0 , 0.125599d0, 3.26588d0 , 104.354d0, 1.91213d0,& ! Y
     17.8765d0, 1.27618d0,  10.948d0,  11.916d0, 5.41732d0 , 0.117622d0, 3.65721d0 , 87.6627d0, 2.06929d0,& ! Zr
     17.6142d0, 1.18865d0, 12.0144d0,  11.766d0, 4.04183d0 , 0.204785d0, 3.53346d0 , 69.7957d0, 3.75591d0,& ! Nb
     3.7025d0,  0.2772d0, 17.2356d0,  1.0958d0, 12.8876d0,  11.004d0,  3.7429d0 , 61.6584d0,  4.3875d0,& ! Mo
     19.1301d0 , 0.864132d0, 11.0948d0, 8.14487d0, 4.64901d0, 21.5707d0, 2.71263d0 , 86.8472d0, 5.40428d0,& ! Tc
     19.2674d0, 0.80852d0, 12.9182d0, 8.43467d0, 4.86337d0, 24.7997d0, 1.56756d0 , 94.2928d0, 5.37874d0,& ! Ru
     19.2957d0 , 0.751536d0, 14.3501d0, 8.21758d0, 4.73425d0, 25.8749d0, 1.28918d0 , 98.6062d0,   5.328d0,& ! Rh
     19.3319d0 , 0.698655d0, 15.5017d0, 7.98929d0, 5.29537d0, 25.2052d0 , 0.605844d0 , 76.8986d0, 5.26593d0,& ! Pd
     19.2808d0,  0.6446d0, 16.6885d0,  7.4726d0,  4.8045d0, 24.6605d0,  1.0463d0 , 99.8156d0,   5.179d0,& ! Ag
     19.2214d0,  0.5946d0, 17.6444d0,  6.9089d0,   4.461d0, 24.7008d0,  1.6029d0 , 87.4825d0,  5.0694d0,& ! Cd
     19.1624d0,  0.5476d0, 18.5596d0,  6.3776d0,  4.2948d0, 25.8499d0,  2.0396d0 , 92.8029d0,  4.9391d0,& ! In
     19.1889d0,  5.8303d0, 19.1005d0,  0.5031d0,  4.4585d0, 26.8909d0,  2.4663d0 , 83.9571d0,  4.7821d0,& ! Sn
     19.6418d0,  5.3034d0, 19.0455d0,  0.4607d0,  5.0371d0, 27.9074d0,  2.6827d0 , 75.2825d0,  4.5909d0,& ! Sb
     19.9644d0, 4.81742d0, 19.0138d0 , 0.420885d0, 6.14487d0, 28.5284d0,  2.5239d0 , 70.8403d0,   4.352d0,& ! Te
     20.1472d0,   4.347d0, 18.9949d0,  0.3814d0,  7.5138d0,  27.766d0,  2.2735d0 , 66.8776d0,  4.0712d0,& ! I
     20.2933d0,  3.9282d0, 19.0298d0,   0.344d0,  8.9767d0, 26.4659d0,    1.99d0 , 64.2658d0,  3.7118d0,& ! Xe
     20.3892d0,   3.569d0, 19.1062d0,  0.3107d0,  10.662d0, 24.3879d0,  1.4953d0 , 213.904d0,  3.3352d0,& ! Cs
     20.3361d0,   3.216d0,  19.297d0,  0.2756d0,  10.888d0, 20.2073d0,  2.6959d0 , 167.202d0,  2.7731d0,& ! Ba
     20.578d0, 2.94817d0,  19.599d0 , 0.244475d0, 11.3727d0, 18.7726d0, 3.28719d0 , 133.124d0, 2.14678d0,& ! La
     21.1671d0, 2.81219d0, 19.7695d0 , 0.226836d0, 11.8513d0, 17.6083d0, 3.33049d0 , 127.113d0, 1.86264d0,& ! Ce
     22.044d0, 2.77393d0, 19.6697d0 , 0.222087d0, 12.3856d0, 16.7669d0, 2.82428d0 , 143.644d0,  2.0583d0,& ! Pr
     22.6845d0, 2.66248d0, 19.6847d0 , 0.210628d0,  12.774d0,  15.885d0, 2.85137d0 , 137.903d0, 1.98486d0,& ! Nd
     23.3405d0,  2.5627d0, 19.6095d0 , 0.202088d0, 13.1235d0, 15.1009d0, 2.87516d0 , 132.721d0, 2.02876d0,& ! Pm
     24.0042d0, 2.47274d0, 19.4258d0 , 0.196451d0, 13.4396d0, 14.3996d0, 2.89604d0 , 128.007d0, 2.20963d0,& ! Sm
     24.6274d0,  2.3879d0, 19.0886d0,  0.1942d0, 13.7603d0, 13.7546d0,  2.9227d0 , 123.174d0,  2.5745d0,& ! Eu
     25.0709d0, 2.25341d0, 19.0798d0 , 0.181951d0, 13.8518d0, 12.9331d0, 3.54545d0 , 101.398d0,  2.4196d0,& ! Gd
     25.8976d0, 2.24256d0, 18.2185d0 , 0.196143d0, 14.3167d0, 12.6648d0, 2.95354d0 , 115.362d0, 3.58324d0,& ! Tb
     26.507d0,  2.1802d0, 17.6383d0 , 0.202172d0, 14.5596d0, 12.1899d0, 2.96577d0 , 111.874d0, 4.29728d0,& ! Dy
     26.9049d0, 2.07051d0,  17.294d0, 0.19794d0, 14.5583d0, 11.4407d0, 3.63837d0 , 92.6566d0, 4.56796d0,& ! Ho
     27.6563d0, 2.07356d0, 16.4285d0 , 0.223545d0, 14.9779d0, 11.3604d0, 2.98233d0 , 105.703d0, 5.92046d0,& ! Er
     28.1819d0, 2.02859d0, 15.8851d0 , 0.238849d0, 15.1542d0, 10.9975d0, 2.98706d0 , 102.961d0, 6.75621d0,& ! Tm
     28.6641d0,  1.9889d0, 15.4345d0 , 0.257119d0, 15.3087d0, 10.6647d0, 2.98963d0 , 100.417d0, 7.56672d0,& ! Yb
     28.9476d0, 1.90182d0, 15.2208d0, 9.98519d0,    15.1d0 , 0.261033d0, 3.71601d0 , 84.3298d0, 7.97628d0,& ! Lu
     29.144d0, 1.83262d0, 15.1726d0,  9.5999d0, 14.7586d0 , 0.275116d0, 4.30013d0, 72.029d0, 8.58154d0,& ! Hf
     29.2024d0, 1.77333d0, 15.2293d0, 9.37046d0, 14.5135d0 , 0.295977d0, 4.76492d0 , 63.3644d0, 9.24354d0,& ! Ta
     29.0818d0, 1.72029d0,   15.43d0,  9.2259d0, 14.4327d0 , 0.321703d0, 5.11982d0, 57.056d0,  9.8875d0,& ! W
     28.7621d0, 1.67191d0, 15.7189d0, 9.09227d0, 14.5564d0,  0.3505d0, 5.44174d0 , 52.0861d0,  10.472d0,& ! Re
     28.1894d0, 1.62903d0,  16.155d0, 8.97948d0, 14.9305d0 , 0.382661d0, 5.67589d0 , 48.1647d0, 11.0005d0,& ! Os
     27.3049d0, 1.59279d0, 16.7296d0, 8.86553d0, 15.6115d0 , 0.417916d0, 5.83377d0 , 45.0011d0, 11.4722d0,& ! Ir
     27.0059d0, 1.51293d0, 17.7639d0, 8.81174d0, 15.7131d0 , 0.424593d0,  5.7837d0 , 38.6103d0, 11.6883d0,& ! Pt
     16.8819d0,  0.4611d0, 18.5913d0,  8.6216d0, 25.5582d0,  1.4826d0,    5.86d0 , 36.3956d0, 12.0658d0,& ! Au
     20.6809d0,   0.545d0, 19.0417d0,  8.4484d0, 21.6575d0,  1.5729d0,  5.9676d0 , 38.3246d0, 12.6089d0,& ! Hg
     27.5446d0, 0.65515d0, 19.1584d0, 8.70751d0,  15.538d0, 1.96347d0, 5.52593d0 , 45.8149d0, 13.1746d0,& ! Tl
     31.0617d0,  0.6902d0, 13.0637d0,  2.3576d0,  18.442d0,   8.618d0,  5.9696d0 , 47.2579d0, 13.4118d0,& ! Pb
     33.3689d0,   0.704d0,  12.951d0,  2.9238d0, 16.5877d0,  8.7937d0,  6.4692d0 , 48.0093d0, 13.5782d0,& ! Bi
     34.6726d0 , 0.700999d0, 15.4733d0, 3.55078d0, 13.1138d0, 9.55642d0, 7.02588d0 , 47.0045d0,  13.677d0,& ! Po
     35.3163d0, 0.68587d0, 19.0211d0, 3.97458d0, 9.49887d0, 11.3824d0, 7.42518d0 , 45.4715d0, 13.7108d0,& ! At
     35.5631d0,  0.6631d0, 21.2816d0,  4.0691d0,  8.0037d0, 14.0422d0,  7.4433d0 , 44.2473d0, 13.6905d0,& ! Rn
     35.9299d0 , 0.646453d0, 23.0547d0, 4.17619d0, 12.1439d0, 23.1052d0, 2.11253d0 , 150.645d0, 13.7247d0,& ! Fr
     35.763d0 , 0.616341d0, 22.9064d0, 3.87135d0, 12.4739d0, 19.9887d0, 3.21097d0 , 142.325d0, 13.6211d0,& ! Ra
     35.6597d0 , 0.589092d0, 23.1032d0, 3.65155d0, 12.5977d0,  18.599d0, 4.08655d0, 117.02d0, 13.5266d0,& ! Ac
     35.5645d0 , 0.563359d0, 23.4219d0, 3.46204d0, 12.7473d0, 17.8309d0, 4.80703d0 , 99.1722d0, 13.4314d0,& ! Th
     35.8847d0 , 0.547751d0, 23.2948d0, 3.41519d0, 14.1891d0, 16.9235d0, 4.17287d0 , 105.251d0, 13.4287d0,& ! Pa
     36.0228d0,  0.5293d0, 23.4128d0,  3.3253d0, 14.9491d0, 16.0927d0,   4.188d0 , 100.613d0, 13.3966d0,& ! U
     36.1874d0 , 0.511929d0, 23.5964d0, 3.25396d0, 15.6402d0, 15.3622d0,  4.1855d0 , 97.4908d0, 13.3573d0,& ! Np
     36.5254d0 , 0.499384d0, 23.8083d0, 3.26371d0, 16.7707d0, 14.9455d0, 3.47947d0, 105.98d0, 13.3812d0& ! Pu
     /),shape(cscatt))

  real*8, parameter :: c2scatt(4,2:94) = reshape((/&
     0.52543d0, -3.43300d0,  4.80070d0, -2.54760d0, & ! He
     0.89463d0, -2.43660d0,  2.32500d0, -0.71949d0, & ! Li
     1.25840d0, -1.94590d0,  1.30460d0, -0.04297d0, & ! Be
     1.66720d0, -1.85560d0,  1.60440d0, -0.65981d0, & ! B
     1.70560d0, -1.56760d0,  1.18930d0, -0.42715d0, & ! C
     1.54940d0, -1.20190d0,  0.51064d0,  0.02472d0, & ! N
     1.30530d0, -0.83742d0, -0.16738d0,  0.47500d0, & ! O
     1.16710d0, -0.63203d0, -0.40207d0,  0.54352d0, & ! F
     1.09310d0, -0.50221d0, -0.53648d0,  0.60957d0, & ! Ne
     0.84558d0, -0.26294d0, -0.87884d0,  0.76974d0, & ! Na
     0.71877d0, -0.13144d0, -1.20900d0,  0.82738d0, & ! Mg
     0.67975d0, -0.08756d0, -0.95431d0,  0.72294d0, & ! Al
     0.70683d0, -0.09888d0, -0.98356d0,  0.55631d0, & ! Si
     0.85532d0, -0.21262d0, -0.37390d0,  0.20731d0, & ! P
     1.10400d0, -0.40325d0,  0.20094d0, -0.26058d0, & ! S
     1.42320d0, -0.63936d0,  0.84722d0, -0.76135d0, & ! Cl
     1.82020d0, -0.92776d0,  1.59220d0, -1.32510d0, & ! Ar
     2.26550d0, -1.24530d0,  2.38330d0, -1.91290d0, & ! K
     2.71740d0, -1.55670d0,  3.13170d0, -2.45670d0, & ! Ca
     3.11730d0, -1.81380d0,  3.71390d0, -2.85330d0, & ! Sc
     3.45360d0, -2.01150d0,  4.13170d0, -3.11710d0, & ! Ti
     3.71270d0, -2.13920d0,  4.34610d0, -3.22040d0, & ! V
     3.87870d0, -2.19000d0,  4.38670d0, -3.17520d0, & ! Cr
     3.98550d0, -2.18850d0,  4.27960d0, -3.02150d0, & ! Mn
     3.99790d0, -2.11080d0,  3.98170d0, -2.71990d0, & ! Fe
     3.95900d0, -1.99650d0,  3.60630d0, -2.37050d0, & ! Co
     3.86070d0, -1.88690d0,  3.12390d0, -1.94290d0, & ! Ni
     3.72510d0, -1.65500d0,  2.60290d0, -1.49760d0, & ! Cu
     3.55950d0, -1.45100d0,  2.03390d0, -1.02160d0, & ! Zn
     3.37560d0, -1.23910d0,  1.46160d0, -0.55471d0, & ! Ga
     3.17800d0, -1.02230d0,  0.89119d0, -0.09884d0, & ! Ge
     2.97740d0, -0.81038d0,  0.34861d0,  0.32231d0, & ! As
     2.78340d0, -0.61110d0, -0.14731d0,  0.69837d0, & ! Se
     2.60610d0, -0.43308d0, -0.57381d0,  1.00950d0, & ! Br
     2.44280d0, -0.27244d0, -0.95570d0,  1.27070d0, & ! Kr
     2.30990d0, -0.14328d0, -1.22600d0,  1.45320d0, & ! Rb
     2.21070d0, -0.04770d0, -1.41100d0,  1.55410d0, & ! Sr
     2.14220d0,  0.01935d0, -1.52240d0,  1.59630d0, & ! Y
     2.12690d0,  0.08618d0, -1.49190d0,  1.51820d0, & ! Zr
     2.12120d0,  0.05381d0, -1.50070d0,  1.50150d0, & ! Nb
     2.18870d0, -0.00655d0, -1.25340d0,  1.24010d0, & ! Mo
     2.25730d0, -0.05737d0, -1.07450d0,  1.06630d0, & ! Tc
     2.37300d0, -0.15040d0, -0.77694d0,  0.79060d0, & ! Ru
     2.50990d0, -0.25906d0, -0.44719d0,  0.49443d0, & ! Rh
     2.67520d0, -0.39137d0, -0.05894d0,  0.15404d0, & ! Pd
     2.88690d0, -0.56119d0,  0.42189d0, -0.25659d0, & ! Ag
     3.08430d0, -0.71450d0,  0.84482d0, -0.60990d0, & ! Cd
     3.31400d0, -0.89697d0,  1.35030d0, -1.03910d0, & ! In
     3.49840d0, -1.02990d0,  1.68990d0, -1.29860d0, & ! Sn
     3.70410d0, -1.18270d0,  2.08920d0, -1.61640d0, & ! Sb
     3.88240d0, -1.30980d0,  2.41170d0, -1.86420d0, & ! Te
     4.08010d0, -1.45080d0,  2.76730d0, -2.13920d0, & ! I
     4.24610d0, -1.56330d0,  3.04200d0, -2.34290d0, & ! Xe
     4.38910d0, -1.65420d0,  3.25450d0, -2.49220d0, & ! Cs
     4.51070d0, -1.72570d0,  3.41320d0, -2.59590d0, & ! Ba
     4.60250d0, -1.77070d0,  3.49970d0, -2.64050d0, & ! La
     4.69060d0, -1.81790d0,  3.60280d0, -2.70670d0, & ! Ce
     4.72150d0, -1.81390d0,  3.56480d0, -2.65180d0, & ! Pr
     4.75090d0, -1.80800d0,  3.51970d0, -2.59010d0, & ! Nd
     4.74070d0, -1.76600d0,  3.37430d0, -2.44210d0, & ! Pm
     4.71700d0, -1.71410d0,  3.20800d0, -2.28170d0, & ! Sm
     4.66940d0, -1.64140d0,  2.98580d0, -2.07460d0, & ! Eu
     4.61010d0, -1.55750d0,  2.73190d0, -1.84040d0, & ! Gd
     4.52550d0, -1.45520d0,  2.43770d0, -1.57950d0, & ! Tb
     4.45230d0, -1.36440d0,  2.17540d0, -1.34550d0, & ! Dy
     4.37660d0, -1.27460d0,  1.92540d0, -1.13090d0, & ! Ho
     4.29460d0, -1.18170d0,  1.67060d0, -0.91467d0, & ! Er
     4.21330d0, -1.09060d0,  1.42390d0, -0.70804d0, & ! Tm
     4.13430d0, -1.00310d0,  1.18810d0, -0.51120d0, & ! Yb
     4.04230d0, -0.90518d0,  0.92889d0, -0.29820d0, & ! Lu
     3.95160d0, -0.80978d0,  0.68951d0, -0.09620d0, & ! Hf
     3.85000d0, -0.70599d0,  0.41103d0,  0.11842d0, & ! Ta
     3.76510d0, -0.61807d0,  0.18568d0,  0.29787d0, & ! W
     3.67600d0, -0.52688d0, -0.04706d0,  0.48180d0, & ! Re
     3.60530d0, -0.45420d0, -0.22529d0,  0.61700d0, & ! Os
     3.53130d0, -0.37856d0, -0.41174d0,  0.75967d0, & ! Ir
     3.47070d0, -0.31534d0, -0.56487d0,  0.87492d0, & ! Pt
     3.41630d0, -0.25987d0, -0.69030d0,  0.96224d0, & ! Au
     3.37350d0, -0.21428d0, -0.79013d0,  1.02850d0, & ! Hg
     3.34590d0, -0.18322d0, -0.84911d0,  1.05970d0, & ! Tl
     3.32330d0, -0.15596d0, -0.89878d0,  1.08380d0, & ! Pb
     3.31880d0, -0.14554d0, -0.90198d0,  1.06850d0, & ! Bi
     3.32030d0, -0.13999d0, -0.89333d0,  1.04380d0, & ! Po
     3.34250d0, -0.15317d0, -0.83350d0,  0.97641d0, & ! At
     3.37780d0, -0.17800d0, -0.74320d0,  0.88510d0, & ! Rn
     3.41990d0, -0.20823d0, -0.64000d0,  0.78354d0, & ! Fr
     3.47530d0, -0.25005d0, -0.50660d0,  0.65836d0, & ! Ra
     3.49020d0, -0.25109d0, -0.49651d0,  0.64340d0, & ! Ac
     3.61060d0, -0.35409d0, -0.18926d0,  0.36849d0, & ! Th
     3.68630d0, -0.41329d0, -0.01192d0,  0.20878d0, & ! Pa
     3.76650d0, -0.47542d0,  0.16850d0,  0.05060d0, & ! U
     3.82870d0, -0.51955d0,  0.29804d0, -0.06566d0, & ! Np
     3.88970d0, -0.56296d0,  0.42597d0, -0.18080d0  & ! Pu
     /),shape(c2scatt))

  ! Contrast-ful colors and point-types for gnuplot
  character*7, parameter :: gplt_rgb(13) = (/&
       "#000000",&
       "#0000FF",&
       "#008B00",&
       "#FF0000",&
       "#FF00FF",&
       "#643700",&
       "#787878",&
       "#00FFFF",&
       "#FFB45A",&
       "#7D26CD",&
       "#CD9B9B",&
       "#CD6D0C",&
       "#00B98B"/) !< 13 high-contrast color strings for gnuplot
  integer, parameter :: gplt_sym(13) = (/&
       4,&
       6,&
       8,&
       10,&
       12,&
       14,&
       3,&
       5,&
       7,&
       9,&
       11,&
       13,&
       1/) !< 13 high-contrast point types for gnuplot

  ! cubic spherical harmonics coefficients for wien2k
  integer, parameter, private :: lmax2 = 14 !< maximum l for LM lattice harmonics expansion of charte inside muffins.
  real*8 :: c_kub(0:lmax2,0:lmax2) !< cubic harmonic coefficients

  ! machine constants
  real*8 :: d1mach_(5)

  character*(20), parameter :: constlist(3) = (/character(len=20)::&
       "pi","e","eps"/)
  character*(20), parameter :: funclist(23) = (/character(len=20)::&
       "abs","exp","sqrt","floor","ceil","ceiling","round","log",&
       "log10","sin","asin","cos","acos","tan","atan","atan2","sinh","cosh",&
       "erf","erfc","min","max","xc"/)

  ! directory separator
#ifdef _WIN32
  character*1, parameter :: dirsep = "\" ! "
#else
  character*1, parameter :: dirsep = "/"
#endif
  character*1, parameter :: ampersand = "&" ! ifort doesn't like the "&&" string

  ! the symbols list
  type(hash) :: vh ! variables hash

  ! module procedure interfaces
  interface
     !xx! proc submodule
     module subroutine param_init()
     end subroutine param_init
  end interface

end module param
