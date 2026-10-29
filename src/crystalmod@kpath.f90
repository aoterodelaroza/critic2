! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Special k-points and band-structure paths in the Brillouin zone, by
! the method of:
!
!   Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, "Band
!   structure diagram paths based on crystallography",
!   Comput. Mater. Sci. 128 (2017) 140-184 (HPKOT).
!   doi:10.1016/j.commatsci.2016.10.015
!
! The tables of points and paths and the choice of the extended
! Bravais lattice are ported from SeeK-path
! (https://github.com/giovannipizzi/seekpath), distributed under the
! MIT license:
!
!   Copyright (c), 2016-2023, Giovanni Pizzi, PAUL SCHERRER INSTITUT
!   (Laboratory for Materials Simulations), ECOLE POLYTECHNIQUE
!   FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS)
!   and National Centre for Computational Design and Discovery of
!   Novel Materials (NCCR MARVEL)), Switzerland. All rights reserved.
!
!   Permission is hereby granted, free of charge, to any person
!   obtaining a copy of this software and associated documentation
!   files (the "Software"), to deal in the Software without
!   restriction, including without limitation the rights to use,
!   copy, modify, merge, publish, distribute, sublicense, and/or sell
!   copies of the Software, and to permit persons to whom the
!   Software is furnished to do so, subject to the following
!   conditions:
!
!   The above copyright notice and this permission notice shall be
!   included in all copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!   OTHER DEALINGS IN THE SOFTWARE.
submodule (crystalmod) kpath
  implicit none

  ! tolerance of the Niggli reduction of the triclinic reciprocal cell
  real*8, parameter :: niggli_eps = 1d-5

  ! coordinates of the points below this are round-off zeros
  real*8, parameter :: xk_zero = 1d-12

contains

  !> Special k-points and suggested band path of crystal c (HPKOT).
  !> Returns the extended Bravais lattice symbol (ext, e.g. cF2), the
  !> space group number (spgnum), the labels of the points in the
  !> HPKOT notation (lbl, e.g. GAMMA, X, A_0) and their coordinates
  !> (xk) in the reciprocal cell of c (crystallographic, no 2*pi) and
  !> in the reciprocal primitive cell of HPKOT (xkp), the vectors of
  !> the latter in the reciprocal crystallographic coordinates of c
  !> (mkp, in columns, so xk = mkp * xkp), and the segments of the path
  !> (ipath, pairs of indices into lbl). The
  !> points refer to the primitive cell of the crystal found by
  !> spglib; if the primitive cell of c (its cell and its centering
  !> vectors) is a supercell of it, issuper is true and the labels do
  !> not correspond to the Brillouin zone of c. If the points cannot be
  !> calculated, errmsg is not empty on output.
  module subroutine special_kpoints(c,ext,spgnum,lbl,xk,xkp,mkp,ipath,issuper,errmsg)
    use spglib, only: SpglibDataset, spg_niggli_reduce
    use tools_math, only: matinv, det3
    use param, only: eye, pi, bohrtoa
    class(crystal), intent(in) :: c
    character(len=3), intent(out) :: ext
    integer, intent(out) :: spgnum
    character(len=8), allocatable, intent(out) :: lbl(:)
    real*8, allocatable, intent(out) :: xk(:,:)
    real*8, allocatable, intent(out) :: xkp(:,:)
    real*8, intent(out) :: mkp(3,3)
    integer, allocatable, intent(out) :: ipath(:,:)
    logical, intent(out) :: issuper
    character(len=:), allocatable, intent(out) :: errmsg

    type(SpglibDataset) :: spg
    character*2 :: brav
    real*8 :: ptr(3,3), rconv(3,3), rprim(3,3), pmat(3,3), kprim(3,3)
    real*8 :: a, b, cc, cosbeta, rlen(3), cosr(3), klen(3), cosk(3), cond(3)
    real*8 :: k(3,3), m2(3,3), m3(3,3)
    integer :: ier

    ! reciprocal lengths in the units of SeeK-path (2*pi/angstrom), for
    ! the Niggli reduction (its tolerance is absolute)
    real*8, parameter :: kunit_seekpath = 2d0 * pi / bohrtoa

    ext = ""
    spgnum = 0
    mkp = 0d0
    issuper = .false.
    errmsg = ""
    if (c%ismolecule) then
       errmsg = "special k-points require a crystal"
       return
    end if

    ! the space group and the transformation to its conventional cell;
    ! always recalculated because c%spg may be a subgroup (e.g. after
    ! an edit with atom indices as types)
    call c%spglib_wrap(spg,.false.,errmsg)
    if (len_trim(errmsg) > 0) return
    spgnum = spg%spacegroup_number
    if (spgnum < 1 .or. spgnum > 230) then
       errmsg = "could not determine the space group"
       return
    end if

    ! The conventional cell, in the Cartesian frame of c (spglib's
    ! standardized cell without the idealization and the rotation):
    ! (a_s b_s c_s) = (a b c) P^-1. The lattice vectors are the columns
    ! of rconv.
    ptr = transpose(spg%transformation_matrix)
    call matinv(ptr,3,ier)
    if (ier /= 0) then
       errmsg = "singular transformation to the conventional cell"
       return
    end if
    rconv = matmul(c%m_x2c,ptr)
    call colpar(rconv,rlen,cosr)
    a = rlen(1)
    b = rlen(2)
    cc = rlen(3)
    cosbeta = cosr(2)

    ! the Bravais lattice and the extended Bravais lattice (HPKOT)
    brav = family(spgnum) // hpkot_centering(spgnum)
    select case (brav)
    case ("cP")
       ext = merge("cP1","cP2",spgnum <= 206)
    case ("cF")
       ext = merge("cF1","cF2",spgnum <= 206)
    case ("cI")
       ext = "cI1"
    case ("tP")
       ext = "tP1"
    case ("tI")
       ext = merge("tI1","tI2",cc <= a)
    case ("oP")
       ext = "oP1"
    case ("oF")
       if (1d0/a**2 > 1d0/b**2 + 1d0/cc**2) then
          ext = "oF1"
       elseif (1d0/cc**2 > 1d0/a**2 + 1d0/b**2) then
          ext = "oF2"
       else
          ext = "oF3"
       end if
    case ("oI")
       ! by the longest vector: c (oI1), a (oI2), b (oI3); ties as
       ! in SeeK-path (b, then a, then c)
       if (b >= a .and. b >= cc) then
          ext = "oI3"
       elseif (a >= cc) then
          ext = "oI2"
       else
          ext = "oI1"
       end if
    case ("oC")
       ext = merge("oC1","oC2",a <= b)
    case ("oA")
       ext = merge("oA1","oA2",b <= cc)
    case ("hP")
       if (any(spgnum == (/143,144,145,146,147,148,149,151,153,157,159,160,161,162,163/))) then
          ext = "hP1"
       else
          ext = "hP2"
       end if
    case ("hR")
       ext = merge("hR1","hR2",sqrt(3d0) * a <= sqrt(2d0) * cc)
    case ("mP")
       ext = "mP1"
    case ("mC")
       if (b < a * sqrt(1d0 - cosbeta**2)) then
          ext = "mC1"
       elseif (-a * cosbeta / cc + a**2 * (1d0 - cosbeta**2) / b**2 <= 1d0) then
          ext = "mC2"
       else
          ext = "mC3"
       end if
    case ("aP")
       ! Niggli-reduce the reciprocal cell (k, vectors in columns), then
       ! permute (M2) and flip (M3) its vectors so it is all-obtuse
       ! (aP2) or all-acute (aP3). M2 and M3 are their own inverse
       ! transposes, so they transform the direct cell in the same way.
       k = recip(rconv)
       k = transpose(k) * kunit_seekpath
       ier = spg_niggli_reduce(k,niggli_eps)
       if (ier == 0) then
          errmsg = "the Niggli reduction of the reciprocal cell failed"
          return
       end if
       k = transpose(k) / kunit_seekpath
       call colpar(k,klen,cosk)
       cond = abs(cosk * (/klen(2)*klen(3),klen(1)*klen(3),klen(1)*klen(2)/))
       select case (minloc(cond,1))
       case (1)
          m2 = reshape((/0d0,1d0,0d0, 0d0,0d0,1d0, 1d0,0d0,0d0/),(/3,3/))
       case (2)
          m2 = reshape((/0d0,0d0,1d0, 1d0,0d0,0d0, 0d0,1d0,0d0/),(/3,3/))
       case default
          m2 = eye
       end select
       k = matmul(k,m2)
       call colpar(k,klen,cosk)
       if ((cosk(1) > 0d0 .eqv. cosk(2) > 0d0) .and. (cosk(1) > 0d0 .eqv. cosk(3) > 0d0)) then
          m3 = eye ! 1a, 1b
       elseif (cosk(2) > 0d0 .eqv. cosk(3) > 0d0) then
          m3 = diag(1d0,-1d0,-1d0) ! 2a, 2b
       elseif (cosk(1) > 0d0 .eqv. cosk(3) > 0d0) then
          m3 = diag(-1d0,1d0,-1d0) ! 3a, 3b
       else
          m3 = diag(-1d0,-1d0,1d0) ! 4a, 4b
       end if
       k = matmul(k,m3)
       call colpar(k,klen,cosk)
       rconv = recip(k)
       if (all(cosk <= 0d0)) then
          ext = "aP2"
       elseif (all(cosk >= 0d0)) then
          ext = "aP3"
       else
          errmsg = "the triclinic reciprocal cell is neither all-obtuse nor all-acute"
          return
       end if
    case default
       errmsg = "unknown Bravais lattice " // brav
       return
    end select

    ! the primitive cell of HPKOT: (a_P b_P c_P) = (a b c) P
    pmat = pmatrix(brav)
    rprim = matmul(rconv,pmat)

    ! the points, from the reciprocal primitive cell of HPKOT to the
    ! reciprocal cell of c through Cartesian coordinates (x_c = A^T k,
    ! A = lattice vectors of c in columns)
    call hpkot_table(ext,a,b,cc,cosbeta,lbl,xkp,ipath)
    if (.not.allocated(lbl)) then
       errmsg = "no table for the extended Bravais lattice " // ext
       return
    end if
    kprim = recip(rprim)
    mkp = matmul(transpose(c%m_x2c),kprim)
    xk = matmul(mkp,xkp)
    where (abs(xk) < xk_zero) xk = 0d0
    where (abs(mkp) < xk_zero) mkp = 0d0

    ! does the primitive cell of c (its cell and centering vectors) miss
    ! lattice translations that spglib found?
    issuper = (nint(abs(c%omega / max(c%ncv,1) / det3(rprim))) > 1)

  contains
    ! Crystal family of space group n (a, m, o, t, h, c)
    function family(n)
      integer, intent(in) :: n
      character*1 :: family

      if (n <= 2) then
         family = "a"
      elseif (n <= 15) then
         family = "m"
      elseif (n <= 74) then
         family = "o"
      elseif (n <= 142) then
         family = "t"
      elseif (n <= 194) then
         family = "h"
      else
         family = "c"
      end if

    end function family

    ! Diagonal 3x3 matrix
    function diag(x1,x2,x3) result(m)
      real*8, intent(in) :: x1, x2, x3
      real*8 :: m(3,3)

      m = 0d0
      m(1,1) = x1
      m(2,2) = x2
      m(3,3) = x3

    end function diag

    ! The matrix P from the conventional to the primitive cell of
    ! Bravais lattice bl (Table 3 in HPKOT)
    function pmatrix(bl) result(p)
      character*2, intent(in) :: bl
      real*8 :: p(3,3)

      select case (bl(2:2))
      case ("F")
         p = reshape((/0d0,1d0,1d0, 1d0,0d0,1d0, 1d0,1d0,0d0/),(/3,3/)) / 2d0
      case ("I")
         p = reshape((/-1d0,1d0,1d0, 1d0,-1d0,1d0, 1d0,1d0,-1d0/),(/3,3/)) / 2d0
      case ("R")
         p = reshape((/2d0,1d0,1d0, -1d0,1d0,1d0, -1d0,-2d0,1d0/),(/3,3/)) / 3d0
      case ("C")
         if (bl(1:1) == "m") then
            p = reshape((/1d0,1d0,0d0, -1d0,1d0,0d0, 0d0,0d0,2d0/),(/3,3/)) / 2d0
         else
            p = reshape((/1d0,-1d0,0d0, 1d0,1d0,0d0, 0d0,0d0,2d0/),(/3,3/)) / 2d0
         end if
      case ("A")
         p = reshape((/0d0,1d0,-1d0, 0d0,1d0,1d0, 2d0,0d0,0d0/),(/3,3/)) / 2d0
      case default
         p = eye
      end select

    end function pmatrix
  end subroutine special_kpoints

  !xx! private procedures

  !> Reciprocal cell (no 2*pi, vectors in columns) of the cell r (vectors
  !> in columns)
  function recip(r) result(k)
    use tools_math, only: matinv
    real*8, intent(in) :: r(3,3)
    real*8 :: k(3,3)

    k = transpose(r)
    call matinv(k,3)

  end function recip

  !> Lengths (rlen) and cosines of the angles (cosang: alpha, beta,
  !> gamma) of the cell r (vectors in columns)
  subroutine colpar(r,rlen,cosang)
    real*8, intent(in) :: r(3,3)
    real*8, intent(out) :: rlen(3), cosang(3)

    integer :: i

    do i = 1, 3
       rlen(i) = norm2(r(:,i))
    end do
    cosang(1) = dot_product(r(:,2),r(:,3)) / rlen(2) / rlen(3)
    cosang(2) = dot_product(r(:,1),r(:,3)) / rlen(1) / rlen(3)
    cosang(3) = dot_product(r(:,1),r(:,2)) / rlen(1) / rlen(2)

  end subroutine colpar

  !xx! begin generated by tools/hpkot-tables.py (do not edit)

  !> Centering letter of the standard setting of each space group
  !> (from SeeK-path's spg_db.py).
  pure function hpkot_centering(n) result(cen)
    integer, intent(in) :: n
    character*1 :: cen

    character(len=230), parameter :: table = &
       "PPPPCPPCCPPCPPCPPPPCCFIIPPPPPPPPPPCCCAAAAFFIIIPPPPPPPPPPPPPP" //&
       "PPCCCCCCFFIIIIPPPPIIPIPPPPIIPPPPPPPPIIPPPPPPPPIIIIPPPPPPPPII" //&
       "IIPPPPPPPPPPPPPPPPIIIIPPPRPRPPPPPPRPPPPRRPPPPRRPPPPPPPPPPPPP" //&
       "PPPPPPPPPPPPPPPFIPIPPFFIPIPPFFIPPIPFIPFIPPPPFFFFII"

    cen = table(n:n)

  end function hpkot_centering

  !> Special k-points (labels lbl and coordinates xk in the reciprocal
  !> primitive cell of HPKOT) and suggested band path (pairs of indices
  !> into the point list) of the extended Bravais lattice ext. The
  !> point coordinates depend on the conventional cell parameters a, b,
  !> c and cos(beta).
  pure subroutine hpkot_table(ext,a,b,c,cosbeta,lbl,xk,ipath)
    character*(*), intent(in) :: ext
    real*8, intent(in) :: a, b, c, cosbeta
    character(len=8), allocatable, intent(out) :: lbl(:)
    real*8, allocatable, intent(out) :: xk(:,:)
    integer, allocatable, intent(out) :: ipath(:,:)

    real*8 :: sinbeta, D, E, F, H, J, K, M, N, P, R, S, U, W, X, Y, Z

    sinbeta = sqrt(1d0 - cosbeta**2)
    select case (ext)
    case ("aP2")
       allocate(lbl(8),xk(3,8),ipath(2,7))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Z"
       xk(:,2) = (/0d0,0d0,1d0/2d0/)
       lbl(3) = "Y"
       xk(:,3) = (/0d0,1d0/2d0,0d0/)
       lbl(4) = "X"
       xk(:,4) = (/1d0/2d0,0d0,0d0/)
       lbl(5) = "V"
       xk(:,5) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(6) = "U"
       xk(:,6) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(7) = "T"
       xk(:,7) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(8) = "R"
       xk(:,8) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       ipath(:,1) = (/1,4/) ! GAMMA-X
       ipath(:,2) = (/3,1/) ! Y-GAMMA
       ipath(:,3) = (/1,2/) ! GAMMA-Z
       ipath(:,4) = (/8,1/) ! R-GAMMA
       ipath(:,5) = (/1,7/) ! GAMMA-T
       ipath(:,6) = (/6,1/) ! U-GAMMA
       ipath(:,7) = (/1,5/) ! GAMMA-V
    case ("aP3")
       allocate(lbl(9),xk(3,9),ipath(2,7))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Z"
       xk(:,2) = (/0d0,0d0,1d0/2d0/)
       lbl(3) = "Y"
       xk(:,3) = (/0d0,1d0/2d0,0d0/)
       lbl(4) = "Y_2"
       xk(:,4) = (/0d0,-1d0/2d0,0d0/)
       lbl(5) = "X"
       xk(:,5) = (/1d0/2d0,0d0,0d0/)
       lbl(6) = "V_2"
       xk(:,6) = (/1d0/2d0,-1d0/2d0,0d0/)
       lbl(7) = "U_2"
       xk(:,7) = (/-1d0/2d0,0d0,1d0/2d0/)
       lbl(8) = "T_2"
       xk(:,8) = (/0d0,-1d0/2d0,1d0/2d0/)
       lbl(9) = "R_2"
       xk(:,9) = (/-1d0/2d0,-1d0/2d0,1d0/2d0/)
       ipath(:,1) = (/1,5/) ! GAMMA-X
       ipath(:,2) = (/3,1/) ! Y-GAMMA
       ipath(:,3) = (/1,2/) ! GAMMA-Z
       ipath(:,4) = (/9,1/) ! R_2-GAMMA
       ipath(:,5) = (/1,8/) ! GAMMA-T_2
       ipath(:,6) = (/7,1/) ! U_2-GAMMA
       ipath(:,7) = (/1,6/) ! GAMMA-V_2
    case ("cF1")
       allocate(lbl(7),xk(3,7),ipath(2,7))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "X"
       xk(:,2) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(3) = "L"
       xk(:,3) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(4) = "W"
       xk(:,4) = (/1d0/2d0,1d0/4d0,3d0/4d0/)
       lbl(5) = "W_2"
       xk(:,5) = (/3d0/4d0,1d0/4d0,1d0/2d0/)
       lbl(6) = "K"
       xk(:,6) = (/3d0/8d0,3d0/8d0,3d0/4d0/)
       lbl(7) = "U"
       xk(:,7) = (/5d0/8d0,1d0/4d0,5d0/8d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-X
       ipath(:,2) = (/2,7/) ! X-U
       ipath(:,3) = (/6,1/) ! K-GAMMA
       ipath(:,4) = (/1,3/) ! GAMMA-L
       ipath(:,5) = (/3,4/) ! L-W
       ipath(:,6) = (/4,2/) ! W-X
       ipath(:,7) = (/2,5/) ! X-W_2
    case ("cF2")
       allocate(lbl(7),xk(3,7),ipath(2,6))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "X"
       xk(:,2) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(3) = "L"
       xk(:,3) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(4) = "W"
       xk(:,4) = (/1d0/2d0,1d0/4d0,3d0/4d0/)
       lbl(5) = "W_2"
       xk(:,5) = (/3d0/4d0,1d0/4d0,1d0/2d0/)
       lbl(6) = "K"
       xk(:,6) = (/3d0/8d0,3d0/8d0,3d0/4d0/)
       lbl(7) = "U"
       xk(:,7) = (/5d0/8d0,1d0/4d0,5d0/8d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-X
       ipath(:,2) = (/2,7/) ! X-U
       ipath(:,3) = (/6,1/) ! K-GAMMA
       ipath(:,4) = (/1,3/) ! GAMMA-L
       ipath(:,5) = (/3,4/) ! L-W
       ipath(:,6) = (/4,2/) ! W-X
    case ("cI1")
       allocate(lbl(4),xk(3,4),ipath(2,6))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "H"
       xk(:,2) = (/1d0/2d0,-1d0/2d0,1d0/2d0/)
       lbl(3) = "P"
       xk(:,3) = (/1d0/4d0,1d0/4d0,1d0/4d0/)
       lbl(4) = "N"
       xk(:,4) = (/0d0,0d0,1d0/2d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-H
       ipath(:,2) = (/2,4/) ! H-N
       ipath(:,3) = (/4,1/) ! N-GAMMA
       ipath(:,4) = (/1,3/) ! GAMMA-P
       ipath(:,5) = (/3,2/) ! P-H
       ipath(:,6) = (/3,4/) ! P-N
    case ("cP1")
       allocate(lbl(5),xk(3,5),ipath(2,7))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "R"
       xk(:,2) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "M"
       xk(:,3) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(4) = "X"
       xk(:,4) = (/0d0,1d0/2d0,0d0/)
       lbl(5) = "X_1"
       xk(:,5) = (/1d0/2d0,0d0,0d0/)
       ipath(:,1) = (/1,4/) ! GAMMA-X
       ipath(:,2) = (/4,3/) ! X-M
       ipath(:,3) = (/3,1/) ! M-GAMMA
       ipath(:,4) = (/1,2/) ! GAMMA-R
       ipath(:,5) = (/2,4/) ! R-X
       ipath(:,6) = (/2,3/) ! R-M
       ipath(:,7) = (/3,5/) ! M-X_1
    case ("cP2")
       allocate(lbl(5),xk(3,5),ipath(2,6))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "R"
       xk(:,2) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "M"
       xk(:,3) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(4) = "X"
       xk(:,4) = (/0d0,1d0/2d0,0d0/)
       lbl(5) = "X_1"
       xk(:,5) = (/1d0/2d0,0d0,0d0/)
       ipath(:,1) = (/1,4/) ! GAMMA-X
       ipath(:,2) = (/4,3/) ! X-M
       ipath(:,3) = (/3,1/) ! M-GAMMA
       ipath(:,4) = (/1,2/) ! GAMMA-R
       ipath(:,5) = (/2,4/) ! R-X
       ipath(:,6) = (/2,3/) ! R-M
    case ("hP1")
       allocate(lbl(7),xk(3,7),ipath(2,10))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "A"
       xk(:,2) = (/0d0,0d0,1d0/2d0/)
       lbl(3) = "K"
       xk(:,3) = (/1d0/3d0,1d0/3d0,0d0/)
       lbl(4) = "H"
       xk(:,4) = (/1d0/3d0,1d0/3d0,1d0/2d0/)
       lbl(5) = "H_2"
       xk(:,5) = (/1d0/3d0,1d0/3d0,-1d0/2d0/)
       lbl(6) = "M"
       xk(:,6) = (/1d0/2d0,0d0,0d0/)
       lbl(7) = "L"
       xk(:,7) = (/1d0/2d0,0d0,1d0/2d0/)
       ipath(:,1) = (/1,6/) ! GAMMA-M
       ipath(:,2) = (/6,3/) ! M-K
       ipath(:,3) = (/3,1/) ! K-GAMMA
       ipath(:,4) = (/1,2/) ! GAMMA-A
       ipath(:,5) = (/2,7/) ! A-L
       ipath(:,6) = (/7,4/) ! L-H
       ipath(:,7) = (/4,2/) ! H-A
       ipath(:,8) = (/7,6/) ! L-M
       ipath(:,9) = (/4,3/) ! H-K
       ipath(:,10) = (/3,5/) ! K-H_2
    case ("hP2")
       allocate(lbl(7),xk(3,7),ipath(2,9))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "A"
       xk(:,2) = (/0d0,0d0,1d0/2d0/)
       lbl(3) = "K"
       xk(:,3) = (/1d0/3d0,1d0/3d0,0d0/)
       lbl(4) = "H"
       xk(:,4) = (/1d0/3d0,1d0/3d0,1d0/2d0/)
       lbl(5) = "H_2"
       xk(:,5) = (/1d0/3d0,1d0/3d0,-1d0/2d0/)
       lbl(6) = "M"
       xk(:,6) = (/1d0/2d0,0d0,0d0/)
       lbl(7) = "L"
       xk(:,7) = (/1d0/2d0,0d0,1d0/2d0/)
       ipath(:,1) = (/1,6/) ! GAMMA-M
       ipath(:,2) = (/6,3/) ! M-K
       ipath(:,3) = (/3,1/) ! K-GAMMA
       ipath(:,4) = (/1,2/) ! GAMMA-A
       ipath(:,5) = (/2,7/) ! A-L
       ipath(:,6) = (/7,4/) ! L-H
       ipath(:,7) = (/4,2/) ! H-A
       ipath(:,8) = (/7,6/) ! L-M
       ipath(:,9) = (/4,3/) ! H-K
    case ("hR1")
       D = a*a/4d0/c/c
       Y = 5d0/6d0-2d0*D
       N = 1d0/3d0+D
       allocate(lbl(20),xk(3,20),ipath(2,7))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "T"
       xk(:,2) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "L"
       xk(:,3) = (/1d0/2d0,0d0,0d0/)
       lbl(4) = "L_2"
       xk(:,4) = (/0d0,-1d0/2d0,0d0/)
       lbl(5) = "L_4"
       xk(:,5) = (/0d0,0d0,-1d0/2d0/)
       lbl(6) = "F"
       xk(:,6) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(7) = "F_2"
       xk(:,7) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(8) = "S_0"
       xk(:,8) = (/N,-N,0d0/)
       lbl(9) = "S_2"
       xk(:,9) = (/1d0-N,0d0,N/)
       lbl(10) = "S_4"
       xk(:,10) = (/N,0d0,-N/)
       lbl(11) = "S_6"
       xk(:,11) = (/1d0-N,N,0d0/)
       lbl(12) = "H_0"
       xk(:,12) = (/1d0/2d0,-1d0+Y,1d0-Y/)
       lbl(13) = "H_2"
       xk(:,13) = (/Y,1d0-Y,1d0/2d0/)
       lbl(14) = "H_4"
       xk(:,14) = (/Y,1d0/2d0,1d0-Y/)
       lbl(15) = "H_6"
       xk(:,15) = (/1d0/2d0,1d0-Y,-1d0+Y/)
       lbl(16) = "M_0"
       xk(:,16) = (/N,-1d0+Y,N/)
       lbl(17) = "M_2"
       xk(:,17) = (/1d0-N,1d0-Y,1d0-N/)
       lbl(18) = "M_4"
       xk(:,18) = (/Y,N,N/)
       lbl(19) = "M_6"
       xk(:,19) = (/1d0-N,1d0-N,1d0-Y/)
       lbl(20) = "M_8"
       xk(:,20) = (/N,N,-1d0+Y/)
       ipath(:,1) = (/1,2/) ! GAMMA-T
       ipath(:,2) = (/2,13/) ! T-H_2
       ipath(:,3) = (/12,3/) ! H_0-L
       ipath(:,4) = (/3,1/) ! L-GAMMA
       ipath(:,5) = (/1,8/) ! GAMMA-S_0
       ipath(:,6) = (/9,6/) ! S_2-F
       ipath(:,7) = (/6,1/) ! F-GAMMA
    case ("hR2")
       Z = 1d0/6d0-c*c/9d0/a/a
       H = 1d0/2d0-2d0*Z
       N = 1d0/2d0+Z
       allocate(lbl(9),xk(3,9),ipath(2,5))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "T"
       xk(:,2) = (/1d0/2d0,-1d0/2d0,1d0/2d0/)
       lbl(3) = "P_0"
       xk(:,3) = (/H,-1d0+H,H/)
       lbl(4) = "P_2"
       xk(:,4) = (/H,H,H/)
       lbl(5) = "R_0"
       xk(:,5) = (/1d0-H,-H,-H/)
       lbl(6) = "M"
       xk(:,6) = (/1d0-N,-N,1d0-N/)
       lbl(7) = "M_2"
       xk(:,7) = (/N,-1d0+N,-1d0+N/)
       lbl(8) = "L"
       xk(:,8) = (/1d0/2d0,0d0,0d0/)
       lbl(9) = "F"
       xk(:,9) = (/1d0/2d0,-1d0/2d0,0d0/)
       ipath(:,1) = (/1,8/) ! GAMMA-L
       ipath(:,2) = (/8,2/) ! L-T
       ipath(:,3) = (/2,3/) ! T-P_0
       ipath(:,4) = (/4,1/) ! P_2-GAMMA
       ipath(:,5) = (/1,9/) ! GAMMA-F
    case ("mC1")
       Z = (2d0+a/c*cosbeta)/4d0/sinbeta/sinbeta
       H = 1d0/2d0-2d0*Z*c*cosbeta/a
       S = 3d0/4d0-b*b/4d0/a/a/sinbeta/sinbeta
       P = S-(3d0/4d0-S)*a*cosbeta/c
       allocate(lbl(16),xk(3,16),ipath(2,9))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Y_2"
       xk(:,2) = (/-1d0/2d0,1d0/2d0,0d0/)
       lbl(3) = "Y_4"
       xk(:,3) = (/1d0/2d0,-1d0/2d0,0d0/)
       lbl(4) = "A"
       xk(:,4) = (/0d0,0d0,1d0/2d0/)
       lbl(5) = "M_2"
       xk(:,5) = (/-1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(6) = "V"
       xk(:,6) = (/1d0/2d0,0d0,0d0/)
       lbl(7) = "V_2"
       xk(:,7) = (/0d0,1d0/2d0,0d0/)
       lbl(8) = "L_2"
       xk(:,8) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(9) = "C"
       xk(:,9) = (/1d0-S,1d0-S,0d0/)
       lbl(10) = "C_2"
       xk(:,10) = (/-1d0+S,S,0d0/)
       lbl(11) = "C_4"
       xk(:,11) = (/S,-1d0+S,0d0/)
       lbl(12) = "D"
       xk(:,12) = (/-1d0+P,P,1d0/2d0/)
       lbl(13) = "D_2"
       xk(:,13) = (/1d0-P,1d0-P,1d0/2d0/)
       lbl(14) = "E"
       xk(:,14) = (/-1d0+Z,1d0-Z,1d0-H/)
       lbl(15) = "E_2"
       xk(:,15) = (/-Z,Z,H/)
       lbl(16) = "E_4"
       xk(:,16) = (/Z,-Z,1d0-H/)
       ipath(:,1) = (/1,9/) ! GAMMA-C
       ipath(:,2) = (/10,2/) ! C_2-Y_2
       ipath(:,3) = (/2,1/) ! Y_2-GAMMA
       ipath(:,4) = (/1,5/) ! GAMMA-M_2
       ipath(:,5) = (/5,12/) ! M_2-D
       ipath(:,6) = (/13,4/) ! D_2-A
       ipath(:,7) = (/4,1/) ! A-GAMMA
       ipath(:,8) = (/8,1/) ! L_2-GAMMA
       ipath(:,9) = (/1,7/) ! GAMMA-V_2
    case ("mC2")
       Z = (a*a/b/b+(1d0+a/c*cosbeta)/sinbeta/sinbeta)/4d0
       M = (1d0+a*a/b/b)/4d0
       D = -a*c*cosbeta/2d0/b/b
       X = 1d0/2d0-2d0*Z*c*cosbeta/a
       P = 1d0+Z-2d0*M
       S = X-2d0*D
       allocate(lbl(16),xk(3,16),ipath(2,6))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Y"
       xk(:,2) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(3) = "A"
       xk(:,3) = (/0d0,0d0,1d0/2d0/)
       lbl(4) = "M"
       xk(:,4) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(5) = "V_2"
       xk(:,5) = (/0d0,1d0/2d0,0d0/)
       lbl(6) = "L_2"
       xk(:,6) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(7) = "F"
       xk(:,7) = (/-1d0+P,1d0-P,1d0-S/)
       lbl(8) = "F_2"
       xk(:,8) = (/1d0-P,P,S/)
       lbl(9) = "F_4"
       xk(:,9) = (/P,1d0-P,1d0-S/)
       lbl(10) = "H"
       xk(:,10) = (/-Z,Z,X/)
       lbl(11) = "H_2"
       xk(:,11) = (/Z,1d0-Z,1d0-X/)
       lbl(12) = "H_4"
       xk(:,12) = (/Z,-Z,1d0-X/)
       lbl(13) = "G"
       xk(:,13) = (/-M,M,D/)
       lbl(14) = "G_2"
       xk(:,14) = (/M,1d0-M,-D/)
       lbl(15) = "G_4"
       xk(:,15) = (/M,-M,-D/)
       lbl(16) = "G_6"
       xk(:,16) = (/1d0-M,M,D/)
       ipath(:,1) = (/1,2/) ! GAMMA-Y
       ipath(:,2) = (/2,4/) ! Y-M
       ipath(:,3) = (/4,3/) ! M-A
       ipath(:,4) = (/3,1/) ! A-GAMMA
       ipath(:,5) = (/6,1/) ! L_2-GAMMA
       ipath(:,6) = (/1,5/) ! GAMMA-V_2
    case ("mC3")
       Z = (a*a/b/b+(1d0+a/c*cosbeta)/sinbeta/sinbeta)/4d0
       R = 1d0-Z*b*b/a/a
       E = 1d0/2d0-2d0*Z*c*cosbeta/a
       F = E/2d0+a*a/4d0/b/b+a*c*cosbeta/2d0/b/b
       U = 2d0*F-Z
       W = c/2d0/a/cosbeta*(1d0-4d0*U+a*a*sinbeta*sinbeta/b/b)
       D = -1d0/4d0+W/2d0-Z*c*cosbeta/a
       allocate(lbl(19),xk(3,19),ipath(2,7))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Y"
       xk(:,2) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(3) = "A"
       xk(:,3) = (/0d0,0d0,1d0/2d0/)
       lbl(4) = "M_2"
       xk(:,4) = (/-1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(5) = "V"
       xk(:,5) = (/1d0/2d0,0d0,0d0/)
       lbl(6) = "V_2"
       xk(:,6) = (/0d0,1d0/2d0,0d0/)
       lbl(7) = "L_2"
       xk(:,7) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(8) = "I"
       xk(:,8) = (/-1d0+R,R,1d0/2d0/)
       lbl(9) = "I_2"
       xk(:,9) = (/1d0-R,1d0-R,1d0/2d0/)
       lbl(10) = "K"
       xk(:,10) = (/-U,U,W/)
       lbl(11) = "K_2"
       xk(:,11) = (/-1d0+U,1d0-U,1d0-W/)
       lbl(12) = "K_4"
       xk(:,12) = (/1d0-U,U,W/)
       lbl(13) = "H"
       xk(:,13) = (/-Z,Z,E/)
       lbl(14) = "H_2"
       xk(:,14) = (/Z,1d0-Z,1d0-E/)
       lbl(15) = "H_4"
       xk(:,15) = (/Z,-Z,1d0-E/)
       lbl(16) = "N"
       xk(:,16) = (/-F,F,D/)
       lbl(17) = "N_2"
       xk(:,17) = (/F,1d0-F,-D/)
       lbl(18) = "N_4"
       xk(:,18) = (/F,-F,-D/)
       lbl(19) = "N_6"
       xk(:,19) = (/1d0-F,F,D/)
       ipath(:,1) = (/1,3/) ! GAMMA-A
       ipath(:,2) = (/3,9/) ! A-I_2
       ipath(:,3) = (/8,4/) ! I-M_2
       ipath(:,4) = (/4,1/) ! M_2-GAMMA
       ipath(:,5) = (/1,2/) ! GAMMA-Y
       ipath(:,6) = (/7,1/) ! L_2-GAMMA
       ipath(:,7) = (/1,6/) ! GAMMA-V_2
    case ("mP1")
       Y = (1d0+a/c*cosbeta)/2d0/sinbeta/sinbeta
       N = 1d0/2d0+Y*c*cosbeta/a
       allocate(lbl(18),xk(3,18),ipath(2,10))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Z"
       xk(:,2) = (/0d0,1d0/2d0,0d0/)
       lbl(3) = "B"
       xk(:,3) = (/0d0,0d0,1d0/2d0/)
       lbl(4) = "B_2"
       xk(:,4) = (/0d0,0d0,-1d0/2d0/)
       lbl(5) = "Y"
       xk(:,5) = (/1d0/2d0,0d0,0d0/)
       lbl(6) = "Y_2"
       xk(:,6) = (/-1d0/2d0,0d0,0d0/)
       lbl(7) = "C"
       xk(:,7) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(8) = "C_2"
       xk(:,8) = (/-1d0/2d0,1d0/2d0,0d0/)
       lbl(9) = "D"
       xk(:,9) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(10) = "D_2"
       xk(:,10) = (/0d0,1d0/2d0,-1d0/2d0/)
       lbl(11) = "A"
       xk(:,11) = (/-1d0/2d0,0d0,1d0/2d0/)
       lbl(12) = "E"
       xk(:,12) = (/-1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(13) = "H"
       xk(:,13) = (/-Y,0d0,1d0-N/)
       lbl(14) = "H_2"
       xk(:,14) = (/-1d0+Y,0d0,N/)
       lbl(15) = "H_4"
       xk(:,15) = (/-Y,0d0,-N/)
       lbl(16) = "M"
       xk(:,16) = (/-Y,1d0/2d0,1d0-N/)
       lbl(17) = "M_2"
       xk(:,17) = (/-1d0+Y,1d0/2d0,N/)
       lbl(18) = "M_4"
       xk(:,18) = (/-Y,1d0/2d0,-N/)
       ipath(:,1) = (/1,2/) ! GAMMA-Z
       ipath(:,2) = (/2,9/) ! Z-D
       ipath(:,3) = (/9,3/) ! D-B
       ipath(:,4) = (/3,1/) ! B-GAMMA
       ipath(:,5) = (/1,11/) ! GAMMA-A
       ipath(:,6) = (/11,12/) ! A-E
       ipath(:,7) = (/12,2/) ! E-Z
       ipath(:,8) = (/2,8/) ! Z-C_2
       ipath(:,9) = (/8,6/) ! C_2-Y_2
       ipath(:,10) = (/6,1/) ! Y_2-GAMMA
    case ("oA1")
       X = (1d0+b*b/c/c)/4d0
       allocate(lbl(10),xk(3,10),ipath(2,11))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Y"
       xk(:,2) = (/-1d0/2d0,1d0/2d0,0d0/)
       lbl(3) = "T"
       xk(:,3) = (/-1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(4) = "Z"
       xk(:,4) = (/0d0,0d0,1d0/2d0/)
       lbl(5) = "S"
       xk(:,5) = (/0d0,1d0/2d0,0d0/)
       lbl(6) = "R"
       xk(:,6) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(7) = "SIGMA_0"
       xk(:,7) = (/X,X,0d0/)
       lbl(8) = "C_0"
       xk(:,8) = (/-X,1d0-X,0d0/)
       lbl(9) = "A_0"
       xk(:,9) = (/X,X,1d0/2d0/)
       lbl(10) = "E_0"
       xk(:,10) = (/-X,1d0-X,1d0/2d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-Y
       ipath(:,2) = (/2,8/) ! Y-C_0
       ipath(:,3) = (/7,1/) ! SIGMA_0-GAMMA
       ipath(:,4) = (/1,4/) ! GAMMA-Z
       ipath(:,5) = (/4,9/) ! Z-A_0
       ipath(:,6) = (/10,3/) ! E_0-T
       ipath(:,7) = (/3,2/) ! T-Y
       ipath(:,8) = (/1,5/) ! GAMMA-S
       ipath(:,9) = (/5,6/) ! S-R
       ipath(:,10) = (/6,4/) ! R-Z
       ipath(:,11) = (/4,3/) ! Z-T
    case ("oA2")
       X = (1d0+c*c/b/b)/4d0
       allocate(lbl(15),xk(3,15),ipath(2,11))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Y"
       xk(:,2) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(3) = "T"
       xk(:,3) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(4) = "T_2"
       xk(:,4) = (/1d0/2d0,1d0/2d0,-1d0/2d0/)
       lbl(5) = "Z"
       xk(:,5) = (/0d0,0d0,1d0/2d0/)
       lbl(6) = "Z_2"
       xk(:,6) = (/0d0,0d0,-1d0/2d0/)
       lbl(7) = "S"
       xk(:,7) = (/0d0,1d0/2d0,0d0/)
       lbl(8) = "R"
       xk(:,8) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(9) = "R_2"
       xk(:,9) = (/0d0,1d0/2d0,-1d0/2d0/)
       lbl(10) = "DELTA_0"
       xk(:,10) = (/-X,X,0d0/)
       lbl(11) = "F_0"
       xk(:,11) = (/X,1d0-X,0d0/)
       lbl(12) = "B_0"
       xk(:,12) = (/-X,X,1d0/2d0/)
       lbl(13) = "B_2"
       xk(:,13) = (/-X,X,-1d0/2d0/)
       lbl(14) = "G_0"
       xk(:,14) = (/X,1d0-X,1d0/2d0/)
       lbl(15) = "G_2"
       xk(:,15) = (/X,1d0-X,-1d0/2d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-Y
       ipath(:,2) = (/2,11/) ! Y-F_0
       ipath(:,3) = (/10,1/) ! DELTA_0-GAMMA
       ipath(:,4) = (/1,5/) ! GAMMA-Z
       ipath(:,5) = (/5,12/) ! Z-B_0
       ipath(:,6) = (/14,3/) ! G_0-T
       ipath(:,7) = (/3,2/) ! T-Y
       ipath(:,8) = (/1,7/) ! GAMMA-S
       ipath(:,9) = (/7,8/) ! S-R
       ipath(:,10) = (/8,5/) ! R-Z
       ipath(:,11) = (/5,3/) ! Z-T
    case ("oC1")
       X = (1d0+a*a/b/b)/4d0
       allocate(lbl(10),xk(3,10),ipath(2,11))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Y"
       xk(:,2) = (/-1d0/2d0,1d0/2d0,0d0/)
       lbl(3) = "T"
       xk(:,3) = (/-1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(4) = "Z"
       xk(:,4) = (/0d0,0d0,1d0/2d0/)
       lbl(5) = "S"
       xk(:,5) = (/0d0,1d0/2d0,0d0/)
       lbl(6) = "R"
       xk(:,6) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(7) = "SIGMA_0"
       xk(:,7) = (/X,X,0d0/)
       lbl(8) = "C_0"
       xk(:,8) = (/-X,1d0-X,0d0/)
       lbl(9) = "A_0"
       xk(:,9) = (/X,X,1d0/2d0/)
       lbl(10) = "E_0"
       xk(:,10) = (/-X,1d0-X,1d0/2d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-Y
       ipath(:,2) = (/2,8/) ! Y-C_0
       ipath(:,3) = (/7,1/) ! SIGMA_0-GAMMA
       ipath(:,4) = (/1,4/) ! GAMMA-Z
       ipath(:,5) = (/4,9/) ! Z-A_0
       ipath(:,6) = (/10,3/) ! E_0-T
       ipath(:,7) = (/3,2/) ! T-Y
       ipath(:,8) = (/1,5/) ! GAMMA-S
       ipath(:,9) = (/5,6/) ! S-R
       ipath(:,10) = (/6,4/) ! R-Z
       ipath(:,11) = (/4,3/) ! Z-T
    case ("oC2")
       X = (1d0+b*b/a/a)/4d0
       allocate(lbl(15),xk(3,15),ipath(2,11))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Y"
       xk(:,2) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(3) = "T"
       xk(:,3) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(4) = "T_2"
       xk(:,4) = (/1d0/2d0,1d0/2d0,-1d0/2d0/)
       lbl(5) = "Z"
       xk(:,5) = (/0d0,0d0,1d0/2d0/)
       lbl(6) = "Z_2"
       xk(:,6) = (/0d0,0d0,-1d0/2d0/)
       lbl(7) = "S"
       xk(:,7) = (/0d0,1d0/2d0,0d0/)
       lbl(8) = "R"
       xk(:,8) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(9) = "R_2"
       xk(:,9) = (/0d0,1d0/2d0,-1d0/2d0/)
       lbl(10) = "DELTA_0"
       xk(:,10) = (/-X,X,0d0/)
       lbl(11) = "F_0"
       xk(:,11) = (/X,1d0-X,0d0/)
       lbl(12) = "B_0"
       xk(:,12) = (/-X,X,1d0/2d0/)
       lbl(13) = "B_2"
       xk(:,13) = (/-X,X,-1d0/2d0/)
       lbl(14) = "G_0"
       xk(:,14) = (/X,1d0-X,1d0/2d0/)
       lbl(15) = "G_2"
       xk(:,15) = (/X,1d0-X,-1d0/2d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-Y
       ipath(:,2) = (/2,11/) ! Y-F_0
       ipath(:,3) = (/10,1/) ! DELTA_0-GAMMA
       ipath(:,4) = (/1,5/) ! GAMMA-Z
       ipath(:,5) = (/5,12/) ! Z-B_0
       ipath(:,6) = (/14,3/) ! G_0-T
       ipath(:,7) = (/3,2/) ! T-Y
       ipath(:,8) = (/1,7/) ! GAMMA-S
       ipath(:,9) = (/7,8/) ! S-R
       ipath(:,10) = (/8,5/) ! R-Z
       ipath(:,11) = (/5,3/) ! Z-T
    case ("oF1")
       J = (1d0+a*a/b/b-a*a/c/c)/4d0
       H = (1d0+a*a/b/b+a*a/c/c)/4d0
       allocate(lbl(9),xk(3,9),ipath(2,9))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "T"
       xk(:,2) = (/1d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "Z"
       xk(:,3) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(4) = "Y"
       xk(:,4) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(5) = "SIGMA_0"
       xk(:,5) = (/0d0,H,H/)
       lbl(6) = "U_0"
       xk(:,6) = (/1d0,1d0-H,1d0-H/)
       lbl(7) = "A_0"
       xk(:,7) = (/1d0/2d0,1d0/2d0+J,J/)
       lbl(8) = "C_0"
       xk(:,8) = (/1d0/2d0,1d0/2d0-J,1d0-J/)
       lbl(9) = "L"
       xk(:,9) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       ipath(:,1) = (/1,4/) ! GAMMA-Y
       ipath(:,2) = (/4,2/) ! Y-T
       ipath(:,3) = (/2,3/) ! T-Z
       ipath(:,4) = (/3,1/) ! Z-GAMMA
       ipath(:,5) = (/1,5/) ! GAMMA-SIGMA_0
       ipath(:,6) = (/6,2/) ! U_0-T
       ipath(:,7) = (/4,8/) ! Y-C_0
       ipath(:,8) = (/7,3/) ! A_0-Z
       ipath(:,9) = (/1,9/) ! GAMMA-L
    case ("oF2")
       J = (1d0+c*c/a/a-c*c/b/b)/4d0
       K = (1d0+c*c/a/a+c*c/b/b)/4d0
       allocate(lbl(9),xk(3,9),ipath(2,9))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "T"
       xk(:,2) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "Z"
       xk(:,3) = (/1d0/2d0,1d0/2d0,1d0/)
       lbl(4) = "Y"
       xk(:,4) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(5) = "LAMBDA_0"
       xk(:,5) = (/K,K,0d0/)
       lbl(6) = "Q_0"
       xk(:,6) = (/1d0-K,1d0-K,1d0/)
       lbl(7) = "G_0"
       xk(:,7) = (/1d0/2d0-J,1d0-J,1d0/2d0/)
       lbl(8) = "H_0"
       xk(:,8) = (/1d0/2d0+J,J,1d0/2d0/)
       lbl(9) = "L"
       xk(:,9) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-T
       ipath(:,2) = (/2,3/) ! T-Z
       ipath(:,3) = (/3,4/) ! Z-Y
       ipath(:,4) = (/4,1/) ! Y-GAMMA
       ipath(:,5) = (/1,5/) ! GAMMA-LAMBDA_0
       ipath(:,6) = (/6,3/) ! Q_0-Z
       ipath(:,7) = (/2,7/) ! T-G_0
       ipath(:,8) = (/8,4/) ! H_0-Y
       ipath(:,9) = (/1,9/) ! GAMMA-L
    case ("oF3")
       H = (1d0+a*a/b/b-a*a/c/c)/4d0
       K = (1d0+b*b/a/a-b*b/c/c)/4d0
       P = (1d0+c*c/b/b-c*c/a/a)/4d0
       allocate(lbl(11),xk(3,11),ipath(2,10))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "T"
       xk(:,2) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "Z"
       xk(:,3) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(4) = "Y"
       xk(:,4) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(5) = "A_0"
       xk(:,5) = (/1d0/2d0,1d0/2d0+H,H/)
       lbl(6) = "C_0"
       xk(:,6) = (/1d0/2d0,1d0/2d0-H,1d0-H/)
       lbl(7) = "B_0"
       xk(:,7) = (/1d0/2d0+K,1d0/2d0,K/)
       lbl(8) = "D_0"
       xk(:,8) = (/1d0/2d0-K,1d0/2d0,1d0-K/)
       lbl(9) = "G_0"
       xk(:,9) = (/P,1d0/2d0+P,1d0/2d0/)
       lbl(10) = "H_0"
       xk(:,10) = (/1d0-P,1d0/2d0-P,1d0/2d0/)
       lbl(11) = "L"
       xk(:,11) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       ipath(:,1) = (/1,4/) ! GAMMA-Y
       ipath(:,2) = (/4,6/) ! Y-C_0
       ipath(:,3) = (/5,3/) ! A_0-Z
       ipath(:,4) = (/3,7/) ! Z-B_0
       ipath(:,5) = (/8,2/) ! D_0-T
       ipath(:,6) = (/2,9/) ! T-G_0
       ipath(:,7) = (/10,4/) ! H_0-Y
       ipath(:,8) = (/2,1/) ! T-GAMMA
       ipath(:,9) = (/1,3/) ! GAMMA-Z
       ipath(:,10) = (/1,11/) ! GAMMA-L
    case ("oI1")
       Z = (1d0+a*a/c/c)/4d0
       H = (1d0+b*b/c/c)/4d0
       D = (b*b-a*a)/4d0/c/c
       N = (a*a+b*b)/4d0/c/c
       allocate(lbl(13),xk(3,13),ipath(2,11))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "X"
       xk(:,2) = (/1d0/2d0,1d0/2d0,-1d0/2d0/)
       lbl(3) = "S"
       xk(:,3) = (/1d0/2d0,0d0,0d0/)
       lbl(4) = "R"
       xk(:,4) = (/0d0,1d0/2d0,0d0/)
       lbl(5) = "T"
       xk(:,5) = (/0d0,0d0,1d0/2d0/)
       lbl(6) = "W"
       xk(:,6) = (/1d0/4d0,1d0/4d0,1d0/4d0/)
       lbl(7) = "SIGMA_0"
       xk(:,7) = (/-Z,Z,Z/)
       lbl(8) = "F_2"
       xk(:,8) = (/Z,1d0-Z,-Z/)
       lbl(9) = "Y_0"
       xk(:,9) = (/H,-H,H/)
       lbl(10) = "U_0"
       xk(:,10) = (/1d0-H,H,-H/)
       lbl(11) = "L_0"
       xk(:,11) = (/-N,N,1d0/2d0-D/)
       lbl(12) = "M_0"
       xk(:,12) = (/N,-N,1d0/2d0+D/)
       lbl(13) = "J_0"
       xk(:,13) = (/1d0/2d0-D,1d0/2d0+D,-N/)
       ipath(:,1) = (/1,2/) ! GAMMA-X
       ipath(:,2) = (/2,8/) ! X-F_2
       ipath(:,3) = (/7,1/) ! SIGMA_0-GAMMA
       ipath(:,4) = (/1,9/) ! GAMMA-Y_0
       ipath(:,5) = (/10,2/) ! U_0-X
       ipath(:,6) = (/1,4/) ! GAMMA-R
       ipath(:,7) = (/4,6/) ! R-W
       ipath(:,8) = (/6,3/) ! W-S
       ipath(:,9) = (/3,1/) ! S-GAMMA
       ipath(:,10) = (/1,5/) ! GAMMA-T
       ipath(:,11) = (/5,6/) ! T-W
    case ("oI2")
       Z = (1d0+b*b/a/a)/4d0
       H = (1d0+c*c/a/a)/4d0
       D = (c*c-b*b)/4d0/a/a
       N = (b*b+c*c)/4d0/a/a
       allocate(lbl(13),xk(3,13),ipath(2,11))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "X"
       xk(:,2) = (/-1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "S"
       xk(:,3) = (/1d0/2d0,0d0,0d0/)
       lbl(4) = "R"
       xk(:,4) = (/0d0,1d0/2d0,0d0/)
       lbl(5) = "T"
       xk(:,5) = (/0d0,0d0,1d0/2d0/)
       lbl(6) = "W"
       xk(:,6) = (/1d0/4d0,1d0/4d0,1d0/4d0/)
       lbl(7) = "Y_0"
       xk(:,7) = (/Z,-Z,Z/)
       lbl(8) = "U_2"
       xk(:,8) = (/-Z,Z,1d0-Z/)
       lbl(9) = "LAMBDA_0"
       xk(:,9) = (/H,H,-H/)
       lbl(10) = "G_2"
       xk(:,10) = (/-H,1d0-H,H/)
       lbl(11) = "K"
       xk(:,11) = (/1d0/2d0-D,-N,N/)
       lbl(12) = "K_2"
       xk(:,12) = (/1d0/2d0+D,N,-N/)
       lbl(13) = "K_4"
       xk(:,13) = (/-N,1d0/2d0-D,1d0/2d0+D/)
       ipath(:,1) = (/1,2/) ! GAMMA-X
       ipath(:,2) = (/2,8/) ! X-U_2
       ipath(:,3) = (/7,1/) ! Y_0-GAMMA
       ipath(:,4) = (/1,9/) ! GAMMA-LAMBDA_0
       ipath(:,5) = (/10,2/) ! G_2-X
       ipath(:,6) = (/1,4/) ! GAMMA-R
       ipath(:,7) = (/4,6/) ! R-W
       ipath(:,8) = (/6,3/) ! W-S
       ipath(:,9) = (/3,1/) ! S-GAMMA
       ipath(:,10) = (/1,5/) ! GAMMA-T
       ipath(:,11) = (/5,6/) ! T-W
    case ("oI3")
       Z = (1d0+c*c/b/b)/4d0
       Y = (1d0+a*a/b/b)/4d0
       D = (a*a-c*c)/4d0/b/b
       M = (c*c+a*a)/4d0/b/b
       allocate(lbl(13),xk(3,13),ipath(2,11))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "X"
       xk(:,2) = (/1d0/2d0,-1d0/2d0,1d0/2d0/)
       lbl(3) = "S"
       xk(:,3) = (/1d0/2d0,0d0,0d0/)
       lbl(4) = "R"
       xk(:,4) = (/0d0,1d0/2d0,0d0/)
       lbl(5) = "T"
       xk(:,5) = (/0d0,0d0,1d0/2d0/)
       lbl(6) = "W"
       xk(:,6) = (/1d0/4d0,1d0/4d0,1d0/4d0/)
       lbl(7) = "SIGMA_0"
       xk(:,7) = (/-Y,Y,Y/)
       lbl(8) = "F_0"
       xk(:,8) = (/Y,-Y,1d0-Y/)
       lbl(9) = "LAMBDA_0"
       xk(:,9) = (/Z,Z,-Z/)
       lbl(10) = "G_0"
       xk(:,10) = (/1d0-Z,-Z,Z/)
       lbl(11) = "V_0"
       xk(:,11) = (/M,1d0/2d0-D,-M/)
       lbl(12) = "H_0"
       xk(:,12) = (/-M,1d0/2d0+D,M/)
       lbl(13) = "H_2"
       xk(:,13) = (/1d0/2d0+D,-M,1d0/2d0-D/)
       ipath(:,1) = (/1,2/) ! GAMMA-X
       ipath(:,2) = (/2,8/) ! X-F_0
       ipath(:,3) = (/7,1/) ! SIGMA_0-GAMMA
       ipath(:,4) = (/1,9/) ! GAMMA-LAMBDA_0
       ipath(:,5) = (/10,2/) ! G_0-X
       ipath(:,6) = (/1,4/) ! GAMMA-R
       ipath(:,7) = (/4,6/) ! R-W
       ipath(:,8) = (/6,3/) ! W-S
       ipath(:,9) = (/3,1/) ! S-GAMMA
       ipath(:,10) = (/1,5/) ! GAMMA-T
       ipath(:,11) = (/5,6/) ! T-W
    case ("oP1")
       allocate(lbl(8),xk(3,8),ipath(2,12))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "X"
       xk(:,2) = (/1d0/2d0,0d0,0d0/)
       lbl(3) = "Z"
       xk(:,3) = (/0d0,0d0,1d0/2d0/)
       lbl(4) = "U"
       xk(:,4) = (/1d0/2d0,0d0,1d0/2d0/)
       lbl(5) = "Y"
       xk(:,5) = (/0d0,1d0/2d0,0d0/)
       lbl(6) = "S"
       xk(:,6) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(7) = "T"
       xk(:,7) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(8) = "R"
       xk(:,8) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       ipath(:,1) = (/1,2/) ! GAMMA-X
       ipath(:,2) = (/2,6/) ! X-S
       ipath(:,3) = (/6,5/) ! S-Y
       ipath(:,4) = (/5,1/) ! Y-GAMMA
       ipath(:,5) = (/1,3/) ! GAMMA-Z
       ipath(:,6) = (/3,4/) ! Z-U
       ipath(:,7) = (/4,8/) ! U-R
       ipath(:,8) = (/8,7/) ! R-T
       ipath(:,9) = (/7,3/) ! T-Z
       ipath(:,10) = (/2,4/) ! X-U
       ipath(:,11) = (/5,7/) ! Y-T
       ipath(:,12) = (/6,8/) ! S-R
    case ("tI1")
       H = (1d0+c*c/a/a)/4d0
       allocate(lbl(7),xk(3,7),ipath(2,8))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "M"
       xk(:,2) = (/-1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(3) = "X"
       xk(:,3) = (/0d0,0d0,1d0/2d0/)
       lbl(4) = "P"
       xk(:,4) = (/1d0/4d0,1d0/4d0,1d0/4d0/)
       lbl(5) = "Z"
       xk(:,5) = (/H,H,-H/)
       lbl(6) = "Z_0"
       xk(:,6) = (/-H,1d0-H,H/)
       lbl(7) = "N"
       xk(:,7) = (/0d0,1d0/2d0,0d0/)
       ipath(:,1) = (/1,3/) ! GAMMA-X
       ipath(:,2) = (/3,2/) ! X-M
       ipath(:,3) = (/2,1/) ! M-GAMMA
       ipath(:,4) = (/1,5/) ! GAMMA-Z
       ipath(:,5) = (/6,2/) ! Z_0-M
       ipath(:,6) = (/3,4/) ! X-P
       ipath(:,7) = (/4,7/) ! P-N
       ipath(:,8) = (/7,1/) ! N-GAMMA
    case ("tI2")
       H = (1d0+a*a/c/c)/4d0
       Z = a*a/2d0/c/c
       allocate(lbl(9),xk(3,9),ipath(2,9))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "M"
       xk(:,2) = (/1d0/2d0,1d0/2d0,-1d0/2d0/)
       lbl(3) = "X"
       xk(:,3) = (/0d0,0d0,1d0/2d0/)
       lbl(4) = "P"
       xk(:,4) = (/1d0/4d0,1d0/4d0,1d0/4d0/)
       lbl(5) = "N"
       xk(:,5) = (/0d0,1d0/2d0,0d0/)
       lbl(6) = "S_0"
       xk(:,6) = (/-H,H,H/)
       lbl(7) = "S"
       xk(:,7) = (/H,1d0-H,-H/)
       lbl(8) = "R"
       xk(:,8) = (/-Z,Z,1d0/2d0/)
       lbl(9) = "G"
       xk(:,9) = (/1d0/2d0,1d0/2d0,-Z/)
       ipath(:,1) = (/1,3/) ! GAMMA-X
       ipath(:,2) = (/3,4/) ! X-P
       ipath(:,3) = (/4,5/) ! P-N
       ipath(:,4) = (/5,1/) ! N-GAMMA
       ipath(:,5) = (/1,2/) ! GAMMA-M
       ipath(:,6) = (/2,7/) ! M-S
       ipath(:,7) = (/6,1/) ! S_0-GAMMA
       ipath(:,8) = (/3,8/) ! X-R
       ipath(:,9) = (/9,2/) ! G-M
    case ("tP1")
       allocate(lbl(6),xk(3,6),ipath(2,9))
       lbl(1) = "GAMMA"
       xk(:,1) = (/0d0,0d0,0d0/)
       lbl(2) = "Z"
       xk(:,2) = (/0d0,0d0,1d0/2d0/)
       lbl(3) = "M"
       xk(:,3) = (/1d0/2d0,1d0/2d0,0d0/)
       lbl(4) = "A"
       xk(:,4) = (/1d0/2d0,1d0/2d0,1d0/2d0/)
       lbl(5) = "R"
       xk(:,5) = (/0d0,1d0/2d0,1d0/2d0/)
       lbl(6) = "X"
       xk(:,6) = (/0d0,1d0/2d0,0d0/)
       ipath(:,1) = (/1,6/) ! GAMMA-X
       ipath(:,2) = (/6,3/) ! X-M
       ipath(:,3) = (/3,1/) ! M-GAMMA
       ipath(:,4) = (/1,2/) ! GAMMA-Z
       ipath(:,5) = (/2,5/) ! Z-R
       ipath(:,6) = (/5,4/) ! R-A
       ipath(:,7) = (/4,2/) ! A-Z
       ipath(:,8) = (/6,5/) ! X-R
       ipath(:,9) = (/3,4/) ! M-A
    end select

  end subroutine hpkot_table

  !xx! end generated

end submodule kpath
