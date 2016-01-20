! Copyright (c) 2015 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
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

!> Space group information from space group symbol.
!> Based on:
!> * Uri Shmueli, Theories and techniques of crystal structure determination,
!>   ed. IUCr, Oxford University Press, 2007
!> * Shmueli, U., Space-group algorithms. I. The space group and its symmetry
!>   elements. Acta Crystallogr. A40 (1984) 559-567.
module spgs
  implicit none

  private

  ! public interface variables, run spgs_driver to calculate these
  integer :: spgs_id !< space group number id (1-306)
  character*(24) :: spgs_shortspg !< short spg identifier
  character*(24) :: spgs_longspg !< long spg identifier
  integer :: spgs_n !< number of spg operations (quotient translation subgroup)
  integer, dimension(3,4,48) :: spgs_m !< spg ops. (quotient translation subgroup). spgs_m(1:3,1:3,i) is the 
                                       !< (integer) rotation matrix. spgs_m(:,4,i) is the translation vector
                                       !< for the operation in units of 1/12th of unit cell.
  integer :: spgs_ncv !< number of centering vectors
  integer, dimension(3,4) :: spgs_cen !< centering vectors (integer)
  integer :: spgs_sys !< crystal system 
  integer, parameter :: spgs_tri = 1 !< Crystal system: trigonal
  integer, parameter :: spgs_mon = 2 !< Crystal system: monoclinic
  integer, parameter :: spgs_ort = 3 !< Crystal system: orthorhombic
  integer, parameter :: spgs_tet = 4 !< Crystal system: tetragonal
  integer, parameter :: spgs_rom = 6 !< Crystal system: rhombohedral
  integer, parameter :: spgs_hex = 7 !< Crystal system: hexagonal
  integer, parameter :: spgs_cub = 8 !< Crystal system: cubic
  integer :: spgs_lcent !< lattice centering type
  logical :: spgs_inv !< true if centrosymmetric
  integer :: spgs_laue !< laue class. One of 1=1bar, 2=2/m, 3=mmm, 4=4/m, 5=4/mm, 6=3bar, 7=3bar/m, 8=6/m, 9=6/mmm, 10=m3bar, 11=m3barm

  ! private to spgr
  integer :: npol !< internal spgr
  integer :: jrt(3,4,24) !< internal spgr
  real*8  :: cen(3,16) !< centering vectors
  real*8  :: rotm(3,4,48) !< symmetry operations
  integer :: laue !< laue group (1=1bar,2=2/m,3=mmm,4=4/m,5=4/mm,6=r3r,7=r3mr,
                  !< 8=3, 9=3m1, 10=31m, 11=6/m, 12=6/mmm, 13=m3, 14=m3m)
  integer :: lcent !< lattice centering 1=P, 2=A, 3=B, 4=C, 5=I, 6=F, 7=R.
  integer :: lsys !< crystal system 1=tri,2=mon,3=ort,4=tet,5=triR,6=triH,7=hex,8=cub.
  integer :: naxis !< unique axis for monoclinic groups (1=x,2=y,3=z)
  integer :: ncent !< 1 for centrosymmetric, 0 for non-cs
  integer :: ncv !< number of centering vectors
  integer :: neqv !< number of symmetry operations

  ! public/protected interface
  public :: spgs_id, spgs_shortspg, spgs_longspg, spgs_n, spgs_m, spgs_ncv
  public :: spgs_cen, spgs_sys, spgs_tri, spgs_mon, spgs_ort, spgs_tet
  public :: spgs_rom, spgs_hex, spgs_cub, spgs_lcent
  public :: spgs_inv, spgs_laue
  protected :: spgs_id, spgs_shortspg, spgs_longspg, spgs_n, spgs_m, spgs_ncv
  protected :: spgs_cen, spgs_sys, spgs_lcent
  protected :: spgs_inv, spgs_laue
  
  ! public routines
  public :: spgs_init
  public :: spgs_driver

  ! long string identifiers
  character*24, parameter, dimension(306) :: spgs_longstr = (/character(len=24)::&
       'PAN$P1A000',&
       'PAC$I1A000',&
       'PMN$P2B000',&              
       'PMN$P2C000',&              
       'PMN$P2B060',&              
       'PMN$P2C006',&              
       'CMN$P2B000',&              
       'CMN$P2B000',&              
       'CMN$P2B000',&              
       'CMN$P2C000',&              
       'CMN$P2C000',&              
       'CMN$P2C000',&              
       'PMN$I2B000',&              
       'PMN$I2C000',&              
       'PMN$I2B006',&              
       'PMN$I2B606',&              
       'PMN$I2B600',&              
       'PMN$I2C600',&              
       'PMN$I2C660',&              
       'PMN$I2C060',&              
       'CMN$I2B000',&              
       'CMN$I2B000',&              
       'CMN$I2B000',&              
       'CMN$I2C000',&              
       'CMN$I2C000',&              
       'CMN$I2C000',&              
       'CMN$I2B006',&              
       'CMN$I2B606',&              
       'CMN$I2B600',&              
       'CMN$I2C600',&              
       'CMN$I2C660',&              
       'CMN$I2C060',&              
       'PMC$I1A000$P2B000',&       
       'PMC$I1A000$P2C000',&       
       'PMC$I1A000$P2B060',&       
       'PMC$I1A000$P2C006',&       
       'CMC$I1A000$P2B000',&       
       'CMC$I1A000$P2B000',&       
       'CMC$I1A000$P2B000',&       
       'CMC$I1A000$P2C000',&       
       'CMC$I1A000$P2C000',&       
       'CMC$I1A000$P2C000',&       
       'PMC$I1A000$P2B006',&       
       'PMC$I1A000$P2B606',&       
       'PMC$I1A000$P2B600',&       
       'PMC$I1A000$P2C600',&       
       'PMC$I1A000$P2C660',&       
       'PMC$I1A000$P2C060',&       
       'PMC$I1A000$P2B066',&       
       'PMC$I1A000$P2B666',&       
       'PMC$I1A000$P2B660',&       
       'PMC$I1A000$P2C606',&       
       'PMC$I1A000$P2C666',&       
       'PMC$I1A000$P2C066',&       
       'CMC$I1A000$P2B006',&       
       'CMC$I1A000$P2B606',&       
       'CMC$I1A000$P2B600',&       
       'CMC$I1A000$P2C600',&       
       'CMC$I1A000$P2C660',&       
       'CMC$I1A000$P2C060',&       
       'PON$P2C000$P2A000',&       
       'PON$P2C006$P2A000',&       
       'PON$P2C000$P2A660',&       
       'PON$P2C606$P2A660',&       
       'CON$P2C006$P2A000',&       
       'CON$P2C000$P2A000',&       
       'FON$P2C000$P2A000',&       
       'ION$P2C000$P2A000',&       
       'ION$P2C606$P2A660',&       
       'PON$P2C000$I2A000',&       
       'PON$P2C006$I2A000',&       
       'PON$P2C000$I2A006',&       
       'PON$P2C000$I2A600',&       
       'PON$P2C006$I2A606',&       
       'PON$P2C000$I2A066',&       
       'PON$P2C606$I2A000',&       
       'PON$P2C000$I2A660',&       
       'PON$P2C006$I2A666',&       
       'PON$P2C000$I2A666',&       
       'CON$P2C000$I2A000',&       
       'CON$P2C006$I2A000',&       
       'CON$P2C000$I2A006',&       
       'AON$P2C000$I2A000',&       
       'AON$P2C000$I2A060',&       
       'AON$P2C000$I2A600',&       
       'AON$P2C000$I2A660',&       
       'FON$P2C000$I2A000',&       
       'FON$P2C000$I2A333',&       
       'ION$P2C000$I2A000',&       
       'ION$P2C000$I2A660',&       
       'ION$P2C000$I2A600',&       
       'POC$I1A000$P2C000$P2A000',&
       'POC$I1A666$P2C000$P2A000',&
       'POC$I1A000$P2C660$P2A066',&
       'POC$I1A000$P2C000$P2A006',&
       'POC$I1A660$P2C000$P2A000',&
       'POC$I1A000$P2C660$P2A060',&
       'POC$I1A000$P2C600$P2A600',&
       'POC$I1A000$P2C600$P2A066',&
       'POC$I1A000$P2C606$P2A000',&
       'POC$I1A000$P2C600$P2A606',&
       'POC$I1A000$P2C000$P2A660',&
       'POC$I1A000$P2C660$P2A606',&
       'POC$I1A000$P2C006$P2A060',&
       'POC$I1A000$P2C000$P2A666',&
       'POC$I1A660$P2C000$P2A660',&
       'POC$I1A000$P2C660$P2A600',&
       'POC$I1A000$P2C666$P2A660',&
       'POC$I1A000$P2C606$P2A660',&
       'POC$I1A000$P2C606$P2A666',&
       'COC$I1A000$P2C006$P2A000',&
       'COC$I1A000$P2C066$P2A000',&
       'COC$I1A000$P2C000$P2A000',&
       'COC$I1A000$P2C000$P2A006',&
       'COC$I1A000$P2C060$P2A000',&
       'COC$I1A066$P2C660$P2A660',&
       'COC$I1A000$P2C600$P2A606',&
       'FOC$I1A000$P2C000$P2A000',&
       'FOC$I1A333$P2C000$P2A000',&
       'FOC$I1A000$P2C990$P2A099',&
       'IOC$I1A000$P2C000$P2A000',&
       'IOC$I1A000$P2C000$P2A660',&
       'IOC$I1A000$P2C606$P2A660',&
       'IOC$I1A000$P2C060$P2A000',&
       'PTN$P4C000',&              
       'PTN$P4C003',&              
       'PTN$P4C006',&              
       'PTN$P4C009',&              
       'ITN$P4C000',&              
       'ITN$P4C063',&              
       'PTN$I4C000',&              
       'ITN$I4C000',&              
       'PTC$I1A000$P4C000',&       
       'PTC$I1A000$P4C006',&       
       'PTC$I1A660$P4C660',&       
       'PTC$I1A000$P4C600',&       
       'PTC$I1A666$P4C666',&       
       'PTC$I1A000$P4C066',&       
       'ITC$I1A000$P4C000',&       
       'ITC$I1A063$P4C063',&       
       'ITC$I1A000$P4C933',&       
       'PTN$P4C000$P2A000',&       
       'PTN$P4C660$P2A660',&       
       'PTN$P4C003$P2A006',&       
       'PTN$P4C663$P2A669',&       
       'PTN$P4C006$P2A000',&       
       'PTN$P4C666$P2A666',&       
       'PTN$P4C009$P2A006',&       
       'PTN$P4C669$P2A663',&       
       'ITN$P4C000$P2A000',&       
       'ITN$P4C063$P2A063',&       
       'PTN$P4C000$I2A000',&       
       'PTN$P4C000$I2A660',&       
       'PTN$P4C006$I2A006',&       
       'PTN$P4C666$I2A666',&       
       'PTN$P4C000$I2A006',&       
       'PTN$P4C000$I2A666',&       
       'PTN$P4C006$I2A000',&       
       'PTN$P4C006$I2A660',&       
       'ITN$P4C000$I2A000',&       
       'ITN$P4C000$I2A006',&       
       'ITN$P4C063$I2A666',&       
       'ITN$P4C063$I2A660',&       
       'PTN$I4C000$P2A000',&       
       'PTN$I4C000$P2A006',&       
       'PTN$I4C000$P2A660',&       
       'PTN$I4C000$P2A666',&       
       'PTN$I4C000$P2D000',&       
       'PTN$I4C000$P2D006',&       
       'PTN$I4C000$P2D660',&       
       'PTN$I4C000$P2D666',&       
       'ITN$I4C000$P2D000',&       
       'ITN$I4C000$P2D006',&       
       'ITN$I4C000$P2A000',&       
       'ITN$I4C000$P2A609',&       
       'PTC$I1A000$P4C000$P2A000',&
       'PTC$I1A000$P4C000$P2A006',&
       'PTC$I1A660$P4C000$P2A000',&
       'PTC$I1A000$P4C600$P2A060',&
       'PTC$I1A666$P4C000$P2A000',&
       'PTC$I1A000$P4C600$P2A066',&
       'PTC$I1A000$P4C000$P2A660',&
       'PTC$I1A000$P4C000$P2A666',&
       'PTC$I1A660$P4C660$P2A660',&
       'PTC$I1A000$P4C600$P2A600',&
       'PTC$I1A660$P4C660$P2A666',&
       'PTC$I1A000$P4C600$P2A606',&
       'PTC$I1A000$P4C006$P2A000',&
       'PTC$I1A000$P4C006$P2A006',&
       'PTC$I1A666$P4C666$P2A006',&
       'PTC$I1A000$P4C606$P2A060',&
       'PTC$I1A666$P4C666$P2A000',&
       'PTC$I1A000$P4C606$P2A066',&
       'PTC$I1A000$P4C006$P2A660',&
       'PTC$I1A000$P4C666$P2A666',&
       'PTC$I1A666$P4C666$P2A666',&
       'PTC$I1A000$P4C606$P2A600',&
       'PTC$I1A666$P4C666$P2A660',&
       'PTC$I1A000$P4C606$P2A606',&
       'ITC$I1A000$P4C000$P2A000',&
       'ITC$I1A000$P4C000$P2A006',&
       'ITC$I1A063$P4C063$P2A063',&
       'ITC$I1A000$P4C393$P2A000',&
       'ITC$I1A063$P4C063$P2A069',&
       'ITC$I1A000$P4C393$P2A006',&
       'PRN$P3C000',&              
       'PRN$P3C004',&              
       'PRN$P3C008',&              
       'RRN$P3C000',&              
       'PRN$P3Q000',&              
       'PRC$I3C000',&              
       'RRC$I3C000',&              
       'PRC$I3Q000',&              
       'PRN$P3C000$P2G000',&       
       'PRN$P3C000$P2F000',&       
       'PRN$P3C004$P2G000',&       
       'PRN$P3C004$P2F008',&       
       'PRN$P3C008$P2G000',&       
       'PRN$P3C008$P2F004',&       
       'RRN$P3C000$P2F000',&       
       'PRN$P3Q000$P2E000',&       
       'PRN$P3C000$I2F000',&       
       'PRN$P3C000$I2G000',&       
       'PRN$P3C000$I2F006',&       
       'PRN$P3C000$I2G006',&       
       'RRN$P3C000$I2F000',&       
       'PRN$P3Q000$I2E000',&       
       'RRN$P3C000$I2F006',&       
       'PRN$P3Q000$I2E666',&       
       'PRC$I3C000$P2G000',&       
       'PRC$I3C000$P2G006',&       
       'PRC$I3C000$P2F000',&       
       'PRC$I3C000$P2F006',&       
       'RRC$I3C000$P2F000',&       
       'PRC$I3Q000$P2E000',&       
       'RRC$I3C000$P2F006',&       
       'PRC$I3Q000$P2E666',&       
       'PHN$P6C000',&              
       'PHN$P6C002',&              
       'PHN$P6C005',&              
       'PHN$P6C004',&              
       'PHN$P6C008',&              
       'PHN$P6C006',&              
       'PHN$I6C000',&              
       'PHC$I1A000$P6C000',&       
       'PHC$I1A000$P6C006',&       
       'PHN$P6C000$P2F000',&       
       'PHN$P6C002$P2F000',&       
       'PHN$P6C005$P2F000',&       
       'PHN$P6C004$P2F000',&       
       'PHN$P6C008$P2F000',&       
       'PHN$P6C006$P2F000',&       
       'PHN$P6C000$I2F000',&       
       'PHN$P6C000$I2F006',&       
       'PHN$P6C006$I2F006',&       
       'PHN$P6C006$I2F000',&       
       'PHN$I6C000$P2G000',&       
       'PHN$I6C006$P2G000',&       
       'PHN$I6C000$P2F000',&       
       'PHN$I6C006$P2F000',&       
       'PHC$I1A000$P6C000$P2F000',&
       'PHC$I1A000$P6C000$P2F006',&
       'PHC$I1A000$P6C006$P2F006',&
       'PHC$I1A000$P6C006$P2F000',&
       'PCN$P3Q000$P2C000$P2A000',&
       'FCN$P3Q000$P2C000$P2A000',&
       'ICN$P3Q000$P2C000$P2A000',&
       'PCN$P3Q000$P2C606$P2A660',&
       'ICN$P3Q000$P2C606$P2A660',&
       'PCC$I3Q000$P2C000$P2A000',&
       'PCC$I3Q666$P2C000$P2A000',&
       'PCC$I3Q000$P2C660$P2A066',&
       'FCC$I3Q000$P2C000$P2A000',&
       'FCC$I3Q333$P2C000$P2A000',&
       'FCC$I3Q000$P2C330$P2A033',&
       'ICC$I3Q000$P2C000$P2A000',&
       'PCC$I3Q000$P2C606$P2A660',&
       'ICC$I3Q000$P2C606$P2A660',&
       'PCN$P3Q000$P4C000$P2D000',&
       'PCN$P3Q000$P4C666$P2D666',&
       'FCN$P3Q000$P4C000$P2D000',&
       'FCN$P3Q000$P4C993$P2D939',&
       'ICN$P3Q000$P4C000$P2D000',&
       'PCN$P3Q000$P4C939$P2D399',&
       'PCN$P3Q000$P4C393$P2D933',&
       'ICN$P3Q000$P4C393$P2D933',&
       'PCN$P3Q000$I4C000$I2D000',&
       'FCN$P3Q000$I4C000$I2D000',&
       'ICN$P3Q000$I4C000$I2D000',&
       'PCN$P3Q000$I4C666$I2D666',&
       'FCN$P3Q000$I4C666$I2D666',&
       'ICN$P3Q000$I4C939$I2D399',&
       'PCC$I3Q000$P4C000$P2D000',&
       'PCC$I3Q666$P4C000$P2D000',&
       'PCC$I3Q000$P4C600$P2D006',&
       'PCC$I3Q000$P4C666$P2D666',&
       'PCC$I3Q666$P4C666$P2D666',&
       'PCC$I3Q000$P4C066$P2D660',&
       'FCC$I3Q000$P4C000$P2D000',&
       'FCC$I3Q000$P4C666$P2D666',&
       'FCC$I3Q333$P4C993$P2D939',&
       'FCC$I3Q000$P4C693$P2D936',&
       'FCC$I3Q999$P4C993$P2D939',&
       'FCC$I3Q000$P4C093$P2D930',&
       'ICC$I3Q000$P4C000$P2D000',&
       'ICC$I3Q000$P4C393$P2D933'/)

  ! short string identifiers
  character*24, parameter, dimension(306) :: spgs_shortstr = (/character(len=24)::&
       'p 1',&           ! 1. P 1
       'p -1',&          ! 2. P -1
       'p 1 2 1',&       ! 3. P 2 (unique axis b)
       'p 1 1 2',&       ! 3. P 2 (unique axis c)
       'p 1 21 1',&      ! 4. P 21 (unique axis b)
       'p 1 1 21',&      ! 4. P 21 (unique axis c)
       'c 1 2 1',&       ! 5. C 2 (unique axis b)
       'a 1 2 1',&       ! 5. A 1 2 1 (alternative ITA setting)
       'i 1 2 1',&       ! 5. I 1 2 1 (alternative ITA setting)
       'a 1 1 2',&       ! 5. A 1 1 2 (alternative ITA setting)
       'b 1 1 2',&       ! 5. B 1 1 2 (alternative ITA setting)
       'i 1 1 2',&       ! 5. I 1 1 2 (alternative ITA setting)
       'p 1 m 1',&       ! 6. P m (unique axis b)
       'p 1 1 m',&       ! 6. P m (unique axis c)
       'p 1 c 1',&       ! 7. P c (unique axis b)
       'p 1 n 1',&       ! 7. P 1 n 1 (alternative ITA setting)       
       'p 1 a 1',&       ! 7. P 1 a 1 (alternative ITA setting)       
       'p 1 1 a',&       ! 7. P 1 1 a (alternative ITA setting)       
       'p 1 1 n',&       ! 7. P 1 1 n (alternative ITA setting)
       'p 1 1 b',&       ! 7. P 1 1 b (alternative ITA setting)
       'c 1 m 1',&       ! 8. C m (unique axis b)
       'a 1 m 1',&       ! 8. A 1 m 1 (alternative ITA setting)       
       'i 1 m 1',&       ! 8. I 1 m 1 (alternative ITA setting)       
       'a 1 1 m',&       ! 8. A 1 1 m (alternative ITA setting)       
       'b 1 1 m',&       ! 8. B 1 1 m (alternative ITA setting)       
       'i 1 1 m',&       ! 8. I 1 1 m (alternative ITA setting)
       'c 1 c 1',&       ! 9. C c (unique axis b)       
       'a 1 n 1',&       ! 9. A 1 n 1 (alternative ITA setting)       
       'i 1 a 1',&       ! 9. I 1 a 1 (alternative ITA setting)       
       'a 1 1 a',&       ! 9. A 1 1 a (alternative ITA setting)       
       'b 1 1 n',&       ! 9. B 1 1 n (alternative ITA setting)       
       'i 1 1 b',&       ! 9. I 1 1 b (alternative ITA setting)       
       'p 1 2/m 1',&     ! 10. P 2/m (unique axis b)  
       'p 1 1 2/m',&     ! 10. P 2/m (unique axis c)  
       'p 1 21/m 1',&    ! 11. P 21/m (unique axis b)   
       'p 1 1 21/m',&    ! 11. P 21/m (unique axis c)
       'c 1 2/m 1',&     ! 12. C 2/m (unique axis b)    
       'a 1 2/m 1',&     ! 12. A 1 2/m 1 (alternative ITA setting)    
       'i 1 2/m 1',&     ! 12. I 1 2/m 1 (alternative ITA setting)    
       'a 1 1 2/m',&     ! 12. A 1 1 2/m (alternative ITA setting)    
       'b 1 1 2/m',&     ! 12. B 1 1 2/m (alternative ITA setting)    
       'i 1 1 2/m',&     ! 12. I 1 1 2/m (alternative ITA setting)  
       'p 1 2/c 1',&     ! 13. P 2/c (unique axis b)      
       'p 1 2/n 1',&     ! 13. P 1 2/n 1 (alternative ITA setting)      
       'p 1 2/a 1',&     ! 13. P 1 2/a 1 (alternative ITA setting)    
       'p 1 1 2/a',&     ! 13. P 1 1 2/a (alternative ITA setting)    
       'p 1 1 2/n',&     ! 13. P 1 1 2/n (alternative ITA setting)    
       'p 1 1 2/b',&     ! 13. P 1 1 2/b (alternative ITA setting)  
       'p 1 21/c 1',&    ! 14. P 21/c (unique axis b)   
       'p 1 21/n 1',&    ! 14. P 1 21/n 1 (alternative ITA setting)   
       'p 1 21/a 1',&    ! 14. P 1 21/a 1 (alternative ITA setting)   
       'p 1 1 21/a',&    ! 14. P 1 1 21/a (alternative ITA setting)   
       'p 1 1 21/n',&    ! 14. P 1 1 21/n (alternative ITA setting)
       'p 1 1 21/b',&    ! 14. P 1 1 21/b (alternative ITA setting)
       'c 1 2/c 1',&     ! 15. C 2/c (unique axis b)    
       'a 1 2/n 1',&     ! 15. A 1 2/n 1 (alternative ITA setting)    
       'i 1 2/a 1',&     ! 15. I 1 2/a 1 (alternative ITA setting)    
       'a 1 1 2/a',&     ! 15. A 1 1 2/a (alternative ITA setting)    
       'b 1 1 2/n',&     ! 15. B 1 1 2/n (alternative ITA setting)  
       'i 1 1 2/b',&     ! 15. I 1 1 2/b (alternative ITA setting)  
       'p 2 2 2',&       ! 16. P 2 2 2
       'p 2 2 21',&      ! 17. P 2 2 21 
       'p 21 21 2',&     ! 18. P 21 21 2  
       'p 21 21 21',&    ! 19. P 21 21 21   
       'c 2 2 21',&      ! 20. C 2 2 21 
       'c 2 2 2',&       ! 21. C 2 2 2
       'f 2 2 2',&       ! 22. F 2 2 2
       'i 2 2 2',&       ! 23. I 2 2 2
       'i 21 21 21',&    ! 24. I 21 21 21   
       'p m m 2',&       ! 25. P m m 2
       'p m c 21',&      ! 26. P m c 21
       'p c c 2',&       ! 27. P c c 2
       'p m a 2',&       ! 28. P m a 2
       'p c a 21',&      ! 29. P c a 21
       'p n c 2',&       ! 30. P n c 2
       'p m n 21',&      ! 31. P m n 21
       'p b a 2',&       ! 32. P b a 2
       'p n a 21',&      ! 33. P n a 21
       'p n n 2',&       ! 34. P n n 2
       'c m m 2',&       ! 35. C m m 2
       'c m c 21',&      ! 36. C m c 21
       'c c c 2',&       ! 37. C c c 2
       'a m m 2',&       ! 38. A m m 2
       'a b m 2',&       ! 39. A e m 2 (also: A b m 2)
       'a m a 2',&       ! 40. A m a 2
       'a b a 2',&       ! 41. A e a 2 (also: A b a 2)
       'f m m 2',&       ! 42. F m m 2
       'f d d 2',&       ! 43. F d d 2
       'i m m 2',&       ! 44. I m m 2
       'i b a 2',&       ! 45. I b a 2
       'i m a 2',&       ! 46. I m a 2
       'p m m m',&       ! 47. P m m m
       'p n n n 1',&     ! 48. P n n n (origin choice 1)
       'p n n n 2',&     ! 48. P n n n (origin choice 2)
       'p c c m',&       ! 49. P c c m
       'p b a n 1',&     ! 50. P b a n (origin choice 1)
       'p b a n 2',&     ! 50. P b a n (origin choice 2)
       'p m m a',&       ! 51. P m m a
       'p n n a',&       ! 52. P n n a
       'p m n a',&       ! 53. P m n a
       'p c c a',&       ! 54. P c c a
       'p b a m',&       ! 55. P b a m
       'p c c n',&       ! 56. P c c n
       'p b c m',&       ! 57. P b c m
       'p n n m',&       ! 58. P n n m
       'p m m n 1',&     ! 59. P m m n (origin choice 1)
       'p m m n 2',&     ! 59. P m m n (origin choice 2)
       'p b c n',&       ! 60. P b c n
       'p b c a',&       ! 61. P b c a
       'p n m a',&       ! 62. P n m a
       'c m c m',&       ! 63. C m c m
       'c m c a',&       ! 64. C m c a
       'c m m m',&       ! 65. C m m m
       'c c c m',&       ! 66. C c c m
       'c m m a',&       ! 67. C m m a
       'c c c a 1',&     ! 68. C c c a (origin choice 1)
       'c c c a 2',&     ! 68. C c c a (origin choice 2)
       'f m m m',&       ! 69. F m m m
       'f d d d 1',&     ! 70. F d d d (origin choice 1)
       'f d d d 2',&     ! 70. F d d d (origin choice 2)
       'i m m m',&       ! 71. I m m m
       'i b a m',&       ! 72. I b a m
       'i b c a',&       ! 73. I b c a
       'i m m a',&       ! 74. I m m a
       'p 4',&           ! 75. P 4  
       'p 41',&          ! 76. P 41   
       'p 42',&          ! 77. P 42
       'p 43',&          ! 78. P 43   
       'i 4',&           ! 79. I 4  
       'i 41',&          ! 80. I 41   
       'p -4',&          ! 81. P -4   
       'i -4',&          ! 82. I -4   
       'p 4/m',&         ! 83. P 4/m
       'p 42/m',&        ! 84. P 42/m
       'p 4/n 1',&       ! 85. P 4/n (origin choice 1)
       'p 4/n 2',&       ! 85. P 4/n (origin choice 2)
       'p 42/n 1',&      ! 86. P 42/n (origin choice 1)
       'p 42/n 2',&      ! 86. P 42/n (origin choice 2)
       'i 4/m',&         ! 87. I 4/m
       'i 41/a 1',&      ! 88. I 41/a (origin choice 1)
       'i 41/a 2',&      ! 88. I 41/a (origin choice 2)
       'p 4 2 2',&       ! 89. P 4 2 2
       'p 4 21 2',&      ! 90. P 4 21 2
       'p 41 2 2',&      ! 91. P 41 2 2
       'p 41 21 2',&     ! 92. P 41 21 2 
       'p 42 2 2',&      ! 93. P 42 2 2
       'p 42 21 2',&     ! 94. P 42 21 2 
       'p 43 2 2',&      ! 95. P 43 2 2
       'p 43 21 2',&     ! 96. P 43 21 2 
       'i 4 2 2',&       ! 97. I 4 2 2
       'i 41 2 2',&      ! 98. I 41 2 2
       'p 4 m m',&       ! 99. P 4 m m
       'p 4 b m',&       ! 100. P 4 b m
       'p 42 c m',&      ! 101. P 42 c m
       'p 42 n m',&      ! 102. P 42 n m
       'p 4 c c',&       ! 103. P 4 c c
       'p 4 n c',&       ! 104. P 4 n c
       'p 42 m c',&      ! 105. P 42 m c
       'p 42 b c',&      ! 106. P 42 b c
       'i 4 m m',&       ! 107. I 4 m m 
       'i 4 c m',&       ! 108. I 4 c m
       'i 41 m d',&      ! 109. I 41 m d
       'i 41 c d',&      ! 110. I 41 c d
       'p -4 2 m',&      ! 111. P -4 2 m
       'p -4 2 c',&      ! 112. P -4 2 c
       'p -4 21 m',&     ! 113. P -4 21 m 
       'p -4 21 c',&     ! 114. P -4 21 c 
       'p -4 m 2',&      ! 115. P -4 m 2
       'p -4 c 2',&      ! 116. P -4 c 2
       'p -4 b 2',&      ! 117. P -4 b 2
       'p -4 n 2',&      ! 118. P -4 n 2
       'i -4 m 2',&      ! 119. I -4 m 2
       'i -4 c 2',&      ! 120. I -4 c 2
       'i -4 2 m',&      ! 121. I -4 2 m
       'i -4 2 d',&      ! 122. I -4 2 d
       'p 4/m m m',&     ! 123. P 4/m m m
       'p 4/m c c',&     ! 124. P 4/m c c
       'p 4/n b m 1',&   ! 125. P 4/n b m (origin choice 1)
       'p 4/n b m 2',&   ! 125. P 4/n b m (origin choice 2)
       'p 4/n n c 1',&   ! 126. P 4/n n c (origin choice 1)
       'p 4/n n c 2',&   ! 126. P 4/n n c (origin choice 2)
       'p 4/m b m',&     ! 127. P 4/m b m
       'p 4/m n c',&     ! 128. P 4/m n c
       'p 4/n m m 1',&   ! 129. P 4/n m m (origin choice 1)
       'p 4/n m m 2',&   ! 129. P 4/n m m (origin choice 2)
       'p 4/n c c 1',&   ! 130. P 4/n c c (origin choice 1)
       'p 4/n c c 2',&   ! 130. P 4/n c c (origin choice 2)
       'p 42/m m c',&    ! 131. P 42/m m c
       'p 42/m c m',&    ! 132. P 42/m c m
       'p 42/n b c 1',&  ! 133. P 42/n b c (origin choice 1)
       'p 42/n b c 2',&  ! 133. P 42/n b c (origin choice 2)
       'p 42/n n m 1',&  ! 134. P 42/n n m (origin choice 1)
       'p 42/n n m 2',&  ! 134. P 42/n n m (origin choice 2)
       'p 42/m b c',&    ! 135. P 42/m b c
       'p 42/m n m',&    ! 136. P 42/m n m
       'p 42/n m c 1',&  ! 137. P 42/n m c (origin choice 1)
       'p 42/n m c 2',&  ! 137. P 42/n m c (origin choice 2)
       'p 42/n c m 1',&  ! 138. P 42/n c m (origin choice 1)
       'p 42/n m c 2',&  ! 138. P 42/n c m (origin choice 2)
       'i 4/m m m',&     ! 139. I 4/m m m
       'i 4/m c m',&     ! 140. I 4/m c m
       'i 41/a m d 1',&  ! 141. I 41/a m d (origin choice 1)
       'i 41/a m d 2',&  ! 141. I 41/a m d (origin choice 2)
       'i 41/a c d 1',&  ! 142. I 41/a c d (origin choice 1)
       'i 41/a c d 2',&  ! 142. I 41/a c d (origin choice 2)
       'p 3',&           ! 143. P 3    
       'p 31',&          ! 144. P 31     
       'p 32',&          ! 145. P 32     
       'r 3 h',&         ! 146. R 3 h (hexagonal axes)      
       'r 3 r',&         ! 146. R 3 r (rhombohedral axes)      
       'p -3',&          ! 147. P -3     
       'r -3 h',&        ! 148. R -3 (hexagonal axes)       
       'r -3 r',&        ! 148. R -3 (rhomobhedral axes)       
       'p 3 1 2',&       ! 149. P 3 1 2 
       'p 3 2 1',&       ! 150. P 3 2 1 
       'p 31 1 2',&      ! 151. P 31 1 2  
       'p 31 2 1',&      ! 152. P 31 2 1  
       'p 32 1 2',&      ! 153. P 32 1 2  
       'p 32 2 1',&      ! 154. P 32 2 1  
       'r 3 2 h',&       ! 155. R 3 2 (hexagonal axes) 
       'r 3 2 r',&       ! 155. R 3 2 (rhombohedral axes) 
       'p 3 m 1',&       ! 156. P 3 m 1 
       'p 3 1 m',&       ! 157. P 3 1 m 
       'p 3 c 1',&       ! 158. P 3 c 1 
       'p 3 1 c',&       ! 159. P 3 1 c 
       'r 3 m h',&       ! 160. R 3 m (hexagonal axes) 
       'r 3 m r',&       ! 160. R 3 m (rhombohedral axes) 
       'r 3 c h',&       ! 161. R 3 c (hexagonal axes) 
       'r 3 c r',&       ! 161. R 3 c (rhombohedral axes)
       'p -3 1 m',&      ! 162. P -3 1 m  
       'p -3 1 c',&      ! 163. P -3 1 c  
       'p -3 m 1',&      ! 164. P -3 m 1  
       'p -3 c 1',&      ! 165. P -3 c 1  
       'r -3 m h',&      ! 166. R -3 m (hexagonal axes)  
       'r -3 m r',&      ! 166. R -3 m (rhomohedral axes)  
       'r -3 c h',&      ! 167. R -3 c (hexagonal axes)  
       'r -3 c r',&      ! 167. R -3 c (rhombohedral axes)  
       'p 6',&           ! 168. P 6    
       'p 61',&          ! 169. P 61     
       'p 65',&          ! 170. P 65     
       'p 62',&          ! 171. P 62     
       'p 64',&          ! 172. P 64     
       'p 63',&          ! 173. P 63     
       'p -6',&          ! 174. P -6     
       'p 6/m',&         ! 175. P6/m
       'p 63/m',&        ! 176. P 63/m
       'p 6 2 2',&       ! 177. P 6 2 2 
       'p 61 2 2',&      ! 178. P 61 2 2  
       'p 65 2 2',&      ! 179. P 65 2 2  
       'p 62 2 2',&      ! 180. P 62 2 2  
       'p 64 2 2',&      ! 181. P 64 2 2  
       'p 63 2 2',&      ! 182. P 63 2 2  
       'p 6 m m',&       ! 183. P 6 m m 
       'p 6 c c',&       ! 184. P 6 c c 
       'p 63 c m',&      ! 185. P 63 c m  
       'p 63 m c',&      ! 186. P 63 m c  
       'p -6 m 2',&      ! 187. P -6 m 2  
       'p -6 c 2',&      ! 188. P -6 c 2  
       'p -6 2 m',&      ! 189. P -6 2 m  
       'p -6 2 c',&      ! 190. P -6 2 c  
       'p 6/m m m',&     ! 191. P 6/m m m
       'p 6/m c c',&     ! 192. P 6/m c c
       'p 63/m c m',&    ! 193. P 63/m c m
       'p 63/m m c',&    ! 194. P 63/m m c
       'p 2 3',&         ! 195. P 2 3
       'f 2 3',&         ! 196. F 2 3
       'i 2 3',&         ! 197. I 2 3
       'p 21 3',&        ! 198. P 21 3
       'i 21 3',&        ! 199. I 21 3
       'p m -3',&        ! 200. P m -3
       'p n -3 1',&      ! 201. P n -3 (origin choice 1)
       'p n -3 2',&      ! 201. P n -3 (origin choice 2)
       'f m -3',&        ! 202. F m -3
       'f d -3 1',&      ! 203. F d -3 (origin choice 1)
       'f d -3 2',&      ! 203. F d -3 (origin choice 2)
       'i m -3',&        ! 204. I m -3
       'p a -3',&        ! 205. P a -3
       'i a -3',&        ! 206. I a -3
       'p 4 3 2',&       ! 207. P 4 3 2
       'p 42 3 2',&      ! 208. P 42 3 2
       'f 4 3 2',&       ! 209. F 4 3 2
       'f 41 3 2',&      ! 210. F 41 3 2
       'i 4 3 2',&       ! 211. I 4 3 2
       'p 43 3 2',&      ! 212. P 43 3 2
       'p 41 3 2',&      ! 213. P 41 3 2
       'i 41 3 2',&      ! 214. I 41 3 2
       'p -4 3 m',&      ! 215. P -4 3 m
       'f -4 3 m',&      ! 216. F -4 3 m
       'i -4 3 m',&      ! 217. I -4 3 m
       'p -4 3 n',&      ! 218. P -4 3 n
       'f -4 3 c',&      ! 219. F -4 3 c
       'i -4 3 d',&      ! 220. I -4 3 d
       'p m -3 m',&      ! 221. P m -3 m
       'p n -3 n 1',&    ! 222. P n -3 n (origin choice 1)
       'p n -3 n 2',&    ! 222. P n -3 n (origin choice 2)
       'p m -3 n',&      ! 223. P m -3 n 
       'p n -3 m 1',&    ! 224. P n -3 m (origin choice 1)
       'p n -3 m 2',&    ! 224. P n -3 m (origin choice 2)
       'f m -3 m',&      ! 225. F m -3 m
       'f m -3 c',&      ! 226. F m -3 c
       'f d -3 m 1',&    ! 227. F d -3 m (origin choice 1)
       'f d -3 m 2',&    ! 227. F d -3 m (origin choice 2)
       'f d -3 c 1',&    ! 228. F d -3 c (origin choice 1)
       'f d -3 c 2',&    ! 228. F d -3 c (origin choice 2)
       'i m -3 m',&      ! 229. I m -3 m
       'i a -3 d'/)      ! 230. I a -3 d

  integer, parameter :: malias = 1000
  character*24 :: spgalias(malias)
  integer :: ialias(malias), nalias

  ! spgs internals
  integer, dimension(4) :: matorder

contains

  !> Generate all symmetry information from the space group label
  subroutine spgs_driver(spgin,usespgr)
    use tools_io

    character(len=*), intent(in) :: spgin
    logical, intent(in) :: usespgr

    integer :: i, n, lpt, lptx, ier
    character(len(spgin)) :: spgin0
    logical :: isblank

    ! clean, compact, and lowercase
    spgin0 = string(spgin)
    spgs_shortspg = ""
    isblank = .false.
    n = 0
    do i = 1, len_trim(spgin0)
       if (spgin0(i:i) == " ") then
          if (.not.isblank) then
             n = n+1
             spgs_shortspg(n:n) = " "
             isblank = .true.
          endif
       else
          isblank = .false.
          n = n+1
          spgs_shortspg(n:n) = spgin0(i:i)
       endif
    end do
    spgs_shortspg = lower(adjustl(spgs_shortspg))

    spgs_id = 0
    if (.not.usespgr) then
       ! check the master list
       do i = 1, 306
          if (trim(spgs_shortstr(i)) == trim(spgs_shortspg)) then
             spgs_id = i
             spgs_longspg = spgs_longstr(i)
             exit
          end if
       end do

       ! check the alias list
       if (spgs_id == 0) then
          do i = 1, nalias
             if (trim(spgalias(i)) == trim(spgs_shortspg)) then
                spgs_id = ialias(i)
                spgs_longspg = spgs_longstr(ialias(i))
                exit
             end if
          end do
       endif
    end if

    ! argh
    if (spgs_id == 0) then
       ! well, let's pass it to SPGR and hope for the best
       rotm = 0d0
       call spgr(spgs_shortspg,-1,-1,ier)

       ! I give up
       if (ier /= 0) &
          call ferror('spgs_driver','Unknown space group symbol',faterr,spgin)

       ! translate to spgs notation
       spgs_longspg = spgs_shortspg
       spgs_n = neqv
       spgs_m(1:3,1:3,1:24) = nint(rotm(1:3,1:3,1:24))
       spgs_m(1:3,4,1:24) = nint(rotm(1:3,4,1:24) * 12d0)
       if (ncent == 1) then
          spgs_n = 2 * neqv
          do i = 1, neqv
             spgs_m(1:3,1:3,neqv+i) = -spgs_m(1:3,1:3,i)
             spgs_m(1:3,4,neqv+i) = spgs_m(1:3,4,i)
          end do
       endif
       spgs_ncv = ncv
       spgs_cen = nint(cen(:,1:4) * 12d0)
       if (lsys==1.or.lsys==2.or.lsys==3.or.lsys==4.or.lsys==7.or.lsys==8) then
          spgs_sys = lsys
       elseif (lsys==5.or.lsys==6) then
          spgs_sys = 6
       endif
       call spgs_getlaue()
       spgs_lcent = lcent
       spgs_inv = (ncent==1)
    else
       ! do the stuff
       call spgs_parse()
       call spgs_generate()
       call spgs_getlaue()
    end if

  end subroutine spgs_driver

  subroutine spgs_parse()
    use tools_io

    integer :: i,l
    logical :: proper

    l = len_trim(spgs_longspg)

    ! lattice symbol
    if (spgs_longspg(1:1) == 'P') then
       spgs_lcent = 1
       spgs_ncv = 1
       spgs_cen(:,1) = (/ 0, 0, 0 /)
    else if (spgs_longspg(1:1) == 'A') then
       spgs_lcent = 2
       spgs_ncv = 2
       spgs_cen(:,1) = (/ 0, 0, 0 /)
       spgs_cen(:,2) = (/ 0, 6, 6 /)
    else if (spgs_longspg(1:1) == 'B') then
       spgs_lcent = 3
       spgs_ncv = 2
       spgs_cen(:,1) = (/ 0, 0, 0 /)
       spgs_cen(:,2) = (/ 6, 0, 6 /)
    else if (spgs_longspg(1:1) == 'C') then
       spgs_lcent = 4
       spgs_ncv = 2
       spgs_cen(:,1) = (/ 0, 0, 0 /)
       spgs_cen(:,2) = (/ 6, 6, 0 /)
    else if (spgs_longspg(1:1) == 'I') then
       spgs_lcent = 5
       spgs_ncv = 2
       spgs_cen(:,1) = (/ 0, 0, 0 /)
       spgs_cen(:,2) = (/ 6, 6, 6 /)
    else if (spgs_longspg(1:1) == 'F') then
       spgs_lcent = 6
       spgs_ncv = 4
       spgs_cen(:,1) = (/ 0, 0, 0 /)
       spgs_cen(:,2) = (/ 0, 6, 6 /)
       spgs_cen(:,3) = (/ 6, 0, 6 /)
       spgs_cen(:,4) = (/ 6, 6, 0 /)
    else if (spgs_longspg(1:1) == 'R') then
       spgs_lcent = 7
       spgs_ncv = 3
       spgs_cen(:,1) = (/ 0, 0, 0 /)
       spgs_cen(:,2) = (/ 4, 8, 8 /)
       spgs_cen(:,3) = (/ 8, 4, 4 /)
    else
       call ferror('spgs','Unrecognized lattice letter',faterr)
    end if

    ! crystal system
    if (spgs_longspg(2:2) == 'A') then
       spgs_sys = spgs_tri
    else if (spgs_longspg(2:2) == 'M') then
       spgs_sys = spgs_mon
    else if (spgs_longspg(2:2) == 'O') then
       spgs_sys = spgs_ort
    else if (spgs_longspg(2:2) == 'T') then
       spgs_sys = spgs_tet
    else if (spgs_longspg(2:2) == 'R') then
       spgs_sys = spgs_rom
    else if (spgs_longspg(2:2) == 'H') then
       spgs_sys = spgs_hex
    else if (spgs_longspg(2:2) == 'C') then
       spgs_sys = spgs_cub
    else
       call ferror('spgs','Unrecognized crystal system letter',faterr)
    end if

    ! centrosymmetry
    spgs_inv = (spgs_longspg(3:3) == 'C')

    ! the identity goes first
    spgs_n = 1
    spgs_m(1,:,1) = (/ 1, 0, 0, 0 /)
    spgs_m(2,:,1) = (/ 0, 1, 0, 0 /)
    spgs_m(3,:,1) = (/ 0, 0, 1, 0 /)
    i = 4
    do while (spgs_longspg(i:i) == '$')
       ! new generator
       spgs_n = spgs_n + 1
       ! proper or improper rotation
       i = i + 1
       proper = (spgs_longspg(i:i) == 'P')
       ! read rotation matrix
       i = i + 1 
       read(spgs_longspg(i:i),*) matorder(spgs_n)
       if (.not.proper .and. mod(matorder(spgs_n),2) == 1) then
          matorder(spgs_n) = matorder(spgs_n) * 2
       end if
       select case (spgs_longspg(i:i+1))
       case ("1A")
          if (.not.proper) then
             spgs_m(1,1:3,spgs_n) = (/ 1, 0, 0 /)
             spgs_m(2,1:3,spgs_n) = (/ 0, 1, 0 /)
             spgs_m(3,1:3,spgs_n) = (/ 0, 0, 1 /)
          else
             ! I already have the identity
             spgs_n = spgs_n - 1
             i = i + 5
             cycle
          end if
       case ("2A")
          spgs_m(1,1:3,spgs_n) = (/ 1, 0, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 0,-1, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0,-1 /)
       case ("2B")
          spgs_m(1,1:3,spgs_n) = (/-1, 0, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 0, 1, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0,-1 /)
       case ("2C")
          spgs_m(1,1:3,spgs_n) = (/-1, 0, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 0,-1, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0, 1 /)
       case ("2D")
          spgs_m(1,1:3,spgs_n) = (/ 0, 1, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 1, 0, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0,-1 /)
       case ("2E")
          spgs_m(1,1:3,spgs_n) = (/ 0,-1, 0 /)
          spgs_m(2,1:3,spgs_n) = (/-1, 0, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0,-1 /)
       case ("2F")
          spgs_m(1,1:3,spgs_n) = (/ 1,-1, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 0,-1, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0,-1 /)
       case ("2G")
          spgs_m(1,1:3,spgs_n) = (/ 1, 0, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 1,-1, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0,-1 /)
       case ("3Q")
          spgs_m(1,1:3,spgs_n) = (/ 0, 0, 1 /)
          spgs_m(2,1:3,spgs_n) = (/ 1, 0, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 1, 0 /)
       case ("3C")
          spgs_m(1,1:3,spgs_n) = (/ 0,-1, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 1,-1, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0, 1 /)
       case ("4C")
          spgs_m(1,1:3,spgs_n) = (/ 0,-1, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 1, 0, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0, 1 /)
       case ("6C")
          spgs_m(1,1:3,spgs_n) = (/ 1,-1, 0 /)
          spgs_m(2,1:3,spgs_n) = (/ 1, 0, 0 /)
          spgs_m(3,1:3,spgs_n) = (/ 0, 0, 1 /)
       case default
          call ferror('spgs','Unrecognized generator rotation matrix',faterr)
       end select
       if (.not.proper) then
          spgs_m(:,1:3,spgs_n) = -spgs_m(:,1:3,spgs_n)
       end if
       ! read translation vector
       i = i + 2
       read (spgs_longspg(i:i),*) spgs_m(1,4,spgs_n)
       i = i + 1
       read (spgs_longspg(i:i),*) spgs_m(2,4,spgs_n)
       i = i + 1
       read (spgs_longspg(i:i),*) spgs_m(3,4,spgs_n)
       if (spgs_m(3,4,spgs_n) == 5) then
          spgs_m(3,4,spgs_n) = 10
       end if
       i = i + 1
       if (i > len(spgs_longspg)) exit
    end do

  end subroutine spgs_parse

  subroutine spgs_generate()

    integer :: i, j, k
    integer, dimension(3,4) :: amat, bmat, cmat
    integer :: ngen

    ! save number of generators
    ngen = spgs_n

    ! put amat to identity and add
    amat(1,:) = (/ 1, 0, 0, 0 /)
    amat(2,:) = (/ 0, 1, 0, 0 /)
    amat(3,:) = (/ 0, 0, 1, 0 /)

    ! generate full list
    if (ngen > 1) then
       do i = 0, matorder(2)
          if (ngen > 2) then
             bmat = amat
             do j = 0, matorder(3)
                if (ngen > 3) then
                   cmat = bmat
                   do k = 0, matorder(4)
                      cmat = spgs_mult(cmat,spgs_m(:,:,4))
                      call spgs_add(cmat)
                   end do
                end if
                bmat = spgs_mult(bmat,spgs_m(:,:,3))
                call spgs_add(bmat)
             end do
          end if
          amat = spgs_mult(amat,spgs_m(:,:,2))
          call spgs_add(amat)
       end do
    end if

  end subroutine spgs_generate

  subroutine spgs_getlaue()
    use tools_io

    integer :: order
    ! laue: 1=1bar, 2=2/m, 3=mmm, 4=4/m, 5=4/mmm, 6=3bar, 7=3bar/m, 8=6/m, 9=6/mmm
    !       10=m3bar, 11=m3barm

    order = spgs_laueorder()

    if (spgs_sys == spgs_tri) then
       spgs_laue = 1
    else if (spgs_sys == spgs_mon) then
       spgs_laue = 2
    else if (spgs_sys == spgs_ort) then
       spgs_laue = 3
    else if (spgs_sys == spgs_tet) then
       if (order == 8) then
          spgs_laue = 4
       else if (order == 16) then
          spgs_laue = 5
       else
          call ferror('spgs_getlaue','Incorrect laue class order',faterr)
       end if
    else if (spgs_sys == spgs_rom) then
       if (order == 6) then
          spgs_laue = 6
       else if (order == 12) then
          spgs_laue = 7
       else
          call ferror('spgs_getlaue','Incorrect laue class order',faterr)
       end if
    else if (spgs_sys == spgs_hex) then
       if (order == 12) then
          spgs_laue = 8
       else if (order == 24) then
          spgs_laue = 9
       else
          call ferror('spgs_getlaue','Incorrect laue class order',faterr)
       end if
    else if (spgs_sys == spgs_cub) then
       if (order == 24) then
          spgs_laue = 10
       else if (order == 48) then
          spgs_laue = 11
       else
          call ferror('spgs_getlaue','Incorrect laue class order',faterr)
       end if
    end if

  end subroutine spgs_getlaue

  function spgs_laueorder()

    integer :: spgs_laueorder
    integer :: i, j
    integer :: repeated

    repeated = 0
    do i = 1, spgs_n
       do j = i+1, spgs_n
          if (all(spgs_m(:,1:3,i) == spgs_m(:,1:3,j))) then
             repeated = repeated + 1
             exit
          end if
       end do
    end do

    spgs_laueorder = spgs_n - repeated
    if (.not.spgs_inv) then
       spgs_laueorder = spgs_laueorder * 2
    end if

  end function spgs_laueorder

  subroutine spgs_add(mat)
    integer, dimension(3,4), intent(in) :: mat

    integer :: i

    do i = 1, spgs_n
       if (all(mat == spgs_m(:,:,i))) then
          return
       end if
    end do

    spgs_n = spgs_n + 1
    spgs_m(:,:,spgs_n) = mat

  end subroutine spgs_add

  function spgs_mult(mat1,mat2)
    integer, dimension(3,4), intent(in) :: mat1, mat2
    integer, dimension(3,4) :: spgs_mult

    spgs_mult(:,1:3) = matmul(mat1(:,1:3),mat2(:,1:3))
    spgs_mult(:,4) = matmul(mat1(:,1:3),mat2(:,4)) + mat1(:,4)

    where (spgs_mult(:,4) >= 12)
       spgs_mult(:,4) = mod(spgs_mult(:,4),12)
    elsewhere (spgs_mult(:,4) < 0)
       spgs_mult(:,4) = mod(spgs_mult(:,4),12) + 12
    end where

  end function spgs_mult

  !> Initialize all the alias names in spgs. If not provided, the default
  !> origin is choice 1, the axes are hexagonal and the unique axis is b.
  !> The space group number and the labels without spaces can be used as well.  
  subroutine spgs_init()

    integer :: n

    ! write the alias list
    n = 0
    ! no origin given, assume origin 1
    n = n + 1
    ialias(n) =  3
    spgalias(n)="p 2"
    n = n + 1
    ialias(n) =  3
    spgalias(n)="p 2 b"
    n = n + 1
    ialias(n) =  4
    spgalias(n)="p 2 c"
    n = n + 1
    ialias(n) =  5
    spgalias(n)="p 21"
    n = n + 1
    ialias(n) =  5
    spgalias(n)="p 21 b"
    n = n + 1
    ialias(n) =  6
    spgalias(n)="p 21 c"
    n = n + 1
    ialias(n) =  7
    spgalias(n)="c 2"
    n = n + 1
    ialias(n) =  7
    spgalias(n)="c 2 b"
    n = n + 1
    ialias(n) = 13
    spgalias(n)="p m"
    n = n + 1
    ialias(n) = 13
    spgalias(n)="p m b"
    n = n + 1
    ialias(n) = 14
    spgalias(n)="p m c"
    n = n + 1
    ialias(n) = 15
    spgalias(n)="p c"
    n = n + 1
    ialias(n) = 15
    spgalias(n)="p c b"
    n = n + 1
    ialias(n) = 21
    spgalias(n)="c m"
    n = n + 1
    ialias(n) = 21
    spgalias(n)="c m b"
    n = n + 1
    ialias(n) = 27
    spgalias(n)="c c"
    n = n + 1
    ialias(n) = 27
    spgalias(n)="c c b"
    n = n + 1
    ialias(n) = 33
    spgalias(n)="p 2/m"
    n = n + 1
    ialias(n) = 33
    spgalias(n)="p 2/m b"
    n = n + 1
    ialias(n) = 34
    spgalias(n)="p 2/m c"
    n = n + 1
    ialias(n) = 35
    spgalias(n)="p 21/m"
    n = n + 1
    ialias(n) = 35
    spgalias(n)="p 21/m b"
    n = n + 1
    ialias(n) = 36
    spgalias(n)="p 21/m c"
    n = n + 1
    ialias(n) = 37
    spgalias(n)="c 2/m"
    n = n + 1
    ialias(n) = 37
    spgalias(n)="c 2/m b"
    n = n + 1
    ialias(n) = 43
    spgalias(n)="p 2/c"
    n = n + 1
    ialias(n) = 43
    spgalias(n)="p 2/c b"
    n = n + 1
    ialias(n) = 49
    spgalias(n)="p 21/c"
    n = n + 1
    ialias(n) = 49
    spgalias(n)="p 21/c b"
    n = n + 1
    ialias(n) = 55
    spgalias(n)="c 2/c"
    n = n + 1
    ialias(n) = 55
    spgalias(n)="c 2/c b"
    n = n + 1
    ialias(n) = 84
    spgalias(n)="a e m 2"
    n = n + 1
    ialias(n) = 86
    spgalias(n)="a e a 2"
    n = n + 1
    ialias(n) = 93
    spgalias(n)="p n n n"
    n = n + 1
    ialias(n) = 95
    spgalias(n)="p c c m"
    n = n + 1
    ialias(n) = 96
    spgalias(n)="p b a n"
    n = n + 1
    ialias(n) =106
    spgalias(n)="p m m n"
    n = n + 1
    ialias(n) =116
    spgalias(n)="c c c a"
    n = n + 1
    ialias(n) =119
    spgalias(n)="f d d d"
    n = n + 1
    ialias(n) =135
    spgalias(n)="p 4/n"
    n = n + 1
    ialias(n) =137
    spgalias(n)="p 42/n"
    n = n + 1
    ialias(n) =140
    spgalias(n)="i 41/a"
    n = n + 1
    ialias(n) =178
    spgalias(n)="p 4/n b m"
    n = n + 1
    ialias(n) =180
    spgalias(n)="p 4/n n c"
    n = n + 1
    ialias(n) =184
    spgalias(n)="p 4/n m m"
    n = n + 1
    ialias(n) =186
    spgalias(n)="p 4/n c c"
    n = n + 1
    ialias(n) =190
    spgalias(n)="p 42/n b c"
    n = n + 1
    ialias(n) =192
    spgalias(n)="p 42/n n m"
    n = n + 1
    ialias(n) =196
    spgalias(n)="p 42/n m c"
    n = n + 1
    ialias(n) =198
    spgalias(n)="p 42/n c m"
    n = n + 1
    ialias(n) =202
    spgalias(n)="i 41/a m d"
    n = n + 1
    ialias(n) =204
    spgalias(n)="i 41/a c d"
    n = n + 1
    ialias(n) =209
    spgalias(n)="r 3 h"
    n = n + 1
    ialias(n) =212
    spgalias(n)="r -3"
    n = n + 1
    ialias(n) =220
    spgalias(n)="r 3 2"
    n = n + 1
    ialias(n) =226
    spgalias(n)="r 3 m"
    n = n + 1
    ialias(n) =228
    spgalias(n)="r 3 c"
    n = n + 1
    ialias(n) =234
    spgalias(n)="r -3 m"
    n = n + 1
    ialias(n) =236
    spgalias(n)="r -3 c"
    n = n + 1
    ialias(n) =271
    spgalias(n)="p n -3"
    n = n + 1
    ialias(n) =274
    spgalias(n)="f d -3"
    n = n + 1
    ialias(n) =294
    spgalias(n)="p n -3 n"
    n = n + 1
    ialias(n) =297
    spgalias(n)="p n -3 m"
    n = n + 1
    ialias(n) =301
    spgalias(n)="f d -3 m"
    n = n + 1
    ialias(n) =303
    spgalias(n)="f d -3 c"

    ! space group number
    n = n + 1
    ialias(n) =  1
    spgalias(n)="1"
    n = n + 1
    ialias(n) =  2
    spgalias(n)="2"
    n = n + 1
    ialias(n) =  3
    spgalias(n)="3"
    n = n + 1
    ialias(n) =  5
    spgalias(n)="4"
    n = n + 1
    ialias(n) =  7
    spgalias(n)="5"
    n = n + 1
    ialias(n) = 13
    spgalias(n)="6"
    n = n + 1
    ialias(n) = 15
    spgalias(n)="7"
    n = n + 1
    ialias(n) = 21
    spgalias(n)="8"
    n = n + 1
    ialias(n) = 27
    spgalias(n)="9"
    n = n + 1
    ialias(n) = 33
    spgalias(n)="10"
    n = n + 1
    ialias(n) = 35
    spgalias(n)="11"
    n = n + 1
    ialias(n) = 37
    spgalias(n)="12"
    n = n + 1
    ialias(n) = 43
    spgalias(n)="13"
    n = n + 1
    ialias(n) = 49
    spgalias(n)="14"
    n = n + 1
    ialias(n) = 55
    spgalias(n)="15"
    n = n + 1
    ialias(n) = 61
    spgalias(n)="16"
    n = n + 1
    ialias(n) = 62
    spgalias(n)="17"
    n = n + 1
    ialias(n) = 63
    spgalias(n)="18"
    n = n + 1
    ialias(n) = 64
    spgalias(n)="19"
    n = n + 1
    ialias(n) = 65
    spgalias(n)="20"
    n = n + 1
    ialias(n) = 66
    spgalias(n)="21"
    n = n + 1
    ialias(n) = 67
    spgalias(n)="22"
    n = n + 1
    ialias(n) = 68
    spgalias(n)="23"
    n = n + 1
    ialias(n) = 69
    spgalias(n)="24"
    n = n + 1
    ialias(n) = 70
    spgalias(n)="25"
    n = n + 1
    ialias(n) = 71
    spgalias(n)="26"
    n = n + 1
    ialias(n) = 72
    spgalias(n)="27"
    n = n + 1
    ialias(n) = 73
    spgalias(n)="28"
    n = n + 1
    ialias(n) = 74
    spgalias(n)="29"
    n = n + 1
    ialias(n) = 75
    spgalias(n)="30"
    n = n + 1
    ialias(n) = 76
    spgalias(n)="31"
    n = n + 1
    ialias(n) = 77
    spgalias(n)="32"
    n = n + 1
    ialias(n) = 78
    spgalias(n)="33"
    n = n + 1
    ialias(n) = 79
    spgalias(n)="34"
    n = n + 1
    ialias(n) = 80
    spgalias(n)="35"
    n = n + 1
    ialias(n) = 81
    spgalias(n)="36"
    n = n + 1
    ialias(n) = 82
    spgalias(n)="37"
    n = n + 1
    ialias(n) = 83
    spgalias(n)="38"
    n = n + 1
    ialias(n) = 84
    spgalias(n)="39"
    n = n + 1
    ialias(n) = 85
    spgalias(n)="40"
    n = n + 1
    ialias(n) = 86
    spgalias(n)="41"
    n = n + 1
    ialias(n) = 87
    spgalias(n)="42"
    n = n + 1
    ialias(n) = 88
    spgalias(n)="43"
    n = n + 1
    ialias(n) = 89
    spgalias(n)="44"
    n = n + 1
    ialias(n) = 90
    spgalias(n)="45"
    n = n + 1
    ialias(n) = 91
    spgalias(n)="46"
    n = n + 1
    ialias(n) = 92
    spgalias(n)="47"
    n = n + 1
    ialias(n) = 93
    spgalias(n)="48"
    n = n + 1
    ialias(n) = 95
    spgalias(n)="49"
    n = n + 1
    ialias(n) = 96
    spgalias(n)="50"
    n = n + 1
    ialias(n) = 98
    spgalias(n)="51"
    n = n + 1
    ialias(n) = 99
    spgalias(n)="52"
    n = n + 1
    ialias(n) =100
    spgalias(n)="53"
    n = n + 1
    ialias(n) =101
    spgalias(n)="54"
    n = n + 1
    ialias(n) =102
    spgalias(n)="55"
    n = n + 1
    ialias(n) =103
    spgalias(n)="56"
    n = n + 1
    ialias(n) =104
    spgalias(n)="57"
    n = n + 1
    ialias(n) =105
    spgalias(n)="58"
    n = n + 1
    ialias(n) =106
    spgalias(n)="59"
    n = n + 1
    ialias(n) =108
    spgalias(n)="60"
    n = n + 1
    ialias(n) =109
    spgalias(n)="61"
    n = n + 1
    ialias(n) =110
    spgalias(n)="62"
    n = n + 1
    ialias(n) =111
    spgalias(n)="63"
    n = n + 1
    ialias(n) =112
    spgalias(n)="64"
    n = n + 1
    ialias(n) =113
    spgalias(n)="65"
    n = n + 1
    ialias(n) =114
    spgalias(n)="66"
    n = n + 1
    ialias(n) =115
    spgalias(n)="67"
    n = n + 1
    ialias(n) =116
    spgalias(n)="68"
    n = n + 1
    ialias(n) =118
    spgalias(n)="69"
    n = n + 1
    ialias(n) =119
    spgalias(n)="70"
    n = n + 1
    ialias(n) =121
    spgalias(n)="71"
    n = n + 1
    ialias(n) =122
    spgalias(n)="72"
    n = n + 1
    ialias(n) =123
    spgalias(n)="73"
    n = n + 1
    ialias(n) =124
    spgalias(n)="74"
    n = n + 1
    ialias(n) =125
    spgalias(n)="75"
    n = n + 1
    ialias(n) =126
    spgalias(n)="76"
    n = n + 1
    ialias(n) =127
    spgalias(n)="77"
    n = n + 1
    ialias(n) =128
    spgalias(n)="78"
    n = n + 1
    ialias(n) =129
    spgalias(n)="79"
    n = n + 1
    ialias(n) =130
    spgalias(n)="80"
    n = n + 1
    ialias(n) =131
    spgalias(n)="81"
    n = n + 1
    ialias(n) =132
    spgalias(n)="82"
    n = n + 1
    ialias(n) =133
    spgalias(n)="83"
    n = n + 1
    ialias(n) =134
    spgalias(n)="84"
    n = n + 1
    ialias(n) =135
    spgalias(n)="85"
    n = n + 1
    ialias(n) =137
    spgalias(n)="86"
    n = n + 1
    ialias(n) =139
    spgalias(n)="87"
    n = n + 1
    ialias(n) =140
    spgalias(n)="88"
    n = n + 1
    ialias(n) =142
    spgalias(n)="89"
    n = n + 1
    ialias(n) =143
    spgalias(n)="90"
    n = n + 1
    ialias(n) =144
    spgalias(n)="91"
    n = n + 1
    ialias(n) =145
    spgalias(n)="92"
    n = n + 1
    ialias(n) =146
    spgalias(n)="93"
    n = n + 1
    ialias(n) =147
    spgalias(n)="94"
    n = n + 1
    ialias(n) =148
    spgalias(n)="95"
    n = n + 1
    ialias(n) =149
    spgalias(n)="96"
    n = n + 1
    ialias(n) =150
    spgalias(n)="97"
    n = n + 1
    ialias(n) =151
    spgalias(n)="98"
    n = n + 1
    ialias(n) =152
    spgalias(n)="99"
    n = n + 1
    ialias(n) =153
    spgalias(n)="100"
    n = n + 1
    ialias(n) =154
    spgalias(n)="101"
    n = n + 1
    ialias(n) =155
    spgalias(n)="102"
    n = n + 1
    ialias(n) =156
    spgalias(n)="103"
    n = n + 1
    ialias(n) =157
    spgalias(n)="104"
    n = n + 1
    ialias(n) =158
    spgalias(n)="105"
    n = n + 1
    ialias(n) =159
    spgalias(n)="106"
    n = n + 1
    ialias(n) =160
    spgalias(n)="107"
    n = n + 1
    ialias(n) =161
    spgalias(n)="108"
    n = n + 1
    ialias(n) =162
    spgalias(n)="109"
    n = n + 1
    ialias(n) =163
    spgalias(n)="110"
    n = n + 1
    ialias(n) =164
    spgalias(n)="111"
    n = n + 1
    ialias(n) =165
    spgalias(n)="112"
    n = n + 1
    ialias(n) =166
    spgalias(n)="113"
    n = n + 1
    ialias(n) =167
    spgalias(n)="114"
    n = n + 1
    ialias(n) =168
    spgalias(n)="115"
    n = n + 1
    ialias(n) =169
    spgalias(n)="116"
    n = n + 1
    ialias(n) =170
    spgalias(n)="117"
    n = n + 1
    ialias(n) =171
    spgalias(n)="118"
    n = n + 1
    ialias(n) =172
    spgalias(n)="119"
    n = n + 1
    ialias(n) =173
    spgalias(n)="120"
    n = n + 1
    ialias(n) =174
    spgalias(n)="121"
    n = n + 1
    ialias(n) =175
    spgalias(n)="122"
    n = n + 1
    ialias(n) =176
    spgalias(n)="123"
    n = n + 1
    ialias(n) =177
    spgalias(n)="124"
    n = n + 1
    ialias(n) =178
    spgalias(n)="125"
    n = n + 1
    ialias(n) =180
    spgalias(n)="126"
    n = n + 1
    ialias(n) =182
    spgalias(n)="127"
    n = n + 1
    ialias(n) =183
    spgalias(n)="128"
    n = n + 1
    ialias(n) =184
    spgalias(n)="129"
    n = n + 1
    ialias(n) =186
    spgalias(n)="130"
    n = n + 1
    ialias(n) =188
    spgalias(n)="131"
    n = n + 1
    ialias(n) =189
    spgalias(n)="132"
    n = n + 1
    ialias(n) =190
    spgalias(n)="133"
    n = n + 1
    ialias(n) =192
    spgalias(n)="134"
    n = n + 1
    ialias(n) =194
    spgalias(n)="135"
    n = n + 1
    ialias(n) =195
    spgalias(n)="136"
    n = n + 1
    ialias(n) =196
    spgalias(n)="137"
    n = n + 1
    ialias(n) =198
    spgalias(n)="138"
    n = n + 1
    ialias(n) =200
    spgalias(n)="139"
    n = n + 1
    ialias(n) =201
    spgalias(n)="140"
    n = n + 1
    ialias(n) =202
    spgalias(n)="141"
    n = n + 1
    ialias(n) =204
    spgalias(n)="142"
    n = n + 1
    ialias(n) =206
    spgalias(n)="143"
    n = n + 1
    ialias(n) =207
    spgalias(n)="144"
    n = n + 1
    ialias(n) =208
    spgalias(n)="145"
    n = n + 1
    ialias(n) =209
    spgalias(n)="146"
    n = n + 1
    ialias(n) =211
    spgalias(n)="147"
    n = n + 1
    ialias(n) =212
    spgalias(n)="148"
    n = n + 1
    ialias(n) =214
    spgalias(n)="149"
    n = n + 1
    ialias(n) =215
    spgalias(n)="150"
    n = n + 1
    ialias(n) =216
    spgalias(n)="151"
    n = n + 1
    ialias(n) =217
    spgalias(n)="152"
    n = n + 1
    ialias(n) =218
    spgalias(n)="153"
    n = n + 1
    ialias(n) =219
    spgalias(n)="154"
    n = n + 1
    ialias(n) =220
    spgalias(n)="155"
    n = n + 1
    ialias(n) =222
    spgalias(n)="156"
    n = n + 1
    ialias(n) =223
    spgalias(n)="157"
    n = n + 1
    ialias(n) =224
    spgalias(n)="158"
    n = n + 1
    ialias(n) =225
    spgalias(n)="159"
    n = n + 1
    ialias(n) =226
    spgalias(n)="160"
    n = n + 1
    ialias(n) =228
    spgalias(n)="161"
    n = n + 1
    ialias(n) =230
    spgalias(n)="162"
    n = n + 1
    ialias(n) =231
    spgalias(n)="163"
    n = n + 1
    ialias(n) =232
    spgalias(n)="164"
    n = n + 1
    ialias(n) =233
    spgalias(n)="165"
    n = n + 1
    ialias(n) =234
    spgalias(n)="166"
    n = n + 1
    ialias(n) =236
    spgalias(n)="167"
    n = n + 1
    ialias(n) =238
    spgalias(n)="168"
    n = n + 1
    ialias(n) =239
    spgalias(n)="169"
    n = n + 1
    ialias(n) =240
    spgalias(n)="170"
    n = n + 1
    ialias(n) =241
    spgalias(n)="171"
    n = n + 1
    ialias(n) =242
    spgalias(n)="172"
    n = n + 1
    ialias(n) =243
    spgalias(n)="173"
    n = n + 1
    ialias(n) =244
    spgalias(n)="174"
    n = n + 1
    ialias(n) =245
    spgalias(n)="175"
    n = n + 1
    ialias(n) =246
    spgalias(n)="176"
    n = n + 1
    ialias(n) =247
    spgalias(n)="177"
    n = n + 1
    ialias(n) =248
    spgalias(n)="178"
    n = n + 1
    ialias(n) =249
    spgalias(n)="179"
    n = n + 1
    ialias(n) =250
    spgalias(n)="180"
    n = n + 1
    ialias(n) =251
    spgalias(n)="181"
    n = n + 1
    ialias(n) =252
    spgalias(n)="182"
    n = n + 1
    ialias(n) =253
    spgalias(n)="183"
    n = n + 1
    ialias(n) =254
    spgalias(n)="184"
    n = n + 1
    ialias(n) =255
    spgalias(n)="185"
    n = n + 1
    ialias(n) =256
    spgalias(n)="186"
    n = n + 1
    ialias(n) =257
    spgalias(n)="187"
    n = n + 1
    ialias(n) =258
    spgalias(n)="188"
    n = n + 1
    ialias(n) =259
    spgalias(n)="189"
    n = n + 1
    ialias(n) =260
    spgalias(n)="190"
    n = n + 1
    ialias(n) =261
    spgalias(n)="191"
    n = n + 1
    ialias(n) =262
    spgalias(n)="192"
    n = n + 1
    ialias(n) =263
    spgalias(n)="193"
    n = n + 1
    ialias(n) =264
    spgalias(n)="194"
    n = n + 1
    ialias(n) =265
    spgalias(n)="195"
    n = n + 1
    ialias(n) =266
    spgalias(n)="196"
    n = n + 1
    ialias(n) =267
    spgalias(n)="197"
    n = n + 1
    ialias(n) =268
    spgalias(n)="198"
    n = n + 1
    ialias(n) =269
    spgalias(n)="199"
    n = n + 1
    ialias(n) =270
    spgalias(n)="200"
    n = n + 1
    ialias(n) =271
    spgalias(n)="201"
    n = n + 1
    ialias(n) =273
    spgalias(n)="202"
    n = n + 1
    ialias(n) =274
    spgalias(n)="203"
    n = n + 1
    ialias(n) =276
    spgalias(n)="204"
    n = n + 1
    ialias(n) =277
    spgalias(n)="205"
    n = n + 1
    ialias(n) =278
    spgalias(n)="206"
    n = n + 1
    ialias(n) =279
    spgalias(n)="207"
    n = n + 1
    ialias(n) =280
    spgalias(n)="208"
    n = n + 1
    ialias(n) =281
    spgalias(n)="209"
    n = n + 1
    ialias(n) =282
    spgalias(n)="210"
    n = n + 1
    ialias(n) =283
    spgalias(n)="211"
    n = n + 1
    ialias(n) =284
    spgalias(n)="212"
    n = n + 1
    ialias(n) =285
    spgalias(n)="213"
    n = n + 1
    ialias(n) =286
    spgalias(n)="214"
    n = n + 1
    ialias(n) =287
    spgalias(n)="215"
    n = n + 1
    ialias(n) =288
    spgalias(n)="216"
    n = n + 1
    ialias(n) =289
    spgalias(n)="217"
    n = n + 1
    ialias(n) =290
    spgalias(n)="218"
    n = n + 1
    ialias(n) =291
    spgalias(n)="219"
    n = n + 1
    ialias(n) =292
    spgalias(n)="220"
    n = n + 1
    ialias(n) =293
    spgalias(n)="221"
    n = n + 1
    ialias(n) =294
    spgalias(n)="222"
    n = n + 1
    ialias(n) =296
    spgalias(n)="223"
    n = n + 1
    ialias(n) =297
    spgalias(n)="224"
    n = n + 1
    ialias(n) =299
    spgalias(n)="225"
    n = n + 1
    ialias(n) =300
    spgalias(n)="226"
    n = n + 1
    ialias(n) =301
    spgalias(n)="227"
    n = n + 1
    ialias(n) =303
    spgalias(n)="228"
    n = n + 1
    ialias(n) =305
    spgalias(n)="229"
    n = n + 1
    ialias(n) =306
    spgalias(n)="230"

    ! no spaces
    n = n + 1
    ialias(n) = 1
    spgalias(n) = 'p1'
    n = n + 1
    ialias(n) = 2
    spgalias(n) = 'p-1'
    n = n + 1
    ialias(n) = 3
    spgalias(n) = 'p121'
    n = n + 1
    ialias(n) = 4
    spgalias(n) = 'p112'
    n = n + 1
    ialias(n) = 5
    spgalias(n) = 'p1211'
    n = n + 1
    ialias(n) = 6
    spgalias(n) = 'p1121'
    n = n + 1
    ialias(n) = 7
    spgalias(n) = 'c121'
    n = n + 1
    ialias(n) = 8
    spgalias(n) = 'a121'
    n = n + 1
    ialias(n) = 9
    spgalias(n) = 'i121'
    n = n + 1
    ialias(n) = 10
    spgalias(n) = 'a112'
    n = n + 1
    ialias(n) = 11
    spgalias(n) = 'b112'
    n = n + 1
    ialias(n) = 12
    spgalias(n) = 'i112'
    n = n + 1
    ialias(n) = 13
    spgalias(n) = 'p1m1'
    n = n + 1
    ialias(n) = 14
    spgalias(n) = 'p11m'
    n = n + 1
    ialias(n) = 15
    spgalias(n) = 'p1c1'
    n = n + 1
    ialias(n) = 16
    spgalias(n) = 'p1n1'
    n = n + 1
    ialias(n) = 17
    spgalias(n) = 'p1a1'
    n = n + 1
    ialias(n) = 18
    spgalias(n) = 'p11a'
    n = n + 1
    ialias(n) = 19
    spgalias(n) = 'p11n'
    n = n + 1
    ialias(n) = 20
    spgalias(n) = 'p11b'
    n = n + 1
    ialias(n) = 21
    spgalias(n) = 'c1m1'
    n = n + 1
    ialias(n) = 22
    spgalias(n) = 'a1m1'
    n = n + 1
    ialias(n) = 23
    spgalias(n) = 'i1m1'
    n = n + 1
    ialias(n) = 24
    spgalias(n) = 'a11m'
    n = n + 1
    ialias(n) = 25
    spgalias(n) = 'b11m'
    n = n + 1
    ialias(n) = 26
    spgalias(n) = 'i11m'
    n = n + 1
    ialias(n) = 27
    spgalias(n) = 'c1c1'
    n = n + 1
    ialias(n) = 28
    spgalias(n) = 'a1n1'
    n = n + 1
    ialias(n) = 29
    spgalias(n) = 'i1a1'
    n = n + 1
    ialias(n) = 30
    spgalias(n) = 'a11a'
    n = n + 1
    ialias(n) = 31
    spgalias(n) = 'b11n'
    n = n + 1
    ialias(n) = 32
    spgalias(n) = 'i11b'
    n = n + 1
    ialias(n) = 33
    spgalias(n) = 'p12/m1'
    n = n + 1
    ialias(n) = 34
    spgalias(n) = 'p112/m'
    n = n + 1
    ialias(n) = 35
    spgalias(n) = 'p121/m1'
    n = n + 1
    ialias(n) = 36
    spgalias(n) = 'p1121/m'
    n = n + 1
    ialias(n) = 37
    spgalias(n) = 'c12/m1'
    n = n + 1
    ialias(n) = 38
    spgalias(n) = 'a12/m1'
    n = n + 1
    ialias(n) = 39
    spgalias(n) = 'i12/m1'
    n = n + 1
    ialias(n) = 40
    spgalias(n) = 'a112/m'
    n = n + 1
    ialias(n) = 41
    spgalias(n) = 'b112/m'
    n = n + 1
    ialias(n) = 42
    spgalias(n) = 'i112/m'
    n = n + 1
    ialias(n) = 43
    spgalias(n) = 'p12/c1'
    n = n + 1
    ialias(n) = 44
    spgalias(n) = 'p12/n1'
    n = n + 1
    ialias(n) = 45
    spgalias(n) = 'p12/a1'
    n = n + 1
    ialias(n) = 46
    spgalias(n) = 'p112/a'
    n = n + 1
    ialias(n) = 47
    spgalias(n) = 'p112/n'
    n = n + 1
    ialias(n) = 48
    spgalias(n) = 'p112/b'
    n = n + 1
    ialias(n) = 49
    spgalias(n) = 'p121/c1'
    n = n + 1
    ialias(n) = 50
    spgalias(n) = 'p121/n1'
    n = n + 1
    ialias(n) = 51
    spgalias(n) = 'p121/a1'
    n = n + 1
    ialias(n) = 52
    spgalias(n) = 'p1121/a'
    n = n + 1
    ialias(n) = 53
    spgalias(n) = 'p1121/n'
    n = n + 1
    ialias(n) = 54
    spgalias(n) = 'p1121/b'
    n = n + 1
    ialias(n) = 55
    spgalias(n) = 'c12/c1'
    n = n + 1
    ialias(n) = 56
    spgalias(n) = 'a12/n1'
    n = n + 1
    ialias(n) = 57
    spgalias(n) = 'i12/a1'
    n = n + 1
    ialias(n) = 58
    spgalias(n) = 'a112/a'
    n = n + 1
    ialias(n) = 59
    spgalias(n) = 'b112/n'
    n = n + 1
    ialias(n) = 60
    spgalias(n) = 'i112/b'
    n = n + 1
    ialias(n) = 61
    spgalias(n) = 'p222'
    n = n + 1
    ialias(n) = 62
    spgalias(n) = 'p2221'
    n = n + 1
    ialias(n) = 63
    spgalias(n) = 'p21212'
    n = n + 1
    ialias(n) = 64
    spgalias(n) = 'p212121'
    n = n + 1
    ialias(n) = 65
    spgalias(n) = 'c2221'
    n = n + 1
    ialias(n) = 66
    spgalias(n) = 'c222'
    n = n + 1
    ialias(n) = 67
    spgalias(n) = 'f222'
    n = n + 1
    ialias(n) = 68
    spgalias(n) = 'i222'
    n = n + 1
    ialias(n) = 69
    spgalias(n) = 'i212121'
    n = n + 1
    ialias(n) = 70
    spgalias(n) = 'pmm2'
    n = n + 1
    ialias(n) = 71
    spgalias(n) = 'pmc21'
    n = n + 1
    ialias(n) = 72
    spgalias(n) = 'pcc2'
    n = n + 1
    ialias(n) = 73
    spgalias(n) = 'pma2'
    n = n + 1
    ialias(n) = 74
    spgalias(n) = 'pca21'
    n = n + 1
    ialias(n) = 75
    spgalias(n) = 'pnc2'
    n = n + 1
    ialias(n) = 76
    spgalias(n) = 'pmn21'
    n = n + 1
    ialias(n) = 77
    spgalias(n) = 'pba2'
    n = n + 1
    ialias(n) = 78
    spgalias(n) = 'pna21'
    n = n + 1
    ialias(n) = 79
    spgalias(n) = 'pnn2'
    n = n + 1
    ialias(n) = 80
    spgalias(n) = 'cmm2'
    n = n + 1
    ialias(n) = 81
    spgalias(n) = 'cmc21'
    n = n + 1
    ialias(n) = 82
    spgalias(n) = 'ccc2'
    n = n + 1
    ialias(n) = 83
    spgalias(n) = 'amm2'
    n = n + 1
    ialias(n) = 84
    spgalias(n) = 'abm2'
    n = n + 1
    ialias(n) = 85
    spgalias(n) = 'ama2'
    n = n + 1
    ialias(n) = 86
    spgalias(n) = 'aba2'
    n = n + 1
    ialias(n) = 87
    spgalias(n) = 'fmm2'
    n = n + 1
    ialias(n) = 88
    spgalias(n) = 'fdd2'
    n = n + 1
    ialias(n) = 89
    spgalias(n) = 'imm2'
    n = n + 1
    ialias(n) = 90
    spgalias(n) = 'iba2'
    n = n + 1
    ialias(n) = 91
    spgalias(n) = 'ima2'
    n = n + 1
    ialias(n) = 92
    spgalias(n) = 'pmmm'
    n = n + 1
    ialias(n) = 93
    spgalias(n) = 'pnnn1'
    n = n + 1
    ialias(n) = 94
    spgalias(n) = 'pnnn2'
    n = n + 1
    ialias(n) = 95
    spgalias(n) = 'pccm'
    n = n + 1
    ialias(n) = 96
    spgalias(n) = 'pban1'
    n = n + 1
    ialias(n) = 97
    spgalias(n) = 'pban2'
    n = n + 1
    ialias(n) = 98
    spgalias(n) = 'pmma'
    n = n + 1
    ialias(n) = 99
    spgalias(n) = 'pnna'
    n = n + 1
    ialias(n) = 100
    spgalias(n) = 'pmna'
    n = n + 1
    ialias(n) = 101
    spgalias(n) = 'pcca'
    n = n + 1
    ialias(n) = 102
    spgalias(n) = 'pbam'
    n = n + 1
    ialias(n) = 103
    spgalias(n) = 'pccn'
    n = n + 1
    ialias(n) = 104
    spgalias(n) = 'pbcm'
    n = n + 1
    ialias(n) = 105
    spgalias(n) = 'pnnm'
    n = n + 1
    ialias(n) = 106
    spgalias(n) = 'pmmn1'
    n = n + 1
    ialias(n) = 107
    spgalias(n) = 'pmmn2'
    n = n + 1
    ialias(n) = 108
    spgalias(n) = 'pbcn'
    n = n + 1
    ialias(n) = 109
    spgalias(n) = 'pbca'
    n = n + 1
    ialias(n) = 110
    spgalias(n) = 'pnma'
    n = n + 1
    ialias(n) = 111
    spgalias(n) = 'cmcm'
    n = n + 1
    ialias(n) = 112
    spgalias(n) = 'cmca'
    n = n + 1
    ialias(n) = 113
    spgalias(n) = 'cmmm'
    n = n + 1
    ialias(n) = 114
    spgalias(n) = 'cccm'
    n = n + 1
    ialias(n) = 115
    spgalias(n) = 'cmma'
    n = n + 1
    ialias(n) = 116
    spgalias(n) = 'ccca1'
    n = n + 1
    ialias(n) = 117
    spgalias(n) = 'ccca2'
    n = n + 1
    ialias(n) = 118
    spgalias(n) = 'fmmm'
    n = n + 1
    ialias(n) = 119
    spgalias(n) = 'fddd1'
    n = n + 1
    ialias(n) = 120
    spgalias(n) = 'fddd2'
    n = n + 1
    ialias(n) = 121
    spgalias(n) = 'immm'
    n = n + 1
    ialias(n) = 122
    spgalias(n) = 'ibam'
    n = n + 1
    ialias(n) = 123
    spgalias(n) = 'ibca'
    n = n + 1
    ialias(n) = 124
    spgalias(n) = 'imma'
    n = n + 1
    ialias(n) = 125
    spgalias(n) = 'p4'
    n = n + 1
    ialias(n) = 126
    spgalias(n) = 'p41'
    n = n + 1
    ialias(n) = 127
    spgalias(n) = 'p42'
    n = n + 1
    ialias(n) = 128
    spgalias(n) = 'p43'
    n = n + 1
    ialias(n) = 129
    spgalias(n) = 'i4'
    n = n + 1
    ialias(n) = 130
    spgalias(n) = 'i41'
    n = n + 1
    ialias(n) = 131
    spgalias(n) = 'p-4'
    n = n + 1
    ialias(n) = 132
    spgalias(n) = 'i-4'
    n = n + 1
    ialias(n) = 133
    spgalias(n) = 'p4/m'
    n = n + 1
    ialias(n) = 134
    spgalias(n) = 'p42/m'
    n = n + 1
    ialias(n) = 135
    spgalias(n) = 'p4/n1'
    n = n + 1
    ialias(n) = 136
    spgalias(n) = 'p4/n2'
    n = n + 1
    ialias(n) = 137
    spgalias(n) = 'p42/n1'
    n = n + 1
    ialias(n) = 138
    spgalias(n) = 'p42/n2'
    n = n + 1
    ialias(n) = 139
    spgalias(n) = 'i4/m'
    n = n + 1
    ialias(n) = 140
    spgalias(n) = 'i41/a1'
    n = n + 1
    ialias(n) = 141
    spgalias(n) = 'i41/a2'
    n = n + 1
    ialias(n) = 142
    spgalias(n) = 'p422'
    n = n + 1
    ialias(n) = 143
    spgalias(n) = 'p4212'
    n = n + 1
    ialias(n) = 144
    spgalias(n) = 'p4122'
    n = n + 1
    ialias(n) = 145
    spgalias(n) = 'p41212'
    n = n + 1
    ialias(n) = 146
    spgalias(n) = 'p4222'
    n = n + 1
    ialias(n) = 147
    spgalias(n) = 'p42212'
    n = n + 1
    ialias(n) = 148
    spgalias(n) = 'p4322'
    n = n + 1
    ialias(n) = 149
    spgalias(n) = 'p43212'
    n = n + 1
    ialias(n) = 150
    spgalias(n) = 'i422'
    n = n + 1
    ialias(n) = 151
    spgalias(n) = 'i4122'
    n = n + 1
    ialias(n) = 152
    spgalias(n) = 'p4mm'
    n = n + 1
    ialias(n) = 153
    spgalias(n) = 'p4bm'
    n = n + 1
    ialias(n) = 154
    spgalias(n) = 'p42cm'
    n = n + 1
    ialias(n) = 155
    spgalias(n) = 'p42nm'
    n = n + 1
    ialias(n) = 156
    spgalias(n) = 'p4cc'
    n = n + 1
    ialias(n) = 157
    spgalias(n) = 'p4nc'
    n = n + 1
    ialias(n) = 158
    spgalias(n) = 'p42mc'
    n = n + 1
    ialias(n) = 159
    spgalias(n) = 'p42bc'
    n = n + 1
    ialias(n) = 160
    spgalias(n) = 'i4mm'
    n = n + 1
    ialias(n) = 161
    spgalias(n) = 'i4cm'
    n = n + 1
    ialias(n) = 162
    spgalias(n) = 'i41md'
    n = n + 1
    ialias(n) = 163
    spgalias(n) = 'i41cd'
    n = n + 1
    ialias(n) = 164
    spgalias(n) = 'p-42m'
    n = n + 1
    ialias(n) = 165
    spgalias(n) = 'p-42c'
    n = n + 1
    ialias(n) = 166
    spgalias(n) = 'p-421m'
    n = n + 1
    ialias(n) = 167
    spgalias(n) = 'p-421c'
    n = n + 1
    ialias(n) = 168
    spgalias(n) = 'p-4m2'
    n = n + 1
    ialias(n) = 169
    spgalias(n) = 'p-4c2'
    n = n + 1
    ialias(n) = 170
    spgalias(n) = 'p-4b2'
    n = n + 1
    ialias(n) = 171
    spgalias(n) = 'p-4n2'
    n = n + 1
    ialias(n) = 172
    spgalias(n) = 'i-4m2'
    n = n + 1
    ialias(n) = 173
    spgalias(n) = 'i-4c2'
    n = n + 1
    ialias(n) = 174
    spgalias(n) = 'i-42m'
    n = n + 1
    ialias(n) = 175
    spgalias(n) = 'i-42d'
    n = n + 1
    ialias(n) = 176
    spgalias(n) = 'p4/mmm'
    n = n + 1
    ialias(n) = 177
    spgalias(n) = 'p4/mcc'
    n = n + 1
    ialias(n) = 178
    spgalias(n) = 'p4/nbm1'
    n = n + 1
    ialias(n) = 179
    spgalias(n) = 'p4/nbm2'
    n = n + 1
    ialias(n) = 180
    spgalias(n) = 'p4/nnc1'
    n = n + 1
    ialias(n) = 181
    spgalias(n) = 'p4/nnc2'
    n = n + 1
    ialias(n) = 182
    spgalias(n) = 'p4/mbm'
    n = n + 1
    ialias(n) = 183
    spgalias(n) = 'p4/mnc'
    n = n + 1
    ialias(n) = 184
    spgalias(n) = 'p4/nmm1'
    n = n + 1
    ialias(n) = 185
    spgalias(n) = 'p4/nmm2'
    n = n + 1
    ialias(n) = 186
    spgalias(n) = 'p4/ncc1'
    n = n + 1
    ialias(n) = 187
    spgalias(n) = 'p4/ncc2'
    n = n + 1
    ialias(n) = 188
    spgalias(n) = 'p42/mmc'
    n = n + 1
    ialias(n) = 189
    spgalias(n) = 'p42/mcm'
    n = n + 1
    ialias(n) = 190
    spgalias(n) = 'p42/nbc1'
    n = n + 1
    ialias(n) = 191
    spgalias(n) = 'p42/nbc2'
    n = n + 1
    ialias(n) = 192
    spgalias(n) = 'p42/nnm1'
    n = n + 1
    ialias(n) = 193
    spgalias(n) = 'p42/nnm2'
    n = n + 1
    ialias(n) = 194
    spgalias(n) = 'p42/mbc'
    n = n + 1
    ialias(n) = 195
    spgalias(n) = 'p42/mnm'
    n = n + 1
    ialias(n) = 196
    spgalias(n) = 'p42/nmc1'
    n = n + 1
    ialias(n) = 197
    spgalias(n) = 'p42/nmc2'
    n = n + 1
    ialias(n) = 198
    spgalias(n) = 'p42/ncm1'
    n = n + 1
    ialias(n) = 199
    spgalias(n) = 'p42/nmc2'
    n = n + 1
    ialias(n) = 200
    spgalias(n) = 'i4/mmm'
    n = n + 1
    ialias(n) = 201
    spgalias(n) = 'i4/mcm'
    n = n + 1
    ialias(n) = 202
    spgalias(n) = 'i41/amd1'
    n = n + 1
    ialias(n) = 203
    spgalias(n) = 'i41/amd2'
    n = n + 1
    ialias(n) = 204
    spgalias(n) = 'i41/acd1'
    n = n + 1
    ialias(n) = 205
    spgalias(n) = 'i41/acd2'
    n = n + 1
    ialias(n) = 206
    spgalias(n) = 'p3'
    n = n + 1
    ialias(n) = 207
    spgalias(n) = 'p31'
    n = n + 1
    ialias(n) = 208
    spgalias(n) = 'p32'
    n = n + 1
    ialias(n) = 209
    spgalias(n) = 'r3h'
    n = n + 1
    ialias(n) = 210
    spgalias(n) = 'r3r'
    n = n + 1
    ialias(n) = 211
    spgalias(n) = 'p-3'
    n = n + 1
    ialias(n) = 212
    spgalias(n) = 'r-3h'
    n = n + 1
    ialias(n) = 213
    spgalias(n) = 'r-3r'
    n = n + 1
    ialias(n) = 214
    spgalias(n) = 'p312'
    n = n + 1
    ialias(n) = 215
    spgalias(n) = 'p321'
    n = n + 1
    ialias(n) = 216
    spgalias(n) = 'p3112'
    n = n + 1
    ialias(n) = 217
    spgalias(n) = 'p3121'
    n = n + 1
    ialias(n) = 218
    spgalias(n) = 'p3212'
    n = n + 1
    ialias(n) = 219
    spgalias(n) = 'p3221'
    n = n + 1
    ialias(n) = 220
    spgalias(n) = 'r32h'
    n = n + 1
    ialias(n) = 221
    spgalias(n) = 'r32r'
    n = n + 1
    ialias(n) = 222
    spgalias(n) = 'p3m1'
    n = n + 1
    ialias(n) = 223
    spgalias(n) = 'p31m'
    n = n + 1
    ialias(n) = 224
    spgalias(n) = 'p3c1'
    n = n + 1
    ialias(n) = 225
    spgalias(n) = 'p31c'
    n = n + 1
    ialias(n) = 226
    spgalias(n) = 'r3mh'
    n = n + 1
    ialias(n) = 227
    spgalias(n) = 'r3mr'
    n = n + 1
    ialias(n) = 228
    spgalias(n) = 'r3ch'
    n = n + 1
    ialias(n) = 229
    spgalias(n) = 'r3cr'
    n = n + 1
    ialias(n) = 230
    spgalias(n) = 'p-31m'
    n = n + 1
    ialias(n) = 231
    spgalias(n) = 'p-31c'
    n = n + 1
    ialias(n) = 232
    spgalias(n) = 'p-3m1'
    n = n + 1
    ialias(n) = 233
    spgalias(n) = 'p-3c1'
    n = n + 1
    ialias(n) = 234
    spgalias(n) = 'r-3mh'
    n = n + 1
    ialias(n) = 235
    spgalias(n) = 'r-3mr'
    n = n + 1
    ialias(n) = 236
    spgalias(n) = 'r-3ch'
    n = n + 1
    ialias(n) = 237
    spgalias(n) = 'r-3cr'
    n = n + 1
    ialias(n) = 238
    spgalias(n) = 'p6'
    n = n + 1
    ialias(n) = 239
    spgalias(n) = 'p61'
    n = n + 1
    ialias(n) = 240
    spgalias(n) = 'p65'
    n = n + 1
    ialias(n) = 241
    spgalias(n) = 'p62'
    n = n + 1
    ialias(n) = 242
    spgalias(n) = 'p64'
    n = n + 1
    ialias(n) = 243
    spgalias(n) = 'p63'
    n = n + 1
    ialias(n) = 244
    spgalias(n) = 'p-6'
    n = n + 1
    ialias(n) = 245
    spgalias(n) = 'p6/m'
    n = n + 1
    ialias(n) = 246
    spgalias(n) = 'p63/m'
    n = n + 1
    ialias(n) = 247
    spgalias(n) = 'p622'
    n = n + 1
    ialias(n) = 248
    spgalias(n) = 'p6122'
    n = n + 1
    ialias(n) = 249
    spgalias(n) = 'p6522'
    n = n + 1
    ialias(n) = 250
    spgalias(n) = 'p6222'
    n = n + 1
    ialias(n) = 251
    spgalias(n) = 'p6422'
    n = n + 1
    ialias(n) = 252
    spgalias(n) = 'p6322'
    n = n + 1
    ialias(n) = 253
    spgalias(n) = 'p6mm'
    n = n + 1
    ialias(n) = 254
    spgalias(n) = 'p6cc'
    n = n + 1
    ialias(n) = 255
    spgalias(n) = 'p63cm'
    n = n + 1
    ialias(n) = 256
    spgalias(n) = 'p63mc'
    n = n + 1
    ialias(n) = 257
    spgalias(n) = 'p-6m2'
    n = n + 1
    ialias(n) = 258
    spgalias(n) = 'p-6c2'
    n = n + 1
    ialias(n) = 259
    spgalias(n) = 'p-62m'
    n = n + 1
    ialias(n) = 260
    spgalias(n) = 'p-62c'
    n = n + 1
    ialias(n) = 261
    spgalias(n) = 'p6/mmm'
    n = n + 1
    ialias(n) = 262
    spgalias(n) = 'p6/mcc'
    n = n + 1
    ialias(n) = 263
    spgalias(n) = 'p63/mcm'
    n = n + 1
    ialias(n) = 264
    spgalias(n) = 'p63/mmc'
    n = n + 1
    ialias(n) = 265
    spgalias(n) = 'p23'
    n = n + 1
    ialias(n) = 266
    spgalias(n) = 'f23'
    n = n + 1
    ialias(n) = 267
    spgalias(n) = 'i23'
    n = n + 1
    ialias(n) = 268
    spgalias(n) = 'p213'
    n = n + 1
    ialias(n) = 269
    spgalias(n) = 'i213'
    n = n + 1
    ialias(n) = 270
    spgalias(n) = 'pm-3'
    n = n + 1
    ialias(n) = 271
    spgalias(n) = 'pn-31'
    n = n + 1
    ialias(n) = 272
    spgalias(n) = 'pn-32'
    n = n + 1
    ialias(n) = 273
    spgalias(n) = 'fm-3'
    n = n + 1
    ialias(n) = 274
    spgalias(n) = 'fd-31'
    n = n + 1
    ialias(n) = 275
    spgalias(n) = 'fd-32'
    n = n + 1
    ialias(n) = 276
    spgalias(n) = 'im-3'
    n = n + 1
    ialias(n) = 277
    spgalias(n) = 'pa-3'
    n = n + 1
    ialias(n) = 278
    spgalias(n) = 'ia-3'
    n = n + 1
    ialias(n) = 279
    spgalias(n) = 'p432'
    n = n + 1
    ialias(n) = 280
    spgalias(n) = 'p4232'
    n = n + 1
    ialias(n) = 281
    spgalias(n) = 'f432'
    n = n + 1
    ialias(n) = 282
    spgalias(n) = 'f4132'
    n = n + 1
    ialias(n) = 283
    spgalias(n) = 'i432'
    n = n + 1
    ialias(n) = 284
    spgalias(n) = 'p4332'
    n = n + 1
    ialias(n) = 285
    spgalias(n) = 'p4132'
    n = n + 1
    ialias(n) = 286
    spgalias(n) = 'i4132'
    n = n + 1
    ialias(n) = 287
    spgalias(n) = 'p-43m'
    n = n + 1
    ialias(n) = 288
    spgalias(n) = 'f-43m'
    n = n + 1
    ialias(n) = 289
    spgalias(n) = 'i-43m'
    n = n + 1
    ialias(n) = 290
    spgalias(n) = 'p-43n'
    n = n + 1
    ialias(n) = 291
    spgalias(n) = 'f-43c'
    n = n + 1
    ialias(n) = 292
    spgalias(n) = 'i-43d'
    n = n + 1
    ialias(n) = 293
    spgalias(n) = 'pm-3m'
    n = n + 1
    ialias(n) = 294
    spgalias(n) = 'pn-3n1'
    n = n + 1
    ialias(n) = 295
    spgalias(n) = 'pn-3n2'
    n = n + 1
    ialias(n) = 296
    spgalias(n) = 'pm-3n'
    n = n + 1
    ialias(n) = 297
    spgalias(n) = 'pn-3m1'
    n = n + 1
    ialias(n) = 298
    spgalias(n) = 'pn-3m2'
    n = n + 1
    ialias(n) = 299
    spgalias(n) = 'fm-3m'
    n = n + 1
    ialias(n) = 300
    spgalias(n) = 'fm-3c'
    n = n + 1
    ialias(n) = 301
    spgalias(n) = 'fd-3m1'
    n = n + 1
    ialias(n) = 302
    spgalias(n) = 'fd-3m2'
    n = n + 1
    ialias(n) = 303
    spgalias(n) = 'fd-3c1'
    n = n + 1
    ialias(n) = 304
    spgalias(n) = 'fd-3c2'
    n = n + 1
    ialias(n) = 305
    spgalias(n) = 'im-3m'
    n = n + 1
    ialias(n) = 306
    spgalias(n) = 'ia-3d'
    n = n + 1
    ialias(n) =  3
    spgalias(n)="p2"
    n = n + 1
    ialias(n) =  3
    spgalias(n)="p2b"
    n = n + 1
    ialias(n) =  4
    spgalias(n)="p2c"
    n = n + 1
    ialias(n) =  5
    spgalias(n)="p21"
    n = n + 1
    ialias(n) =  5
    spgalias(n)="p21b"
    n = n + 1
    ialias(n) =  6
    spgalias(n)="p21c"
    n = n + 1
    ialias(n) =  7
    spgalias(n)="c2"
    n = n + 1
    ialias(n) =  7
    spgalias(n)="c2b"
    n = n + 1
    ialias(n) = 13
    spgalias(n)="pm"
    n = n + 1
    ialias(n) = 13
    spgalias(n)="pmb"
    n = n + 1
    ialias(n) = 14
    spgalias(n)="pmc"
    n = n + 1
    ialias(n) = 15
    spgalias(n)="pc"
    n = n + 1
    ialias(n) = 15
    spgalias(n)="pcb"
    n = n + 1
    ialias(n) = 21
    spgalias(n)="cm"
    n = n + 1
    ialias(n) = 21
    spgalias(n)="cmb"
    n = n + 1
    ialias(n) = 27
    spgalias(n)="cc"
    n = n + 1
    ialias(n) = 27
    spgalias(n)="ccb"
    n = n + 1
    ialias(n) = 33
    spgalias(n)="p2/m"
    n = n + 1
    ialias(n) = 33
    spgalias(n)="p2/mb"
    n = n + 1
    ialias(n) = 34
    spgalias(n)="p2/mc"
    n = n + 1
    ialias(n) = 35
    spgalias(n)="p21/m"
    n = n + 1
    ialias(n) = 35
    spgalias(n)="p21/mb"
    n = n + 1
    ialias(n) = 36
    spgalias(n)="p21/mc"
    n = n + 1
    ialias(n) = 37
    spgalias(n)="c2/m"
    n = n + 1
    ialias(n) = 37
    spgalias(n)="c2/mb"
    n = n + 1
    ialias(n) = 43
    spgalias(n)="p2/c"
    n = n + 1
    ialias(n) = 43
    spgalias(n)="p2/cb"
    n = n + 1
    ialias(n) = 49
    spgalias(n)="p21/c"
    n = n + 1
    ialias(n) = 49
    spgalias(n)="p21/cb"
    n = n + 1
    ialias(n) = 55
    spgalias(n)="c2/c"
    n = n + 1
    ialias(n) = 55
    spgalias(n)="c2/cb"
    n = n + 1
    ialias(n) = 84
    spgalias(n)="aem2"
    n = n + 1
    ialias(n) = 86
    spgalias(n)="aea2"
    n = n + 1
    ialias(n) = 93
    spgalias(n)="pnnn"
    n = n + 1
    ialias(n) = 95
    spgalias(n)="pccm"
    n = n + 1
    ialias(n) = 96
    spgalias(n)="pban"
    n = n + 1
    ialias(n) =106
    spgalias(n)="pmmn"
    n = n + 1
    ialias(n) =116
    spgalias(n)="ccca"
    n = n + 1
    ialias(n) =119
    spgalias(n)="fddd"
    n = n + 1
    ialias(n) =135
    spgalias(n)="p4/n"
    n = n + 1
    ialias(n) =137
    spgalias(n)="p42/n"
    n = n + 1
    ialias(n) =140
    spgalias(n)="i41/a"
    n = n + 1
    ialias(n) =178
    spgalias(n)="p4/nbm"
    n = n + 1
    ialias(n) =180
    spgalias(n)="p4/nnc"
    n = n + 1
    ialias(n) =184
    spgalias(n)="p4/nmm"
    n = n + 1
    ialias(n) =186
    spgalias(n)="p4/ncc"
    n = n + 1
    ialias(n) =190
    spgalias(n)="p42/nbc"
    n = n + 1
    ialias(n) =192
    spgalias(n)="p42/nnm"
    n = n + 1
    ialias(n) =196
    spgalias(n)="p42/nmc"
    n = n + 1
    ialias(n) =198
    spgalias(n)="p42/ncm"
    n = n + 1
    ialias(n) =202
    spgalias(n)="i41/amd"
    n = n + 1
    ialias(n) =204
    spgalias(n)="i41/acd"
    n = n + 1
    ialias(n) =209
    spgalias(n)="r3h"
    n = n + 1
    ialias(n) =212
    spgalias(n)="r-3"
    n = n + 1
    ialias(n) =220
    spgalias(n)="r32"
    n = n + 1
    ialias(n) =226
    spgalias(n)="r3m"
    n = n + 1
    ialias(n) =228
    spgalias(n)="r3c"
    n = n + 1
    ialias(n) =234
    spgalias(n)="r-3m"
    n = n + 1
    ialias(n) =236
    spgalias(n)="r-3c"
    n = n + 1
    ialias(n) =271
    spgalias(n)="pn-3"
    n = n + 1
    ialias(n) =274
    spgalias(n)="fd-3"
    n = n + 1
    ialias(n) =294
    spgalias(n)="pn-3n"
    n = n + 1
    ialias(n) =297
    spgalias(n)="pn-3m"
    n = n + 1
    ialias(n) =301
    spgalias(n)="fd-3m"
    n = n + 1
    ialias(n) =303
    spgalias(n)="fd-3c"

    ! additional alias
    n = n + 1
    ialias(n) =301
    spgalias(n)="f d 3 m"
    n = n + 1
    ialias(n) =301
    spgalias(n)="fd3m"
    n = n + 1
    ialias(n) =299
    spgalias(n)="f m 3 m"
    n = n + 1
    ialias(n) =299
    spgalias(n)="fm3m"
    nalias = n

  end subroutine spgs_init

  subroutine spgr(spg, lpt, lptx, ier)
    use struct
    use tools_io
    use param
!-----------------------------------------------------------------------
!
!     SSSSSSSS      PPPPPPPPPP        GGGGGGGG      RRRRRRRRRR
!    SSSSSSSSSS     PPPPPPPPPPP      GGGGGGGGGG     RRRRRRRRRRR
!   SSS             PPP      PPP    GGG      GGG    RRR      RRR
!   SSS             PPP      PPP    GGG      GGG    RRR      RRR
!    SSSSSSSSS      PPP      PPP    GGG             RRR      RRR
!     SSSSSSSSS     PPP      PPP    GGG             RRR      RRR
!            SSS    PPPPPPPPPPP     GGG   GGGGGG    RRRRRRRRRRR
!            SSS    PPPPPPPPPP      GGG   GGGGGG    RRRRRRRRRR
!   SSS      SSS    PPP             GGG      GGG    RRR    RRR
!   SSS      SSS    PPP             GGG      GGG    RRR     RRR
!    SSSSSSSSSS     PPP              GGGGGGGGGG     RRR      RRR
!     SSSSSSSS      PPP               GGGGGGGG      RRR       RRR
!
!.....Space group symbol recognition library.
!
!.....Based upon work by Allen C. Larson (the National Research
!     Council of Canada) in the sixties.
!
!.....Rewritten in Fortran 77, cleaned up and corrected by Angel Martin
!     Pendas (Universidad de Oviedo).
!     Adapted to REAL*8 by Miguel Alvarez Blanco (Univ. de Oviedo)
!
!----------------------------------------------------------------------
!
!
!.....spgr - interprets the space group symbol
!
!     data in the calling sequence are
!
!     spg ...... space group symbol
!     laue ..... output the laue group number.
!                1=1bar, 2=2/m, 3=mmm, 4=4/m, 5=4/mm, 6=r3r, 7=r3mr,
!                8=3, 9=3m1, 10=31m, 11=6/m, 12=6/mmm, 13=m3 AND 14=m3m
!     naxis .... unique axis in monoclinic space groups
!     ncent .... 1bar flag  (0/1) for (acentric/centric)
!     lcent .... lattice centering number.
!                1=P, 2=A, 3=B, 4=C, 5=I, 6=F and 7=R
!     neqv ..... Number of symmetry matrices generated.
!     jrt ...... neqv (3,4,neqv) matrices
!     cen ...... lattice centering vectors
!     ncv ...... Number of lattice centering vectors.
!     lpt ...... stdout
!                if (lpt.lt.0) no listing will be produced
!     lptx ..... stderr
!                if (lptx.lt.0) no error report will be produced
!     ier ...... if (ier.ne.0). An error has been found.
!
!
      implicit none

      character*(24)    spg
      integer           lpt, lptx, ier
!
      character*(1)     chr(25)
      real*8            rt(5,4,25), d(3,3)
      integer           l(4,4)
      integer           lcen(7)
      integer           length
      integer           i, j, k, m, ierx, n, nromb, ncubic, ifound
      integer           nmenos, id, nxi, njs, m2, localident, nxl, nx
      integer           m1
!
!.....data in lcen for  C B A P F I R  respectively
!
      data&
       lcen(1),lcen(2),lcen(3),lcen(4),lcen(5) /4,3,2,1,6/,&
       lcen(6),lcen(7)                         /5,7/
!
!                  1    2    3    4    5    6    7    8    9    10
      data chr  / ' ', 'C', 'B', 'A', 'P', 'F', 'I', 'R', 'M', 'N',&
                  'D', '1', '2', '3', '4', '5', '6', '-', '/', 'H',&
                  '.', ' ', ' ', ' ', ' ' /
!                  11   12   13   14   15   16   17   18   19   20
!
      do i = 1, 4
         do j = 1, 4
            l(i,j) = 0
         enddo
      enddo
      ier    = 0
      ncent  = 0
      laue   = 0
      naxis  = 0
      ierx   = 0
      n      = 0
      nromb  = 0
      ncubic = 0
!
!.....transform symbol to uppercase
!
      spg = upper(adjustl(spg))
      length = len_trim(spg)
!
!.....break the space group symbol into the 4 fields and code
!     as numerical values for manipulation
      k = 0
      m = 0
      do 110 j=1,length
         ifound=0
         do 100 i=1,21
!...........Exit if found
            if ( spg(j:j).eq.chr(i) ) then 
                ifound=i
                go to 101
            endif
 100     continue
 101     if (ifound.ne.0) then 
            if (k+m+ifound.ne.1) then
                if ( ifound.ne.1) then 
                   if ( m.eq.0 ) k=k+1
                   m = m+1
                   l(m,k) = ifound
                else
                   m=0
!..................exit loop if more than 3 fields read
                   if (k.gt.3) goto 200
                endif
                if ( ifound.lt.12 .or. m.ge.4) then
                   m=0
!..................exit loop if more than 3 fields read
                   if (k.gt.3) goto 200
                endif
            endif
         endif      
 110  continue
!
!
  200 if ( k.le.1 .or. l(1,1).gt.8 ) then
         ier=1
         if (l(1,1).gt.8) ier=2
         if (lptx.ge.0) call sgerrs (spg(1:24), ier, lptx)
         return
      endif
!
!.....convert the -N notation to the nb(ar) notation
!
      nmenos=0
      if ( l(1,2).eq.18 ) then
         nmenos=nmenos+1
         call sglpak(l(1,2),ier)
         if (ier .gt. 0 ) then
            if (lptx.ge.0) call sgerrs (spg(1:24), ier, lptx)
            return
         endif
      endif
      if ( l(1,3).eq.18 ) then
         nmenos=nmenos+1
         call sglpak(l(1,3),ier)
         if (ier .gt. 0 ) then
            if (lptx.ge.0) call sgerrs (spg(1:24), ier, lptx)
            return
         endif
      endif
      if ( l(1,4).eq.18 ) then
         nmenos=nmenos+1
         call sglpak(l(1,4),ier)
         if (ier .gt. 0 ) then
            if (lptx.ge.0) call sgerrs (spg(1:24), ier, lptx)
            return
         endif
      endif
!
!.....only one bar allowed
!
      if ( nmenos.gt.1 ) then
         ier=24
         if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
         return
      endif
!
!.....set the matrix count n to 2
!
      n = 2
!
!.....set the translation flags
!
      d(1,1) = zero
      d(1,2) = zero
      d(1,3) = zero
      d(2,1) = zero
      d(2,2) = zero
      d(2,3) = zero
      d(3,1) = zero
      d(3,2) = zero
      d(3,3) = zero
!
!.....set the lattice centering flag. 1=P, 2=A, 3=B, 4=C, 5=I, 6=F, 7=R
!
      lcent = l(1,1)-1
      lcent = lcen(lcent)
      if ( lcent.eq.7 ) then
!
!........rhombohedral lattice.
!        make sure that there is a 3-axis.
!
         if ( l(1,2).ne.14 ) then
            ier=3
            if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
            return
         endif
         if ( l(1,k).ne.8 ) then
!
!...........hexagonal axes. retain r centering and set laue to 8 or 9
!
            if ( l(1,k).eq.20 ) k=k-1
            laue = k+6
         else
!
!...........rhombohedral axes. delete R centering.set laue to 6 or 7
!
            lcent = 1
            k = k-1
            laue = k+4
            nromb=1
         endif
      else
!
!........determine laue and some preliminary data
!
         ier = 0
         ncubic = 0
         call sglatc(k, l, d, lcent, laue, naxis, ier, ncubic, id)
         if ( ier.gt.0 ) then
            if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
            return
         endif
      endif
!
      if ( ncubic.eq.1 .or. nromb.eq.1 ) then
!
!.........cubic or rhombohedral cell. insert the 3-fold axis
!
          call sgrmat(rt(1,1,2),0,1,0,0,0,1,1,0,0)
          call sgrmat(rt(1,1,3),0,0,1,1,0,0,0,1,0)
          n = 4
      endif
!
!.....insert the unit matrix
!
      call sgrmat(rt,1,0,0,0,1,0,0,0,1)
!
!.....decode the last 3 fields of the symbol
!
      do 20 m=2,k
         if ( l(1,m).eq.0 ) then
            ier=6
            if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
            return
         endif
         i = iabs(l(1,m)-5)
  10     if ( i.le.0.or.i.gt.15 ) then
            ier=7
            if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
            return
         endif
         nxi = n
!
!       A   B   C   M   N   D   1   2   3   4   5   6   -   /   H
!  i=   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
!
         if (i.ge.1 .and. i.le.5) then
!
!...........a mirror is needed
!
            if (m .eq. 2) then
!
!         A   B   C  AXIS
!  m=     1   2   3
!
               if ( laue.gt.3 )  then
                  if (i.eq.3 ) then
                      ier=10
                      if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
                      return
                  endif
!
!.................a C-axis mirror
!
                  call sgrmat(rt(1,1,n),1,0,0,0,1,0,0,0,-1)
                  rt(3,4,n) = d(3,3)
                  if ( i.eq.1.or.i.eq.5 ) rt(1,4,n) = half 
                  if ( i.eq.2.or.i.eq.5 ) rt(2,4,n) = half
                  if ( m.eq. 2.and.l(1,2).eq.17) then
!
!....................if this is a 63-axis the mirror is at 1/4
!
                     if ( l(2,2).eq.14 ) rt(3,4,n)=half
                  endif
               else if (k.eq.2) then
                  if (i.eq.2 ) then
                     ier=9
                     if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
                     return
                  endif
!
!.................a B-axis mirror
!
                  call sgrmat(rt(1,1,n),1,0,0,0,-1,0,0,0,1)
                  rt(2,4,n) = d(2,2)
                  if ( i.eq.1.or.i.eq.5 ) rt(1,4,n) = half
                  if ( i.eq.3.or.i.eq.5 ) rt(3,4,n) = half
               else if (i.eq.1) then
                  ier=8
                  if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
                  return
               else
!
!.................a A-axis mirror
!
                  call sgrmat(rt(1,1,n),-1,0,0,0,1,0,0,0,1)
                  rt(1,4,n) = d(1,1)
                  if ( i.eq.2.or.i.eq.5 ) rt(2,4,n)=half
                  if ( i.eq.3.or.i.eq.5 ) rt(3,4,n)=half
               endif
            else if (m.eq.3) then
               if ( l(1,2).eq.14.or.l(1,2).eq.17 ) then
!
                  if ( laue.eq.7 ) then
!
!....................a diagonal mirrror perpendicular to -110
!
                     call sgrmat(rt(1,1,n),0,1,0,1,0,0,0,0,1)
                     rt(1,4,n) = d(2,2)
                     rt(2,4,n) = -d(2,2)
                     if ( i.eq.3.or.i.eq.5 ) rt(3,4,n) = half
                     if (( laue.ne.7 .or. i.ne.3) .and. (  i.ne.3 .and. i.ne.4  )) then
                        if ( lcent.eq.6.or.lcent.eq.4 ) then
!
!.......................F or C-centered tetragonal. glides are 1/4,1/4
!
                           rt(1,4,n) = fourth+rt(1,4,n)
                           rt(2,4,n) = fourth+rt(2,4,n)
                        else
                           rt(1,4,n) = half+rt(1,4,n)
                           rt(2,4,n) = half+rt(2,4,n)
                        endif
                     endif
                  else
!
!....................mirror normal to (1000) in hex cell
!
                     call sgrmat(rt(1,1,n),-1,1,0,0,1,0,0,0,1)
                     if ( i.eq.3 ) rt(3,4,n)=half
                  endif
               else if (l(1,2).eq.15) then
                  if (i.eq.1) then
                     ier=8
                     if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
                     return
                  else
!
!...................a A-axis mirror
!
                    call sgrmat(rt(1,1,n),-1,0,0,0,1,0,0,0,1)
                    rt(1,4,n) = d(1,1)
                    if ( i.eq.2.or.i.eq.5 ) rt(2,4,n)=half
                    if ( i.eq.3.or.i.eq.5 ) rt(3,4,n)=half
                  endif
               else
                  if (i.eq.2 ) then
                     ier=9
                     if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
                     return
                  endif
!
!.................a B-axis mirror
!
                  call sgrmat(rt(1,1,n),1,0,0,0,-1,0,0,0,1)
                  rt(2,4,n) = d(2,2)
                  if ( i.eq.1.or.i.eq.5 ) rt(1,4,n) = half
                  if ( i.eq.3.or.i.eq.5 ) rt(3,4,n) = half
               endif
            else if (m.eq.4) then
!
!..............it is not cubic or tetragonal nor trigonal or hexagonal
!
               if (( l(1,3).eq.14.or.l(1,2).eq.15 ) .or. ( l(1,2).eq.14.or.l(1,2).eq.17 )) then
!
!.................a diagonal mirrror perpendicular to -110
!
                  call sgrmat(rt(1,1,n),0,1,0,1,0,0,0,0,1)
                  rt(1,4,n) = d(2,2)
                  rt(2,4,n) = -d(2,2)
                  if ( i.eq.3.or.i.eq.5 ) rt(3,4,n) = half
                  if (( laue.ne.7 .or. i.ne.3) .and. (  i.ne.3 .and. i.ne.4  )) then
                     if ( lcent.eq.6.or.lcent.eq.4 ) then
!
!....................F or C-centered tetragonal. glides are 1/4,1/4
!
                        rt(1,4,n) = fourth+rt(1,4,n)
                        rt(2,4,n) = fourth+rt(2,4,n)
                     else
                        rt(1,4,n) = half+rt(1,4,n)
                        rt(2,4,n) = half+rt(2,4,n)
                     endif
                  endif
               else
                  if ( i.eq.3 ) then
                     ier=10
                     if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
                     return
                  endif
!
!.................a c-axis mirror
!
                  call sgrmat(rt(1,1,n),1,0,0,0,1,0,0,0,-1)
                  rt(3,4,n) = d(3,3)
                  if ( i.eq.1.or.i.eq.5 ) rt(1,4,n) = half
                  if ( i.eq.2.or.i.eq.5 ) rt(2,4,n) = half
                  if ( m.eq. 2.and.l(1,2).eq.17) then
!
!....................if this is a 63-axis the mirror is at 1/4
!
                     if ( l(2,2).eq.14 ) rt(3,4,n)=half
                  endif
               endif
            endif
         else if (i.eq.6 ) then
!
!...........d type mirror
!
            if ( lcent.le.1 ) then
               ier=11
               if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
               return
            endif
!
!           go to (500,246,247,248),m
!
            if (m.eq.2) then
               if ( laue.gt.3 ) then
                  call sgrmat(rt(1,1,n),1,0,0,0,1,0,0,0,-1)
                  rt(1,4,n) = fourth
                  rt(2,4,n) = fourth
                  if ( id.eq.2 ) rt(3,4,n)=fourth
               else if (k.eq.2) then
                  call sgrmat(rt(1,1,n),1,0,0,0,-1,0,0,0,1)
                  rt(1,4,n) = fourth
                  if ( id.eq.2 ) rt(2,4,n)=fourth
                  if ( laue.eq.5 ) rt(2,4,n) = d(2,1)
                  rt(3,4,n) = fourth
               else
                  call sgrmat(rt(1,1,n),-1,0,0,0,1,0,0,0,1)
                  if ( id.eq.2 ) rt(1,4,n)=fourth
                  rt(2,4,n) = fourth
                  rt(3,4,n) = fourth
               endif
            else if (m.eq.3) then
               call sgrmat(rt(1,1,n),1,0,0,0,-1,0,0,0,1)
               rt(1,4,n) = fourth
               if ( id.eq.2 ) rt(2,4,n)=fourth
               if ( laue.eq.5 ) rt(2,4,n) = d(2,1)
               rt(3,4,n) = fourth
            else if (m.eq.4) then
               if ( l(1,2).eq.15.or.l(1,3).eq.14 ) then
!
!.................cubic or tetragonal. d-glide along diagonal
!
                  call sgrmat(rt(1,1,n),0,1,0,1,0,0,0,0,1)
                  rt(1,4,n) = fourth
                  rt(2,4,n) = fourth
                  rt(3,4,n) = fourth
                  if (l(1,3).eq.13) then
                     rt(1,4,n) = zero
                     rt(2,4,n) = half
                  endif
               else
!
!.................it is not tetragonal or cubic
!
                  call sgrmat(rt(1,1,n),1,0,0,0,1,0,0,0,-1)
                  rt(1,4,n) = fourth
                  rt(2,4,n) = fourth
                  if ( id.eq.2 ) rt(3,4,n)=fourth
               endif
            endif
         else if (i.eq.7) then
!
!...........1 fold rotation
!
            if ( l(2,m).eq.3 ) then
!
!..............we have a center of symmetry
!
               ncent = 1
            endif
!
!...........Loop exit
!
            goto 20
!
         else if (i.eq.8) then
!
!...........2 fold rotation axis
!...........we will not allow a -2 axis.
!
            if ( l(2,m).eq.3 )then
               ier=19
               if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
               return
            endif
!
!           go to (500,256,257,258),m
!
            if (m.eq.2) then
               if ( k.eq.2 )then
!
!.................rotation about the b-axis
!
                  call sgrmat(rt(1,1,n),-1,0,0,0,1,0,0,0,-1)
                  rt(1,4,n) = d(1,2)
                  rt(3,4,n) = d(3,2)
                  if ( l(2,m).eq.12 ) rt(2,4,n)=half
               else
!
!..................rotation about the a-axis. (orthogonal cell)
!
                   call sgrmat(rt(1,1,n),1,0,0,0,-1,0,0,0,-1)
                   rt(2,4,n) = d(2,1)
                   rt(3,4,n) = d(3,1)
                   if ( iabs(l(2,m)-13).eq.1 ) rt(1,4,n) = half
               endif
            else if (m.eq.3) then
               if ( l(1,2).eq.14 .or. l(1,2).eq.17 ) then
                  if (l(1,2).eq.17 .and. l(1,4).ne.12) then
!
!.....................2-axis normal to (-2110)
!.....................used for the p 6n22 groups
!
                      call sgrmat(rt(1,1,n),1,-1,0,0,-1,0,0,0,-1)
                  else
                     if ( laue.eq.7 ) then
!
!........................2-axis normal to (110)
!
                         call sgrmat(rt(1,1,n),0,-1,0,-1,0,0,0,0,-1)
                     else
!
!........................2-axis along to (11-20) trig  and (110) tetrag
!........................used for the p 3n21 groups
!
                         call sgrmat(rt(1,1,n),0,1,0,1,0,0,0,0,-1)
                         rt(1,4,n) = d(2,1)
                         if ( l(2,m).eq.12 ) rt(1,4,n)=rt(1,4,n)+half
                         rt(2,4,n) = -d(2,1)
                         rt(3,4,n) = d(3,1)
                     endif
                  endif
               else
!
!.................it is not a hexagonal or trigonal group
!.................rotation about the b-axis
!
                  call sgrmat(rt(1,1,n),-1,0,0,0,1,0,0,0,-1)
                  rt(1,4,n) = d(1,2)
                  rt(3,4,n) = d(3,2)
                  if ( l(2,m).eq.12 ) rt(2,4,n)=half
               endif
            else if (m.eq.4) then
!
               if ( l(1,2).ge.14 .or. l(1,3).eq.14) then
                  if ( l(1,2).eq.15 ) then
!
!.....................2-axis along to (11-20) trig  and (110) tetrag
!.....................used for the p 3n21 groups
!
                      call sgrmat(rt(1,1,n),0,1,0,1,0,0,0,0,-1)
                      rt(1,4,n) = d(2,1)
                      if ( l(2,m).eq.12 ) rt(1,4,n)=rt(1,4,n)+half
                      rt(2,4,n) = -d(2,1)
                      rt(3,4,n) = d(3,1)
                  else
!cc                  if ( l(1,2).eq.17.and.l(1,3).ne.12 )
!
!.....................2-axis normal to (10-10)
!
                     call sgrmat(rt(1,1,n),1,0,0,1,-1,0,0,0,-1)
                  endif
               else
                  call sgrmat(rt(1,1,n),-1,0,0,0,-1,0,0,0,1)
                  rt(1,4,n) = d(1,3)
                  rt(2,4,n) = d(2,3)
                  if ( iabs(l(2,m)-13).eq.1 ) rt(3,4,n) = half
                  if (l(2,m).eq.16) rt(3,4,n) = half
               endif
            endif
         else if (i.eq.9) then
            if (m.eq.2 .or. m.eq.3) then
               if (m.eq.2 .and. laue.gt.7) then
                  call sgrmat(rt(1,1,n),0,-1,0,1,-1,0,0,0,1)
                  if ( l(2,m).eq.12 ) rt(3,4,n)=third
                  if ( l(2,m).eq.13 ) rt(3,4,n)=2d0/3d0
                  if ( l(2,2).eq.3 ) ncent=1
               else
!
!.................1 fold rotation
!
                  if ( l(2,m).eq.3 ) then
!
!.....................we have a center of symmetry
!
                      ncent = 1
                   endif
!
!..................Loop exit
!
                   goto 20
               endif
            else if (m.eq.4) then
               ier=25
               if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
               return
            endif
         else if (i.eq.10 ) then
!
!...........four fold axis
!
            if ( m.ne.2 ) then
               ier=12
               if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
               return
            endif
            if ( l(2,2).eq. 3 ) then
               call sgrmat(rt(1,1,n),0,1,0,-1,0,0,0,0,-1)
               rt(1,4,n) = d(1,3)
               rt(2,4,n) = d(2,3)
               rt(3,4,n) = d(3,3)
            else
               call sgrmat(rt(1,1,n),0,-1,0,1,0,0,0,0,1)
               rt(1,4,n) = d(1,3)
               rt(2,4,n) = d(2,3)
               if ( l(2,2).eq.12 ) rt(3,4,n) = fourth
               if ( l(2,2).eq.13 ) rt(3,4,n) = half
               if ( l(2,2).eq.14 ) rt(3,4,n) = 0.75d0
            endif
         else if (i.eq.11) then
            ier=25
            if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
            return
         else if (i.eq.12) then
!
!...........6-axis
!
            if ( m.ne.2 )then
               ier=13
               if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
               return
            endif
            if ( l(2,2).eq.3 ) then
               call sgrmat(rt(1,1,n),-1,1,0,-1,0,0,0,0,-1)
               if ( l(1,3).eq.2.or.l(1,4).eq.2 ) rt(3,4,n)=half
            else
               call sgrmat(rt(1,1,n),1,-1,0,1,0,0,0,0,1)
            if ( l(2,2).gt.11.and.l(2,2).lt.18 ) rt(3,4,n)=(l(2,2)-11)/6.0d0
            endif
!
!        else if (i.eq.13 .or. i.eq.14 .or. i.eq.15) then
!
         endif
!
         rt(1,4,n) = mod(rt(1,4,n)+5d0,1d0)
         rt(2,4,n) = mod(rt(2,4,n)+5d0,1d0)
         rt(3,4,n) = mod(rt(3,4,n)+5d0,1d0)
         rt(5,2,n) = 1728*rt(1,4,n)+144*rt(2,4,n)+12*rt(3,4,n)
         njs=n-1
         do 69 m2=1,njs
            if (rt(5,1,m2).eq.rt(5,1,n) .or. rt(5,1,m2).eq.-rt(5,1,n))  then
               if (rt(5,1,m2).eq.-rt(5,1,n)) ncent=1
               if (rt(5,2,n).ne.rt(5,2,m2)) then
                  call sgtrcf(m,rt,n,m2,lcent,laue,ier,lptx)
                  if (ier.gt.0) ierx=ier
                  ier=0
               endif
!
!..............Loop exit
!
               goto 20
            endif
 69      continue
         n = n+1
         if ( n.gt.25 ) then
            ier=14
            if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
            return
         endif
! 59      continue
 59      localident=0
         nxl = n-1
!
!........exit
!
         if ( nxl.lt.nxi )  goto 49
         do 195 nx=nxi,nxl
            do 196 m1=2,nx
               call sgmtml(rt,m1,rt,nx,rt,n)
               njs=n-1
               do 194 m2=1,njs
                  if (rt(5,1,m2).eq.rt(5,1,n) .or. rt(5,1,m2).eq.-rt(5,1,n))  then
                     if (rt(5,1,m2).eq.-rt(5,1,n)) ncent=1
!
!....................Loop exit
!
                     go to 196
                  endif
194            continue
               n = n+1
               if ( n.gt.25 ) then
                  ier=15
                  if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
                  return
               endif
196         continue
!
!............Loop exit
!
            if ( n-1.eq.nxl ) goto 49
195      continue
         nxi = nxl+1
!
!........Repeat whole cycle
!
         go to 59
!
!........After completing loop cycle
!
 49      continue
!
!........search for a / . a mirror perpendicular to this axis
!........Next m
!
         if ( l(1,m).lt.12 .or. l(2,m).eq.3 )  goto 20
         do 39 i=2,3
!
!...........Next m
!
            if ( l(i,m).eq.0 ) goto 20
!
!...........exit loop 39
!
            if ( l(i,m).eq.19 ) goto 29
            if ( l(i,m).lt.12 ) then
               ier=16
               if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
               return
            endif
 39      continue
!
!........Test for another matrix element
!
 29      if ( l(i+1,m).le.1 ) then
            ier=17
            if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
            return
         endif
         i = iabs(l(i+1,m)-5)
!
!........Next matrix element
!        Non ortodox control flow
!
         goto 10
!
!.....end of loop
!
 20   continue
      neqv = n-1
      do 11 i=1,3
         do 21 k=1,neqv
            do 31 j=1,3
               jrt(i,j,k) = int(rt(i,j,k))
 31         continue
            jrt(i,4,k) = int(12*rt(i,4,k)+144.1)
            jrt(i,4,k) = jrt(i,4,k)-12*(jrt(i,4,k)/12)
 21      continue
 11   continue
!
!.....print results and construct interface symmetry matrix
!
      call sgprnt(spg(1:24),lpt)
      if ( ierx.eq.0 ) return
      ier = ierx
!
!.....Error handling
!
      if (lptx.ge.0) call sgerrs(spg(1:24),ier,lptx)
      return
  end subroutine spgr

  !> Determines parametrs for I-F centered lattices
  subroutine ifcent (l, d, lcent, ier)
      implicit none
      real*8            d(3,3)
      integer           l(4,4)
      integer           lcent, ier
!
!.....is the lattice I or F-centered?
      if ( lcent.lt.5 ) then
!
!........is the lattice C-centered?
         if ( lcent.ne.4 ) then
            ier=23
            return
         endif
!
!........c-centred.  an A normal to C
         if ( l(3,2).eq.4.or.l(4,2).eq.4 )  then
!
!............C 4*/A * *
             d(1,3) = 0.25
             d(2,3) = 0.25
             if ( l(1,4).eq.3.or.l(1,4).eq.10 ) d(1,1)=0.5
             return
         endif
!
!........if ( l(3,2).eq.0 ) d(1,1)=2.0*d(2,2)+d(1,1)
         if ( d(1,1).eq.0.0 ) d(1,1)=2.0*d(2,2)
!
!........is there a 21 on the diagonal?
         if ( l(1,4).eq.13.and.l(2,4).eq.12 )  then
!
!...........P 4N/N * *
            d(1,3) = 0.5
            return
         endif
         if ( l(2,2).le.0 ) return
!
!........is there a N-glide normal to (110)?
         if ( l(1,4).ne.10 ) return
         if ( l(2,2).gt.15 ) return
         d(1,1) = d(1,1)-2.0*d(2,2)
!
!........P 4N/N * *
         d(1,3) = 0.5
         return
      endif
!
!.....yes.
!
!.....if there is a C along (110) place the D at Y=1/4
      if ( l(1,4).eq.2 ) d(2,1)=d(2,1)+0.5
!
!.....is this I 41/A * * or F 41/D * * ?
      if ( l(4,2).ne.4.and.l(4,2).ne.11 ) then
!
!........is there a 41 present?
         if ( l(2,2).ne.12 ) return
!
!........yes.
!........F-centered?
         if ( lcent.eq.6 ) then
!
!...........F 41 * *
            d(1,3) = 0.25
            d(2,3) = 0.75
            return
         endif
         d(2,3) = 0.5
!
!........set the B-axis translation flags for I 41 2 2
         if ( lcent.eq.5 ) d(1,2) = 0.5
         d(3,2) = 0.75
         return
      endif
!
!.....yes
!
      d(1,3) = 0.25
      if ( lcent.eq.5 ) d(2,3) = 0.75
      return
  end subroutine ifcent

  subroutine ortho (d, l, laue, lcent, id)
      implicit none
      integer           laue, lcent, id
      real*8            d(3,3)
      integer           l(4,4)

      integer           im, ir, ia, ib, ic, i21
!
!.....it is orthorhombic
      laue = 3
!
!.....set up counts of the various types of mirrors.
      im = 0
      ir = 0
      ia = 0
      ib = 0
      ic = 0
      id = 0
      i21 = 0
!
!.....do we have a 2-axis along A
      if ( l(1,2).ne.13 ) then
!
         ir = 1
         if ( l(1,2).eq.9 ) im=4
!        if ( l(1,2).eq.4 ) ier =5  return
         if ( l(1,2).eq.3 ) ib=1
         if ( l(1,2).eq.2 ) ic=1
         if ( l(1,2).eq.11 ) id=1
         if ( l(1,3).eq.4.or.l(1,3).eq.10 ) d(1,1)=0.5
         if ( l(1,4).eq.4.or.l(1,4).eq.10 ) d(1,1)=d(1,1)+0.5
!
!.....yes is it a 21?
      else if(l(2,2).eq.12) then
         d(1,2) = 0.5
         d(1,3) = 0.5
         i21 = 4
      endif
!
!.....do we have a 2-axis along B
      if ( l(1,3).ne.13 ) then
!
         ir = ir+1
         if ( l(1,3).eq.9 ) im=im+2
         if ( l(1,3).eq.4 ) ia=1
!        if ( l(1,3).eq.3 ) ier =5 return
         if ( l(1,3).eq.2 ) ic=ic+1
         if ( l(1,3).eq.11 ) id=id+1
         if ( l(1,2).eq.3.or.l(1,2).eq.10 ) d(2,2)=0.5
         if ( l(1,4).eq.3.or.l(1,4).eq.10 ) d(2,2)=d(2,2)+0.5
!
!.....yes, is it a 21?
      else if (l(2,3).eq.12) then
         d(2,1) = 0.5
         d(2,3) = 0.5
         i21 = i21+2
      endif
!
!.....do we have a 2-axis along C
      if ( l(1,4).ne.13 ) then
         ir = ir+1
         if ( l(1,4).eq.9 ) im=im+1
         if ( l(1,4).eq.4 ) ia=ia+1
         if ( l(1,4).eq.3 ) ib=ib+1
!        if ( l(1,4).eq.2 ) ier=5 return
         if ( l(1,4).eq.11 ) id=id+1
         if ( l(1,2).eq.2.or.l(1,2).eq.10 ) d(3,3)=0.5
         if ( l(1,3).eq.2.or.l(1,3).eq.10 ) d(3,3)=d(3,3)+0.5
!
!.....yes , is it a 21
      else if (l(2,4).eq.12) then
         d(3,1) = 0.5
         d(3,2) = 0.5
         i21 = i21+1
      endif
!.....if there are 3 mirrors check for centering
!     which may alter the origin location
!
      if ( ir.eq.3 )  then
!........3 mirrors present.  is the lattice centered?
         if ( lcent.eq.1 ) return
!........yes . is it A-centered?
         if ( lcent.eq.2 ) then
!...........an A-centered lattice.
!           only one B or C glide present.relocate the mirrors by A
            if ( ib+ic.ne.1 ) return
            if ( ia.eq.2 ) return
            d(2,2) = d(2,2)+0.5
            d(3,3) = d(3,3)+0.5
            return
         endif
!
!........no.   is it B-centered?
         if ( lcent.eq.3 ) then
!
!...........a B-centered lattice
            if ( ia+ic.ne.1 ) return
            if ( ib.eq.2 ) return
            d(1,1) = d(1,1)+0.5
            d(3,3) = d(3,3)+0.5
            return
         endif
!
!........no.  is it C-centered?
         if ( lcent.eq.4 ) then
!
!...........a C-centered lattice
            if ( ia+ib.ne.1 ) return
            if ( ic.eq.2 ) return
            d(1,1) = d(1,1)+0.5
            d(2,2) = d(2,2)+0.5
            return
         endif
!
!........no.   is it I-centered?
         if ( lcent.ne.5 ) return
!
!........yes.  if only 1 glide plane shift the mirrors by I
         if ( ia+ib+ic.ne.1 ) return
         d(1,1) = d(1,1)+0.5
         d(2,2) = d(2,2)+0.5
         d(3,3) = d(3,3)+0.5
         return
      endif
!
!.....less than 3 mirrors. set up the 2-axes locations
      if ( i21.eq.4.or.i21.eq.5.or.i21.eq.7 ) d(1,2)=0.0
      if ( i21.eq.6.or.i21.eq.7 ) d(1,3)=0.0
      if ( i21.eq.3 ) d(2,1)=0.0
      if ( i21.eq.2.or.i21.eq.6.or.i21.eq.7 ) d(2,3)=0.0
      if ( i21.eq.1.or.i21.eq.3.or.i21.eq.7 ) d(3,1)=0.0
      if ( i21.eq.5 ) d(3,2)=0.0
      if ( im.le.0 ) return
      if ( im.eq.1.and.(i21.eq.4.or.i21.eq.2)) then
         if ( d(3,3).eq.0.0 ) return
         d(3,3)=0.0
         d(3,2) = d(3,2)+0.5
         return
      endif
      if ( im.eq.2.and.(i21.eq.4.or.i21.eq.1)) then
         if ( d(2,2).eq.0.0 ) return
         d(2,2)=0.0
         d(2,1) = d(2,1)+0.5
         return
      endif
      if ( im.eq.4.and.(i21.eq.2.or.i21.eq.1)) then
         if ( d(1,1).eq.0.0 ) return
         d(1,1)=0.0
         d(1,3) = d(1,3)+0.5
         return
      endif
!     if ( im.eq.3.amd.i21.eq.4 ) ier=5 return
!     if ( im.eq.5.and.i21.eq.2 ) ier=5 return
!     if ( im.eq.6.and.i21.eq.1 ) ier=5 return
      return
  end subroutine ortho

  subroutine sgerrs (sgp, ierr, lpt)
      implicit none
      character*(30)    sgp
      integer           ierr, lpt

      integer           ier
!
      ier=ierr+1
      write (lpt,999) ier,sgp
!
!.....gigant if
!
      if (ier.eq.1) then
         write (lpt,1)
      else if (ier.eq.2) then
         write (lpt,2)
      else if (ier.eq.3) then
         write (lpt,3)
      else if (ier.eq.4) then
         write (lpt,4)
      else if (ier.eq.5) then
         write (lpt,5)
      else if (ier.eq.6) then
         write (lpt,6)
      else if (ier.eq.7) then
         write (lpt,7)
      else if (ier.eq.8) then
         write (lpt,8)
      else if (ier.eq.9) then
         write (lpt,9)
      else if (ier.eq.10) then
         write (lpt,10)
      else if (ier.eq.11) then
         write (lpt,11)
      else if (ier.eq.12) then
         write (lpt,12)
      else if (ier.eq.13) then
         write (lpt,13)
      else if (ier.eq.14) then
         write (lpt,14)
      else if (ier.eq.15) then
         write (lpt,15)
      else if (ier.eq.16) then
         write (lpt,16)
      else if (ier.eq.17) then
         write (lpt,17)
      else if (ier.eq.18) then
         write (lpt,18)
      else if (ier.eq.19) then
         write (lpt,19)
      else if (ier.eq.20) then
         write (lpt,20)
      else if (ier.eq.21) then
         write (lpt,21)
      else if (ier.eq.22) then
         write (lpt,22)
      else if (ier.eq.23) then
         write (lpt,23)
      else if (ier.eq.24) then
         write (lpt,24)
      else if (ier.eq.25) then
         write (lpt,25)
      else if (ier.eq.26) then
         write (lpt,26)
      endif
!
      return
!
!.....formats
!
    1 format ('either a 5-axis anywhere or a 3-axis in field 4')
    2 format ('less than 2 operator fields were found')
    3 format ('lattice operator was not a p, a, b, c, i, f or r')
    4 format ('rhombohedral lattice without a 3-axis')
    5 format ('minus sign does not preceed 1, 2, 3, 4 or 6')
    6 format ('lattice s. r. found an error')
    7 format ('1st operator in a field was a space. impossible')
    8 format ('i for computed go to out of range.')
    9 format ('an a-glide mirror normal to a')
   10 format ('a b-glide mirror normal to b')
   11 format ('a c-glide mirror normal to c')
   12 format ('d-glide in a primitive lattice')
   13 format ('a 4-axis not in the 2nd operator field')
   14 format ('a 6-axis not in the 2nd operator field')
   15 format ('more than 24 matrices needed to define the group')
   16 format ('more than 24 matrices needed to define the group')
   17 format ('improper construction of a rotation operator')
   18 format ('no mirror following a /')
   19 format ('a translation conflict between operators')
   20 format ('the 2bar operator is not allowed')
   21 format ('3 fields are legal only in r-lattices and m3 cubic')
   22 format ('syntax error. expected i-43d at this point.')
   23 format ('error unknown')
   24 format ('a or b centered tetragonal?  impossible!!!!')
   25 format ('only one bar allowed')
   26 format ('either a 5-axis anywhere or a 3-axis in field 4')
!
  999 format ('error no. ',i3,' in processing space group symbol ',a30)
!
  end subroutine sgerrs

  subroutine sglatc (k, l, d, lcent, laue, naxis, ier, ncubic, id)
      implicit none
      integer           l(4,4), k, lcent, laue, naxis, ier, ncubic, id
      real*8            d(3,3)
!
      integer           im, ia, ic, na, nb, nc
!
      if (k .lt.3 ) then
         if(l(1,2).eq.17) then
!
!...........6/m
!
            laue=11
            return
         else if(l(1,2).eq.14) then
!
!...........2bar
!
            laue=8
            return
         else if(l(1,2).eq.12) then
!
!...........1bar
!
            laue=1
            return
         else if(l(1,2).ne.15) then
!
!...........2/m, B-axis unique
!
            im=2
            laue = 2
            naxis = 2
            ia = 4
            ic = 2
            na = 1
            nb = 2
            nc = 3
            if ( l(2,im).eq.12 ) d(nb,naxis)=0.5d0
            if ( l(3,im).eq.ia.or.l(3,im).eq.10 ) d(na,naxis)=0.5d0
            if ( l(3,im).eq.ic.or.l(3,im).eq.10 ) d(nc,naxis)=0.5d0
            if ( l(4,im).eq.ia.or.l(4,im).eq.10 ) d(na,naxis)=0.5d0
            if ( l(4,im).eq.ic.or.l(4,im).eq.10 ) d(nc,naxis)=0.5d0
            return
         else
!
!...........4/m
!
            laue=4
!
!...........is there a n-glide normal to c?
!
            if ( lcent.le.3 ) then
               if (l(3,2).eq.10) then
!
!.................P 4N/N * *
!
                  d(1,3) = 0.5d0
                  return
                endif
                if (l(4,2).eq.10) d(2,3)=0.5d0
                return
!
!...........is it c-centered?
!
            else if (lcent.eq.4) then
!
!..............if there is no A-glide normal to C we are through
!
               if ( l(3,2).ne.4.and.l(4,2).ne.4 ) return
               d(1,3) = 0.25d0
               d(2,3) = 0.25d0
               if ( l(4,2).eq.4 ) d(2,3)=0.75d0
               return
!
!...........is it I-centered or F-centered
!
            else if (lcent.ge.5) then
!
!..............is there an A-glide or a D-glide normal to C?
!
               if ( l(4,2).ne.4.and.l(4,2).ne.11 ) then
!
!.................no
!
                  call ifcent (l,d,lcent,ier)
                  return
               else
!
!.................yes
!
                  d(1,3) = 0.75d0
                  if ( lcent.eq.5 ) d(2,3) = 0.25d0
                  return
               endif
            endif
         endif
!
      else if (k .eq. 3) then
!
!........3 fields.  M3 cubic. (R3R has been taken care of)
!
         if ( l(1,3).ne.14 ) then
            ier=20
            return
         endif
         laue = 13
!
!........set the B-axis translation flag if A 21 along A
!
         if ( l(2,2).eq.12 ) d(2,1)=0.5d0
!
!........set the C-axis translation flag if an A-glide normal to C
!
         if ( l(1,2).eq.3.or.l(1,2).eq.4 ) d(3,3)=0.5d0
         ncubic=1
         return
!
      else if (k .gt. 3) then
!
!.........four fields read.
!
          if ( l(1,3).eq.14 ) then
!
!...........it is m3m cubic
!
            laue = 14
!
!...........set the C-axis translation flag if an A-glide normal to C
!
            if ( l(1,2).eq.3.or.l(1,2).eq.4 ) d(3,3)=0.5d0
!
!...........is a 4N-axis specified?
!
            if ( l(1,2).ne.15 ) then
               ncubic=1
               return
            endif
!
!.......... yes.  is it 4bar?
!
            if ( l(2,2).eq.3 ) then
!
!..............4B.  is it 4B 3 M?
!
              if ( l(1,4).eq.9 )  then
                 ncubic=1
                 return
              endif
!
!.............no.  is it 4B 3 D?
!
              if ( l(1,4).eq.11 ) then
!
!................I 4B 3 D we hope
!
                 if ( lcent.ne.5 ) then
                    ier=21
                    return
                 endif
                 d(1,3) = 0.75d0
                 d(2,3) = 0.25d0
                 d(3,3) = 0.75d0
                 ncubic=1
                 return
              endif
!
!.............no.
!
              d(1,3) = 0.5d0
              d(2,3) = 0.5d0
              d(3,3) = 0.5d0
              ncubic=1
              return
           endif
!
!..........no.  is it a 4?
!
           if ( l(2,2).lt.12 .or. l(2,2).gt.14 ) then
              ncubic=1
              return
           endif
!
!..........no. is it a 41?
!
           if ( l(2,2).eq.12 ) then
!
!.............41-axis. is it F 41 3 2?
!
              if ( lcent.eq.6 ) then
!
!................F 41 3 2
!
                 d(1,3) = 0.25d0
                 d(2,3) = 0.25d0
                 ncubic=1
                 return
              endif
!
!.............no. it is either P 41 3 2 or I 41 3 2
!
              d(1,3) = 0.25d0
              d(2,3) = -0.25d0
              ncubic=1
              return
           endif
!
!..........no . is it a 42?
!
           if ( l(2,2).eq.13 ) then
!
!.............P 42 3 2
!
              d(1,3) = 0.5d0
              d(2,3) = 0.5d0
              ncubic=1
              return
           endif
!
!..........no.  it must be a 43
!
!..........P 43 3 2
!
           if ( lcent.eq.6 ) then
!
!.............F 41 3 2
!
              d(1,3) = 0.25d0
              d(2,3) = 0.25d0
              ncubic=1
              return
           endif
           d(1,3) = 0.75d0
           d(2,3) = 0.25d0
           ncubic=1
           return
         else if (l(1,2).eq.17) then
            if ( l(1,3).eq.12 ) then
!
!..............we have something like P 6N 1 *
!
               if ( l(1,4).ne.12 ) then
!
!.................it is hexagonal 6/mmm
!
                  laue = 12
                  return
               endif
!
!..............6/M
!
               laue = 11
               return
            endif
!
!...........it is hexagonal 6/mmm
!
            laue = 12
            return
!
         else if (l(1,2).eq.14) then
!
!...........it is trigonal P3**
!
            if ( l(1,3).eq.12 ) then
!
!..............it is trigonal 31*
!
               if ( l(1,4).eq.12 ) then
!
!.................3BAR
!
                  laue = 8
                  return
               endif
!
!..............it is trigonal 31M
!
               laue = 10
               return
            endif
            if ( l(1,4).ne.12 ) then
!
!..............it is hexagonal 6/mmm
!
               laue = 12
               return
            endif
!
!...........it is trigonal 3M1
!
            laue = 9
            return
         else if (l(1,2).eq.15) then
!
!...........it is tetragonal 4/mmm
!
            laue = 5
!
!...........if there is an N-glide normal to C place any
!
!           mirror normal to A at 1/4
            if ( l(3,2).eq.10.or.l(4,2).eq.10 ) d(1,1)=0.5d0
!
!...........if there is an A-glide normal to C place any
!
!           mirror normal to (110) at 1/4
            if ( l(3,2).eq.4.or.l(4,2).eq.4 ) d(2,2)=0.25d0
!
!...........if there is a 21 along B move place it at X=1/4
!
            if ( l(1,3).eq.13.and.l(2,3).eq.12 ) d(1,2)=0.5d0
!
!...........if there is a 21 along (110) move place it at X=1/4
!
            if ( l(1,4).eq.13.and.l(2,4).eq.12 ) d(2,1)=0.5d0
!
!...........if there is either a B- or N-glide normal to the A-axis
!           shift the mirror by 1/4 along the A-axis
!
            if ( l(1,3).eq.3.or.l(1,3).eq.10 ) d(1,1)=d(1,1)+0.5d0
!
!...........if there is either a B- OR N-glide normal to (110)
!           shift the mirror by 1/4 along the A-axis
!
            if ( l(1,4).eq.3.or.l(1,4).eq.10 ) d(2,2)=d(2,2)+0.25d0
!
!...........set the Z location for 2-axes along (110)
!
            if ( l(2,2).gt.11.and.l(2,2).lt.15.and.l(2,3).ne.12 ) d(3,1)=-(l(2,2)-11)/4.0d0
!
!...........set the Z-translation for 21-axes along (110)
!
            if ( l(1,4).eq.13.and.l(2,4).ne.12 ) then
            else
             if ( l(2,2).gt.11.and.l(2,2).lt.15 ) d(3,1)=(l(2,2)-11)/4.0d0
            endif
!
!...........set the Z-translation for 21-axes along B
!
            if ( l(1,3).eq.13.and.l(2,3).ne.12 ) then
            else
            if ( l(2,2).gt.11.and.l(2,2).lt.15 ) d(3,2)=(l(2,2)-11)/4.0d0
            endif
!
!...........place the D in F 4* D * at Y=7/8
!
            if ( l(1,3)+l(3,2).eq.11.and.lcent.eq.6 ) d(2,1)=0.75d0
!
!...........M in F 4** * * at X=1/8 if there is a C along (110)
!
            if ( l(1,4).eq.2.and.lcent.eq.6 ) d(1,1)=0.5d0
!
!...........is this a 4bar?
!
            if ( l(2,2).eq.3 ) then
!
!..............translations due to diagonal symmetry operation
!
               if ( l(1,3).eq.11.and.lcent.eq.6 ) d(3,1)=0.25d0
!
!..............if * 4B * 21 we want the mirror at X=1/4
!
               if ( l(1,4).eq.13.and.l(2,4).eq.12 ) d(1,1)=0.5d0
!
!..............C- or a N-glide along (110) set the 2-axis at Z=1/4
!
               if ( l(1,4).eq.2.or.l(1,4).eq.10 ) d(3,2)=0.5d0
!
!..............a B- or a N-glide along (110) set the 2 at X=1/4
!
               if ( l(1,4).eq.3.or.l(1,4).eq.10 ) d(1,2)=0.5d0
               if ( l(1,4).ne.11 ) return
               if ( lcent.eq.5 ) d(1,2) = 0.5d0
               d(3,2) = 0.75d0
               return
            endif
!
!...........is the lattice primitive?
!
            if ( lcent.gt.1 ) then
               call ifcent  (l,d,lcent,ier)
               return
            endif
!
!...........yes.  do we have a N-glide normal to C?
!
            if ( l(3,2).eq.10.or.l(4,2).eq.10 ) then
!
!..............P 4N/N * *
!
               d(1,3) = 0.5d0
               return
            endif
!
!...........no.  do we have a 21 along B?
!
            if ( l(1,3).eq.13.and.l(2,3).eq.12 ) then
               d(1,3) = 0.5d0
               d(2,3) = 0.5d0
               return
            endif
!
!...........no. do we have a N-glide normal to A?
!
            if ( l(1,3).ne.10 ) return
            if ( l(2,2).le.0 )  return
            if ( l(2,2).gt.15)  return
            d(1,3) = 0.5d0
            d(2,3) = 0.5d0
            return
         else if (l(1,2).eq.12) then
            if ( l(1,3).eq.12 ) then
               if ( l(1,4).eq.12 ) then
!
!.................1BAR
!
                  laue = 1
                  return
               endif
!
!..............it is C-axis unique monoclinic
!
               laue = 2
               naxis = 3
               ia = 4
               ic = 3
               na = 1
               nb = 3
               nc = 2
               im = 4
            else
!
!..............it is not C-axis unique monoclinic
!
               if ( l(1,4).ne.12 ) then
                  call ortho (d, l, laue, lcent, id)
                  return
               endif
               im = 3
!
!..............it is B-axis unique monoclinic. (FULL SYMBOL USED)
!
               laue = 2
               naxis = 2
               ia = 4
               ic = 2
               na = 1
               nb = 2
               nc = 3
            endif
            if ( l(2,im).eq.12 ) d(nb,naxis)=0.5d0
            if ( l(3,im).eq.ia.or.l(3,im).eq.10 ) d(na,naxis)=0.5d0
            if ( l(3,im).eq.ic.or.l(3,im).eq.10 ) d(nc,naxis)=0.5d0
            if ( l(4,im).eq.ia.or.l(4,im).eq.10 ) d(na,naxis)=0.5d0
            if ( l(4,im).eq.ic.or.l(4,im).eq.10 ) d(nc,naxis)=0.5d0
            return
         else if (l(1,3).eq.12) then
            if ( l(1,4).ne.12 ) then
               call ortho (d, l, laue, lcent, id)
               return
            endif
!
!...........it is A-axis unique monoclinic
!
            laue = 2
            naxis = 1
            ia = 3
            ic = 2
            na = 2
            nb = 1
            nc = 3
            im = 2
            if ( l(2,im).eq.12 ) d(nb,naxis)=0.5d0
            if ( l(3,im).eq.ia.or.l(3,im).eq.10 ) d(na,naxis)=0.5d0
            if ( l(3,im).eq.ic.or.l(3,im).eq.10 ) d(nc,naxis)=0.5d0
            if ( l(4,im).eq.ia.or.l(4,im).eq.10 ) d(na,naxis)=0.5d0
            if ( l(4,im).eq.ic.or.l(4,im).eq.10 ) d(nc,naxis)=0.5d0
            return
         else
!
!...........It may be orthorhombic
!
            call ortho (d, l, laue, lcent, id)
            return
         endif
      endif
!
      return
  end subroutine sglatc

  !> Converts -N notation into Nbar notation
  subroutine sglpak (l, ier)
      implicit none
      integer           l(4), ier

      if ( l(2).lt.12 ) ier = 4
      if ( l(2).gt.17 ) ier = 4
      l(1) = l(2)
      l(2) = 3
      return
  end subroutine sglpak

  !> Space group matrix multiplication
  subroutine sgmtml (x, i, y, j, z, k)
      implicit none
      real*8         x(5,4,25), y(5,4,25), z(5,4,25)
      integer        i, j, k
!
      integer        l, m, n
!
      do 100 l=1,4
         do 101 m=1,4
            z(l,m,k) = 0.0d0
            do 102 n=1,4
               z(l,m,k) = z(l,m,k)+y(l,n,j)*x(n,m,i)
102         continue
101      continue
100   continue
      z(1,4,k) = mod(5.0d0+z(1,4,k),1.0d0)
      z(2,4,k) = mod(5.0d0+z(2,4,k),1.0d0)
      z(3,4,k) = mod(5.0d0+z(3,4,k),1.0d0)
      z(5,1,k) = 81*(2*z(1,1,k)+3*z(1,2,k)+4*z(1,3,k))+9*(2*z(2,1,k)+3*z(2,2,k)+4*z(2,3,k))+2*z(3,1,k)+3*z(3,2,k)+4*z(3,3,k)
      z(5,2,k) = 1728*z(1,4,k)+144*z(2,4,k)+12*z(3,4,k)
      z(5,3,k) = 0.d0
      z(5,4,k) = 0.d0
!
      return
  end subroutine sgmtml

  !> Prints out spatial symmetry information
  subroutine sgprnt (spg, lpt)
    use struct
    use tools_io
    use param
    implicit none
!
      character*(*)     spg
      integer           lpt
!
      character*(1)     nax(3)
      character*(3)     outl(3,2,24), xyz(12), polar(8), tra(11)
      character*(8)     chlaue(14)
      character*(12)    ltyp(7), syst(8)
!
      integer           ncvt(7), nsys(14)
      real*8            cenv(3,6), tra2(11)
!      dimension         ncenv2(3,6)
      integer           multgen
      integer           i, j, k, ij, ik, i1, npx, npy, npz
      integer           npxyz, npyxz, naxi
!
!.....datas
!
      data&
       tra2( 1),tra2( 2),tra2( 3) /zero     ,zero     ,sixth    /,&
       tra2( 4),tra2( 5),tra2( 6) /fourth   ,third    ,zero     /,&
       tra2( 7),tra2( 8),tra2( 9) /half     ,zero     ,twothird /,&
       tra2(10),tra2(11)          /threefourth        ,fivesixth/
!
!      data
!     & ncenv2(1,1),ncenv2(2,1),ncenv2(3,1) /0,6,6/,
!     & ncenv2(1,2),ncenv2(2,2),ncenv2(3,2) /6,0,6/,
!     & ncenv2(1,3),ncenv2(2,3),ncenv2(3,3) /6,6,0/,
!     & ncenv2(1,4),ncenv2(2,4),ncenv2(3,4) /6,6,6/,
!     & ncenv2(1,5),ncenv2(2,5),ncenv2(3,5) /4,8,8/,
!     & ncenv2(1,6),ncenv2(2,6),ncenv2(3,6) /8,4,4/
!
      data ltyp /   'primitive   ', 'A-centered  ', 'B-centered  ',&
                    'C-centered  ', 'I-centered  ', 'F-centered  ',&
                    'R-centered  ' /
      data syst /   'triclinic   ', 'monoclinic  ', 'orthorhombic',&
                    'tetragonal  ', 'trigonal-R  ', 'trigonal-P  ',&
                    'hexagonal   ', 'cubic       ' /
!
      data chlaue / '1bar    ', '2/m     ', 'mmm     ', '4/m     ',&
                    '4/mmm   ', '3bar    ', '3barm   ', '3bar    ',&
                    '3barm1  ', '3bar1m  ', '6/m     ', '6/mmm   ',&
                    'm3      ', 'm3m     ' /
!
      data polar  / 'x  ',   'y  ',   'x y',   'z  ',   'x z', 'y z',   'xyz',   '111' /
!
      data nax    / 'a',      'b',      'c' /
!
      data xyz   /  '-z ',   '-y ',   '-x ',   'x-y',   'err',&
                    'y-x',   ' x ',   ' y ',   ' z ',   '+x ',&
                    '+y ',   '+z ' /
!
      data tra   /  '   ',   'err',   '1/6',   '1/4',   '1/3',&
                    'err',   '1/2',   'err',   '2/3',   '3/4',&
                    '5/6' /
!
      data&
       nsys( 1),nsys( 2),nsys( 3),nsys( 4),nsys( 5) /1,2,3,4,4/,&
       nsys( 6),nsys( 7),nsys( 8),nsys( 9),nsys(10) /5,5,6,6,6/,&
       nsys(11),nsys(12),nsys(13),nsys(14)          /7,7,8,8/
!
      data&
       ncvt(1),ncvt(2),ncvt(3),ncvt(4),ncvt(5) /1,2,2,2,2/,&
       ncvt(6),ncvt(7)                         /4,3/
!
      data&
       cenv(1,1),cenv(2,1),cenv(3,1) /zero      ,half      ,half      /,&
       cenv(1,2),cenv(2,2),cenv(3,2) /half      ,zero      ,half      /,&
       cenv(1,3),cenv(2,3),cenv(3,3) /half      ,half      ,zero      /,&
       cenv(1,4),cenv(2,4),cenv(3,4) /half      ,half      ,half      /,&
       cenv(1,5),cenv(2,5),cenv(3,5) /third,twothird,twothird/&
       cenv(1,6),cenv(2,6),cenv(3,6) /twothird,third,third/
!
      ncv = ncvt(lcent)
      multgen = ncv*neqv*(ncent+1)
      lsys = nsys(laue)
      cen(1,1) = zero
      cen(2,1) = zero
      cen(3,1) = zero
      if ( ncv.gt.1 )  then
         j = lcent-1
         if ( lcent.eq.6 ) j=1
         if ( lcent.eq.7 ) j=5
         do 100 i=2,ncv
            cen(1,i) = cenv(1,j)
            cen(2,i) = cenv(2,j)
            cen(3,i) = cenv(3,j)
            j = j+1
100      continue
      endif
      npx = 1
      npy = 2
      npz = 4
      npxyz = 0
      npyxz = 1
      do 120 i = 1, neqv
         if ( jrt(1,1,i).le.0 ) npx = 0
         if ( jrt(2,2,i).le.0 ) npy = 0
         if ( jrt(3,3,i).le.0 ) npz = 0
         if ( jrt(1,3,i).gt.0 ) npxyz = 8
         if ( jrt(1,3,i).lt.0 ) npyxz = 0
120   continue
      npol = (npx+npy+npz+npxyz*npyxz)*(1-ncent)
      naxi = naxis
      if (lsys.eq.4.or.lsys.eq.6.or.lsys.eq.7) naxi = 3
!
!.....Create rotm matrix
!
      do 211 i=1,neqv
         do 212 j=1,3
            do 213 k=1,3
               rotm(j,k,i)=jrt(j,k,i)
 213        continue
            rotm(j,4,i)=tra2(jrt(j,4,i)+1)
 212     continue
 211  continue
!
!.....Write the results
!
      if (lpt.ge.0) then
         write (lpt,'("* SPGR: space group symbol recognition library"/)')
         write (lpt,1005) trim(adjustl(spg))
         if (ncent.eq.0) write (lpt,1006)
         if (ncent.eq.1) write (lpt,1007)
         write (lpt,1008) ltyp(lcent)
         write (lpt,1009) syst(lsys)
         write (lpt,1010) chlaue(laue)
         write (lpt,1011) multgen
         if (naxi.gt.0) write (lpt,1012) nax(naxi)
         if (naxi.eq.0) write (lpt,1013)
         if (npol.gt.0) write (lpt,1014) polar(npol)
         if (npol.le.0) write (lpt,1015)
         write (lpt,1016) neqv
      endif
!
!.....redefine npol
!
      npol=npol*10000+lsys*1000+multgen
!
!.....first part concluded, start with symmetry relations
!
      do 200 i=1,neqv
         do 190 j=1,3
            ij = 2*jrt(j,1,i)+3*jrt(j,2,i)+4*jrt(j,3,i)+5
            ik = jrt(j,4,i)+1
            if (ik.gt.1 .and. ij.gt.6) ij=ij+3
            outl(j,2,i) = xyz(ij)
            outl(j,1,i) = tra(ik)
 190     continue
 200  continue
!
!.....write the symmetry matrices and lattice centering vectors.
!
      if (lpt.ge.0) then
         write (lpt,1026) ((outl(i1,1,i),outl(i1,2,i),i1=1,3),i=1,neqv)
         write (lpt, 7001) neqv
         write (lpt, 7020) ncv
         do 7030 i=1,ncv
            write (lpt, 7025) (cen(j,i), j=1,3)
 7030    continue
         write (lpt,*)
      endif
!
!.....Formats.
!
 1026 format (3(2x,'(',2(2a3,1h,),2a3,') '))
 1005 format ("* Space group: ",a)
 1006 format ('  the structure is noncentrosymmetric')
 1007 format ('  the structure is centrosymmetric')
 1008 format ('  the lattice is ',a12)
 1009 format ('  the system is ',a12)
 1010 format ('  the laue group is ',a8)
 1011 format ('  the multiplicity of the general position is ',i3)
 1012 format ('  the unique axis is ',a1)
 1013 format ('  there is no unique axis')
 1014 format ('  the space group is polar along ',a3)
 1015 format ('  the space group is not polar')
 1016 format ('+ Number of generated symmetry relations: ',i2)
 7001 format ('+ Number of symmetry matrices:',i5)
 7020 format ('+ Number of Lattice Centering Vectors is ',i5)
 7025 format ('  Cent. vect.: ', 3f10.7)
!
      return
  end subroutine sgprnt

  !> Builds the rt(,) matrix out of the a..o rotational elements.
  !> rt(,) is a 5x4 matrix, which contains a 3x3 rotation
  !> plus a column for the traslation, a row that could
  !> represent a deformation of the axes (set to 0), and ...
  subroutine sgrmat (rt, a, b, c, d, e, f, g, h, o)
      implicit none
      real*8            rt(5,4)
      integer           a, b, c, d, e, f, g, h, o
!
      rt(1,1) = a
      rt(1,2) = b
      rt(1,3) = c
      rt(1,4) = 0.0d0
      rt(2,1) = d
      rt(2,2) = e
      rt(2,3) = f
      rt(2,4) = 0.0d0
      rt(3,1) = g
      rt(3,2) = h
      rt(3,3) = o
      rt(3,4) = 0.0d0
      rt(4,1) = 0.0d0
      rt(4,2) = 0.0d0
      rt(4,3) = 0.0d0
      rt(4,4) = 1.0d0
      rt(5,1) = 81 * (2 * rt(1,1) + 3 * rt(1,2) + 4 * rt(1,3)) + 9 * (2 * rt(2,1) +&
           3 * rt(2,2) + 4 * rt(2,3)) + 2 * rt(3,1) + 3 * rt(3,2) + 4 * rt(3,3)
      rt(5,2) = 0.0d0
      rt(5,3) = 10.0d0
      rt(5,4) = 20.d0
!
      return
  end subroutine sgrmat

  !> Check rotation/traslation compatibility.
  subroutine sgtrcf (m, rt, n, m2, lcent, laue, ier, lpt)
      implicit none
      integer           m, n, m2, lcent, laue, ier, lpt
      real*8            rt(5,4,25)
!
      integer           icenv(3,5), ncvt(7), jcvt(7)
      integer           irn, irm, irx, iry, irz, ncv, icv, jcv
      integer           irx1, iry1, irz1, m2z, itottr
!
      data&
       icenv(1,1),icenv(2,1),icenv(3,1) /0,0,0/,&
       icenv(1,2),icenv(2,2),icenv(3,2) /0,6,6/,&
       icenv(1,3),icenv(2,3),icenv(3,3) /6,0,6/,&
       icenv(1,4),icenv(2,4),icenv(3,4) /6,6,0/,&
       icenv(1,5),icenv(2,5),icenv(3,5) /6,6,6/
!
      data&
       ncvt(1),ncvt(2),ncvt(3),ncvt(4),ncvt(5) /1,2,3,4,5/,&
       ncvt(6),ncvt(7)                         /4,1/
!
      data&
       jcvt(1),jcvt(2),jcvt(3),jcvt(4),jcvt(5) /1,1,2,3,4/,&
       jcvt(6),jcvt(7)                         /1,1/
!
      ier = 0
      irn = int(rt(5,2,n))
      irm = int(rt(5,2,m2))
      irx = mod((irn/144+irm/144),12)
      iry = mod((irn/12 +irm/12 ),12)
      irz = mod(irn+irm,12)
      ncv = ncvt(lcent)
      jcv = jcvt(lcent)
!
      do 100 icv=1,ncv,jcv
         irx1 = mod(irx+icenv(1,icv),12)
         iry1 = mod(iry+icenv(2,icv),12)
         irz1 = mod(irz+icenv(3,icv),12)
!
!........does this pair make a 1bar?
!
         m2z = m2
!
!........yes
!
         if ( rt(5,1,n)+rt(5,1,m2).eq.0 ) m2z=1
!
!........no
!
         if ( rt(3,3,n)+rt(3,3,m2z).le.0 ) irz1=0
!
!........is this an operator operating along the face diagonal?
!
         if ( laue.le.3.or.m.ne.4 ) then
!
!...........no
!
            if ( rt(1,1,n)+rt(1,1,m2z).le.0 ) irx1=0
            if ( rt(2,2,n)+rt(2,2,m2z).le.0 ) iry1=0
         else
!
!...........yes
!
            irx1 = mod(irx1+iry1,12)
            iry1 = 0
         endif
         itottr = 144*irx1+12*iry1+irz1
         if ( itottr.eq.0 ) return
 100  continue
!
!.....An error has been found
!
      if (lpt.ge.0) write (lpt,2990) rt(5,2,n),rt(5,2,m2),itottr,irx,iry,irz,rt(5,1,n),rt(5,1,m2)
      ier = 18
      if (lpt.ge.0) write (lpt,2991) m,n,m2
!
!.....formats
!
2990  format (2(1x,f9.1),4(1x,i4),2(1x,f9.1))
2991  format('operator ',i2,' generates matrix ',i3,' which has a',' translation conflict',2(1x,i2))
!
      return
  end subroutine sgtrcf
  
end module spgs
