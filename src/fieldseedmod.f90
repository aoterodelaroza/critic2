! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
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

! Field seed class. 
module fieldseedmod
  use wfn_private, only: molden_type_psi4
  use param, only: ifformat_unknown, mlen, mmlen
  implicit none

  private

  !> Minimal amount of information to generate a field (usually, a
  !> pointer to an external field).
  type fieldseed
     integer :: iff = ifformat_unknown !< Field format/type
     character(len=mlen) :: errmsg = "" !< Error message
     integer :: nfile = 0 !< Number of external files provided
     character(len=mlen), allocatable :: file(:) !< External scalar field file names
     character*10, allocatable :: piat(:) !< pi atomic symbols
     integer :: n(3) !< grid size for load as
     logical :: isry = .false. !< use rydberg units in LOAD AS POT
     integer :: clm1, clm2 !< clm fields for add and sub
     character(len=mlen) :: ids = "" !< sizeof for load as promolecular/core; id.s for lap/grad; id1.s in clm
     character(len=mlen) :: ids2 = "" !< id2.s in clm
     character(len=mmlen) :: expr = "" !< expression in load as
     character(len=mmlen) :: elseopt = "" !< options to be handled elsewhere
     logical :: testrmt = .true. !< whether to test rmt (wien/elk)
     logical :: readvirtual = .false. !< Read the virtual orbitals (fchk/molden)
     integer :: vaspblk = 1 !< Which VASP block to read (CHGCAR/CHG/...)
     character(len=mlen) :: fid = "" !< field ID
     logical :: nou = .false. !< wannier option nou
     integer :: molden_type = molden_type_psi4 !< type of molden file
     ! pwc
     integer :: pwcspin = 0 !< read which spin (0=all,1=up,2=dn)
     integer, allocatable :: pwcikpt(:) !< read which kpts
     integer, allocatable :: pwcibnd(:) !< read which bands
     real*8 :: pwcemin = -1d40 !< read minimum of energy range 
     real*8 :: pwcemax = 1d40 !< read maximum of energy range
   contains
     procedure :: end => fieldseed_end !< Terminate a fieldseed; set the type to unknown
     procedure :: parse => fieldseed_parse !< Build a field seed from an external command
     procedure :: parse_options => fieldseed_parse_options !< Parse field options from a command
  end type fieldseed
  public :: fieldseed
  
  interface
     module subroutine fieldseed_end(f)
       class(fieldseed), intent(inout) :: f
     end subroutine fieldseed_end
     module subroutine fieldseed_parse(f,line,withoptions,lp0)
       class(fieldseed), intent(inout) :: f
       character*(*) :: line
       logical, intent(in) :: withoptions
       integer, intent(inout), optional :: lp0
     end subroutine fieldseed_parse
     module subroutine fieldseed_parse_options(f,line,lp0)
       class(fieldseed), intent(inout) :: f
       character*(*) :: line
       integer, intent(inout), optional :: lp0
     end subroutine fieldseed_parse_options
  end interface

end module fieldseedmod
