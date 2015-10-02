! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------
Module internal_types
USE Precision_Model
IMPLICIT NONE
private

!
! Named constants, to enhance readability
!
INTEGER, public, PARAMETER ::             &
                      Simplex = 1 ,       &
                      Hyperrectangle = 2, &
                      EPSTABLENGHT = 52

!
! The following record stores information that determines the
! behaviour of the global adaptive integration procedure.
! This information is passed to the region processor.
!
type, public :: integrator_info
    integer         :: key, nrsub
    REAL(kind=stnd) :: tune
    logical         :: uniform_subdiv 
end type integrator_info

!
! The following record stores information that is specific
! for the region collection and other data structures.
!
type, public :: collection_info
    integer :: dimens, nrvert, niinfo, nrinfo
end type collection_info

!
! The following record stores scalar information supplied by the user
! or meant for the user related to the integration procedure.
!
type, public :: user_info
! Input
    integer         :: numfun,numrgn,minpts,maxpts
    REAL(kind=stnd) :: epsabs,epsrel
    logical         :: restart
! Output
    integer   :: neval,ifail
end type user_info

type, public :: epsalg_mem
    LOGICAL :: HEURISTIC_USED
    INTEGER :: DIVLEVEL
    INTEGER, POINTER, DIMENSION(:) :: NRRCOPY
    REAL(kind=stnd) ::  ERRORMAXPOOL, EPSABS, EPSREL
    REAL(kind=stnd), POINTER, DIMENSION(:) :: ERLARG,  &
                                              RESULT1, &
                                              ABSERR1
    REAL(kind=stnd), POINTER, DIMENSION(:,:) :: RCOPY, &
                                                RESLA
end type epsalg_mem

end module internal_types
