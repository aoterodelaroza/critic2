! This file is F-compatible, except for upper/lower case conventions.
!--------------------------------------------------------------------

MODULE Error_Handling

Implicit NONE

PUBLIC  :: Print_Error, Handle_Error
PRIVATE :: Message
INTEGER, PRIVATE, PARAMETER :: English = 1,        &
                            MESSAGE_LENGHT = 51

CONTAINS

   FUNCTION Message( LANGUAGE, NUMBER )RESULT(text)
  
   INTEGER, INTENT(IN) :: NUMBER, LANGUAGE
   CHARACTER (LEN=MESSAGE_LENGHT) :: text

   SELECT CASE (LANGUAGE)

   CASE(English)
      SELECT CASE (NUMBER)
         
        !
        ! Error messages from the Check_Input Module
        !
        CASE (1)  !  7 < INFORM  < 16384 
           text = "Wrong input parameters for CUBATR.                 "
        CASE (2)  ! INFORM >= 8192
           text = "- This happend during a restart.                   "
        CASE (3)  ! INFORM = 8192 (should never occur)
           text = "-> There is nothing to restart from.               "
        CASE (8)  ! INFORM = 256
           text = "-> Wrong limits on the number of evaluations.      "
        CASE (9)  ! INFORM = 512
           text = "-> Reliability parameter TUNE > 1 or < 0.          "
        CASE (10)  ! INFORM = 1024
           text = "-> Requested accuracy must be nonnegative.         "
        CASE (11)  ! INFORM = 2048
           text = "-> On entry, IFAIL must be -1, 0 or 1.             "
        CASE (12)  ! INFORM = 4096
           text = "-> On entry, |JOB| must be 0, 1, 2 or 11.          "
        !
        ! Error messages when RESTART = false
        !
        CASE (13)   ! INFORM = 8
           text = "-> The DIMension must be at least 1                "
        CASE (14)   ! INFORM = 16
           text = "-> The NUM of FUNctions must be positive.          "
        CASE (15)   ! INFORM = 32
           text = "-> The NUMber of ReGioNs must be positive.         "
        CASE (16)   ! INFORM = 64
           text = "-> No support for this type of region.             "
        CASE (17)   ! INFORM = 128
           text = "-> For DIMension > 3, only JOB = 1 is implemented  "
             !
             ! Error messages when RESTART = true
             !
        CASE (18)  ! INFORM = 8
           text = "-> The DIMension may not change for a restart.     "
        CASE (19)  ! INFORM = 16
           text = "-> The NUM of FUNctions may not be changed.        "
        CASE (20)  ! INFORM = 32
           text = "-> The length of the integer workspace was changed."
        CASE (21)  ! INFORM = 64
           text = "-> The length of the real workspace was changed.   "
        CASE (22)  ! INFORM = 128
           text = "-> Unbelievable low length of integer workspace.   "
             !
             ! Error messages from Global_Adapt
             !
        CASE (23)  ! INFORM = 1
            text = " -> Allowed number of function evaluations reached."
        CASE (24)  ! INFORM = 2
            text = " -> Work space is full.                            "
        CASE (25)  ! INFORM = 3
            text = " -> Can't restart with extrapolation after JOB/=2. "
        CASE (26)  ! INFORM = 4
            text = " -> fail = 4 , no meaning allocated.               "
        CASE (27)  ! INFORM = 5
            text = " -> Rule_: No support for this type of region.     "
        CASE (28)  ! INFORM = 6
            text = " -> Divide: No support for this type of region.    "
        CASE (29)  ! INFORM = 7
            text = " -> Divide: No support for this subdivision.       "
        CASE (30)  ! INFORM >= 16384
            text = " -> DS_Routines: Internal Error => Consult experts."
        CASE DEFAULT
            text = " -> Error_Handling: Error => Consult experts.      "

      END SELECT

   END SELECT

   RETURN
   END FUNCTION Message
!---------------------------------------------------------------------------
   SUBROUTINE Print_Error(INFO)

   INTEGER, INTENT(IN) :: INFO

   INTEGER, PARAMETER :: MAXERR=12
   INTEGER :: I,ERRNUM,MESNR,OFSET,Inform,Language
!
   Language = English
   Inform = INFO
!     IF (NAME == "CUBATR") THEN
   IF (Inform < 2**(MAXERR+1)) THEN
       IF (Inform >= 8) THEN
           write(unit=*,fmt=*) Message(Language,1)
           IF (Inform >= 2** (MAXERR+1)) THEN
               write(unit=*,fmt=*) Message(Language,2)
               Inform = Inform - 2** (MAXERR+1)
               IF (Inform == 0) THEN
                                     write(unit=*,fmt=*) Message(Language,3)
               END IF
               OFSET = 15
           ELSE
               OFSET = 10
           END IF
!
!             Wrong input parameters
!
           ERRNUM = 2**MAXERR
           DO I = MAXERR,3,-1
               IF (Inform >= ERRNUM) THEN
                   Inform = Inform - ERRNUM
                   IF (I >= 8) THEN
                       MESNR = I
                   ELSE
                       MESNR = OFSET + I
                   END IF
                   write(unit=*,fmt=*) Message(Language,MESNR)
               END IF
               ERRNUM = ERRNUM/2
           END DO
       END IF
!
!     Errors from integration procedure
!
       IF (Inform > 0) THEN
           write(unit=*,fmt=*) Message(Language,Inform+22)
       END IF
   ELSE
!
!     Errors from DS_routines
!
       write(unit=*,fmt=*) Message(Language,30)
       write(unit=*,fmt=*) "  Internal error with code ", Inform
   END IF

   RETURN
   END SUBROUTINE Print_Error
!---------------------------------------------------------------------------
   SUBROUTINE Handle_Error(Inform,IFAIL)
   INTEGER, INTENT(IN) :: Inform
   INTEGER, INTENT(IN OUT), OPTIONAL :: IFAIL
   IF (Inform /= 0) THEN
      IF (PRESENT(IFAIL)) THEN
         SELECT CASE (IFAIL)
         CASE (1)
              ! Soft silent error: do nothing
              IFAIL = Inform
              RETURN
         CASE (-1)
              ! Soft noisy error
              CALL Print_Error(Inform)
              IFAIL = Inform
              RETURN
         CASE (0)
              ! Hard noisy error
              CALL Print_Error(Inform)
              write(unit=*,fmt=*) "Hard error: program terminated"
              STOP  ! "Hard error: program terminated"
         END SELECT
      ELSE
         ! Soft noisy error
         CALL Print_Error(Inform)
         RETURN
      END IF
   ELSE
      IF (PRESENT(IFAIL)) THEN
         IFAIL = 0
      END IF
   END IF
   RETURN
   END SUBROUTINE Handle_Error

END MODULE Error_Handling
