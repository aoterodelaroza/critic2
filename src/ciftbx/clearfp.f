      subroutine clearfp
C
C     This is a null version of a
C     subroutine to clear IEEE floating point exceptions
C     for inexact and underflow under SUN OS 4 f77
C     For most other systems, no action is needed.
C
C     character*1 out
C     ii = ieee_flags('clear','exception','underflow',out)
C     ii = ieee_flags('clear','execption','inexact',out)
      return
      end
