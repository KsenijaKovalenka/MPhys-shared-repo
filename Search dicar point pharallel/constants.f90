MODULE CONSTANTS
        IMPLICIT NONE
        REAL*8 KBOLTZ,TWOPI,PI
        END MODULE CONSTANTS
        SUBROUTINE INIT_CONSTANTS
        USE CONSTANTS
        IMPLICIT NONE
        INCLUDE ‘mpif.h’
         KBOLTZ = 8.617D-5
         PI     = 4.0D0*ATAN(1.0D0)
         TWOPI    = PI*2D0
        END SUBROUTINE INIT_CONSTANTS
        