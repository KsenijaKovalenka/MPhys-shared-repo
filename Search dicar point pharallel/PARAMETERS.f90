	MODULE PARAMETERS
	IMPLICIT NONE
!	INCLUDE 'mpif.h'
!   MPI varibales
	INTEGER IERR,MYID,NUMPROCS
!   input file units
	INTEGER IUH, INNKP,IIN
!   k-points & bands
	INTEGER NBAND,NKX,NKY,NKZ,NKPT,NBMIN,NBMAX
	REAL*8  KX_FRAC,KY_FRAC,KZ_FRAC
	REAL*8  KX_OFFSET,KY_OFFSET,KZ_OFFSET
!   real-space sampling
	INTEGER NRPTS
	INTEGER,ALLOCATABLE::NDEG(:),IRVEC(:,:)
!   real-space hamiltonian
	COMPLEX*16,ALLOCATABLE::HAMR(:,:,:)
!   Fermi level
	REAL*8 EFERMI
!   lattice parameters
	REAL*8 AVEC(3,3),BVEC(3,3)
!   prefix of the wannier outputs
	CHARACTER(LEN=100)PREFIX
!   work arrays
        REAL(8),ALLOCATABLE::EIGEN(:,:),VELOC(:,:),MAGNET(:,:,:)
        INTEGER,ALLOCATABLE:: IKMAP(:,:,:)
!   control keys for plotting velocity and magnetization
        LOGICAL PLOT_VELOC, PLOT_MAGNET
	END MODULE PARAMETERS

!============================================================
! initialization
!============================================================
        SUBROUTINE INIT_PARAM
        USE CONSTANT          ,          ONLY: PI2
        USE PARAMETERS         ,          ONLY:  IERR,MYID,NUMPROCS,IIN,             &
                                                 IUH,INNKP,PREFIX,NBAND,NRPTS,       &
                                                 NKX,NKY,NKZ,NKPT,NBMIN,NBMAX,       &
                                                 KX_FRAC,KY_FRAC,KZ_FRAC,            &
                                                 KX_OFFSET,KY_OFFSET,KZ_OFFSET,      &
                                                 EIGEN,IKMAP,PLOT_VELOC,PLOT_MAGNET
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        CHARACTER(LEN=100)LINE,POSLINE,NNKPLINE
        CHARACTER(LEN=300)TOKEN,VALUE
        LOGICAL CASE_HR,CASE_POS,CASE_NNKP,CASE_IN
! initialize input units
        IUH  = 99
        INNKP= 98
        IIN  = 97
! open master input 
        INQUIRE(FILE='INPUT',EXIST=CASE_IN)
        IF(.NOT.CASE_IN) GO TO 149
        OPEN(IIN,FILE='INPUT',STATUS='OLD',ERR=149)
!============================================================
!  READ MANDATORY PARAMETERS
!============================================================
! 1. read prfix and inquire the related files
        TOKEN='PREFIX'
        CALL GETPARAM(TOKEN,VALUE,IIN)
	
        PREFIX=TRIM(ADJUSTL(VALUE))
        WRITE(LINE,'(2A)')TRIM(ADJUSTL(PREFIX)),'_hr.dat'
        INQUIRE(FILE=TRIM(ADJUSTL(LINE)),EXIST=CASE_HR)
        IF(.NOT.CASE_HR) GO TO 150
        WRITE(NNKPLINE,'(2A)')TRIM(ADJUSTL(PREFIX)),'.nnkp'
        INQUIRE(FILE=TRIM(ADJUSTL(NNKPLINE)),EXIST=CASE_NNKP)
        IF(.NOT.CASE_NNKP) GO TO 152

! 2. read NKX, NKY and NKZ and set NKPT
        TOKEN='NKX'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
	
           READ(VALUE,*,ERR=154) NKX
        ELSE
           GO TO 154
        ENDIF
        TOKEN='NKY'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
           READ(VALUE,*,ERR=155) NKY
	
        ELSE
           GO TO 155
        ENDIF
        TOKEN='NKZ'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
           READ(VALUE,*,ERR=156) NKZ
	
        ELSE
           GO TO 156
        ENDIF
! 3. read NBMIN and NBMAX
        TOKEN='NBMIN'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
           READ(VALUE,*,ERR=157) NBMIN
	
        ELSE
           GO TO 157
        ENDIF

        TOKEN='NBMAX'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
           READ(VALUE,*,ERR=158) NBMAX
	
        ELSE
           GO TO 158
        ENDIF
        NKPT=(2*NKX+1)*(2*NKY+1)*(2*NKZ+1)
        ALLOCATE(IKMAP(-NKX:NKX,-NKY:NKY,-NKZ:NKZ))
!============================================================
!  READ OPTIONAL PARAMETERS
!============================================================
! 1. read k-offsets
        TOKEN='KX_OFFSET'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
           READ(VALUE,*,ERR=159) KX_OFFSET
           KX_OFFSET=KX_OFFSET*PI2
	
        ELSE
           KX_OFFSET=0D0
        ENDIF

        TOKEN='KY_OFFSET'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
           READ(VALUE,*,ERR=160) KY_OFFSET
           KY_OFFSET=KY_OFFSET*PI2
	
        ELSE
           KY_OFFSET=0D0
        ENDIF

        TOKEN='KZ_OFFSET'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
           READ(VALUE,*,ERR=161) KZ_OFFSET
           KZ_OFFSET=KZ_OFFSET*PI2
	
        ELSE
           KZ_OFFSET=0D0
        ENDIF
! 2. read k fractions
 
        TOKEN='KX_FRAC'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
	
           READ(VALUE,*,ERR=162) KX_FRAC
        ELSE
           KX_FRAC=0D0
        ENDIF

        TOKEN='KY_FRAC'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
	
           READ(VALUE,*,ERR=163) KY_FRAC
        ELSE
           KY_FRAC=0D0
        ENDIF

        TOKEN='KZ_FRAC'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
	
           READ(VALUE,*,ERR=164) KZ_FRAC
        ELSE
           KZ_FRAC=0D0
        ENDIF
! 3. read control keys plot_veloc and plot_magnet
        TOKEN='PLOT_VELOC'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN
        
           READ(VALUE,*,ERR=165) PLOT_VELOC
        ELSE
           PLOT_VELOC=.FALSE.
        ENDIF

        TOKEN='PLOT_MAGNET'
        CALL GETPARAM(TOKEN,VALUE,IIN)
        IF(LEN_TRIM(VALUE).NE.0) THEN

           READ(VALUE,*,ERR=166) PLOT_MAGNET
        ELSE
           PLOT_MAGNET=.FALSE.
        ENDIF
!============================================================
!  Done with INPUT; close it
!============================================================
        CLOSE(IIN)
!============================================================
!  READ OTHER PARAMETERS in case_hr.dat and case.nnk
!============================================================
        OPEN(IUH  ,FILE=TRIM(ADJUSTL(LINE)),STATUS='OLD',ERR=150)
        OPEN(INNKP,FILE=TRIM(ADJUSTL(NNKPLINE)),STATUS='OLD',ERR=152)
! FIND THE NUMBER OF BANDS and R-vectors
        READ(IUH,*)
        READ(IUH,*,END=153)NBAND
        READ(IUH,*,END=153)NRPTS
        REWIND(IUH)
        RETURN
! Error messages
 149    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '  FATAL ERROR: MASTER FILE "INPUT" NOT FOUND.  ', &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 150    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,3A,/,A,/,A)') &
          '-----------------------------------------------', &
          '  FATAL ERROR: ',TRIM(ADJUSTL(LINE)),' NOT FOUND.', &
          '                 PROGRAM ABORTED               ', & 
          '-----------------------------------------------' 
        ENDIF
        STOP
 152    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,3A,/,A,/,A)') &
          '-----------------------------------------------', &
          '  FATAL ERROR: ',TRIM(ADJUSTL(NNKPLINE)),' NOT FOUND.', &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------' 
        ENDIF
        STOP
 153    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,3A,/,A,/A)') &
          '-----------------------------------------------', &
          '  FATAL ERROR: REACHED END OF ',TRIM(ADJUSTL(LINE)),'.', &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------' 
        ENDIF
        STOP
 154    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '            FATAL ERROR: INVALID NKX.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 155    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '            FATAL ERROR: INVALID NKY.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 156    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '            FATAL ERROR: INVALID NKZ.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 157    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '           FATAL ERROR: INVALID NBMIN.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 158    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '           FATAL ERROR: INVALID NBMAX.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 159    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '         FATAL ERROR: INVALID KX_OFFSET.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 160    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '         FATAL ERROR: INVALID KY_OFFSET.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 161    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '         FATAL ERROR: INVALID KZ_OFFSET.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 162    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '          FATAL ERROR: INVALID KX_FRAC.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 163    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '          FATAL ERROR: INVALID KY_FRAC.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 164    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '          FATAL ERROR: INVALID KZ_FRAC.  ',         &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 165    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '         FATAL ERROR: INVALID PLOT_VELOC.  ',       &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
 166    IF(MYID.eq.0) THEN
          WRITE(*,'(A,/,A,/,A,/,A)') &
          '-----------------------------------------------', &
          '        FATAL ERROR: INVALID PLOT_MAGNET.  ',       &
          '                 PROGRAM ABORTED               ', &
          '-----------------------------------------------'
        ENDIF
        STOP
        END SUBROUTINE INIT_PARAM


