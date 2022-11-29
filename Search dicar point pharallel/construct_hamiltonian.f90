      SUBROUTINE construct_hamiltonian
      USE CONSTANTS          ,          ONLY:  TWOPI !unused but may be usefull later?
      USE PARAMETERS         ,          ONLY:  IERR,MYID,NUMPROCS,               & !for mpi
                                             IUH_TRIVIAL,IUH_TOPOLOGICAL,INNKP,  & !file indicies
                                             PREFIX,NBAND,NRPTS,                 & !basic stuff
                                             NKX,NKY,NKZ,NKPT,NBMIN,NBMAX,       & !bounds for kmesh
                                             KX_HBOX,KY_HBOX,KZ_HBOX,            & !kmesh
                                             KX_ORIGIN,KY_ORIGIN,KZ_ORIGIN,      & !initial k-point in lattice basis
                                             AVEC, BVEC, RVEC                      !basis + rvector array
                                             
                                             
      IMPLICIT NONE
      INCLUDE 'mpif.h'

!-------------------local variables----------------------------
      integer*4 npartitions=10                                                !may move to init
      integer*4 ipart, ikx, iky, ikz, ECOUNTS
      real*8  alpha, gap(npartitions), ktemp(3)
      character(len=100):: line, hamil_file_triv, hamil_file_top
      integer*4 i,j,k,i1,i2
      real*8 phase,r1,r2,r3
      real*8,allocatable:: rwork(:)
      INTEGER*4 lwork,lrwork,info
      COMPLEX*16,ALLOCATABLE:: work(:),Hk_triv(:,:),Hk_top(:,:),Hk(:,:), &
      Hr_triv(:,:,:),Hr_top(:,:,:),Hr_alpha(:,:)
      real*8,allocatable:: ENE(:,:), EIGEN(:,:), klist(:,:)      

!---------------  reciprocal vectors
      open(INNKP,file=trim(adjustl(nnkp)),err=333)
110   read(INNKP,'(a)')line
      if(trim(adjustl(line)).ne."begin real_lattice") goto 110
    
      read(INNKP,*)AVEC
      
111   read(INNKP,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      
      read(INNKP,*)BVEC
      close(INNKP)
      
!  set up work arrays for ZHEEV
      lwork  = MAX(1,2*NBAND-1)
      lrwork = MAX(1,3*NBAND-1)
      ALLOCATE(work(lwork),rwork(lrwork),Hk(NBAND,NBAND))

!  distribute k-points among the avaialble processors
!  first check if NKPT is divisible by NUMPROCS
      IF (MOD(NKPT,NUMPROCS).EQ.0) THEN
         KNUM = NKPT / NUMPROCS  !integer division, number of k-points in the batch
         KMIN = 1 + KNUM * MYID
         KMAX = KNUM * (MYID + 1)
         ALLOCATE(ENE(NBAND,1:KNUM))
         ALLOCATE(EIGEN(NBAND,NKPT))

      ELSE
         KNUM = (NKPT / NUMPROCS)+1          !number of k-points in the batch
         KMIN = 1 + KNUM * MYID
         KMAX = KNUM * (MYID + 1)
         ALLOCATE(ENE(NBAND,1:KNUM))
         ALLOCATE(EIGEN(NBAND,KNUM*NUMPROCS))

      ENDIF 
      
!---------defining k point which previously gave the minimum energy in carthesian co-ords     
      data ktemp / KX_ORIGIN, KY_ORIGIN, KZ_ORIGIN /
      KX_ORIGIN = ktemp(1) * BVEC(1,1) + ktemp(2) * BVEC (1,2) + ktemp(3) * BVEC (1,3)
      KY_ORIGIN = ktemp(1) * BVEC(2,1) + ktemp(2) * BVEC (2,2) + ktemp(3) * BVEC (2,3)
      KZ_ORIGIN = ktemp(1) * BVEC(3,1) + ktemp(2) * BVEC (3,2) + ktemp(3) * BVEC (3,3)
             
          
!  generate a uniform 3D k-mesh
      ALLOCATE(KLIST(3,NKPT))
      K=0
      DO IKX=-NKX,NKX
       DO IKY=-NKY,NKY
        DO IKZ=-NKZ,NKZ
         K=K+1
         klist(1,K) = (float(ikx)/float(nkx))*KX_HBOX + KX_ORIGIN
         klist(2,K) = (float(iky)/float(nky))*KY_HBOX + KY_ORIGIN
         klist(3,K) = (float(ikz)/float(nkz))*KZ_HBOX + KZ_ORIGIN
        ENDDO
       ENDDO
      ENDDO 

!-----define file units for trivial and topological files
      IUH_TRIVIAL = 96
      IUH_TOPOLOGICAL = 95

!--------for now have the two filenames       
      write(hamil_file_triv,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
      write(hamil_file_top,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"

!------read H(R) trivial
      open(IUH_TRIVIAL,file=trim(adjustl(hamil_file_triv)),err=444)
      read(IUH_TRIVIAL,*)
      read(IUH_TRIVIAL,*)
      read(IUH_TRIVIAL,*)
!      read(IUH_TRIVIAL,*)NBANDS,nr
      allocate(RVEC(3,NRPTS),Hk_triv(NBANDS,NBANDS),Hr_triv(NBANDS,NBANDS,NRPTS),ndeg(NRPTS)) 
      ! read the weighting array
      read(IUH_TRIVIAL,*)ndeg
      do k=1,NRPTS
         do i=1,NBANDS
            do j=1,NBANDS
               read(IUH_TRIVIAL,*)r1,r2,r3,i1,i2,a,b
               RVEC(1,k)=r1*AVEC(1,1) + r2*AVEC(1,2) + r3*AVEC(1,3)
               RVEC(2,k)=r1*AVEC(2,1) + r2*AVEC(2,2) + r3*AVEC(2,3)
               RVEC(3,k)=r1*AVEC(3,1) + r2*AVEC(3,2) + r3*AVEC(3,3)
               Hr_triv(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
      close(IUH_TRIVIAL)
      
!------read H(R) topological
      open(IUH_TOPOLOGICAL,file=trim(adjustl(hamil_file_top)),err=445)
      read(IUH_TOPOLOGICAL,*)
      read(IUH_TOPOLOGICAL,*)
      read(IUH_TOPOLOGICAL,*)
!      read(IUH_TOPOLOGICAL,*)NBANDS,nr
      allocate(Hk_top(NBANDS,NBANDS),Hr_top(NBANDS,NBANDS,NRPTS),Hr_alpha(NBANDS,NBANDS)) !ndeg is weight matrix - same for both
      read(IUH_TOPOLOGICAL,*)ndeg
      do k=1,NRPTS
         do i=1,NBANDS
            do j=1,NBANDS
               read(IUH_TOPOLOGICAL,*)r1,r2,r3,i1,i2,a,b
               Hr_top(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
      close(IUH_TOPOLOGICAL)


!---- Fourrier transform H(R) to H(k)

    ENE=0d0 !why this one is initialised here and the other inside the loop?
    do ipart=1,npartitions
      alpha=(float(ipart-1)/float(npartitions-1)*0.15d0) + 0.7d0 !offset and scaling for closer interval
          EIGEN=0D0
          !ENE=0d0 ?
          do K=KMIN,MAX(KMAX, NKPT)
             Hk=(0d0,0d0)
             do j=1,NRPTS
                phase=0.0d0
                do i=1,3
                   phase=phase+klist(i,k)*RVEC(i,j)
                enddo
                Hr_alpha = alpha*Hr_top(:,:,j)+(1d0-alpha)*Hr_triv(:,:,j)
                Hk=Hk+Hr_alpha*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))

             enddo
             !-----------find energies---------------------------------
             call zheev('V','U',NBANDS,Hk,NBANDS,ENE(:,k),work,lwork,rwork,info)
          enddo

          !  gather the information corresponding to the distributed k-points
          ECOUNTS=KNUM*NBAND
          CALL MPI_GATHER(ENE,ECOUNTS,MPI_DOUBLE_PRECISION,   &
                      EIGEN,ECOUNTS,MPI_DOUBLE_PRECISION,     &
                      0,MPI_COMM_WORLD,IERR)
                      
          GAP(IPART) = MINVAL(EIGEN(13,:)) - MAXVAL(EIGEN(12,:))
          !EF(IPART) = (MINVAL(EIGEN(13,:)) + MAXVAL(EIGEN(12,:)))/2D0
         
    enddo
    IF(MYID.NE.0) DEALLOCATE(EIGEN)
    DEALLOCATE(klist,Hr_triv,Hr_top,Hr_alpha,Hk_triv,Hk_top,Hk,RVEC) !ndeg?
    DEALLOCATE(work,rwork,ENE)
    RETURN    
    
      
!---------error traces--------------------------------------- 

444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_triv)),' not found'
      stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_top)),' not found'
      stop

      END SUBROUTINE construct_hamiltonian
      


