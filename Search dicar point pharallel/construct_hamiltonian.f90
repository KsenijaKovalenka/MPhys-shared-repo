      SUBROUTINE construct_hamiltonian
      USE CONSTANTS          ,          ONLY:  TWOPI
      USE PARAMETERS         ,          ONLY:  IERR,MYID,NUMPROCS,IIN,           &
                                             IUH_TRIVIAL,IUH_TOPOLOGICAL,INNKP,  & 
                                             PREFIX,NBAND,NRPTS,NPARTITIONS,     &
                                             NKX,NKY,NKZ,NKPT,NBMIN,NBMAX,       &
                                             KX_HBOX,KY_HBOX,KZ_HBOX,            &
                                             KX_ORIGIN,KY_ORIGIN,KZ_ORIGIN,      &
                                             avec, bvec, RVEC
                                             
                                             
      IMPLICIT NONE
      INCLUDE 'mpif.h'

!------------------------------------------------------
      integer*4 ik,ikmax,ipart,ikx, iky, ikz, gaploc, ECOUNTS
      real*8  alpha, kx, ky, kz, ktemp(3)
      character(len=30):: partnumber
      character(len=100):: line
      integer*4 i,j,k,nr,i1,i2,nb
      real*8 phase,jk,a,b,dk(3), r1, r2, r3
      real*8 korigin(3)
      real*8,allocatable:: rwork(:)
      INTEGER*4 LWORK,LRWORK,INFO
      COMPLEX*16,ALLOCATABLE:: WORK(:),Hk_triv(:,:),Hk_top(:,:),Hk(:,:),HR_triv(:,:,:),HR_top(:,:,:)
      real*8,allocatable:: ene(:,:), EIGEN(:,:)      

!---------------  reciprocal vectors
      open(INNKP,file=trim(adjustl(nnkp)),err=333)
110   read(INNKP,'(a)')line
      if(trim(adjustl(line)).ne."begin real_lattice") goto 110
    
      read(INNKP,*)avec
      
111   read(INNKP,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      
      read(INNKP,*)bvec
      close(INNKP)
      
!  set up work arrays for ZHEEV
      LWORK  = MAX(1,2*NBAND-1)
      LRWORK = MAX(1,3*NBAND-1)
      ALLOCATE(WORK(LWORK),RWORK(LRWORK),HK(NBAND,NBAND))
!  distribute k-points among the avaialble processors
!  first check if NKPT is divisible by NUMPROCS
      IF (MOD(NKPT,NUMPROCS).EQ.0) THEN
         KNUM = NKPT / NUMPROCS  ! integer division 
         KMIN = 1 + KNUM * MYID
         KMAX = KNUM * (MYID + 1)
         ALLOCATE(ENE(NBAND,1:KNUM))
         ALLOCATE(EIGEN(NBAND,NKPT))

      ELSE
         KNUM = (NKPT / NUMPROCS)+1
         KMIN = 1 + KNUM * MYID
         KMAX = KNUM * (MYID + 1)
         ALLOCATE(ENE(NBAND,1:KNUM))
         ALLOCATE(EIGEN(NBAND,KNUM*NUMPROCS))

      ENDIF 
      
!---------defining k point previously giving the minimum energy      

      ktemp = korigin/twopi
      korigin(1) = ktemp(1) * bvec(1,1) + ktemp(2) * bvec (1,2) + ktemp(3) * bvec (1,3)
      korigin(2) = ktemp(1) * bvec(2,1) + ktemp(2) * bvec (2,2) + ktemp(3) * bvec (2,3)
      korigin(3) = ktemp(1) * bvec(3,1) + ktemp(2) * bvec (3,2) + ktemp(3) * bvec (3,3)
      data khbox / 0.1d0, 0.1d0, 0.1d0 /
             
          
!  generate a uniform 3D k-mesh
      ALLOCATE(KLIST(3,NKPT))
      K=0
      DO IKX=-NKX,NKX
       DO IKY=-NKY,NKY
        DO IKZ=-NKZ,NKZ
         K=K+1
         klist(1,K) = (float(ikx)/float(nkx))*KX_HBOX + korigin(1)
         klist(2,K) = (float(iky)/float(nky))*KY_HBOX + korigin(2)
         klist(3,K) = (float(ikz)/float(nkz))*KZ_HBOX + korigin(3)
        ENDDO
       ENDDO
      ENDDO      


!------read H(R) trivial
      read(IUH_TRIVIAL,*)
      read(IUH_TRIVIAL,*)nb,nr
      allocate(rvec(3,nr),Hk_triv(nb,nb),Hr_triv(nb,nb,nr),ndeg(nr),ene(nb,nkpt), &
      Hk(nb,nb)) !ene just storing( temporary) array - same for two 
      read(IUH_TRIVIAL,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(IUH_TRIVIAL,*)r1,r2,r3,i1,i2,a,b
               rvec(1,k)=r1*avec(1,1) + r2*avec(1,2) + r3*avec(1,3)
               rvec(2,k)=r1*avec(2,1) + r2*avec(2,2) + r3*avec(2,3)
               rvec(3,k)=r1*avec(3,1) + r2*avec(3,2) + r3*avec(3,3)
               Hr_triv(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
      close(IUH_TRIVIAL)
      
!-------check the r point
      open(888,file='rpoint.dat')
      do k=1,nr
          write(888, '(3(x,f12.8))') rvec(1,k), rvec(2,k), rvec(3,k)
      enddo
      close(888)
      
!------read H(R) topological
      open(IUH_TOPOLOGICAL,file=trim(adjustl(hamil_file_top)),err=445)
      read(IUH_TOPOLOGICAL,*)
      read(IUH_TOPOLOGICAL,*)nb,nr
      allocate(Hk_top(nb,nb),Hr_top(nb,nb,nr),Hr_alpha(nb,nb)) !ndeg is weight matrix - same for both
      read(IUH_TOPOLOGICAL,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(IUH_TOPOLOGICAL,*)r1,r2,r3,i1,i2,a,b
               Hr_top(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
      close(IUH_TOPOLOGICAL)


!---- Fourrier transform H(R) to H(k)

    ene=0d0
    do ipart=1,npartitions
      !alpha=(float(ipart-1)/float(npartitions-1)*0.3d0) + 0.6d0 !offset and scaling for closer interval
      alpha=(float(ipart-1)/float(npartitions-1)*0.15d0) + 0.7d0
          EIGEN=0D0
          do K=KMIN,MAX(KMAX, NKPT)
             HK=(0d0,0d0)
             do j=1,nr
                phase=0.0d0
                do i=1,3
                   phase=phase+klist(i,k)*rvec(i,j)
                enddo
                Hr_alpha = alpha*Hr_top(:,:,j)+(1d0-alpha)*Hr_triv(:,:,j)
                !Hr_alpha = Hr_top(:,:,j)
                HK=HK+Hr_alpha*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))

             enddo
             !-----------find energies---------------------------------
             call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)
          enddo
!  gather the infor corresponding to the distributed k-points
          ECOUNTS=KNUM*NBAND
          CALL MPI_GATHER(ENE  ,ECOUNTS,MPI_DOUBLE_PRECISION,   &
                      EIGEN,ECOUNTS,MPI_DOUBLE_PRECISION, &
                      0,MPI_COMM_WORLD,IERR)
                      
          GAP(IPART) = MINVAL(EIGEN(13,:)) - MAXVAL(EIGEN(12,:))
          EF(IPART) = (MINVAL(EIGEN(13,:)) + MAXVAL(EIGEN(12,:)))/2D0
         
    enddo
    IF(MYID.NE.0) DEALLOCATE(EIGEN)
    DEALLOCATE(KLIST,HR_TRIV,HR_TOP,HK_TRIV,HK_TOP,HK,RVEC)
    DEALLOCATE(WORK,RWORK,ENE)
    RETURN    
    
      
!---------error traces--------------------------------------- 
     
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_triv)),' not found'
      stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_top)),' not found'
      stop

      end SUBROUTINE construct_hamiltonian
      


