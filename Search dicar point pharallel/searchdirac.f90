      Program search
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,npartitions=20,nkx=50,nky=50,nkz=50
      integer,parameter::nkpt=(2*nkx+1)*(2*nky+1)*(2*nkz+1)
!------------------------------------------------------
      integer*4 ik,ikmax,ipart,ikx, iky, ikz, gaploc
      real*8 ef(npartitions), alpha, gap(npartitions),kx, ky, kz
      character(len=30)::klabel(nkpath), partnumber
      character(len=80) hamil_file_triv,hamil_file_top,nnkp,line, partition_file
      integer*4 i,j,k,nr,i1,i2,nb,lwork,info
      real*8,parameter::third=1d0/3d0
      real*8 phase,twopi,jk,a,b,dk(3),khbox(3), r1, r2, r3
      real*8 korigin(3), ktemp(3),bvec(3,3),avec(3,3),klist(3,nkpt)
      real*8,allocatable:: rvec(:,:),ene(:,:),rwork(:)
      integer*4,allocatable:: ndeg(:)
      complex*16,allocatable:: Hk(:,:),Hk_triv(:,:),Hk_top(:,:),Hr_triv(:,:,:), &
      Hr_top(:,:,:),work(:),Hr_alpha(:,:)

!------------------------------------------------------
      write(hamil_file_triv,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
      write(hamil_file_top,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

      twopi=4.0d0*atan(1.0d0)*2.0d0
      

!---------------  reciprocal vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
110   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin real_lattice") goto 110
    
      read(98,*)avec
      
111   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      
      read(98,*)bvec
      close(98)
      
!---------defining k point previously giving the minimum energy      
      !data korigin / 0.104720, 0.104720, 3.141593 /
      data korigin / 0.104720, 0.104720, 3.141593 /
      ktemp = korigin/twopi
      korigin(1) = ktemp(1) * bvec(1,1) + ktemp(2) * bvec (1,2) + ktemp(3) * bvec (1,3)
      korigin(2) = ktemp(1) * bvec(2,1) + ktemp(2) * bvec (2,2) + ktemp(3) * bvec (2,3)
      korigin(3) = ktemp(1) * bvec(3,1) + ktemp(2) * bvec (3,2) + ktemp(3) * bvec (3,3)
      data khbox / 0.1d0, 0.1d0, 0.1d0 /
      
!---------------kmesh
      open(777,file='kmesh.dat')
      i = 0
      do ikx=-nkx,nkx
            do iky=-nky,nky
                do ikz=-nkz,nkz
                    i = i + 1
                    klist(1,i) = (float(ikx)/float(nkx))*khbox(1) + korigin(1)
                    klist(2,i) = (float(iky)/float(nky))*khbox(2) + korigin(2)
                    klist(3,i) = (float(ikz)/float(nkz))*khbox(3) + korigin(3)
                    write(777, '(3(x,f12.8))') klist(1,i), klist(2,i), klist(3,i)
                enddo 
            enddo 
      enddo

!------read H(R) trivial
      open(99,file=trim(adjustl(hamil_file_triv)),err=444)
      read(99,*)
      read(99,*)nb,nr
      allocate(rvec(3,nr),Hk_triv(nb,nb),Hr_triv(nb,nb,nr),ndeg(nr),ene(nb,nkpt), &
      Hk(nb,nb)) !ene just storing( temporary) array - same for two 
      read(99,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)r1,r2,r3,i1,i2,a,b
               rvec(1,k)=r1*avec(1,1) + r2*avec(1,2) + r3*avec(1,3)
               rvec(2,k)=r1*avec(2,1) + r2*avec(2,2) + r3*avec(2,3)
               rvec(3,k)=r1*avec(3,1) + r2*avec(3,2) + r3*avec(3,3)
               Hr_triv(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
      close(99)
      
!-------check the r point
      open(888,file='rpoint.dat')
      do k=1,nr
          write(888, '(3(x,f12.8))') rvec(1,k), rvec(2,k), rvec(3,k)
      enddo
      close(888)
      
!------read H(R) topological
      open(99,file=trim(adjustl(hamil_file_top)),err=445)
      read(99,*)
      read(99,*)nb,nr
      allocate(Hk_top(nb,nb),Hr_top(nb,nb,nr),Hr_alpha(nb,nb)) !ndeg is weight matrix - same for both
      read(99,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)r1,r2,r3,i1,i2,a,b
               Hr_top(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
      enddo
      close(99)
      
!-----parameters for LAPACK diagonalisation routine to work
     lwork=max(1,2*nb-1)
     allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

!---- Fourrier transform H(R) to H(k)

    open(101, file='gap.dat')
    ene=0d0
    do ipart=1,npartitions
      !alpha=(float(ipart-1)/float(npartitions-1)*0.3d0) + 0.6d0 !offset and scaling for closer interval
      alpha=(float(ipart-1)/float(npartitions-1)*0.15d0) + 0.7d0
      !alpha=float(ipart-1)/float(npartitions-1) !- full range
      print*, 'itteration',ipart
          do k=1,nkpt
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
          !-----------calculating fermi level and minimum band gap
          gap(ipart)=minval(ene(13,:)-ene(12,:))
          gaploc = minloc(ene(13,:),dim=1)
          if (ipart==1) write(2234, '(f12.6)') ene(13,:)
          if (ipart==npartitions) write(2235, '(f12.6)') ene(13,:)
          !-----------writing band gap results------------------
          write(101, '(2(x,f12.6), x, i6, 3(x,f12.6))') alpha,gap(ipart),gaploc, &
          klist(1, gaploc),klist(2, gaploc),klist(3, gaploc)

    enddo
    stop
      
!---------error traces---------------------------------------      
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_triv)),' not found'
      stop
445   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_top)),' not found'
      stop

      end program search
      


