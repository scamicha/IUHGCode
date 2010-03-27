!-*-f90-*-     
      PROGRAM Ams
      IMPLICIT NONE

      INCLUDE 'mpif.h'
      
      INTEGER,PARAMETER :: DOUBLE = SELECTED_REAL_KIND(15,300)
      INTEGER,PARAMETER :: jmax=512,kmax=64,lmax=128
      INTEGER,PARAMETER :: start = 110200
      INTEGER,PARAMETER :: finish = 250000
      INTEGER,PARAMETER :: skip = 200
      INTEGER,PARAMETER :: size = ((finish-start)/skip)+1
      INTEGER,PARAMETER :: jstart = 90

      INTEGER :: jmax2,kmax2,mmax,amcount,process
      INTEGER :: count,j,k,l,m,i,m_return,sender
      INTEGER :: numranks,mpierr,myrank
      INTEGER :: status(MPI_STATUS_SIZE)
      REAL(DOUBLE),PARAMETER :: tconv = 1605.63
      REAL(DOUBLE),PARAMETER :: pi = 3.14159265358979323846d0
      REAL(DOUBLE),PARAMETER :: twopi = 2.d0*pi
      REAL(DOUBLE),DIMENSION(:,:,:),ALLOCATABLE    :: rho
      REAL(DOUBLE),DIMENSION(:,:),ALLOCATABLE      :: am,bm
      REAL(DOUBLE)           :: a0tot
      REAL(DOUBLE),DIMENSION(:),ALLOCATABLE        :: avgaslice,a0slice
      REAL(DOUBLE),DIMENSION(:,:),ALLOCATABLE      :: avga,avgamid,a0
      REAL(DOUBLE),DIMENSION(size) :: timearr
      REAL(DOUBLE) :: time,phi,dphi
      LOGICAL EXISTSTAT
      CHARACTER outfile*80,indir*80
      CHARACTER rhofile*80,amfile*80,filenum*8,str*80

      type answer_return
         sequence
         REAL(DOUBLE) :: a, amid
      end type answer_return

      type (answer_return) answer

      outfile = "indiram15AU.dat"
      indir = "../WAN_RHO/"

      call MPI_INIT(mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numranks, mpierr)

      mmax   = LMAX/2
      jmax2 = jmax+2
      kmax2 = kmax+2
      dphi  = twopi/lmax

      ALLOCATE(rho(jmax2,kmax2,lmax))
 
      avga(:,:)          = 0.d0
      avgamid(:,:)       = 0.d0
      timearr(:)         = 0.d0
 
      IF(myrank.eq.0) THEN

         ALLOCATE(avga(LMAX/2,size))
         ALLOCATE(avgamid(LMAX/2,size))

         count = 0
         write(filenum,'(I6.6)')start
         rhofile=trim(indir)//'rho3d.'//filenum
         OPEN(UNIT=12,FILE=trim(rhofile),FORM='UNFORMATTED')

         READ(12) rho
         READ(12) time
         CLOSE(12)

         count = count+1
         timearr(count) = time/tconv

         DO i=start,finish,skip

            AMCOUNT = 1
            
            CALL MPI_BCAST(RHO,JMAX2*KMAX2*LMAX,MPI_DOUBLE_PRECISION,&
                 0,MPI_COMM_WORLD,mpierr)

            DO PROCESS=1,min(numranks-1,mmax)
               CALL MPI_SEND(AMCOUNT,1,MPI_INTEGER,PROCESS,AMCOUNT,&
                    MPI_COMM_WORLD,mpierr)
               AMCOUNT = AMCOUNT+1
            ENDDO

911         CONTINUE

            IF(I.lt.finish) THEN
               WRITE(filenum,'(I6.6)')i
               rhofile = indir//"rho3d."//filenum
         
               INQUIRE(FILE=rhofile,EXIST=EXISTSTAT)
         
               IF(.not.EXISTSTAT) THEN
                  print*,"file ",rhofile, "does not exist"
                  STOP
!                  I = I+1
!                  GO TO 911
               ENDIF
         
!         print*," AMS OUT -> OPENING FILE: ", rhofile
               OPEN(UNIT=12,FILE=rhofile,FORM="UNFORMATTED")

               READ(12) rho
               READ(12) time
               CLOSE(12)
         
!         print*,' AMS OUT -> READ FILE: ',rhofile
               count = count+1
               timearr(count) = time/tconv
            ENDIF

            DO PROCESS=1,MMAX+1
               call MPI_RECV(answer,2,MPI_DOUBLE_PRECISION,&
                    MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,&
                    mpierr)
               
               sender = status(MPI_SOURCE)
               m_return = status(MPI_TAG)

               avga(m_return,count) = answer%a
               avgamid(m_return,count) = answer%amid
               
               IF(AMCOUNT.lt.MMAX) THEN
                  call MPI_SEND(AMCOUNT,1,MPI_INTEGER,sender,AMCOUNT,&
                       MPI_COMM_WORLD,mpierr)
                  AMCOUNT = AMCOUNT+1
               ELSE
                  call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER,sender,MMAX+1,&
                       MPI_COMM_WORLD,mpierr)
               ENDIF
            ENDDO
         ENDDO

      ELSE
         do I=start,finish,skip
            CALL MPI_BCAST(RHO,JMAX2*KMAX2*LMAX,MPI_DOUBLE_PRECISION,&
                 0,MPI_COMM_WORLD,mpierr)
            IF (myrank.lt.mmax) THEN

300            CALL MPI_RECV(AMCOUNT,1,MPI_INTEGER,0,MPI_ANY_TAG,&
                    MPI_COMM_WORLD,status,mpierr)
               
               IF (status(MPI_TAG).eq.MMAX+1) go to 400

               ALLOCATE(a0(jmax2,kmax2))
               ALLOCATE(am(jmax2,kmax2))
               ALLOCATE(bm(jmax2,kmax2))
               ALLOCATE(avgaslice(kmax2))
               ALLOCATE(a0slice(kmax2))
               am(:,:)        = 0.d0
               bm(:,:)        = 0.d0
               a0(:,:)        = 0.d0
               avgaslice(:)   = 0.d0
               a0slice(:)     = 0.d0
               a0tot          = 0.d0
               answer%a       = 0.d0
!$OMP PARALLEL DO DEFAULT(SHARED)&
!$OMP PRIVATE(phi,j,k,l)
               DO j=2,jmax+1
                  DO k=2,kmax+1
                     DO l=1,lmax
                        phi = twopi*dble(l)/dble(lmax)
                        am(j,k) = am(j,k)+rho(j,k,l)*&
                             cos(amcount*phi)
                        bm(j,k) = bm(j,k)+rho(j,k,l)*&
                             sin(amcount*phi)
                        a0(j,k) = a0(j,k)+rho(j,k,l)
                     ENDDO
                  ENDDO
               ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED)&
!$OMP PRIVATE(am,bm,j,k) REDUCTION(+:answer%a,a0tot)
               DO k=2,kmax+1
                  DO j=jstart,jmax+1
                     avgaslice(k) = avgaslice(k)+(sqrt((am(j,k))**2+&
                          (bm(j,k))**2)*(dble(J+1)**2-dble(J)**2))
                     a0slice(k) = a0slice(k)+a0(j,k)
                  ENDDO
                  answer%a = answer%a+avgaslice(k)
                  a0tot = a0tot+a0slice(k)
               ENDDO
!$OMP END PARALLEL DO
               answer%amid = avgaslice(2)*2.d0/a0slice(2)
               answer%a = answer%a*2.d0/a0tot
                  


!                DO k=2,kmax+1
!                   DO j=jstart,jmax+1
!                      DO l=1,lmax
!                         phi = twopi*dble(l)/dble(lmax)        
!                         am(j,k) = am(j,k)+rho(j,k,l)*&
!                              cos(amcount*phi)
!                      bm(j,k) = bm(j,k)+rho(j,k,l)*&
!                           sin(amcount*phi)
!                      a0(j,k) = a0(j,k)+rho(j,k,l)
!                   ENDDO
!                   ambmslice(j,k) = sqrt((am(j,k))**2+
!      &                 (bm(j,k))**2)*(dble(J+1)**2-dble(J)**2)
!                   avgaslice(k) = avgaslice(k)
!      &                 +ambmslice(j,k)
!                   a0tot(m,count) = a0tot(m,count)+a0(j,k)
!      &                 *(dble(J+1)**2-dble(J)**2)
!                ENDDO
!                IF(K==2)THEN
!                   avgamid(m,count) = avgaslice(k)
!                   a0mid(m,count) = a0tot(m,count)
!                ENDIF
!                avga(m,count) = avga(m,count)+avgaslice(k)
!             ENDDO
            DEALLOCATE(a0)
            DEALLOCATE(am)
            DEALLOCATE(bm)
            DEALLOCATE(avgaslice)
            DEALLOCATE(a0slice)

            CALL MPI_SEND(answer,2,MPI_DOUBLE_PRECISION,0,AMCOUNT,&
                 MPI_COMM_WORLD,mpierr)
            GO TO 300
         ENDIF
400      CONTINUE
      ENDDO
   ENDIF

!      print*,'FORTRAN COUNT = ',count

   IF(myrank.eq.0)THEN

      OPEN(UNIT=12,FILE=outfile,FORM='UNFORMATTED')
      WRITE(12) mmax,count
      WRITE(12) avga
      WRITE(12) avgamid
      WRITE(12) timearr
      CLOSE(12)
   ENDIF

   DEALLOCATE(rho)
   CALL MPI_FINALIZE(mpierr)
      
      
 END PROGRAM Ams

      
      
