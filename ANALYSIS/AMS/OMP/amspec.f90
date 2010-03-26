!-*-f90-*-     
      PROGRAM Ams
      IMPLICIT NONE
      
      INTEGER,PARAMETER :: DOUBLE = SELECTED_REAL_KIND(15,300)
      INTEGER :: jmax,kmax,lmax
      INTEGER :: start,finish,skip,size,jstart

      INTEGER :: jmax2,kmax2,mmax
      INTEGER :: count,j,k,l,m,i,numargs
      REAL(DOUBLE),PARAMETER :: tconv = 1605.63
      REAL(DOUBLE),PARAMETER :: pi = 3.14159265358979323846d0
      REAL(DOUBLE),PARAMETER :: twopi = 2.d0*pi
      REAL(DOUBLE),DIMENSION(:,:,:),ALLOCATABLE    :: rho
      REAL(DOUBLE),DIMENSION(:,:),ALLOCATABLE      :: am,bm,ambmslice
      REAL(DOUBLE),DIMENSION(:,:),ALLOCATABLE      :: a0tot,a0mid,a0
      REAL(DOUBLE),DIMENSION(:),ALLOCATABLE        :: avgaslice
      REAL(DOUBLE),DIMENSION(:,:),ALLOCATABLE      :: avga,avgamid
      REAL(DOUBLE),DIMENSION(:),ALLOCATABLE        :: timearr
      REAL(DOUBLE) :: time,phi,dphi
      LOGICAL EXISTSTAT
      CHARACTER outfile*80,indir*80
      CHARACTER rhofile*80,amfile*80,filenum*8,str*80
      CHARACTER jmaxin*10,kmaxin*10,lmaxin*10,startin*10,finishin*10
      CHARACTER jstartin*10,skipin*10

      numargs = IARGC()

      IF (numargs.ne.9) THEN
         print*,"Incorrect number of arguments"
         STOP
      ENDIF
      
      call getarg(1,jmaxin)
      call getarg(2,kmaxin)
      call getarg(3,lmaxin)
      call getarg(4,startin)
      call getarg(5,finishin)
      call getarg(6,skipin)
      call getarg(7,outfile)
      call getarg(8,jstartin)
      call getarg(9,indir)
      read(jmaxin,*)jmax
      read(kmaxin,*)kmax
      read(lmaxin,*)lmax
      read(startin,*)start
      read(finishin,*)finish
      read(skipin,*)skip
      read(jstartin,*)jstart

      size  = ((finish-start)/skip)+1
      mmax  = LMAX/2
      jmax2 = jmax+2
      kmax2 = kmax+2
      dphi  = twopi/lmax

      ALLOCATE(rho(jmax2,kmax2,lmax))
      ALLOCATE(avga(LMAX/2,size))
      ALLOCATE(avgamid(LMAX/2,size))
      ALLOCATE(a0tot(mmax,size))
      ALLOCATE(a0mid(mmax,size))
      ALLOCATE(timearr(size))

      avga(:,:)          = 0.d0
      avgamid(:,:)       = 0.d0
      timearr(:)         = 0.d0
      a0tot(:,:)         = 0.d0
      a0mid(:,:)         = 0.d0
      count              = 0

      
      DO i=start,finish,skip
         WRITE(filenum,'(I6.6)')i
         rhofile = TRIM(indir)//"rho3d."//filenum
         
         INQUIRE(FILE=rhofile,EXIST=EXISTSTAT)
         
         IF(.not.EXISTSTAT) THEN
            print*,"file ",rhofile, "does not exist"
            CYCLE
         ENDIF
         
!         print*," AMS OUT -> OPENING FILE: ", rhofile
         OPEN(UNIT=12,FILE=rhofile,FORM="UNFORMATTED")

         READ(12) rho
         READ(12) time
         CLOSE(12)
         
!         print*,' AMS OUT -> READ FILE: ',rhofile
         count = count+1
         timearr(count) = time/tconv
         

!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE(am,bm,a0,avgaslice,ambmslice,phi,j,k,l)

         DO m=1,MMAX
            ALLOCATE(a0(jmax2,kmax2))
            ALLOCATE(am(jmax2,kmax2))
            ALLOCATE(bm(jmax2,kmax2))
            ALLOCATE(ambmslice(jmax2,kmax2))
            ALLOCATE(avgaslice(kmax2))
            am(:,:)        = 0.d0
            bm(:,:)        = 0.d0
            a0(:,:)        = 0.d0
            avgaslice(:)   = 0.d0
            ambmslice(:,:) = 0.d0
            DO k=2,kmax+1
               DO j=jstart,jmax+1
                  DO l=1,lmax
                     phi = twopi*dble(l)/dble(lmax)        
                     am(j,k) = am(j,k)+rho(j,k,l)*&
                          cos(m*phi)
                     bm(j,k) = bm(j,k)+rho(j,k,l)*& 
                          sin(m*phi)
                     a0(j,k) = a0(j,k)+rho(j,k,l)
                  ENDDO
                  ambmslice(j,k) = sqrt((am(j,k))**2+&
                       (bm(j,k))**2)*(dble(J+1)**2-dble(J)**2)
                  avgaslice(k) = avgaslice(k)&
                       +ambmslice(j,k)
                  a0tot(m,count) = a0tot(m,count)+a0(j,k)&
                       *(dble(J+1)**2-dble(J)**2)
               ENDDO
               IF(K==2)THEN
                  avgamid(m,count) = avgaslice(k)
                  a0mid(m,count) = a0tot(m,count)
               ENDIF
               avga(m,count) = avga(m,count)+avgaslice(k)
            ENDDO
            DEALLOCATE(a0)
            DEALLOCATE(am)
            DEALLOCATE(bm)
            DEALLOCATE(ambmslice)
            DEALLOCATE(avgaslice)
         ENDDO
!$OMP END PARALLEL DO
      ENDDO

      avgamid(:,:) = avgamid(:,:)*2.d0/a0mid(:,:)
      avga(:,:)    = avga(:,:)*2.d0/a0tot(:,:)      
!      avga(:,:)    = LOG10(avga(:,:))
!      avgamid(:,:) = LOG10(avgamid(:,:))
!      print*,'FORTRAN COUNT = ',count

!      print*,avga(2,:)
!      print*,mmax, count

      OPEN(UNIT=12,FILE=TRIM(outfile),FORM='UNFORMATTED')
      WRITE(12) mmax,count
      WRITE(12) avga
      WRITE(12) avgamid
      WRITE(12) timearr
      CLOSE(12)

      DEALLOCATE(rho)
      
    END PROGRAM Ams

      
      
