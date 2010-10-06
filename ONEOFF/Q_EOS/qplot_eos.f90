PROGRAM  QPLOT_EOS

  IMPLICIT NONE
      
  INCLUDE 'hydroparam.h'
  INCLUDE 'globals.h'
  
  integer, parameter :: double = selected_real_kind(15,300)
  REAL(DOUBLE), PARAMETER :: torp = 1605.63
  REAL(DOUBLE), PARAMETER :: pi = 3.14159265358979323846d0
  REAL(DOUBLE), PARAMETER :: twopi = 2.d0*pi
  INTEGER I,NUMFILES,COUNTER
  REAL(DOUBLE) :: time1,time2,dr
  REAL(DOUBLE), DIMENSION(:,:),ALLOCATABLE :: qomega,qkappa
  REAL(DOUBLE), DIMENSION(:),ALLOCATABLE :: timearr
  CHARACTER :: savedfile*80,filenum*8
  LOGICAL FILE_EXIST

  time1 = 13.0
  time2 = 21.0

  NUMFILES = ((IEND - ISTART)/ISKIP) + 1
  
  ALLOCATE(qomega(JMAX2,NUMFILES))
  ALLOCATE(qkappa(JMAX2,NUMFILES))
  ALLOCATE(timearr(NUMFILES))

  CALL INIT()
  CALL INITENGTABLE()

  COUNTER = 0

  DO I=ISTART,IEND,ISKIP
     WRITE(filenum,'(I8.8)')I
     savedfile=trim(datadir)//'saved.'//filenum
     INQUIRE(file=savedfile,exist=FILE_EXIST)
     IF(FILE_EXIST) THEN
        OPEN(UNIT=8, FILE=trim(savedfile),FORM='UNFORMATTED',  $
        STATUS="OLD")
      
        READ(8) S
        READ(8) T
        READ(8) A
        READ(8) RHO
        READ(8) EPS
        READ(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
        read(8) tmassini
        CLOSE(8)
               
        print*, 'Opened file: ',trim(savedfile)
        print*, '   Your time is ', time/torp
     ELSE
        CYCLE
     ENDIF
     dr = ROF3N
     dz = ZOF3N
     
     IF ((time < time1*torp).or.(time > time2*torp)) CYCLE 

     IF (counter == 0) then
        DO j = 1,jmax2 
           r(j) = (FLOAT(j-1)*dr+0.5*dr)*rdiskau/(FLOAT(JREQ-1)*dr  $
           +0.5*dr)
           rhf(J) = r(J)*(FLOAT(JREQ-1)*dr+0.5*dr)/rdiskau
        ENDDO
     ENDIF

     COUNTER = COUNTER + 1
     timearr(COUNTER) = time/torp

     CALL TEMPFIND()
!$OMP PARALLEL DEFAULT(SHARED)

!$OMP DO SCHEDULE(STATIC)
         DO J=1,JMAX2
            DO L=1,LMAX
                  p(J,2,L) = bkmpcgs*tempk(J,2,L)*rhoconv*rho(J,2,L)   $
                  / (muc*pconv)
            ENDDO
         ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
     DO J=2,JMAX1
        omegavg(J)  = 0.d0
        csavg(J)    = 0.d0
        sigmavg(J)  = 0.d0
        DO L=1,LMAX
           DO K=2,KMAX1
              sigmavg(J) = sigmavg(J) + rho(J,K,L)*dz
           ENDDO
           omegavg(J) = omegavg(J)+ a(j,k,l)/rho(j,k,l)/rhf(j)**2
           csavg(J) = csavg(J) + sqrt(gamma1(J,2,L)*p(J,2,L)/rho(J,2,L))
        ENDDO
        sigmaavg(J) = 2.d0*sigmavg(J)/DBLE(LMAX)
        omegavg(J)  = omegavg(J)/DBLE(LMAX)
        csavg(J)    = csavg(J)/DBLE(LMAX)
     ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
     DO J=3,JMAX-1
        kappa(j)=sqrt(ABS(rhf(j)*omegavg(j)*(omegavg(j+1)-   $
        omegavg(j-1))/dr + 4.d0*omegavg(j)**2))
        qkappa(j,counter)=kappa(j)*csavg(j)/(pi*sigmavg(j))
        qomega(j,counter)=omegavg(j)*csavg(j)/(pi*sigmavg(j))
     ENDDO
!$OMP END DO

!$OMP END PARALLEL

     OPEN(UNIT=15,FILE=qfile,FORM='UNFORMATTED')
     WRITE(15)JMAX2,NUMFILES,COUNTER
     WRITE(15)ISTART,IEND,ISKIP
     WRITE(15)timearr
     WRITE(15)r
     WRITE(15)qomega
     WRITE(15)qkappa
     CLOSE(15)

     DEALLOCATE(timearr)
     DEALLOCATE(qomega)
     DEALLOCATE(qkappa)
