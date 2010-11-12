PROGRAM  QPLOT_EOS

  IMPLICIT NONE

  INCLUDE 'hydroparam.h'
  INCLUDE 'globals.h'
  INCLUDE 'mpif.h'

  integer, parameter :: double = selected_real_kind(15,300)
  REAL(DOUBLE), PARAMETER :: torp = 1605.63
  REAL(DOUBLE), PARAMETER :: pi = 3.14159265358979323846d0
  REAL(DOUBLE), PARAMETER :: twopi = 2.d0*pi
  INTEGER I,NUMFILES,COUNTER,J,K,L,JREQ,IERR
  REAL(DOUBLE) :: time_begin,time_end,engtmp,dummy,limit
  REAL(DOUBLE) :: dr,dz,elost,sound,ommax,time0,totengtmp,cooltmp
  REAL(DOUBLE), DIMENSION(:,:),ALLOCATABLE :: qomega,qkappa,colcool
  REAL(DOUBLE), DIMENSION(:),ALLOCATABLE :: timearr,omegavg,csavg,sigmavg
  REAL(DOUBLE), DIMENSION(:),ALLOCATABLE :: kappa,toteng,vol
  CHARACTER :: savedfile*80,filenum*8,qfile*80,coolfile*80
  LOGICAL FILE_EXIST

  CALL MPI_INIT(IERR)

  time_begin = MPI_Wtime()

  NUMFILES = ((IEND - ISTART)/ISKIP) + 1
  
  ALLOCATE(qomega(JMAX2,NUMFILES))
  ALLOCATE(qkappa(JMAX2,NUMFILES))
  ALLOCATE(colcool(JMAX2,NUMFILES))
  ALLOCATE(timearr(NUMFILES))
  ALLOCATE(toteng(NUMFILES))
  ALLOCATE(omegavg(JMAX2))
  ALLOCATE(csavg(JMAX2))
  ALLOCATE(sigmavg(JMAX2))
  ALLOCATE(kappa(JMAX2))
  ALLOCATE(vol(JMAX2))

  COUNTER = 0

  DO I=ISTART,IEND,ISKIP
     WRITE(filenum,'(I8.8)')I
     savedfile=trim(datadir)//'saved.'//filenum
     coolfile=trim(datadir)//'coolheat_full.'//filenum
     INQUIRE(file=savedfile,exist=FILE_EXIST)
     IF(FILE_EXIST) THEN
        INQUIRE(file=coolfile,exist=FILE_EXIST)
     ELSE
        CYCLE
     ENDIF
        
     IF(FILE_EXIST) THEN
        OPEN(UNIT=8, FILE=trim(savedfile),FORM='UNFORMATTED',    &
        STATUS="OLD")
      
        READ(8) S
        READ(8) T
        READ(8) A
        READ(8) RHO
        READ(8) EPS
        READ(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
        read(8) tmassini
        CLOSE(8)

        OPEN(UNIT=9,FILE=trim(coolfile),FORM='UNFORMATTED',      &
             STATUS='OLD')

        READ(9) divflux
        READ(9) lambda
        READ(9) hgamma 
        READ(9) igamma
        READ(9) dummy
        READ(9) TempK
        READ(9) TeffK
        READ(9) TphK
        READ(9) time0
        CLOSE(9)
        
        IF (ABS((time/time0)-1).gt.0.01) STOP

        print*, 'Opened file: ',trim(savedfile),' and ',trim(coolfile)
        print*, '   Your time is ', time/torp
     ELSE
        CYCLE
     ENDIF
     dr = ROF3N
     dz = ZOF3N
     limit = phylim*den
     
     IF ((time < time1*torp).or.(time > time2*torp)) CYCLE 

     IF (counter == 0) then
        DO j = 1,jmax2 
           r(j) = (FLOAT(j-2)*dr)
           rhf(J) = r(J)+0.5*dr
        ENDDO
        CALL INIT()
        CALL INITENGTABLE()
     ENDIF

     COUNTER = COUNTER + 1
     timearr(COUNTER) = time/torp

     CALL TEMPFIND()

     DO J=2,JMAX1
        vol(J) = pi*(r(J+1)**2-r(J)**2)*dz
     ENDDO
!$OMP PARALLEL DEFAULT(SHARED)

     totengtmp = 0.d0
!$OMP DO REDUCTION(+:totengtmp) FIRSTPRIVATE(vol)
     DO L=1,LMAX
        DO K=2,KMAX1
           DO J=2,JMAX1
              totengtmp = totengtmp + eps(J,K,L)*vol(J)
           ENDDO
        ENDDO
     ENDDO
!$OMP END DO
     
     engtmp  = 0.d0
     cooltmp = 0.d0
     DO J=2,JMAX1
        engtmp  = 0.d0
        cooltmp = 0.d0
        DO L=1,LMAX
           DO K=2,KMAX1
              IF(rho(J,K,L).ge.limit)THEN
                 engtmp = engtmp + eps(J,K,L)
                 cooltmp = cooltmp + divflux(J,K,L)
              ENDIF
           ENDDO
        ENDDO
        colcool(J,COUNTER) = engtmp/cooltmp
     ENDDO

!$OMP DO SCHEDULE(STATIC)
         DO J=1,JMAX2
            DO L=1,LMAX
                  p(J,2,L) = bkmpcgs*tempk(J,2,L)*rhoconv*rho(J,2,L)   &
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
           omegavg(J) = omegavg(J)+ a(j,2,l)/rho(j,2,l)/rhf(j)**2
           csavg(J) = csavg(J) + sqrt(gamma1(J,2,L)*p(J,2,L)/rho(J,2,L))
        ENDDO
        sigmavg(J) = 2.d0*sigmavg(J)/DBLE(LMAX)
        omegavg(J)  = omegavg(J)/DBLE(LMAX)
        csavg(J)    = csavg(J)/DBLE(LMAX)
     ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC) FIRSTPRIVATE(omegavg)
     DO J=3,JMAX-1
        kappa(j)=sqrt(ABS(rhf(j)*omegavg(j)*(omegavg(j+1)-   &
        omegavg(j-1))/dr + 4.d0*omegavg(j)**2))
        qkappa(j,counter)=kappa(j)*csavg(j)/(pi*sigmavg(j))
        qomega(j,counter)=omegavg(j)*csavg(j)/(pi*sigmavg(j))
     ENDDO
!$OMP END DO

!$OMP END PARALLEL
     toteng(COUNTER) = totengtmp

  ENDDO

!$OMP PARALLEL DO

  DO J=1,JMAX2
     r(J) = rhf(j)*rdiskau/rhf(JREQ)
  ENDDO
!$OMP END PARALLEL DO

  time_end = MPI_Wtime()
  
  CALL MPI_Finalize(IERR)

  print*,'Execution time ',time_end-time_begin,' seconds.'
  
  write (filenum,'(I8.8)')IEND
  qfile=trim(outdir)//trim(outfile)//filenum

  OPEN(UNIT=15,FILE=qfile,FORM='UNFORMATTED')
  WRITE(15)JMAX2,NUMFILES,COUNTER
  WRITE(15)ISTART,IEND,ISKIP
  WRITE(15)timearr
  WRITE(15)toteng
  WRITE(15)r
  WRITE(15)qomega
  WRITE(15)qkappa
  WRITE(15)colcool
  CLOSE(15)

  DEALLOCATE(timearr)
  DEALLOCATE(qomega)
  DEALLOCATE(qkappa)
  DEALLOCATE(colcool)
  DEALLOCATE(toteng)
  DEALLOCATE(omegavg)
  DEALLOCATE(csavg)
  DEALLOCATE(sigmavg)
  DEALLOCATE(kappa)
  DEALLOCATE(vol)
END PROGRAM QPLOT_EOS
