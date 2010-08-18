!-*-f90-*-
PROGRAM PERIODOGRAM
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)
  integer, parameter :: OFAC = 2, HIFAC = 2
  integer :: jmax, kmax, lmax, istart, iend, iskip, numargs
  integer :: jmax2, kmax2, lmax2,modes,numfiles,fileiter
  integer :: j,k,l,i,m,tstartind,tendind,nout,nsub
  integer :: jmax1,kmax1,jreq,ik
  REAL(DOUBLE), PARAMETER :: torp = 1605.63
  REAL(DOUBLE), PARAMETER :: pi = 3.14159265358979323846d0
  REAL(DOUBLE), PARAMETER :: twopi = 2.d0*pi
  real(double) :: ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,OMMAX,tmassini
  real(double) :: aujreq,tstart,tend,prob,phi,dr,omin,omax,delta,massum
  real(DOUBLE) :: rconv
  real(double),dimension(:,:,:),allocatable :: rho,s,t,a,eps,omega
  real(double),dimension(:,:,:),allocatable :: angle,angsub,results,frequencies
  real(double),dimension(:),allocatable :: omegavg,kappa,timearr,am,bm,tsub
  real(double),dimension(:),allocatable :: X1,Y1,oneang,rhf
  character :: rhodir*80, savedir*80,rhofile*94,savedfile*94,outfile*80
  character :: jmaxin*8,kmaxin*8,lmaxin*8,istartin*8,iendin*8,iskipin*8,modein*8
  character :: filenum*8,aujreqin*10,tstartin*8,tendin*8

  numargs = IARGC()

  if(numargs.ne. 13) then
     print*,"Incorrect number of arguments"
     STOP
  ENDIF
  
  call getarg(1,jmaxin)
  call getarg(2,kmaxin)
  call getarg(3,lmaxin)
  call getarg(4,modein)
  call getarg(5,istartin)
  call getarg(6,iendin)
  call getarg(7,iskipin)
  call getarg(8,tstartin)
  call getarg(9,tendin)
  call getarg(10,aujreqin)
  call getarg(11,outfile)
  call getarg(12,rhodir)
  call getarg(13,savedir)
  read(jmaxin,*)jmax
  read(kmaxin,*)kmax
  read(lmaxin,*)lmax
  read(modein,*)modes
  read(istartin,*)istart
  read(iendin,*)iend
  read(iskipin,*)iskip
  read(aujreqin,*)aujreq
  read(tstartin,*)tstart
  read(tendin,*)tend

  numfiles = ((iend-istart)/iskip)+1
  jmax1 = jmax+1
  kmax1 = kmax+1
  jmax2 = jmax+2
  kmax2 = kmax+2
  lmax2 = lmax/2

  allocate(timearr(numfiles))
  allocate(rho(jmax2,kmax2,lmax))
  allocate(s(jmax2,kmax2,lmax))
  allocate(t(jmax2,kmax2,lmax))
  allocate(a(jmax2,kmax2,lmax))
  allocate(eps(jmax2,kmax2,lmax))
  allocate(omega(jmax,kmax,lmax))
  allocate(omegavg(jmax))
  allocate(kappa(jmax))
  allocate(angle(jmax,modes,numfiles))
  allocate(rhf(jmax))

  write (filenum,'(I6.6)')iend
  savedfile=trim(savedir)//'saved.'//filenum
  OPEN(UNIT=8, FILE=trim(savedfile),FORM='UNFORMATTED',  &
       STATUS="OLD")      
  READ(8) S
  READ(8) T
  READ(8) A
  READ(8) RHO
  READ(8) EPS
  READ(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
  read(8) tmassini
  CLOSE(8)

  dr = ROF3N

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO j = 1,jmax
     rhf(j) = (DBLE(J)-1.d0)*dr+(0.5d0*dr)
  ENDDO
!$OMP END DO
  rconv = aujreq/rhf(JREQ)

!$OMP DO PRIVATE(J,K)
  DO l=1,LMAX
     DO j=1,jmax
        DO k=1,kmax
           omega(j,k,l)=a(j+1,k+1,l)/(rhf(j)**2)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  DO J=1,jmax
     omegavg(j) = 0.d0
     massum = 0.d0
     DO k=1,kmax
        DO l=1,lmax
           omegavg(j) = omegavg(j) + omega(j,k,l)
           massum = massum + rho(j+1,k+1,l)
        ENDDO
     ENDDO
     omegavg(j) = omegavg(j)/massum
  ENDDO

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO J=2,jmax-1
     kappa(j) = sqrt(ABS(rhf(j)*omegavg(j)*(omegavg(j+1)- &
          omegavg(j-1))/dr + 4.d0*omegavg(j)**2))*(torp/twopi)
     rhf(J) = rhf(J)*rconv     
  ENDDO
!$OMP END DO

!$OMP DO
  DO J=2,jmax-1
     omegavg(j) = omegavg(j)*(torp/twopi)
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  omegavg(1) = omegavg(1)*(torp/twopi)
  kappa(1) = omegavg(1)
  rhf(1) = rhf(1)*rconv
  omegavg(jmax) = omegavg(jmax)*(torp/twopi)
  kappa(jmax) = omegavg(jmax)
  rhf(jmax) = rhf(jmax)*rconv

  DO I = 1,NUMFILES
     fileiter = ((I-1)*iskip)+istart
     write (filenum,'(I6.6)')fileiter
     rhofile = trim(rhodir)//'rho3d.'//filenum
     
     print*," PERIOD OUT -> OPENING FILE: ", rhofile
     
     OPEN(UNIT=8, FILE=trim(rhofile),FORM='UNFORMATTED',  &
          STATUS='OLD')
     READ(8) RHO
     READ(8) time
     CLOSE(8)
     
     print*,' PERIOD OUT -> READ FILE: ',rhofile
     print*,' PERIOD OUT -> AT TIME: ',time/torp,' ORPS'
     
     timearr(I)=time/torp

!...   COMPUTE PHASE ANGLES
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(am,bm,phi,j,l)

     DO m = 1,modes
        allocate(am(jmax))
        allocate(bm(jmax))
        am(:)   = 0.d0
        bm(:)   = 0.d0
        DO j = 1,jmax
           DO l=1,lmax
              phi = twopi*dble(l-1)/dble(lmax)
              am(j) = am(j)+rho(j+1,2,l)*cos(dble(m)*phi)
              bm(j) = bm(j)+rho(j+1,2,l)*sin(dble(m)*phi)
           ENDDO

           IF (am(J) > 0 ) THEN
              angle(J,M,I) = atan(bm(J)/am(J))*180./pi
           ELSE IF (am(J) < 0 ) THEN
              angle(J,M,I) = (atan(bm(J)/am(J)) + pi)*180./pi
           ELSE 
              IF (bm(J) > 0 ) THEN
                 angle(J,M,I) = 90.d0
              ELSE IF (bm(J) < 0 ) THEN
                 angle(J,M,I) = -90.d0
              else
                 angle(J,M,I) = 720.
              endif
           endif

!           angle(J,M,I) = atan2(bm(J),am(J))+pi
        ENDDO
        deallocate(am)
        deallocate(bm)
     ENDDO
!$OMP END PARALLEL DO
  ENDDO

!...    Extract time interval

  I = 0
  DO WHILE ((timearr(I).lt.tstart).and.(I.lt.numfiles))
     I = I+1
  ENDDO
  
  tstartind = I
  IF(tstartind.eq.numfiles) THEN
     print*,"Beginning of interval is not in this data set"
     STOP
  ENDIF

  DO WHILE ((timearr(I).lt.tend).and.(I.lt.numfiles))
     I = I+1
  ENDDO

  tendind = I
  IF(tendind.eq.numfiles) THEN
     print*,"!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!"
     print*,"!!!  The end time specified in the   !!!"
     print*,"!!!  interval is beyond the data set !!!"
     print*,"!!!  This may adversely affect the   !!!"
     print*,"!!!  periodogram results             !!!"
     print*,"!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!"
  ENDIF

  nsub = tendind - tstartind + 1

  print*,"NUMFILES = ",numfiles
  print*,"NSUB = ",nsub
  print*,"TSTARTIND = ",tstartind
  print*,"TENDIND = ",tendind
  print*,"modes = ",modes
  print*,"jmax2 = ",jmax2
  
  allocate(tsub(nsub))
  allocate(angsub(jmax,modes,nsub))

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE(J,I)
  DO M=1,modes
     DO J=1,JMAX
        DO I=1,nsub
           angsub(J,M,I)=angle(J,M,I+tstartind-1)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO 
  DO I = 1,nsub
     tsub(I)=timearr(I+tstartind-1)*torp
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  nout = 0.5*OFAC*HIFAC*nsub

  print*,"about to allocate"


  allocate(X1(nout))
  allocate(results(JMAX,modes,2*nsub))
  allocate(frequencies(JMAX,modes,2*nsub))

  print*,"about to call periodogram"
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(oneang,J,I,Y1,omin,delta,omax,ik)
  DO M=1,modes
     allocate(oneang(nsub))
     allocate(Y1(2*nsub))
     DO J=1,JMAX
        DO I = 1,nsub
           oneang(I) = cos(angsub(J,M,I)*pi/180.)
        ENDDO
        DO I=1,2*nsub+1
           Y1(I) = 0.d0
        ENDDO
        
        omin = 0.d0
        delta = 0.d0

        call PERIOD_ORIG(oneang,tsub,nsub,omin,omax,delta,torp,M,ik,Y1)


        omin = omin*(torp/(twopi*M))
        omax = omax*(torp/(twopi*M))
        delta = delta*(torp/(twopi*M))
        
        DO I = 1,ik
           results(J,M,I) = Y1(I)
           frequencies(J,M,I) = omin + (delta)*(I-1)
!           frequencies(J,M,I) = X1(I)*(torp/(twopi*M))
        ENDDO
        
     ENDDO
     deallocate(oneang)
     deallocate(Y1)
  ENDDO
!$OMP END PARALLEL DO

  OPEN(UNIT=12,FILE=TRIM(outfile),FORM='UNFORMATTED')
  WRITE(12) JMAX,modes,2*nsub,tstart,tend,aujreq
  WRITE(12) results
  WRITE(12) frequencies
  WRITE(12) rhf
  WRITE(12) omegavg
  WRITE(12) kappa
  CLOSE(12)

           
        

  
!...    Subroutine uses periodogram technique for generating
!...    power spectra (see references below).  The Power Spectra are
!...    plotted in the frequency-radius.  The spectra can be added into
!...    one curve for OMPAT extraction.  Note that OMPAT=2pi/PPAT.

END PROGRAM PERIODOGRAM

SUBROUTINE PERIOD_ORIG(x,ti,nmax,omin,omax,delta,TORP,mode,ik,power)

  IMPLICIT NONE
  INTEGER, parameter :: double = selected_real_kind(15,30)


  real(double),dimension(nmax)          :: x,ti
  real(double),dimension(2*nmax+1)      :: power
  REAL(double)          :: SMEGA,OMEGA,T,TMIN,omin,omax,delta,torp
  real(double)          :: TWOPI,FNI,SUM,AVE,PMAX,PMIN,TOP,BOT,THETA1
  real(double)          :: THETA2,XC,THETA,TAU,A,B,AA,BB,C,S,VAR,pi
  INTEGER               :: N0,NMAX,mode,ARRAYSIZE,IK,M,L,j,K


  TWOPI = 2.d0*ACOS(-1.d0)
  PI = ACOS(-1.d0)
  
  SUM = 0.d0
  DO M=1,nmax
     SUM = SUM+X(M)
  ENDDO
  AVE=SUM/nmax

  SUM = 0.d0
  DO M=1,nmax
     SUM = SUM+(X(M)-AVE)**2
  ENDDO
  VAR = SUM/nmax

  DO M=1,nmax
     X(M)=X(M)-AVE
  ENDDO

  DO M=2,nmax
     TI(M)=TI(M)-TI(1)
  ENDDO
  TI(1)=0.d0

  T = TI(nmax)

  OMIN = TWOPI/T
  OMAX = 4.d0*mode*TWOPI/torp
  IK = 2.d0*nmax

  DELTA = (omax-omin)/IK

  OMEGA = OMIN

  DO J=1,IK+1

     SMEGA=OMEGA-0.0001*DELTA
     CALL TTT(TI,SMEGA,nmax,TOP,BOT)
     THETA1=ATAN(TOP/BOT)
     CALL TTT(TI,OMEGA,nmax,TOP,BOT)
     THETA2=ATAN(TOP/BOT)

     XC=TOP/BOT
     IF((XC.GE.0.d0).AND.(THETA2.GT.THETA1)) THETA=THETA2
     IF((XC.GE.0.d0).AND.(THETA2.LT.THETA1)) THETA=PI-THETA2
     IF((XC.LE.0.d0).AND.(THETA2.LT.THETA1)) THETA=PI+THETA2
     IF((XC.LE.0.d0).AND.(THETA2.GT.THETA1)) THETA=TWOPI-THETA2

     TAU=THETA/(2.d0*OMEGA)

     A=0.d0
     B=0.d0
     AA=0.d0
     BB=0.d0
     DO K=1,nmax
        A=A+(COS(OMEGA*(TI(K)-TAU)))**2
        B=B+(SIN(OMEGA*(TI(K)-TAU)))**2
        AA=AA+X(K)*COS(OMEGA*(TI(K)-TAU))
        BB=BB+X(K)*SIN(OMEGA*(TI(K)-TAU))
     ENDDO

     C=AA*1.0/SQRT(A)
     S=BB*1.0/SQRT(B)
     POWER(j)=0.5*(C**2+S**2)/VAR
     IF(Mode.eq.2)THEN
     ENDIF
     OMEGA=OMEGA+DELTA
  ENDDO

  RETURN
END SUBROUTINE PERIOD_ORIG

SUBROUTINE TTT(TI,OMEGA,N0,TOP,BOT)

!...Companion to Periodogram
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)

  REAL(double),dimension(N0)          :: TI
  REAL(DOUBLE)                        :: OMEGA,E,F,TOP,BOT
  INTEGER                             :: K,N0

  TOP=0.d0
  BOT=0.d0
  DO K=1,N0
     E=SIN(2.d0*OMEGA*TI(K))
     F=COS(2.d0*OMEGA*TI(K))
     TOP=TOP+E
     BOT=BOT+F
  ENDDO
  RETURN
END SUBROUTINE TTT


SUBROUTINE PERIOD(X,Y,N,OFAC,HIFAC,PX,PY,NP,PROB)
  IMPLICIT NONE

  integer, parameter :: double = selected_real_kind(15,300)
  INTEGER :: i,j,N,OFAC,HIFAC,NP
  REAL(DOUBLE),parameter :: pi = 3.14159265358979323846d0
  REAL(DOUBLE),parameter :: twopi = 2.d0*pi
  REAL(DOUBLE) :: ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss
  REAL(DOUBLE) :: sumc,sumcy,sums,sumsh,sumsy,swtau,var,wtau
  REAL(DOUBLE) :: xave,xdif,xmax,xmin,yy,arg,wtemp
  REAL(DOUBLE) :: X(*),Y(*),PX(*),PY(*)
  REAL(DOUBLE) :: PROB
  REAL(DOUBLE),dimension(:),allocatable :: wi,wpi,wpr,wr

  allocate(wi(n))
  allocate(wpi(n))
  allocate(wpr(n))
  allocate(wr(n))

  CALL AVEVAR(Y,N,AVE,VAR)

  DO J=2,N
     X(J) = X(J) - X(1)
  ENDDO
  
  X(1) = 0.d0

  XMIN=X(1)
  XMAX=XMIN
  DO J=1,N
     IF(X(J).lt.XMIN)XMIN=X(J)
     IF(X(J).gt.XMAX)XMAX=X(J)
  ENDDO
  XDIF=XMAX-XMIN
  XAVE=0.5*(XMAX+XMIN)
  pymax = 0.d0
  pnow = 1.d0/(XDIF*OFAC)

  DO J=1,N
     arg = twopi*((X(J)-xave)*pnow)
     wpr(J) = -2.d0*(sin(0.5d0*arg))**2
     wpi(J) = sin(arg)
     wr(J)  = cos(arg)
     wi(J)  = wpi(J)
  ENDDO

  DO I=1,NP
     px(i) = pnow
     sumsh = 0.d0
     sumc  = 0.d0
     DO J=1,N
        c = wr(J)
        s = wi(J)
        sumsh = sumsh + s*c
        sumc = sumc + (c-s)*(c+s)
     ENDDO
     
     wtau  = 0.5d0*atan2(2.d0*sumsh,sumc)
     swtau = sin(wtau)
     cwtau = cos(wtau)
     sums  = 0.d0
     sumc  = 0.d0
     sumsy = 0.d0
     sumcy = 0.d0
     DO J=1,N
        s = wi(J)
        c = wr(J)
        ss= s*cwtau - c*swtau
        cc= c*cwtau + s*swtau
        sums = sums + ss**2
        sumc = sumc + cc**2
        yy = y(J) - ave
        sumsy = sumsy + yy*ss
        sumcy = sumcy + yy*cc
        wtemp = wr(J)
        wr(J) = (wtemp*wpr(J)-wi(J)*wpi(J))+wr(J)
        wi(J) = (wi(J)*wpr(J)+wtemp*wpi(J))+wi(J)
     ENDDO
     py(I) = 0.5d0*(sumcy**2/sumc+sumsy**2/sums)/var
     IF (py(I).ge.pymax) pymax = py(I)
     pnow = pnow + 1.d0/(ofac*xdif)
  ENDDO
  
  expy = exp(-pymax)
  effm = 2.d0*np/ofac
  prob = effm*expy
  IF (prob.lt.0.01d0) prob = 1.d0-(1.d0-expy)**effm
  
  deallocate(wr)
  deallocate(wpr)
  deallocate(wpi)
  deallocate(wi)
  
  RETURN
END SUBROUTINE PERIOD



SUBROUTINE FASPER(X,Y,N,OFAC,HIFAC,WK1R,WK2R,NOUT,PROB)
  IMPLICIT NONE

  integer, parameter :: double = selected_real_kind(15,300)
  INTEGER, PARAMETER :: MACC=4
  INTEGER :: NOUT,NFREQT,NFREQ,NDIM,FNDIM,J,K,N
  INTEGER :: OFAC,HIFAC
  REAL(DOUBLE) :: X(*),Y(*),WK1R(*),WK2R(*)
  REAL(DOUBLE),DIMENSION(:),allocatable :: wk1,wk2
  REAL(DOUBLE) :: AVE,VAR,XMIN,XMAX,FAC,CK,CKK,DF,PMAX
  REAL(DOUBLE) :: HYPO,HC2WT,HS2WT,CWT,SWT,DEN,CTERM
  REAL(DOUBLE) :: STERM,EXPY,EFFM,PROB,XDIF

  NFREQT = OFAC*HIFAC*N*MACC
  NFREQ = 64
  DO WHILE (NFREQ.lt.NFREQT)
     NFREQ = NFREQ*2
  ENDDO
  NDIM = NFREQ*2
  allocate(wk1(NDIM))
  allocate(wk2(NDIM))
  CALL AVEVAR(Y,N,AVE,VAR)
! Do some data rescaling. Subtract the mean from the angles.
! Also subtract the start time from the time array, making
! this array go from 0.0 to interval.
  
  DO J=2,N
     X(J) = X(J) - X(1)
!     Y(J) = Y(J) - AVE
  ENDDO
  
  X(1) = 0.d0
!  Y(1) = Y(1) - AVE

  XMIN=X(1)
  XMAX=XMIN
  DO J=1,N
     IF(X(J).lt.XMIN)XMIN=X(J)
     IF(X(J).gt.XMAX)XMAX=X(J)
  ENDDO
  XDIF = XMAX - XMIN
  wk1(:) = 0.d0
  wk2(:) = 0.d0
  FAC = NDIM/(XDIF*OFAC)
  FNDIM = NDIM
  DO J=1,N
     CK = 1.d0+DMOD((X(J)-XMIN)*FAC,dble(FNDIM))
     CKK = 1.d0+DMOD(2.d0*(CK-1.d0),dble(FNDIM))
     CALL SPREAD(Y(J)-AVE,wk1,NDIM,CK,MACC)
     CALL SPREAD(1.d0,wk2,NDIM,CKK,MACC)
  ENDDO
  CALL REALFT(wk1,NFREQ)
  CALL REALFT(wk2,NFREQ)
  DF = 1.d0/(XDIF*OFAC)
  K = 3
  PMAX = -1.d0

  DO J=1,NOUT
     HYPO = sqrt(wk2(K)**2+wk2(k+1)**2)
     HC2WT = 0.5d0*wk2(k)/HYPO
     HS2WT = 0.5d0*wk2(k+1)/HYPO
     CWT = sqrt(0.5d0+HC2WT)
     SWT = sign(sqrt(0.5d0-HC2WT),HS2WT)
     DEN = 0.5d0*N+HC2WT*wk2(k)+HS2WT*wk2(k+1)
     CTERM = (CWT*wk1(k) + SWT*wk1(k+1))**2/DEN
     STERM = (CWT*wk1(K+1)-SWT*wk1(k))**2/(N-DEN)
     wk1(J) = J*DF
     wk2(J) = (CTERM+STERM)/(2.d0*VAR)
     IF (wk2(J).gt.PMAX) PMAX = wk2(J)
     K = K+2
  ENDDO
  EXPY = exp(-PMAX)
  EFFM = 2.d0*NOUT/OFAC
  PROB = EFFM*EXPY
  IF (PROB.gt.0.01d0)PROB = 1.d0-(1.d0-EXPY)**EFFM
  DO J = 1,NOUT
     WK1R(J) = wk1(J)
     WK2R(J) = wk2(J)
  ENDDO
  RETURN
END SUBROUTINE FASPER
  

  

!...    Following routines from Numerical Recipes.
SUBROUTINE FOUR1(data, nn)
  IMPLICIT NONE  

  integer, parameter :: double = selected_real_kind(15,300)
  INTEGER :: n,mmax,m,j,istep,i,nn
  REAL(DOUBLE) :: wtemp,wr,wpr,wpi,wi,theta,tempr,tempi
  REAL(DOUBLE) :: data(*)
  
  n = ISHFT(nn,1)
  j = 1
  DO I=1,N-1,2
     IF(j.gt.i) THEN
        call SWAP(data(j),data(i))
        call SWAP(data(j+1),data(i+1))
     ENDIF
     m = ISHFT(n,-1)
     DO WHILE ((m.ge.2).and.(j.gt.m))
        j = j - m
        m = ISHFT(m,-1)
     ENDDO
     j = j + m
  ENDDO

  mmax = 2
  DO WHILE (n.gt.mmax)
     istep = ISHFT(mmax,1)
     theta = (6.28318530717959/mmax)
     wtemp = sin(0.5d0*theta)
     wpr = -2.d0*wtemp**2
     wpi = sin(theta)
     wr = 1.d0
     wi = 0.d0
     DO M=1,MMAX-1,2
        DO I=M,N,istep
           j = i + mmax
           tempr = wr*data(j)-wi*data(j+1)
           tempi = wr*data(j+1)+wi*data(j)
           data(j) = data(i)-tempr
           data(j+1) = data(i+1)-tempi
           data(i) = data(i)+tempr
           data(i+1) = data(i+1)+tempi
        ENDDO
        wtemp = wr
        wr = wtemp*wpr-wi*wpi+wr
        wi = wi*wpr+wtemp*wpi+wi
     ENDDO
     mmax=istep
  ENDDO
  RETURN
END SUBROUTINE FOUR1

SUBROUTINE REALFT(dataft,n)
  IMPLICIT NONE

  integer, parameter :: double = selected_real_kind(15,300)
  INTEGER :: i, i1,i2,i3,i4,np3,dum,n
  REAL :: c1=0.5,c2,h1r,h1i,h2r,h2i
  REAL(DOUBLE) :: wr,wi,wpr,wpi,wtemp,theta
  REAL(DOUBLE) :: dataft(*)

  dum = ISHFT(n,-1)
  theta = 3.141592653589793/dble(dum)
  c2 = 0.5
  call four1(dataft,dum)
  wtemp = sin(0.5d0*theta)
  wpr = -2.d0*wtemp**2
  wpi = sin(theta)
  wr = 1.d0+wpr
  wi = wpi
  np3 = n+3
  dum = ISHFT(n,-2)
  DO i=2,dum
     i1 = i+i+1
     i2 = 1+i1
     i3 = np3-i2
     i4 = 1+i3
     h1r = c1*(dataft(i1)+dataft(i3))
     h1i = c1*(dataft(i2)-dataft(i4))
     h2r = -c2*(dataft(i2)+dataft(i4))
     h2i = c2*(dataft(i1)-dataft(i3))
     dataft(i1) = h1r+wr*h2r-wi*h2i
     dataft(i2) = h1i+wr*h2i+wi*h2r
     dataft(i3) = h1r-wr*h2r+wi*h2i
     dataft(i1) = -h1i+wr*h2i+wi*h2r
     wtemp = wr
     wr = wtemp*wpr-wi*wpi+wr
     wi = wi*wpr+wtemp*wpi+wi
  ENDDO
  h1r = dataft(1)
  dataft(1) = h1r + dataft(2)
  dataft(2) = h1r - dataft(2)
  RETURN
END SUBROUTINE REALFT
  
SUBROUTINE SPREAD(Y,YY,N,X,M)
  IMPLICIT NONE

  integer, parameter :: double = selected_real_kind(15,300)
  REAL(DOUBLE) :: Y,X,fac,nden
  REAL(DOUBLE) :: YY(*)
  INTEGER :: n,m,ix,ilo,ihi,j
  INTEGER,dimension(10),parameter :: nfac=(/1,1,2,6,24,120,720,5040,40320,362880/)
  
  ix = int(x)
  if (x.eq.dble(ix)) THEN
     yy(ix) = yy(ix) + y
  ELSE
     ilo = MIN(MAX(int(x-0.5*m+1.0),1),n-m+1)
     ihi = ilo+m-1
     nden = dble(nfac(m))
     fac = x-ilo
     DO J=ilo+1,IHI
        fac = fac*(x-j)
     ENDDO
     yy(ihi) = yy(ihi)+y*fac/(nden*(x-ihi))
     DO J=ihi-1,ilo,-1
        nden=(nden/(j+1-ilo))*(j-ihi)
        yy(j) = yy(j) + y*fac/(nden*(x-j))
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE SPREAD
     
SUBROUTINE SWAP(A,B)
  IMPLICIT NONE

  integer, parameter :: double = selected_real_kind(15,300)
  REAL(DOUBLE) :: A,B,TMP
  TMP = A
  a = B
  B = TMP
  RETURN
END SUBROUTINE SWAP

SUBROUTINE AVEVAR(X,NUM,AVG,VAR)
  IMPLICIT NONE
  
  integer, parameter :: double = selected_real_kind(15,300)
  INTEGER :: NUM, J
  REAL(DOUBLE) :: X(*)
  REAL(DOUBLE) :: AVG,VAR,S,EP

  AVG = 0.d0
  DO J=1,NUM
     AVG = AVG + X(J)
  ENDDO
  AVG = AVG/NUM
  VAR = 0.d0
  EP  = 0.d0
  DO J=1,NUM
     S = X(J)-AVG
     EP = EP + S
     VAR = VAR + S**2
  ENDDO
  VAR = (VAR-EP**2/NUM)/(NUM-1)
  RETURN
END SUBROUTINE AVEVAR
