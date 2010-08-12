!-*-f90-*-
PROGRAM PERIODOGRAM
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)
  integer, parameter :: OFAC = 4, HIFAC = 2
  integer :: jmax, kmax, lmax, jreq, istart, iend, iskip, numargs
  integer :: jmax2, kmax2, lmax2,modes,numfiles,fileiter
  integer :: j,k,l,i,m,tstartind,tendind,nout,nsub
  integer :: jmax1,kmax1
  REAL(DOUBLE), PARAMETER :: torp = 1605.63
  REAL(DOUBLE), PARAMETER :: pi = 3.14159265358979323846d0
  REAL(DOUBLE), PARAMETER :: twopi = 2.d0*pi
  real(double) :: ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,OMMAX,tmassini
  real(double) :: dr,aujreq,tstart,tend,prob,phi
  real(double),dimension(:,:,:),allocatable :: rho,s,t,a,eps,omega
  real(double),dimension(:,:,:),allocatable :: angle,angsub,results,frequencies
  real(double),dimension(:),allocatable :: omegavg,kappa,massum,timearr,am,bm,tsub
  real(double),dimension(:),allocatable :: X1,Y1,oneang,r,rhf
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
  allocate(omega(jmax2,kmax2,lmax))
  allocate(omegavg(jmax))
  allocate(kappa(jmax))
  allocate(massum(jmax2))
  allocate(angle(jmax2,modes,numfiles))
  allocate(r(jmax2))
  allocate(rhf(jmax2))

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

!$OMP PARALLEL DEFAULT(SHARED)

!$OMP DO
  DO j = 1,jmax2
     r(j)   = (DBLE(J)-2.d0)*dr
     rhf(j) = r(j)+(0.5d0*dr)
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=1,jmax2
     DO k=1,kmax2
        DO l=1,LMAX
           omega(j,k,l)=a(j,k,l)/(rhf(j)**2)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO REDUCTION(+:omegavg,massum)
  DO J=2,jmax1
     omegavg(j-1) = 0.d0
     massum(j) = 0.d0
     DO k=1,kmax2
        DO l=1,lmax
           omegavg(j-1) = omegavg(j-1) + omega(j,k,l)
           massum(j) = massum(j) + rho(j,k,l)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
     
!$OMP DO
  DO J=1,jmax
     omegavg(j) = omegavg(j)/massum(j)
     kappa(j) = sqrt(ABS(rhf(j)*omegavg(j)*(omegavg(j+1)- &
          omegavg(j-1))/dr + 4.d0*omegavg(j)**2))*(torp/(2.d0*pi))
     omegavg(j) = omegavg(j)*(torp/(2.d0*pi))
!!!!!!!!!!!!!!!!!!!!! If not needed put change rhf to AU here.........
     
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

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
!OMP PARALLEL DO DEFAULT(SHARED) &
!OMP PRIVATE(am,bm,phi,j,l)

     DO m = 1,modes
        allocate(am(jmax2))
        allocate(bm(jmax2))
        am(:)   = 0.d0
        bm(:)   = 0.d0
        DO j = 2,jmax+1
           DO l=1,lmax
              phi = twopi*dble(l)/dble(lmax)
              am(j) = am(j)+rho(j,2,l)*cos(dble(m)*phi)
              bm(j) = bm(j)+rho(j,2,l)*sin(dble(m)*phi)
           ENDDO
           angle(J,M,I) = atan2(am(J),bm(J))+pi
        ENDDO
        deallocate(am)
        deallocate(bm)
     ENDDO
!OMP END PARALLEL DO
  ENDDO

!...    Extract time interval

  I = 1
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
  nsub = tendind - tstartind
  
  allocate(tsub(nsub))
  allocate(angsub(jmax2,modes,nsub))

  DO I = tstartind,tendind
     tsub(I+1-tstartind)=timearr(I)
     DO J=1,JMAX2
        DO M=1,modes
           angsub(J,M,I+1-tstartind)=angle(J,M,I)
        ENDDO
     ENDDO
  ENDDO

  nout = 0.5*OFAC*HIFAC*nsub
  allocate(oneang(nsub))
  allocate(X1(nout))
  allocate(Y1(nout))
  allocate(results(JMAX,modes,nout))
  allocate(frequencies(JMAX,modes,nout))

  print*,"about to call periodogram"
!OMP PARALLEL DO DEFAULT(SHARED) &
!OMP PRIVATE(oneang,J,I,X1,Y1,prob)
  DO M=1,modes
     DO J=2,JMAX+1
        DO I = 1,nsub
           oneang(I) = angsub(J,M,I)
        ENDDO
        
        call FASPER(tsub,oneang,nsub,OFAC,HIFAC,X1,Y1,nout,prob)
        
        DO I = 1,nout
           results(J,M,I) = Y1(I)
           frequencies(J,M,I) = X1(I)
        ENDDO
     ENDDO
  ENDDO
!OMP END PARALLEL DO

  OPEN(UNIT=12,FILE=TRIM(outfile),FORM='UNFORMATTED')
  WRITE(12) JMAX,modes,nout
  WRITE(12) results
  WRITE(12) frequencies
  WRITE(12) omegavg
  WRITE(12) kappa
  CLOSE(12)

           
        

  
!...    Subroutine uses periodogram technique for generating
!...    power spectra (see references below).  The Power Spectra are
!...    plotted in the frequency-radius.  The spectra can be added into
!...    one curve for OMPAT extraction.  Note that OMPAT=2pi/PPAT.

END PROGRAM PERIODOGRAM

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
