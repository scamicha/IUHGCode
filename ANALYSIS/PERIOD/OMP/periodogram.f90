!-*-f90-*-
PROGRAM PERIODOGRAM
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)
  integer :: jmax, kmax, lmax, jreq, istart, iend, iskip, numargs
  integer :: jmax2, kmax2, lmax2,modes,numfiles,jreq,fileiter
  integer :: j,k,l,i,m,tstartind,tendind
  REAL(DOUBLE), PARAMETER :: torp = 1605.63
  REAL(DOUBLE), PARAMETER :: pi = 3.14159265358979323846d0
  REAL(DOUBLE), PARAMETER :: twopi = 2.d0*pi
  real(double) :: ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,OMMAX,tmassini
  real(double) :: dr,pi,torp,aujreq,tstart,tend
  real(double),dimension(:,:,:),allocatable :: rho,s,t,a,eps,omega
  real(double),dimesion(:,:,:),allocatable :: angle,angsub
  real(double),dimension(:),allocatable :: omegavg,kappa,massum,timearr,am,bm,tsub
  real(double),dimension(:,:),allocatable :: rhomid
  character :: rhodir*80, savedir*80,rhofile*94,savedfile*94
  character :: jmaxin*8,kmaxin*8,lmaxin*8,istartin*8,iendin*8,iskipin*8,modein*8
  character :: filenum*8,aujreqin*10,tstartin*8,tendin*8

  numargs = IARGC()

  if(numargs.ne. 11) then
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
  allocate(rhomid(jmax2,lmax))
  allocate(s(jmax2,kmax2,lmax))
  allocate(t(jmax2,kmax2,lmax))
  allocate(a(jmax2,kmax2,lmax))
  allocate(eps(jmax2,kmax2,lmax))
  allocate(omega(jmax2,kmax2,lmax))
  allocate(omegavg(jmax2))
  allocate(kappa(jmax2))
  allocate(massum(jmax2))
  allocate(angle(jmax2,modes,numfiles))

  write (filenum,'(I8.8)')iend
  savedfile=trim(savedir)//'saved.'//filenum
  OPEN(UNIT=8, FILE=trim(savedfile),FORM='UNFORMATTED',  &
       STATUS="OLD",ERR=911)      
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

!$OMP DO DEFAULT(SHARED) REDUCTION(+:omegavg,massum)
  DO J=2,jmax1
     omegavg(j) = 0.d0
     massum(j) = 0.d0
     DO k=1,kmax2
        DO l=1,lmax
           omegavg(j) = omegavg(j) + omega(j,k,l)
           massum(j) = massum(j) + rho(j,k,l)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
     
!$OMP DO DEFAULT(SHARED)
  DO J=2,jmax1
     omegavg(j) = omegavg(j)/massum(j)
     kappa(j) = sqrt(ABS(rhf(j)*omegavg(j)*(omegavg(j+1)- &
          omegavg(j-1))/dr + 4.d0*omegavg(j)**2))*(torp/(2.d0*pi))
     omegavg(j) = omegavg(j)*(torp/(2.d0*pi))
!!!!!!!!!!!!!!!!!!!!! If not needed put change rhf to AU here.........
     
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  DO I = 1,NUMFILES
     fileiter = (I*stepskip)+stepbegin
     write (filenum,'(I8.8)')fileiter
     rhofile = trim(rhodir)//'rho3d.'//filenum
     OPEN(UNIT=8, FILE=trim(rhofile),FORM='UNFORMATTED',  &
          STATUS='OLD',ERR=911)
     READ(8) RHO
     READ(8) time
     CLOSE(8)
     timearr(I)=time
!$OMP PARALLEL DO 
     DO L=1,LMAX
        DO J=1,JMAX2
           rhomid(J,L) = rho(J,2,L)
        ENDDO
     ENDDO
!OMP END PARALLEL DO


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
              am(j) = am(j)+rhomid(j,l)*cos(dble(m)*phi)
              bm(j) = bm(j)+rhomid(j,l)*sin(dble(m)*phi)
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
  DO WHILE (timearr(I).lt.(tstart*torp))
     I = I+1
  ENDDO
  
  tstartind = I

  DO WHILE (timearr(I).lt.(tend*torp))
     I = I+1
  ENDDO

  tendind = I
  
  allocate(tsub(tendind-tstartind))
  allocate(angsub(jmax2,modes,tendind-tstartind))

  DO I = tstartind,tendind
     tsub(I+1-tstartind)=timearr(I)
     DO J=1,JMAX2
        DO M=1,modes
           angsub(J,M,I+1-tstartind)=angle(J,M,I)
        ENDDO
     ENDDO
  ENDDO

  DO M=1,modes
     DO J=2,JMAX+1
        

  
!...    Subroutine uses periodogram technique for generating
!...    power spectra (see references below).  The Power Spectra are
!...    plotted in the frequency-radius.  The spectra can be added into
!...    one curve for OMPAT extraction.  Note that OMPAT=2pi/PPAT.

END PROGRAM PERIODOGRAM

SUBROUTINE PERIOD(

PRO slowperiod

TWOPI = 6.283185307179586476d0


MAXSIZE = 500000L

mintime    = 12.03d0
maxtime    = 19.57d0
str =''
timestring = ''
stepdum   = 0L
starcount = 0L
percount  = 0L
macc      = 4L
ofac      = 4L
hifac     = 0.5d0
xdum      = 0.d0
ydum      = 0.d0
rdum      = 0.d0
timedum   = 0.d0
angdum    = 0.d0
dz        = 0.202180656d0
mratio    = 8.075411422d0
pi        = 3.1415926353892d0
maxx      = 0.d0
minx      = 0.d0
maxy      = 0.d0
miny      = 0.d0
maxr      = 0.d0

stepcom = LONARR(MAXSIZE)
xcom    = DBLARR(MAXSIZE)
ycom    = DBLARR(MAXSIZE)
rcom    = DBLARR(MAXSIZE)
timecom = DBLARR(MAXSIZE)
angcom  = DBLARR(MAXSIZE)

tconv   = 1605.63d0

openr,lun1,"comfile50.dat",/GET_LUN

 WHILE NOT EOF(lun1) DO BEGIN
     READF,lun1,str,FORMAT='(A256)'
     READS,str,stepdum,timedum,xdum,ydum,rdum,angdum
     stepcom(starcount)=stepdum
     rcom(starcount) = rdum
     xcom(starcount) = xdum 
     ycom(starcount) = ydum
     timecom(starcount)=timedum
     angcom(starcount)=angdum
     starcount++
 ENDWHILE

 rcom    = rcom(0:starcount-1)
 xcom    = xcom(0:starcount-1)
 ycom    = ycom(0:starcount-1)
 stepcom = stepcom(0:starcount-1)
 timecom = timecom(0:starcount-1)
 angcom  = angcom(0:starcount-1)
 x       = timecom/tconv
 y       = ycom
 indexsub= WHERE((x gt mintime) and (x lt maxtime))
 x       = x(indexsub)
 y       = y(indexsub)
n       = starcount 
 wi      = DBLARR(n)
 wpi     = DBLARR(n)
 wpr     = DBLARR(n)
 wr      = DBLARR(n)
 nout    = LONG(0.5*ofac*hifac*n) 
 px      = DBLARR(nout)
 py      = DBLARR(nout) 
 ave     = MEAN(y,/DOUBLE)
 var     = VARIANCE(y,/DOUBLE)
 IF(var eq 0.d0) THEN BEGIN
     print,"Error: zero variance in period"
     RETURN
 ENDIF

 xmax   = MAX(x,MIN=xmin)
 xdif   = xmax-xmin
 xave   = 0.5d0*(xmax+xmin)
 pymax  = 0.d0
 jmax   = 0L
 pnow   = 1.d0/(xdif*ofac)
 arg    = TWOPI*((x-xave)*pnow)
 wpr    = -2.d0*(sin(0.5d0*arg))^2.d0
 wpi    = sin(arg)
 wr     = cos(arg)
 wi     = wpi
 FOR i=0L,nout-1L DO BEGIN
     px(i)  = pnow
     sumsh  = TOTAL(wr*wi,/DOUBLE)
     sumc   = TOTAL((wr-wi)*(wr+wi),/DOUBLE)
     wtau   = 0.5d0*ATAN(2.d0*sumsh,sumc)
     swtau  = SIN(wtau)
     cwtau  = COS(wtau)
     sums   = TOTAL(((wi*cwtau)-(wr*swtau))^2.d0,/DOUBLE)
     sumc   = TOTAL(((wr*cwtau)+(wi*swtau))^2.d0,/DOUBLE)
     yy     = y - ave
     sumsy  = TOTAL(((wi*cwtau)-(wr*swtau))*yy,/DOUBLE)
     sumcy  = TOTAL(((wr*cwtau)+(wi*swtau))*yy,/DOUBLE)
     wtemp  = wr
     wr     = ((wr*wpr)-(wi*wpi))+wr
     wi     = ((wi*wpr)+(wtemp*wpi))+wi
     py(i)  = 0.5d0*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var
     IF (py(i) ge pymax) THEN BEGIN
         pymax = py(i)
         jmax  = i 
     ENDIF
     pnow += 1.d0/(ofac*xdif)
 ENDFOR

 expy  = exp(-pymax)
 effm  = 2.d0*nout/ofac
 prob  = effm*expy
 IF (prob gt 0.01d0) THEN prob = 1.d0 - (1.d0-expy)^effm

 print, prob
