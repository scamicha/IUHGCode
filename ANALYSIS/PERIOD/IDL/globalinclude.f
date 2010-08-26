!-*-f90-*-
subroutine Globalcoef(prefix,TORP,JMAX,KMAX,LMAX,DR,TIME,CN)
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)

  real(double), allocatable, dimension(:,:,:) :: rho,rhosm
  real(double), allocatable, dimension(:,:)   :: an,bn
  real(double), allocatable, dimension(:)     :: rhf,rho0
  real(double), dimension(-1:JMAX,0:LMAX/2,0:1) :: cn
  real(double) :: time, angle, pi, torp, dr

  integer :: J, K, L, JMAX, KMAX, LMAX
  integer :: N, I, JMIN,PISYM

  character*(*)prefix 

! begin I/O

  pi = acos(-1.d0)
  JMIN  = 15
  PISYM = 0



! allocate variables
  IF(PISYM.EQ.1)THEN
     allocate( rho( -1:JMAX, -1:KMAX, 0:(2*LMAX)-1)  )
     allocate( rhosm( -1:JMAX, -1:KMAX, 0:LMAX-1)  )
  ELSE
     allocate( rho( -1:JMAX, -1:KMAX, 0:LMAX-1)  )
  ENDIF

  allocate(  an(-1:JMAX,0:LMAX/2)             )
  allocate(  bn(-1:JMAX,0:LMAX/2)             )
  allocate( rhf(-1:JMAX)                      )
  allocate(rho0(-1:JMAX)                      )
  
  do J = -1, JMAX
     rhf(J) = (dble(J)+.5d0)*dr
  enddo

! start loop and read from file
     print*,prefix

     print *," GLOBAL OUT-> opening file ",prefix

     open(unit=10,file=prefix(1:23),form="UNFORMATTED",ERR=10)
     
     IF (PISYM.EQ.1) THEN
        READ(10) rhosm
     ELSE
        read(10)rho
     ENDIF

     read(10)time
    
     close(10)

     IF(PISYM.EQ.1)THEN
        DO J=-1,JMAX
           DO K=-1,KMAX
              DO L=0,LMAX-1
                 RHO(J,K,L)=RHOSM(J,K,L)
              ENDDO
              DO L=LMAX,(2*LMAX)-1
                 RHO(J,K,L)=RHOSM(J,K,L-LMAX)
              ENDDO
           ENDDO
        ENDDO
     ENDIF

157  format(A27,1f9.3,A3,1f9.3)

     write(6,157), "GLOBAL OUT-> Time in ORPS (",torp,"): ",time/torp
     

        

!$OMP PARALLEL

!$OMP DO  
     do J = 0, JMAX-1
        rho0(J) = 0.d0
        do L = 0, LMAX-1
           
           rho0(J) = rho0(J) + rho(J,0,L)

        enddo
     enddo
!$OMP ENDDO

     rho0 = rho0/dble(LMAX)
     
!$OMP DO 
     do J = 0, JMAX
        an(J,0) = rho0(J)*256 ! I have no idea why this is times 256...
        do N = 1, LMAX/2
           
           an(J,N) = 0.d0
           bn(J,N) = 0.d0
           cn(J,N,0) = 0.d0
!           cn(J,N,1) = 0.d0

 
           if ( J >= JMIN) then
              do L = 0, LMAX-1
                 
                 angle = (dble(L))/(LMAX)*2.d0*pi
                 an(J,N) = an(J,N) + rho(J,0,L)*cos(dble(N)*angle)
                 bn(J,N) = bn(J,N) + rho(J,0,L)*sin(dble(N)*angle)            
              enddo
           endif

           an(J,N) = an(J,N)/(dble(LMAX))
           bn(J,N) = bn(J,N)/(dble(LMAX))
           
! compute the amplitude
           cn(J,N,0) = 2.d0*sqrt(an(J,N)**2+bn(J,N)**2)/rho0(J)

! compute the phase angle
           if (an(J,N) > 0 ) then
              cn(J,N,1) = atan(bn(J,N)/an(J,N))*180./pi
           else if (an(J,N) < 0 ) then
              cn(J,N,1) = (atan(bn(J,N)/an(J,N)) + pi)*180./pi
           else 
              if (bn(J,N) > 0 ) then
                 cn(J,N,1) = 90.d0
              else if (bn(J,N) < 0 ) then
                 cn(J,N,1) = -90.d0
              else
                 cn(J,N,1) = 720.
              endif
           endif

        enddo
        cn(J,LMAX/2,0) = cn(J,LMAX/2,0)/2.d0 
        cn(J,LMAX/2,1) = an(J,0) 
        
     enddo
!$OMP ENDDO

!$OMP END PARALLEL

     
10   continue

  deallocate(rho)
  deallocate(rho0)
  deallocate(an)
  deallocate(bn)
  deallocate(rhf)

  RETURN
END subroutine Globalcoef

SUBROUTINE SURFPRGM(nmax,JTOP,mode,cirp,nn,RMIN,RMAX,    & 
                    IK,aa,ti,CX,NPX)

!...    Subroutine uses periodogram technique for generating
!...    power spectra (see references below).  The Power Spectra are
!...    plotted in the frequency-radius.  The spectra can be added into
!...    one curve for OMPAT extraction.  Note that OMPAT=2pi/PPAT.
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)
  
  real(double),dimension(JTOP,nmax)        :: aa
  real(double),dimension(JTOP,2*nmax)      :: c
  real(double),dimension(2*nmax,JTOP)      :: cx
  real(double),dimension(nmax)             :: ti,bp,bpp,npp,x
  real(double),dimension(2*nmax+1)         :: power
  real(double),dimension(2*nmax)           :: npx
  real(double)                             :: omin,omax,delta,nnmax,cirp
  real(double)              ::  omegamin,omegamax,delplot,dnpr,dj,dome,pi
  real(double) :: minb,maxb,dox,doy,bigom,bigb,dpx,dpy,dpx1,dpy1,dpx2,dpy2
  integer      :: j,n,np,rmin,rmax,mode,nmax,k,p,addmin,addmax,ik,JTOP
  integer      :: nn,rmintemp,rmaxtemp
  character       mm*2,rmi*3,rma*3,omi*7,oma*7
  character       filename*30,bg*6,psfile*30,amin*2,amax*3

  pi = ACOS(-1.d0)

  write(mm,'(i2)')mode
  write(rmi,'(i3)')rmin
  write(rma,'(i3)')rmax
  

  rmintemp=rmin-1
  rmaxtemp=rmax-1

  do n=1,nmax
     do j=1,JTOP
        aa(j,n)=cos(aa(j,n)*pi/180.)	
     end do
  end do
        
  do j=rmintemp,rmaxtemp
     do n=1,nmax
        x(n)=aa(j,n)
        power(n)=0.
     ENDDO
     print*,"Omin = ", omin
     call periodogram(x,ti,nmax,omin,omax,delta,cirp,mode,ik,power)
!...Where does this power() comes from? 

     do np=1,ik
        c(j,np)=power(np)
        cx(np,j)=c(j,np) ! ? Why flip? 
     ENDDO
     print*,'Finished j= ',j+1
  ENDDO

  do n=1,nmax
     npp(n)=(omin + delta*(n-1))*(cirp/(2.*pi*mode))
  end do
  omin = omin*(cirp/(2.*pi*mode))
  omax = omax*(cirp/(2.*pi*mode))
  delta = delta*(cirp/(2.*pi*mode))


  
  if(nn.eq.2) then
     do n=1,ik   
        npx(n)= omin + (delta)*(n-1)
     end do
  end if

  RETURN
END SUBROUTINE SURFPRGM
 
SUBROUTINE PERIODOGRAM(x,ti,nmax,omin,omax,delta,cirp,mode,ik,power)
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)

!
!...Periodogram...
!

!  THIS PROGRAM DOES A PERIOD SEARCH BY THE TECHNIQUE OF 
!   SCARGLE (APJ 263,835,1982) AS MODIFIED BY HORNE AND BALIUNAS
!   (APJ 302,757, 1986). 

  real(double),dimension(nmax)          :: x,ti
  real(double),dimension(2*nmax+1)      :: power
  REAL(double)          :: SMEGA,OMEGA,T,TMIN,omin,omax,delta,cirp
  real(double)          :: TWOPI,FNI,SUM,AVE,PMAX,PMIN,TOP,BOT,THETA1
  real(double)          :: THETA2,XC,THETA,TAU,A,B,AA,BB,C,S,VAR
  INTEGER               :: N0,NMAX,mode,ARRAYSIZE,IK,M,L,j,K

  
  TWOPI=2.0*ACOS(-1.d0)


!  COMPUTE THE NUMBER OF INDEPENDENT FREQUENCIES

  FNI=-6.362+1.193*N0+0.00098*N0**2


!  COMPUTE THE VARIANCE

  n0=nmax
  SUM=0.0
  DO M=1,N0
     SUM=SUM+X(M)
  ENDDO
  AVE=SUM/N0
  SUM=0.
  DO M=1,N0
     SUM=SUM+(X(M)-AVE)**2
  ENDDO
  VAR=SUM/N0

!  RESCALE THE DATA BY SUBTRACTING THE MEAN

  DO M=1,N0
     X(M)=X(M)-AVE
  ENDDO

!  RESCALE THE TIME DATA BY SUBTRACTING THE FIRST TIME

  DO M=2,N0
     TI(M)=TI(M)-TI(1)
  ENDDO
  TI(1)=0.0

!  COMPUTE THE TIME INTERVAL OF OBSERVATION

  T=TI(N0)
!
!  COMPUTE THE FREQUENCY LIMITS AND INCREMENTS
!
!     FIND THE  MIN TIME INTERVAL BETWEEN OBSERVATIONS
  TMIN=TI(2)-TI(1)
  DO L=3,N0
     IF(TI(L)-TI(L-1) .LT. TMIN) TMIN=TI(L)-TI(L-1)
     if(tmin .eq. 0.) THEN
        write(6,1244) l
1244    format(3x,'time interval at data point ',i4,' is zero.')
        stop
     endif
  ENDDO
  OMIN=TWOPI/T

! ...  For MIRP case

  OMAX=4.0*mode*twopi/(cirp)

! ...  For using ORP to analyze Tcool = const. runs

  PMAX=TWOPI*(1.0/OMIN)
  PMIN=TWOPI*(1.0/OMAX)
  IK=2.*N0 
!       
!  COMPUTE THE SIZE OF THE FREQUENCY STEPS
!
  print*," Omin per = ",Omin
  print*," Omax per = ",Omax
  DELTA=(OMAX-OMIN)/IK
  print*," Delta = ",delta

!  INITIALIZE THE FREQUENCY

  OMEGA=OMIN

!  COMPUTE TAU

!      LOOP OVER THE FREQ RANGE

  DO J=1,IK + 1

!      TAKE A VERY SMALL STEP BACK IN FREQ.

     SMEGA=OMEGA-0.0001*DELTA
     CALL TTT(TI,SMEGA,N0,TOP,BOT,nmax)
     THETA1=ATAN(TOP/BOT)
     CALL TTT(TI,OMEGA,N0,TOP,BOT,nmax)
     THETA2=ATAN(TOP/BOT)

!  CHECK FOR QUADRANT AMBIGUITY

     XC=TOP/BOT
     IF((XC.GE.0.0).AND.(THETA2.GT.THETA1)) THETA=THETA2
     IF((XC.GE.0.0).AND.(THETA2.LT.THETA1)) THETA=3.14159265-THETA2
     IF((XC.LE.0.0).AND.(THETA2.LT.THETA1)) THETA=3.14159265+THETA2
     IF((XC.LE.0.0).AND.(THETA2.GT.THETA1)) THETA=6.28318531-THETA2
     
     TAU=THETA/(2.0*OMEGA)

!  COMPUTE SUMS

     A=0.
     B=0.
     AA=0.
     BB=0.
     DO K=1,N0
        A=A+(COS(OMEGA*(TI(K)-TAU)))**2
        B=B+(SIN(OMEGA*(TI(K)-TAU)))**2
        AA=AA+X(K)*COS(OMEGA*(TI(K)-TAU))
        BB=BB+X(K)*SIN(OMEGA*(TI(K)-TAU))
     ENDDO

!  COMPUTE POWER AND PROBABILITY

     C=AA*1.0/SQRT(A)
     S=BB*1.0/SQRT(B)
     POWER(j)=0.5*(C**2+S**2)/VAR
     OMEGA=OMEGA+DELTA
  ENDDO
  RETURN
END SUBROUTINE PERIODOGRAM

SUBROUTINE TTT(TI,OMEGA,N0,TOP,BOT,ARRAYSIZE)

!...Companion to Periodogram
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)

  REAL(double),dimension(ARRAYSIZE)   :: TI
  REAL                                :: OMEGA,E,F,TOP,BOT
  INTEGER                             :: ARRAYSIZE,K,N0

  TOP=0.0
  BOT=0.0
  DO K=1,N0
     E=SIN(2.0*OMEGA*TI(K))
     F=COS(2.0*OMEGA*TI(K))
     TOP=TOP+E
     BOT=BOT+F
  ENDDO
  RETURN
END SUBROUTINE TTT


!SUBROUTINE PERIOD(b,JTOP,COUNT,c)
!
!...    Auto-Correlation Function Subroutine
!
!  implicit none
!  integer, parameter :: double = selected_real_kind(15,300)

!  real(double), dimension(JTOP,COUNT) :: aa,c,b
!  real(double)                        :: pi,pioneighty
!  integer                             :: count,j,n,JTOP,I
  
!  pi = ACOS(-1.d0)
!  pioneighty = pi/180.d0
	
!  DO j=1,JTOP
!     DO n=1,count
!        aa(j,n)=cos(b(j,n)*pioneighty)
!     ENDDO
!  ENDDO

!  do n=1,count/2
!     do j=1,JTOP
!        c(j,n)=0.0
!        do i=1,count-n
!           c(j,n)=c(j,n)+aa(j,i)*aa(j,i+n)	
!        end do
!     end do
!  end do

!  do n=1,count/2
!     do j=1,JTOP
!        c(j,n)=c(j,n)/(count-n)
!     end do
!  end do
	
!  RETURN
!END SUBROUTINE PERIOD

