      PROGRAM  SAPOT3
C*******************************************************************************
C  Revision notes     (Don Berry)
C
C  v-0: Original code.
C
C  v-1: o Changed DO loops to use DO...ENDDO construct.
C       o Moved most index arithmetic directly into array subscript lists.
C       o Completely reworked calculation of philag to use comprehensible
C           IF..ELSE..ENDIF constructs.
C       o Inserted DOACROSS directives.
C
C*******************************************************************************


     
      IMPLICIT real*8 (a-h,o-z)      

      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'
      INCLUDE 'globals_mpi.h'
      INCLUDE 'mpif.h'

      integer COUNTER
      INTEGER AMCOUNT,answer_width

      COMMON /COEFS/COEF(pot3JMAX2,pot3KMAX2,LMAX2,2)

      COMMON /BLOK6/DTHETA,PI,GRAV
      real*8 DTHETA, PI, GRAV
      real*8 dummy,dummy1,phimode
      real(KIND=8),DIMENSION(:,:,:),ALLOCATABLE  :: gt,
     &     alpha_grav,alpha_sum
      real(KIND=8),DIMENSION(:,:),ALLOCATABLE :: rt,alpha_reyn
      real(KIND=8),DIMENSION(:,:),ALLOCATABLE :: omega_full
      real*8, dimension(JMAX2,KMAX2) :: axi
      real*8, dimension(JMAX2,KMAX2) :: aofm, bofm
      real*8 timer(2)
      real*8,dimension(:),ALLOCATABLE   :: timearr
      real*8,dimension(4) :: startup_array

C$OMP THREADPRIVATE(/BLOK6/)

C A1 and B1  should be dimensioned lmax.  They're used in fft.
      DIMENSION A1(LMAX),B1(LMAX)

C These next arrays are used in blktri. The size of wfw(max) is given by
C max. ge. (2*(kmax3+2)*(log2(kmax3+1)-1) + 6*jmax3+2).
      PARAMETER (JMAX3=JMAX-1,KMAX3=KMAX-1,LMAX1=LMAX2+1)
      DIMENSION AN(KMAX3),BN(KMAX3),CN(KMAX3),
     &          AM(JMAX3),BM(JMAX3),CM(JMAX3),
     &          Y(JMAX3,KMAX3),WFW(2*(KMAX3+2)*4+6*JMAX3+2)

C The following arrays store quantities that are used to calculate blktri
C coefficients.  C(jmax-1) stores laplacian's angular operator.
      DIMENSION C(JMAX3),DENOMR(JMAX),RD2(JMAX),RD3(JMAX),DENOMZ(KMAX),
     &     ZD2(KMAX),XLAMM(LMAX1)

C These common blocks are used in blktri and its subroutines. (I should not
C need to declare them in pot3, but I am not sure OpenMP compiler will do
C the right thing if I do not.)
      COMMON /CBLKT/ NPP,Kx,EPSLON,CNV,Lx,NCMPLX,IK,IZ,DUM(1)
      COMMON /BLOKJ1/IWCN,IW1,IW2,IW3,IWD,IWW,IWU
C$OMP THREADPRIVATE(/CBLKT/,/BLOKJ1/)


      LOGICAL  INIT_BLKTRI

      CHARACTER savedfile*80,phifile*80,filenum*6
      INTEGER I, NUMFILES,PROCESS,m_return,sender
      INTEGER status(MPI_STATUS_SIZE)
CSAM....Read a saved file to set up r and z grid

      type answer_return
      sequence
      real*8, dimension(JMAX2) :: rt,gt,alpha_r,alpha_g
      real*8 :: omega_full(JMAX2)
      end type answer_return

      type (answer_return) answer

      call MPI_INIT(mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numranks, mpierr)

      ISYM   = 2
      MAXTRM = 10
      REDGE  = 0.d0
      GRAV = 1.0
      PI = ACOS(-1.0)
      DTHETA = 2.0*PI/dble(LMAX)

      IF(myrank.eq.0) THEN

         write (filenum,'(I6.6)')IEND
         phifile=trim(outdir)//trim(outfile)//filenum      


         COUNTER=0
         NUMFILES = ((IEND - ISTART)/ISKIP) + 1
      
         ALLOCATE(timearr(NUMFILES))
         ALLOCATE(rt(JMAX2,NUMFILES))
         ALLOCATE(alpha_reyn(JMAX2,NUMFILES))
         ALLOCATE(omega_full(JMAX2,NUMFILES))
         ALLOCATE(gt(JMAX2,LMAX2+1,NUMFILES))
         ALLOCATE(alpha_grav(JMAX2,LMAX2+1,NUMFILES))
         ALLOCATE(alpha_sum(JMAX2,LMAX2+1,NUMFILES))

         write (filenum,'(I6.6)')ISTART
         savedfile=trim(datadir)//'saved.'//filenum
         print *, savedfile
         OPEN(UNIT=8, FILE=trim(savedfile),FORM='UNFORMATTED',
     &        STATUS="OLD",ERR=911)
      
         READ(8) S
         READ(8) T
         READ(8) A
         READ(8) RHOSAVE
         READ(8) EPS
         READ(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
         read(8) tmassini
         CLOSE(8)

         do I=ISTART,IEND,ISKIP
      
            first_pass = 1

            startup_array(1) = ROF3N
            startup_array(2) = ZOF3N
            startup_array(3) = DEN
            startup_array(4) = tmassini            
         
            COUNTER = COUNTER + 1

            timearr(COUNTER) = time/1605.63      
         
            print *, ' Your time is ', time/1605.63

            AMCOUNT = 0

            call MPI_BCAST(RHOSAVE,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(EPS,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(S,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(A,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(startup_array,4,MPI_DOUBLE_PRECISION,0,
     &           MPI_COMM_WORLD,mpierr)
            call MPI_BCAST(JREQ,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(first_pass,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           mpierr)
            
            DO PROCESS=1,min(numranks-1,LMAX2)
               call MPI_SEND(AMCOUNT,1,MPI_INTEGER,PROCESS,AMCOUNT,
     &              MPI_COMM_WORLD,mpierr)
               AMCOUNT = AMCOUNT+1
            ENDDO

            AMCOUNT = AMCOUNT-1
            answer_width = 5*JMAX2
            IF (I.lt.IEND) THEN
               write (filenum,'(I6.6)')I+ISKIP
               savedfile=trim(datadir)//'saved.'//filenum
               print *, savedfile
               OPEN(UNIT=8, FILE=trim(savedfile),FORM='UNFORMATTED',
     &              STATUS="OLD",ERR=911)
      
               READ(8) S
               READ(8) T
               READ(8) A
               READ(8) RHOSAVE
               READ(8) EPS
               READ(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
               read(8) tmassini
               CLOSE(8)
            ENDIF

            DO PROCESS=1,LMAX/2+1
               call MPI_RECV(answer,answer_width,MPI_DOUBLE_PRECISION,
     &              MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,
     &              mpierr)
               sender = status(MPI_SOURCE)
               m_return = status(MPI_TAG)
               
               DO J = 1,JMAX2
                  gt(J,m_return+1,COUNTER) = answer%gt(J)
                  alpha_grav(J,m_return+1,COUNTER) = answer%alpha_g(J)
               ENDDO

               IF(m_return.eq.0) THEN
                  DO J = 1,JMAX2
                     alpha_reyn(J,COUNTER) = answer%alpha_r(J)
                     rt(J,COUNTER) = answer%rt(J)
                     omega_full(J,COUNTER) = answer%omega_full(J)
                  ENDDO
               ENDIF
               
               IF(AMCOUNT.lt.LMAX2) THEN
                  call MPI_SEND(AMCOUNT+1,1,MPI_INTEGER,sender,
     &                 AMCOUNT+1,MPI_COMM_WORLD,mpierr)
                  AMCOUNT = AMCOUNT+1
               ELSE
                  call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER,sender,LMAX2+1,
     &                 MPI_COMM_WORLD,mpierr)

               ENDIF
            ENDDO
         ENDDO
               
      ELSE
         do I=ISTART,IEND,ISKIP
            
            call MPI_BCAST(RHOSAVE,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(EPS,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(S,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(A,JMAX2*KMAX2*LMAX,
     &           MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(startup_array,4,MPI_DOUBLE_PRECISION,0,
     &           MPI_COMM_WORLD,mpierr)
            call MPI_BCAST(JREQ,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           mpierr)
            call MPI_BCAST(first_pass,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           mpierr)

            ROF3N = startup_array(1)
            ZOF3N = startup_array(2)
            DEN =   startup_array(3)
            tmassini = startup_array(4)

            IF (myrank.lt.LMAX/2+1) THEN
CSAM....Set up grid    
               DO j=1, JMAX2
                  r(j)=(j-2)*ROF3N
                  rhf(j) =((j-2)*ROF3N)+(ROF3N/2)
               ENDDO
            
               DO k=1,KMAX2
                  z(k) = (k-2)*ZOF3N
                  zhf(k) = ((k-2)*ZOF3N)+(ZOF3N/2)
               ENDDO
               rholmt = den*1.d-12

               CALL SETBDY(0,ISYM)

 300           call MPI_RECV(AMCOUNT,1,MPI_INTEGER,0,MPI_ANY_TAG,
     &              MPI_COMM_WORLD,status,mpierr)

!!!!!!!!!!!!!!!!!!!!!! work is done jump to end of loop
               IF(status(MPI_TAG).eq.LMAX2+1) go to 400

!$OMP PARALLEL DO         
               DO J=1,JMAX2
                  DO K=1,KMAX2
                     axi(J,K) = 0.d0
                     aofm(J,K) = 0.d0
                     bofm(J,K) = 0.d0
                     DO L=1,LMAX
                        rho(J,K,L) = 0.d0
                     ENDDO
                  ENDDO
               ENDDO
!$OMP END PARALLEL DO
                  
C$OMP PARALLEL DEFAULT(PRIVATE) SHARED(RHOSAVE,axi,aofm,bofm,AMCOUNT)
C$OMP&  PRIVATE(J,K,L)
C$OMP&  COPYIN(/BLOK6/,/GRID/)

!$OMP DO REDUCTION(+:axi)

               do L = 1, LMAX
                  do K = 2, KMAX
                     do J = 2, JMAX
                  
                        axi(J,K) = axi(J,K) + rhosave(J,K,L)*dtheta/
     &                       pi
                     
                     enddo
                  enddo
               enddo

!$OMP END DO

!$OMP DO REDUCTION(+:aofm)
               DO L = 1,LMAX
                  DO K = 2,KMAX
                     DO J = 2,JMAX
                        aofm(J,K) = aofm(J,K) + rhosave(J,K,L) 
     &                       * cos(AMCOUNT*L*dtheta)*dtheta/
     &                       pi
                     ENDDO
                  ENDDO
               ENDDO
!$OMP END DO

!$OMP DO REDUCTION(+:bofm)
               DO L = 1,LMAX
                  DO K = 2,KMAX
                     DO J = 2,JMAX
                        bofm(J,K) = bofm(J,K) + rhosave(J,K,L) 
     &                       * sin(AMCOUNT*L*dtheta)*dtheta/
     &                       pi
                     ENDDO
                  ENDDO
               ENDDO
!$OMP END DO

!$OMP END PARALLEL

C$OMP PARALLEL DEFAULT(PRIVATE) PRIVATE(J,K,L) 
C$OMP&  SHARED(PHI,RHO,COEF,axi,aofm,bofm,AMCOUNT,rholmt,RHOSAVE)
C$OMP&  COPYIN(/BLOK6/,/GRID/) 

!$OMP DO FIRSTPRIVATE(AMCOUNT,rholmt,RHOSAVE)
               do L = 1, LMAX
                  do K = 2, KMAX
                     do J = 2, JMAX
                        IF (AMCOUNT.EQ.0) THEN
                           rho(J,K,L)=RHOSAVE(J,K,L)
                        ELSE
                           rho(J,K,L) = axi(J,K)*0.5d0
     &                          +aofm(J,K)*cos(AMCOUNT*L*dtheta)
     &                          +bofm(J,K)*sin(AMCOUNT*L*dtheta)
                        ENDIF

                        if (rho(J,K,L) < 0.d0 ) rho(J,K,L) = rholmt

                     enddo
                  enddo
               enddo
!$OMP END DO            
!$OMP END PARALLEL

               CALL BDYGEN(MAXTRM,ISYM,REDGE)

C$OMP PARALLEL DEFAULT(PRIVATE) SHARED(PHI,RHO,COEF)
C$OMP&  PRIVATE(J,K,L,M)
C$OMP&  COPYIN(/BLOK6/,/GRID/)

!                  dbgf(dbgl)=DBG_POT3
!                  dbgl=dbgl+1

               N=LMAX/2
               SF1=0.5/DBLE(N)
               SF2=0.5

               DO L=1,LMAX                                                     
                  A1(L)=0.0
                  B1(L)=0.0
               ENDDO
      
C For convenience, put all rho's into phi array -- leave boundary phi's
C alone.
C$OMP DO
               DO L=1,LMAX                                                     
                  DO K=2,KMAX
                     DO J=2,JMAX
                        PHI(J,K,L)=RHO(J,K,L)
                     ENDDO
                  ENDDO
               ENDDO
C$OMP ENDDO


C-----------------------------------------------------------------------
C
C  If lmax>1 then for all j,k, calculate fourier transform of rho's
C  and boundary phi's.

               IF(LMAX.gt.1) THEN
         
C$OMP   DO
                  DO K=2,KMAX1
                     DO J=2,JMAX1
                        DO L=1,N
                           A1(L) = PHI(J,K,2*L-1)
                           B1(L) = PHI(J,K,2*L)
                        ENDDO 
                        CALL FFT(A1,B1,N,N,N,1)
                        CALL REALTR(A1,B1,N,1)
C     Cosine coefficients are in a1? Sine coef's are in b1.
C     Put cosine coef's in phi(j,k,1:n+1)?  Put sine coef's in      
C     phi(j,k,n+2:lmax).  Normalize computed values with 0.5/n.     
C     Compute amplitudes of modes.
                        X1=ABS(A1(1))
                        IF(X1.NE.0.0) THEN
                           DO L=2,N
                              COEF(J,K,L-1,1) = 2.0*SQRT( A1(L)**2 + 
     &                             B1(L)**2 )/X1
                           ENDDO
                           COEF(J,K,N,1) = ABS(A1(N+1))/X1
                        ENDIF
                           
C     Compute phase angles of modes. The following code simply calculates
C     philag=arctan(b/a) in the range -pi/2 <= philag < 3*pi/2. It sets
C     philag=4*pi if a and b are both 0.
                        DO L=2,N
                           IF(A1(L).GT.0.0) THEN
                              PHILAG = ATAN(B1(L)/A1(L))
                           ELSE IF(A1(L).LT.0.0) THEN
                              PHILAG = ATAN(B1(L)/A1(L))+PI
                           ELSE
                              IF(B1(L).GT.0.0) THEN
                                 PHILAG = 0.5*PI
                              ELSE IF(B1(L).LT.0.0) THEN
                                 PHILAG = -0.5*PI
                              ELSE
                                 PHILAG = 4.0*PI
                              ENDIF
                           ENDIF
                           COEF(J,K,L-1,2)=PHILAG*180./PI
                        ENDDO
                        COEF(J,K,N,2)=A1(1)
                           
C     Copy fourier transform back into phi array.
                        DO L=1,N+1
                           PHI(J,K,L)=A1(L)*SF1
                           IF(L.LT.N) THEN
                              PHI(J,K,L+N+1)=B1(L+1)*SF1
                           ENDIF
                        ENDDO
                        DO L=N+1,LMAX
                           A1(L) = 0.0
                           B1(L) = 0.0
                        ENDDO
                     ENDDO
                  ENDDO                                                           
C$OMP   ENDDO
         
               ENDIF
            

C-----------------------------------------------------------------------
C  Now solve the transformed problem. The fourier transform in the L
C  direction has turned the 3-D Poisson problem into a set of LMAX
C  uncoupled 2-D problems which may be solved concurrently.
C     

               DTHET2 = 1.0/(DTHETA*DTHETA)
               PIG4   = 4.0*PI*GRAV
               RD2(1) = 0.5/RHF(2)
               RD3(1) = 1.0/RHF(1)
               DO J=2,JMAX                                                    
                  DENOMR(J)=1.0/(RHF(J+1)-RHF(J-1))
                  RD2(J) = 1.0/(RHF(J+1)-RHF(J))
                  RD3(J) = 1.0/RHF(J)
               ENDDO
               ZD2(1) = 0.5/ZHF(2)
               DO K=2,KMAX                                                    
                  DENOMZ(K) = 1.0/(ZHF(K+1)-ZHF(K-1))
                  ZD2(K) = 1.0/(ZHF(K+1)-ZHF(K))
               ENDDO                                                           
               DO J=2,JMAX                                                    
                  CM(J-1) = 2.*(RD2(J)+0.5*RD3(J))*DENOMR(J)
                  AM(J-1) = 2.*(RD2(J-1)-0.5*RD3(J))*DENOMR(J)
                  C(J-1)  = RD3(J)*RD3(J)*DTHET2
               ENDDO                                                           
               DO K=2,KMAX                                                    
                  CN(K-1) = 2.*ZD2(K)*DENOMZ(K)
                  AN(K-1) = 2.*ZD2(K-1)*DENOMZ(K)
                  BN(K-1) = -CN(K-1)-AN(K-1)
               ENDDO
                  
C  Special conditions:
               CMMAX = CM(JMAX-1)
               CM(JMAX-1) = 0.0
               AMMIN = AM(1)
               AM(1) = 0.0
               CNMAX = CN(KMAX-1)
               CN(KMAX-1) = 0.0
               BN(1) = -CN(1)
               AN(1) = 0.0
               DO M=1,N+1                                                     
                  XLAMM(M)=COS((M-1)*DTHETA)
               ENDDO
                  

C  Coefficients bm and y will vary with m, so they haven't been calculated
C  yet. Now for each value of l, thru lmax, calculate phi's in transformed
C  space.
      
               INIT_BLKTRI=.TRUE.
C$OMP DO
               DO L=1,LMAX                                                
                  IF(L.LE.N+1) THEN
                     M=L-1                                                  
                  ELSE
                     M=L-(N+1)                                                 
                  ENDIF
                  POWRM=(-1.0)**M
                     
C (remember, densities are in array phi right now.)
                  DO K=2,KMAX                                                    
                     DO J=2,JMAX
                        Y(J-1,K-1)=PIG4*PHI(J,K,L)
                     ENDDO
                  ENDDO                                                           

C  Treat y on boundaries where phi is known.
                  DO J=2,JMAX                                                    
                     Y(J-1,KMAX-1)=Y(J-1,KMAX-1)-CNMAX*PHI(J,KMAX+1,L)
                  ENDDO
                  DO K=2,KMAX                                                    
                     Y(JMAX-1,K-1)=Y(JMAX-1,K-1)-CMMAX*PHI(JMAX+1,K,L)
                  ENDDO
                     
                  BM(1) = -CM(1)+(POWRM-1.0)*AMMIN+2.0*(XLAMM(M+1)-1.0)*
     &                 C(1)
                  DO J=2,JMAX-2
                     BM(J) = -CM(J)-AM(J)+2.0*(XLAMM(M+1)-1.0)*C(J)
                  ENDDO
                  BM(JMAX-1) = -CMMAX-AM(JMAX-1)+2.0*(XLAMM(M+1)-1.0)*
     &                 C(JMAX-1)
C  This completes the set up of coefficients for given m.
C     
                  IF(INIT_BLKTRI) THEN
                     CALL BLKTRI(0,1,KMAX-1,AN,BN,CN,1,JMAX-1,
     &                    AM,BM,CM,JMAX-1,Y,IERROR,WFW)
                     INIT_BLKTRI=.FALSE.
                  ENDIF
                  CALL BLKTRI(1,1,KMAX-1,AN,BN,CN,1,JMAX-1,
     &                 AM,BM,CM,JMAX-1,Y,IERROR,WFW)         
C     Solution on 2-D grid is complete.

C     Put transformed phi's from y into phi array.
                  DO K=2,KMAX                                                    
                     DO J=2,JMAX
                        PHI(J,K,L)=Y(J-1,K-1)
                     ENDDO
                  ENDDO
               ENDDO
C$OMP END DO

C-----------------------------------------------------------------------
C
C     All transformed phi's have been calculated.  Now obtain real
C     phi's by inverse fourier transform.

               IF(LMAX.NE.1) THEN
C$OMP   DO
                  DO K=2,KMAX1                                                   
                     DO J=2,JMAX1
                        DO L=1,LMAX
                           A1(L)=0.0
                           B1(L)=0.0
                        ENDDO                                                           
                        DO L=1,N+1
                           A1(L)=PHI(J,K,L)
                           IF(L.LT.N) THEN
                              B1(L+1)=PHI(J,K,L+N+1)
                           ENDIF
                        ENDDO                                                 
                        CALL REALTR(A1,B1,N,-1)
                        CALL FFT(A1,B1,N,N,N,-1)
C     Now real phi's are in first n places of both a1 and b1.  With
C     increasing l, phi's alternate between a1 and b1.  Put phi's into
C     array phi, normalizing them by 0.5.
C     
                        DO L=1,N
                           PHI(J,K,2*L-1) = SF2*A1(L)
                           PHI(J,K,2*L)   = SF2*B1(L)
                        ENDDO
                     ENDDO
                  ENDDO
C$OMP   END DO
               ENDIF                                                        

C     Calculation of phi's is now finished. 

C$OMP DO
               DO L=1,N                                                       
                  DO K=2,KMAX1                                                   
                     PHI(1,K,L)   = PHI(2,K,N+L)
                     PHI(1,K,N+L) = PHI(2,K,L)
                  ENDDO                                                           
                  DO J=1,JMAX1                                                   
                     PHI(J,1,L)   = PHI(J,2,L)
                     PHI(J,1,N+L) = PHI(J,2,N+L)
                  ENDDO
               ENDDO                                                           
C$OMP ENDDO

C$OMP END PARALLEL

               CALL ZAXPHI(10,0)
               
               CALL torqueout(AMCOUNT)
            
!     dbgl=dbgl-1

               do J = 3, JMAX
                  gt1(J) = gt1(J-1) + gt1(J) * (r(J)**2 - r(J-1)**2)*
     &                 econv * dtheta * zof3n
                  answer%gt(J) = gt1(J)
                  answer%rt(J) = treyn(J)*2.0*pi*r(J)**2*econv
               enddo

               tgrav(3:JMAX) = answer%gt(3:JMAX)/
     &              (2.d0*pi*r(3:JMAX)**2*econv)  

               do J = 3, JMAX

                  dummy = 2.d0/3.d0
                  if ( omega_fc(J+1) > 0.d0 .and. omega_fc(J) > 0.d0)
     &                 dummy = 1.d0/abs(LOG(omega_fc(J+1)/omega_fc(J))/
     &                 LOG(r(J+1)/r(J)))


                  answer%alpha_g(J) = dummy * tgrav(J)/
     &                 cs_sig(J)
               ENDDO
               IF(AMCOUNT.EQ.0)THEN
                  DO J = 3,JMAX
                     answer%alpha_r(J) = -dummy * treyn(J)/
     &                    cs_sig(J) ! should be negative
                     answer%omega_full(J) = omega_fc(J)
!                  treyn(J) = 0.d0
!                  cs_sig(J) = 0.d0
!                  omega_fc(J) = 0.d0
                  ENDDO   
               
!               DO J = 1,2
!                  treyn(J) = 0.d0
!                  treyn(JMAX+J) = 0.d0
!                  cs_sig(J) = 0.d0
!                  cs_sig(JMAX+J) = 0.d0
!                  omega_fc(J) = 0.d0
!                  omega_fc(JMAX+J) = 0.d0
!               ENDDO
               ENDIF

            
               call MPI_SEND(answer,5*JMAX2,MPI_DOUBLE_PRECISION,0,
     &              AMCOUNT,MPI_COMM_WORLD,mpierr)
               first_pass = 0

               GO TO 300
            ENDIF
               
 400        CONTINUE

         ENDDO
      ENDIF
 911  CONTINUE      


      IF(myrank.eq.0) THEN
         print *, "The number of files used is ", COUNTER

!      gt = gt/dble(COUNTER)
!      treyn = treyn/dble(COUNTER)

!      omega_fc = omega_fc/dble(COUNTER)
!      cs_sig = cs_sig/dble(COUNTER)






       !.......Writing out the phi data.
         OPEN(UNIT=15,FILE=phifile, FORM='UNFORMATTED')
         WRITE(15)JMAX2,LMAX2+1,NUMFILES,COUNTER
         WRITE(15)ISTART,IEND,ISKIP
         WRITE(15)timearr
         WRITE(15)gt
         WRITE(15)rt
         WRITE(15)alpha_grav
         WRITE(15)alpha_reyn
         WRITE(15)alpha_sum
         WRITE(15)omega_full
         CLOSE(15)

         DEALLOCATE(timearr)
         DEALLOCATE(gt)
         DEALLOCATE(rt)
         DEALLOCATE(alpha_grav)
         DEALLOCATE(alpha_reyn)
         DEALLOCATE(alpha_sum)
         DEALLOCATE(omega_full)
       
      ENDIF
      
      call MPI_FINALIZE(mpierr)

      END                                                                


c*******************************************************************************

      SUBROUTINE ZAXPHI(NPOINT,IPRINT)
      IMPLICIT real*8 (a-h,o-z)      

      include 'hydroparam.h'
      include 'globals.h'

      COMMON /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
      COMMON /COEFS/COEF(POT3JMAX2,POT3KMAX2,LMAX2,2)

      DIMENSION X(10),PHI2(10),PH(LMAX),PHAC(KMAX2)
C     PHAC IS DIMENSIONED KMAX2? PH IS DIMENSIONED LMAX/2.
C But it would not hurt to dimension PH to LMAX.

 
!      dbgf(dbgl)=DBG_ZAXPHI
!      dbgl=dbgl+1

C With X and PHI2 dimensioned 10, 10-point interpolation is the
C maximum you can do.
      IF(NPOINT.GT.10) THEN
        WRITE(3,10000) NPOINT
10000 FORMAT(///,' YOU ARE LIMITED TO A 10-POINT INTERPOLATION HERE IN',
     &       ' SUBROUTINE ZAXPHI.',/,' BUT YOU PUT NPOINT =',I4,///)
        NPOINT=10
      ENDIF
      XPOINT=0.0
      LM=LMAX/2
      N=NPOINT
      ICHK=MOD(N,2)
      N=N-ICHK
      NUM=N/2
      ISPECL=NUM+NUM+1
      JSP=NUM+1
      DO J=2,JSP
        I1=NUM+(J-1)
        I2=NUM-(J-2)
        X(I1)=RHF(J)
        X(I2)=-RHF(J)
      ENDDO
      IF(ICHK.EQ.1) X(ISPECL)=RHF(NUM+2)

      DO 50 K=2,KMAX1
        PHMAX=0.0
C       For given k, calculate phi on z-axis at all angles.  store
C       results in PH(L) array and in PHI(JMAX2,K,L).
        DO 20 L=1,LM
          LP=L+LM
          DO J=2,JSP
            I1=NUM+(J-1)
            I2=NUM-(J-2)
            PHI2(I1)=PHI(J,K,L)
            PHI2(I2)=PHI(J,K,LP)
          ENDDO
          IF(ICHK.EQ.1) PHI2(ISPECL)=PHI(NUM+2,K,L)
          IM=NPOINT-1
          DO I=1,IM
            P1=PHI2(I)
            XINV=1.0/X(I)
            XR=XPOINT*XINV
            IST=I+1
            DO J=IST,NPOINT
              XRATIO=X(J)*XINV
              P2=PHI2(J)
              PHI2(J)=(P1*(XRATIO-XR)+P2*(XR-1.0))/(XRATIO-1.0)
            ENDDO
          ENDDO
          PH(L)=PHI2(NPOINT)
          IF(ABS(PH(L)).GT.ABS(PHMAX))PHMAX=PH(L)
          PHI(JMAX2,K,L)=PH(L)
          PHI(JMAX2,K,LP)=PH(L)
   20   CONTINUE
C       At this k, find maximum deviation in ph(l)'s? Put result in phac(k).
        ERR=0.0
        DO L=1,LM
          ER=1.0-PH(L)/PHMAX
          IF(ABS(ER).GT.ABS(ERR))ERR=ER
        ENDDO
        PHAC(K)=ERR
   50 CONTINUE

C Find largest of all deviations in PHI's on z-axis and put value in PHICHK
C and k-value in KLOCAT.
      ERR=0.0
      DO K=2,KMAX1
        XX=ABS(PHAC(K))
        IF(XX.GT.ABS(ERR)) THEN
          ERR=PHAC(K)
          KLOCAT=K
        ENDIF
      ENDDO
      PHICHK=ERR
 
C Write the amplitude and phase to the results file.
      IF(IPRINT.EQ.1)THEN
        WRITE(3,10100) 'PHAC(1:KMAX1):'
        WRITE(3,10040)(PHAC(K),K=2,KMAX1)
        WRITE(3,10100) 'Fourier amplitudes and phases:'
10040   FORMAT(1P11E11.3)
10100   FORMAT(a)
        IF(LMAX.GT.8) THEN
          MR=1
          MP=8
C         LQ=LM/8
C         DO K=1,LQ
            WRITE(3,10110)(J,((COEF(J,2,M,I),M=MR,MP),I=1,2),J=2,JMAX)
            MR=MR+8
            MP=MP+8
C         ENDDO
        END IF
      END IF
10110 FORMAT(20X,I5,' COEFS=',1P8E12.3,/,25X,' PHASE=',0P8F12.2)
10112 FORMAT(20X,I5,' COEFS=',1P2E12.3,/,25X,' PHASE=',0P2F12.2)
10114 FORMAT(20X,I5,' COEFS=',1P4E12.3,/,25X,' PHASE=',0P4F12.2)

!      dbgl=dbgl-1

      RETURN
      END

