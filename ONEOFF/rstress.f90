! calculates the Reynold's stress
PROGRAM rstress
 IMPLICIT NONE
!
! creating the map.  dr converts from cells to AU.
! dx is zof3n
!
 INTEGER :: JMAX, KMAX, LMAX, YMAX, XMAX, JREQ
 INTEGER :: ISTART, IEND, ISKIP,JS,LS,LSL,LSH,JSH
 INTEGER :: II, JJ, IRR, IAN, IAN2, I, L, J, K
 INTEGER :: AVGRNUM,AVGPHINUM

 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: rho,vr,vphi,vz,eps
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:)::stress,ringmass,avr,avphi
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:)::r,rhf

 REAL(KIND=8) :: dz,dr,torp,sconv,avgvphi,avgvr,mass,elost,den,sound,ommax
 REAL(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time, bg
 REAL(KIND=8) :: tmassini,mstar,rdiskau,kconst,delt,mtot,rconv,vcent
 CHARACTER :: filein*72, fileout*72,filenum*6,dum(13)*72
 LOGICAL :: FILE_EXIST

  CALL GetArg(1,dum(1))
  CALL GetArg(2,dum(2))
  CALL GetArg(3,dum(3))
  CALL GetArg(4,dum(4))
  CALL GetArg(5,dum(5))
  CALL GetArg(6,dum(6))
  CALL GetArg(7,dum(7))
  CALL GetArg(8,dum(8))
  CALL GetArg(9,dum(9))
  CALL GetArg(10,dum(10))
  CALL GetArg(11,dum(11))
  CALL GetArg(12,dum(12))
  CALL GetArg(13,dum(13))

  IF ( LEN(TRIM(dum(1))) == 0 .OR. &
       LEN(TRIM(dum(2))) == 0 .OR. &
       LEN(TRIM(dum(3))) == 0 .OR. &
       LEN(TRIM(dum(4))) == 0 .OR. &
       LEN(TRIM(dum(5))) == 0 .OR. &
       LEN(TRIM(dum(6))) == 0 .OR. &
       LEN(TRIM(dum(7))) == 0 .OR. &
       LEN(TRIM(dum(8))) == 0 .OR. &
       LEN(TRIM(dum(9))) == 0 .OR. &
       LEN(TRIM(dum(10)))== 0 .OR. &
       LEN(TRIM(dum(11)))== 0 .OR. &
       LEN(TRIM(dum(12)))== 0 .OR. &
       LEN(TRIM(dum(13)))== 0        ) THEN
    PRINT *, "A CHYMTOOL SPONSORED EVENT"
    PRINT *, " "
    PRINT *, " "
    PRINT *, "rstress v 1. 17.06.2011. A. C. Boley"
    PRINT *, " "
    PRINT *, " "
    PRINT *, "TORQUE USAGE: rstress {stepbegin} {stepend} {stepskip}"
    PRINT *, "                       {jmax} {kmax} {lmax} {torp (not used at the moment)} {dz}"
    PRINT *, "                       {cell-to-size unit conversion  }"
    PRINT *, "                       {constant torque (usually just set to 0.)}"
    PRINT *, "                       {code to cgs torque conversion}"
    PRINT *, "                       {filein prefix} {fileout prefix}" 
    PRINT *, " "
    PRINT *, " "
    PRINT *, "Calling Cthulu...(run)"
    STOP
  ENDIF

  PRINT *, "A CHYMTOOL SPONSORED EVENT"
  PRINT *, " "
  PRINT *, " "
  PRINT *, "rstress v 1. 17.06.2011. A. C. Boley"
  PRINT *, " "
  PRINT *, " "
 
  dum(1) = TRIM(dum(1)) 
  dum(2) = TRIM(dum(2)) 
  dum(3) = TRIM(dum(3)) 
  dum(4) = TRIM(dum(4)) 
  dum(5) = TRIM(dum(5)) 
  dum(6) = TRIM(dum(6)) 
  dum(7) = TRIM(dum(7)) 
  dum(8) = TRIM(dum(8)) 
  dum(9) = TRIM(dum(9)) 
  dum(10) = TRIM(dum(10)) 
  dum(11) = TRIM(dum(11)) 

  READ(dum(1),"(i8)")ISTART
  READ(dum(2),"(i8)")IEND
  READ(dum(3),"(i8)")ISKIP
  READ(dum(4),"(i8)")JMAX
  READ(dum(5),"(i8)")KMAX
  READ(dum(6),"(i8)")LMAX
  READ(dum(7),"(f15.8)")torp
  READ(dum(8),"(f15.8)")dz
  READ(dum(9),"(f15.8)")dr
  READ(dum(10),"(f15.8)")bg
  READ(dum(11),"(f15.8)")sconv

  PRINT *, " TORQUE OUT: ISTART -> ", ISTART
  PRINT *, " TORQUE OUT: IEND -> ", IEND
  PRINT *, " TORQUE OUT: ISKIP -> ", ISKIP
  PRINT *, " TORQUE OUT: JMAX -> ", JMAX
  PRINT *, " TORQUE OUT: KMAX -> ", KMAX
  PRINT *, " TORQUE OUT: LMAX -> ", LMAX
  PRINT *, " TORQUE OUT: torp -> ", torp
  PRINT *, " TORQUE OUT: dz -> ", dz
  PRINT *, " TORQUE OUT: dr conversion -> ", dr
  PRINT *, " TORQUE OUT: background -> ", bg
  PRINT *, " TORQUE OUT: torque conversion -> ", sconv
  PRINT *, " TORQUE OUT: filein prefix -> ", TRIM(dum(12))
  PRINT *, " TORQUE OUT: fileout prefix-> ", TRIM(dum(13))

  ALLOCATE(rho(-1:JMAX,-1:KMAX,0:LMAX-1))
  ALLOCATE(vr(-1:JMAX,-1:KMAX,0:LMAX-1))
  ALLOCATE(vphi(-1:JMAX,-1:KMAX,0:LMAX-1))
  ALLOCATE(vz(-1:JMAX,-1:KMAX,0:LMAX-1))
  ALLOCATE(eps(-1:JMAX,-1:KMAX,0:LMAX-1))
  ALLOCATE(stress(-1:JMAX))
  ALLOCATE(ringmass(-1:JMAX))
  ALLOCATE(avphi(-1:JMAX))
  ALLOCATE(avr(-1:JMAX))
  ALLOCATE(r(-1:JMAX))
  ALLOCATE(rhf(-1:JMAX))
  
  LSH   = LMAX/32
  JSH   = 2
  AVGRNUM = 2*JSH+1
  AVGPHINUM = 2*LSH+1
  mstar = 1.d0
  rdiskau = 40.d0

 DO I = ISTART,IEND,ISKIP 
    WRITE (filenum,'(I6.6)')I
    filein=TRIM(dum(12))//filenum//' ' 
    fileout=TRIM(dum(13))//filenum//' ' 
    INQUIRE(file=TRIM(filein),exist=FILE_EXIST)
    IF(.NOT.FILE_EXIST) THEN
       PRINT*, "File ",TRIM(filein)," does not exist. Skipping."
       CYCLE
    ENDIF
    PRINT "(a,1x,a)", TRIM(filein), TRIM(fileout) 
    OPEN(UNIT=12, FILE=TRIM(filein),FORM='UNFORMATTED', STATUS="OLD")
    OPEN(unit=13,file=TRIM(fileout))

    READ(12)vr
    READ(12)vz
    READ(12)vphi
    READ(12)rho
    READ(12)eps
    READ(12)dr,dz,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
    READ(12)tmassini
    CLOSE(12)

    DO j=-1, JMAX
       r(j)=DBLE(j)*dr
       rhf(j) =(DBLE(j)*dr)+(dr/2.d0)
    ENDDO

    mtot   = mstar/(1.d0-tmassini)
    kconst = rdiskau*mtot**(0.333333333333333)/(r(JREQ)*7.93d-3)
    sconv = 2.24d48*mtot**(7.d0/3.d0)/kconst
    rconv = rdiskau/r(JREQ)

    PRINT *, " File read."

 
    stress=0d0;avphi=0d0;avr=0d0
    dphi=2d0*ACOS(-1d0)/DBLE(LMAX)
    DO L = 0, LMAX-1
       DO K = 0, KMAX-1
          DO J = JSH, JMAX-JSH-1
             avgvr=0d0
             avgvphi=0d0
             DO JS=J-JSH,J+JSH
                DO LS=L-LSH,L+LSH
                   LSL=LS
                   IF(LS<0)LSL=LS+LMAX
                   IF(LS>LMAX-1)LSL=LS-LMAX
                   vcent = (vr(JS,K,LSL)/(rho(JS-1,K,LSL)+rho(JS,K,LSL)))+  &
                        (vr(JS+1,K,LSL)/(rho(JS,K,LSL)+rho(JS+1,K,LSL)))
                   avgvr=avgvr+vcent
                   
                   avgvphi=avgvphi+vphi(JS,K,LSL)/(rho(JS,K,LSL)*dr*(DBLE(J)+0.5d0))
                ENDDO
             ENDDO
             avgvr=avgvr/(AVGRNUM*AVGPHINUM)
             avgvphi = avgvphi/(AVGRNUM*AVGPHINUM)

             IF(L==0.AND.K==0)avphi(J)=avgvphi
             IF(L==0.AND.K==0)avr(J)=avgvr 

             vcent=(vr(J,K,L)/(rho(J-1,K,L)+rho(J,K,L)))+  &
                        (vr(J+1,K,L)/(rho(J,K,L)+rho(J+1,K,L)))
             stress(J)=stress(J)+ &
                  rho(J,K,L)*(avgvphi-vphi(J,K,L)/(dr*(DBLE(J)+0.5d0)*rho(J,K,L)))&
                  *(avgvr-vcent)*2d0

             ringmass(J)=ringmass(J)+rho(J,K,L)*dphi*dr*dz*dr*(DBLE(J)+0.5d0)
          ENDDO
       ENDDO
    ENDDO

    DO J=0,JMAX-1
       WRITE(13,'(8(1pe15.8,1x))') time/torp,rconv*dr*( DBLE(J)+0.5d0 ),stress(J)*sconv+bg,ringmass(J),avphi(J)-&
            vphi(J,0,0)/(dr*(DBLE(J)+0.5)*rho(J,0,0)),avr(J)-vr(J,0,0)/rho(J,0,0),rho(J,0,0),&
            rho(J,0,0)*(avphi(J)-vphi(J,0,0)/(dr*(DBLE(J)+0.5)*rho(J,0,0)))&
            *(avr(J)-vr(J,0,0)/rho(J,0,0))*dr*dr*dz*dphi*(DBLE(J)+0.5)
    ENDDO


 ENDDO ! end while loop


 STOP
END    
