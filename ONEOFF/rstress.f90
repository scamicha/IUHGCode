! calculates the Reynold's stress
program rstress
 implicit none
!
! creating the map.  dr converts from cells to AU.
! dx is zof3n
!
 integer :: JMAX, KMAX, LMAX, YMAX, XMAX, JREQ
 integer :: ISTART, IEND, ISKIP,JS,LS,LSL,LSH=2,JSH=2
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K

 real(KIND=8), allocatable, dimension(:,:,:) :: rho,vr,vphi,vz,eps
 real(KIND=8), allocatable, dimension(:)::stress,ringmass,avr,avphi
 real(KIND=8), allocatable, dimension(:)::r,rhf

 real(KIND=8) :: dz,dr,torp,sconv,avgvphi,avgvr,mass,elost,den,sound,ommax
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time, bg
 real(KIND=8) :: tmassini,mstar,rdiskau,kconst,delt,mtot
 character :: filein*72, fileout*72,filenum*6,dum(13)*72

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

  if ( len(trim(dum(1))) == 0 .or. &
       len(trim(dum(2))) == 0 .or. &
       len(trim(dum(3))) == 0 .or. &
       len(trim(dum(4))) == 0 .or. &
       len(trim(dum(5))) == 0 .or. &
       len(trim(dum(6))) == 0 .or. &
       len(trim(dum(7))) == 0 .or. &
       len(trim(dum(8))) == 0 .or. &
       len(trim(dum(9))) == 0 .or. &
       len(trim(dum(10)))== 0 .or. &
       len(trim(dum(11)))== 0 .or. &
       len(trim(dum(12)))== 0 .or. &
       len(trim(dum(13)))== 0        ) then
    print *, "A CHYMTOOL SPONSORED EVENT"
    print *, " "
    print *, " "
    print *, "rstress v 1. 17.06.2011. A. C. Boley"
    print *, " "
    print *, " "
    print *, "TORQUE USAGE: rstress {stepbegin} {stepend} {stepskip}"
    print *, "                       {jmax} {kmax} {lmax} {torp (not used at the moment)} {dz}"
    print *, "                       {cell-to-size unit conversion  }"
    print *, "                       {constant torque (usually just set to 0.)}"
    print *, "                       {code to cgs torque conversion}"
    print *, "                       {filein prefix} {fileout prefix}" 
    print *, " "
    print *, " "
    print *, "Calling Cthulu...(run)"
    stop
  endif

  print *, "A CHYMTOOL SPONSORED EVENT"
  print *, " "
  print *, " "
  print *, "rstress v 1. 17.06.2011. A. C. Boley"
  print *, " "
  print *, " "
 
  dum(1) = trim(dum(1)) 
  dum(2) = trim(dum(2)) 
  dum(3) = trim(dum(3)) 
  dum(4) = trim(dum(4)) 
  dum(5) = trim(dum(5)) 
  dum(6) = trim(dum(6)) 
  dum(7) = trim(dum(7)) 
  dum(8) = trim(dum(8)) 
  dum(9) = trim(dum(9)) 
  dum(10) = trim(dum(10)) 
  dum(11) = trim(dum(11)) 

  read(dum(1),"(i8)")ISTART
  read(dum(2),"(i8)")IEND
  read(dum(3),"(i8)")ISKIP
  read(dum(4),"(i8)")JMAX
  read(dum(5),"(i8)")KMAX
  read(dum(6),"(i8)")LMAX
  read(dum(7),"(f15.8)")torp
  read(dum(8),"(f15.8)")dz
  read(dum(9),"(f15.8)")dr
  read(dum(10),"(f15.8)")bg
  read(dum(11),"(f15.8)")sconv

  print *, " TORQUE OUT: ISTART -> ", ISTART
  print *, " TORQUE OUT: IEND -> ", IEND
  print *, " TORQUE OUT: ISKIP -> ", ISKIP
  print *, " TORQUE OUT: JMAX -> ", JMAX
  print *, " TORQUE OUT: KMAX -> ", KMAX
  print *, " TORQUE OUT: LMAX -> ", LMAX
  print *, " TORQUE OUT: torp -> ", torp
  print *, " TORQUE OUT: dz -> ", dz
  print *, " TORQUE OUT: dr conversion -> ", dr
  print *, " TORQUE OUT: background -> ", bg
  print *, " TORQUE OUT: torque conversion -> ", sconv
  print *, " TORQUE OUT: filein prefix -> ", trim(dum(12))
  print *, " TORQUE OUT: fileout prefix-> ", trim(dum(13))

  allocate(rho(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(vr(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(vphi(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(vz(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(eps(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(stress(-1:JMAX))
  allocate(ringmass(-1:JMAX))
  allocate(avphi(-1:JMAX))
  allocate(avr(-1:JMAX))
  allocate(r(-1:JMAX))
  allocate(rhf(-1:JMAX))
  
  mstar = 1.d0
  rdiskau = 40.d0

 I = ISTART
 do while ( I <= IEND )
    write (filenum,'(I6.6)')I
    filein=trim(dum(12))//filenum//' ' 
    fileout=trim(dum(13))//filenum//' ' 
    print "(a,1x,a)", trim(filein), trim(fileout) 
    OPEN(UNIT=12, FILE=trim(filein),FORM='UNFORMATTED', STATUS="OLD")
    open(unit=13,file=trim(fileout))

 read(12)vr
 read(12)vz
 read(12)vphi
 read(12)rho
 read(12)eps
 READ(12)dr,dz,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
 read(12)tmassini
 close(12)

 DO j=-1, JMAX
    r(j)=DBLE(j)*ROF3N
    rhf(j) =(DBLE(j)*ROF3N)+(ROF3N/2.d0)
 ENDDO

 mtot   = mstar/(1.d0-tmassini)
 kconst = rdiskau*mtot**(0.333333333333333)/(r(JREQ)*7.93d-3)
 sconv = 2.24d48*mtot**(7.d0/3.d0)/kconst

 print *, " File read."

 
 stress=0d0;avphi=0d0;avr=0d0
 dphi=2d0*acos(-1d0)/dble(LMAX)
 do L = 0, LMAX-1
  do K = 0, KMAX-1
    do J = JSH, JMAX-JSH
      mass=0d0
      avgvr=0d0
      avgvphi=0d0
      do JS=J-JSH,J+JSH
         avgvr=avgvr+vr(JS,K,L)*(dble(J)+0.5d0)! vr is momentum density
         mass=mass+rho(JS,K,L)*(dble(J)+0.5d0)
      enddo
      avgvr=avgvr/mass
      mass=0d0
      do LS=L-LSH,L+LSH
        LSL=LS
        if(LS<0)LSL=LS+LMAX
        if(LS>LMAX-1)LSL=LS-LMAX
        avgvphi=avgvphi+vphi(J,K,LSL) ! vphi is angular momentum density
        mass=mass+rho(J,K,LSL)
      enddo
      avgvphi=avgvphi/( (dr*(dble(J)+0.5d0))  * mass )
      if(L==0.and.K==0)avphi(J)=avgvphi
      if(L==0.and.K==0)avr(J)=avgvr 
      stress(J)=stress(J)+ &
         rho(J,K,L)*(avgvphi-vphi(J,K,L)/(dr*(dble(J)+0.5d0)*rho(J,K,L)))&
                   *(avgvr-vr(J,K,L)/rho(J,K,L))*dz*dphi*2d0*(dr*(dble(J)+0.5d0))**2
      ringmass(J)=ringmass(J)+rho(J,K,L)*dphi*dr*dz*dr*(dble(J)+0.5d0)
    enddo
  enddo
 enddo

 do J=0,JMAX-1
  write(13,'(7(1pe15.8,1x))') dr*( dble(J)+0.5d0 ),stress(J)*sconv+bg,ringmass(J),avphi(J)-&
       vphi(J,0,0)/(dr*(dble(J)+0.5)*rho(J,0,0)),avr(J)-vr(J,0,0)/rho(J,0,0),rho(J,0,0),&
       rho(J,0,0)*(avphi(J)-vphi(J,0,0)/(dr*(dble(J)+0.5)*rho(J,0,0)))&
          *(avr(J)-vr(J,0,0)/rho(J,0,0))*dr*dr*dz*dphi*(dble(J)+0.5)
 enddo


 I = I+ISKIP
 enddo ! end while loop


 stop
end
  

    
