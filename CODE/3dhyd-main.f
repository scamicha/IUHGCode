      program IU3DHYDRO
C************************************************************************
C************************************************************************
C...THIS IS A NEW VERSION OF THE HYDRO CODE WHICH IS SECOND-ORDER IN BOTH
C...SPACE AND TIME.  S, T ARE FACE-CENTERED; A,RHO,& EPS ARE CELL-
C...CENTERED.   AN ENERGY EQUATION IS INCLUDED FOR ENERGY TRANSPORT.
C...The grid can start at j!=2, as long as the potential of the star is
C...read from the cthuphi_axi.dat file. acm2000.
!
C...Also, opacities added such that cooling can be calculated. Need to 
C...define real dimensions of the disk in the source routine.  acm 2002.
C...ARTIFICAL VISCOSITY IS USED TO TREAT SHOCKS AND HEAT THE MODEL.
!
! Substantial changes have been made to the hydrodynamics code.  These
! changes include fluxing eps directly, and not eps^(1/gamma).  The 
! internal energy of H2 is directly calculated at the beginning of the
! simulations, and the temperature is determined based on the eps, rho
! and the assumed ortho:para hydrogen statistics. The pressure is 
! calculated by p~rho T.  The mean molecular weight for the gas
! is calculated based on what is used for calculating the specific 
! internal energy.  This does create a slight difference when compared
! with D'Alessio's opacities, but it should be small.  In addition,
! a new radiative transfer routine is used, which uses vertical rays.
! A. C. Boley
!
!
C...THE FOLLOWING SUBROUTINES ARE CALLED:
C
c
C     *SETUP   :  Read starting models, impose perts, etc. (io.f)
!     *initengtable : Creates specific energy table for calculating temp. (initengtable.f) ACB
!     *TempFind:  Finds temperature for each cell. (source.f) ACB
!     *TempFindSpec: Finds temperature based on input eps, rho. (source.f) ACB
!     *VLIMIT  :  Limit maximum velocity (usually < 2*SOUND). (housekeeping.f)
!     *DELTA   :  Calculate maximum safe delta time. (housekeeping.f)
!     *RITE    :  Write output information. (io.f)
!     *SLOPE   :  Calculate van Leer Slope. (flux.f)
!     *VLI     :  Calculate van Leer Monotonic Interpolation. (flux.f)
!     *FLUX    :  Advect S,T,A,RHO, and EPS. (flux.f)
!     *CLEANUP :  Fix velocities, densities and energy on boundaries.
!                 (housekeeping.f)
!     *SOURCE   :  Source S,T,A, and EPS. (source.f)
!     *VELOCITY :  From momentum densities, calculate velocities.
!     *CENTMASS :  Calculate Center of Mass. (housekeeping.f)
!     *STATE    :  Equation of State. (state.f)
!     *RAD     : Old routine for "radiative physics." (rad.f) 
!     *REALTR  :  -|Together These Perform a  (fft.f)
!     *FFT     :  -|Fast Fourier Transform.   (fft.f)
!     *POT3    :  Potential Solver. (pot3.f)
!     *ZAXPHI  :  -|                (pot3.f)
!     *BLKTRI  :  -|                (pot3.f)
!     *BLKTR1  :  -|Various Functions used in/with the (pot3.f)
!     *PRDCT   :  -|Potential Solvers. (pot3.f)
!     *COMBP   :  -|                   (pot3.f)
!     *TQLRAT  :  -|                   (pot3.f)
!     *SETBDY  :  Initialization before BDYGEN. (boundary.f)
!     *BDYGEN  :  Boundary Potential Solver.    (boundary.f)
!     *SORT    :  Sort. (housekeeping.f)
!     *CLEARUP :  Clearup after SETMODE or DAMP call. (housekeeping.f)
!     *TORQUE  :  Calculate instaneous torques.
!     *dumpphi :  dump the potential grid for comparison (io.f)
!     *cyl3d   :  dump density grid for use in gen. 3d images (io.f)
!     *lowresrandom: low resolution random pert. (io.f)
!     *hybrid  : Vertical ray radiative physics routine. (hybrid.f) ACB
!     *initengtable: initialize energy table. (initengtable.f) ACB
!     *tempfind: find temp based on energy table. (source.f) ACB
!     *tempfindspec: find the temp of just one cell based on the energy
!                    table. (source.f) ACB
!***********************************************************************
!  
      IMPLICIT real*8 (a-h,o-z)
      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'
      INCLUDE 'units.h'

      COMMON /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
      COMMON /TIMEST/INDX,ISOADI,ALLOW,DMAX,CHGMAX
      COMMON /COEFS/COEF(POT3JMAX2,POT3KMAX2,LMAX2,2)
      COMMON /LIMIT/SOUND
      COMMON /ITS/ITSTRT,ITSTOP,ITSTEP,ISTOR

      real*8 epsjr,rhojr,ommax,mirp
      common /misc/ epsjr,rhojr,ommax

C The following arrays are local to the main program, and should not need 
C to be placed in common. However, some systems may try to put them on the
C stack, which may result in stack overflow. Putting them in a common block
C guarantees they will go on the heap. ! COMMON BLOCK REMOVED. ACB
      REAL*8 SS(JMAX2,KMAX2,LMAX),
     &       TT(JMAX2,KMAX2,LMAX),
     &       AA(JMAX2,KMAX2,LMAX),
     &       RRHO(JMAX2,KMAX2,LMAX),
     &       EEPS(JMAX2,KMAX2,LMAX),
     &       Tauross(jmax2,kmax2,lmax)

      CHARACTER  np*2,tw*3
      CHARACTER  tim*8,rhofile*80,tempfile*80,starfile*80,
     &           restart_star*80
      CHARACTER  coeffile*80,index*8,modefile*80,centfile*80
      CHARACTER  tcoeffile*80,tcoeftenfile*80,massfile*80,epsfile*80
      CHARACTER  resfile*80,savedfile*80,epsfull*80,planetfile*80
      CHARACTER  planeten*80
      real*8     massc(jmax2), massf(jmax2),limmf(10),masst,mdot
      real*8     sumeven(8),sumodd(8),theta(lmax),tcool,theat,tflux,tirr
      real*8     volume(jmax2),totcoef(8)
      REAL*8     mass1,moneav,totcoefcym(10,8)
      DATA ISTRDN,CHANGD/25,2.00/,NCONS,DELCON/10,2.0/
      integer    OMP_GET_MAX_THREADS
      real*8     dwalltime00, dummy, cooltime
      real*8     iphix,iphiy,x,y,disp,alf  ! these are tallies for the indirect potential
      REAL*8     ipot(jmax2,kmax2,lmax)
      REAL*8     planetx,planety,planetvx,planetvy,planetax,planetay  !for planet integration
      REAL*8 planetr,planetphi
      REAL*8 orbitenergy,disken,r_hill,smooth
      REAL*8 initial_x,initial_y,initial_vx,initial_vy,initial_r    
      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),absfrac
      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum
      COMMON /PLANETINIT/ initial_x,initial_y,initial_vx,initial_vy
      COMMON /PLANETWRITE/ planetx,planety,planetvx,planetvy

!$    write(6,*) 'MP enabled.  Nproc = ', OMP_GET_MAX_THREADS()
!
! FIRST TOUCH: INITIALIZE ALL ARRAYS.  THIS IS NOT ONLY GOOD PRACTICE, BUT
! IT CAN SIGNICANTLY SPEED UP THE OPENMP CODE.  THE ARRAYS ARE MORE 
! EFFICIENTLY DISTRIBUTED AMONG PROCESSORS IF TOUCHED ALL AT ONCE.

      call first_touch()  

!$OMP PARALLEL DO SCHEDULE(STATIC)
      DO L=1,LMAX
         DO K=1,KMAX2
            DO J=1,JMAX2
               SS(J,K,L)=0.0
               TT(J,K,L)=0.0
               AA(J,K,L)=0.0
               RRHO(J,K,L)=0.0
               EEPS(J,K,L)=0.0
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO    


C....START EVOLUTION

      TEMPCD=100.0*CHANGD
      MSTORE=TEMPCD
      ITSTEP=0
      MAXTRM=10
      ISYM=2
      ENON=1.d0
      CONSTP=0
      BDYTEM=1.E-3
      CORMAS=1.E-2
      
!bkp..cs=0 for no shear viscosity (not implemented yet, so don't diddle!)
!rpl..xxxtodo: make this parameter a read-in.      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ DATA AND SET UP GRID.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL SETUP(ITSTRT,ITSTOP,IDIAG,ISOADI,ISTOR,ITSTEP,
     &           ISYM,MAXTRM)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OPEN RECORDING FILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(index,'(i8.8)')ITSTRT+ISTOR-MOD(ITSTRT,ISTOR)
      coeffile=trim(outpath)//'coefs.'//index
      modefile=trim(outpath)//'modes.'//index
      centfile=trim(outpath)//'c_o_m.'//index
      tcoeffile=trim(outpath)//'tcoef.'//index
      tcoeftenfile=trim(outpath)//'tctot.'//index
      massfile=trim(outpath)//'massflow.'//index
      epsfile=trim(outpath)//'coolheat.'//index
      planetfile=trim(outpath)//'planet_pos.'//index
      planeten  =trim(outpath)//'planet_en.'//index
      
      OPEN(unit=9,file=coeffile)
      OPEN(unit=19,file=centfile)
      open(unit=13,file=tcoeffile)
      open(unit=15,file=tcoeftenfile)
      open(unit=17,file=massfile,form='formatted')
      open(unit=21,file=epsfile,form='formatted')

      SELECT CASE (star_state)
      CASE("wiggle")
         starfile='starrysim.'//index//' '
         restart_star='starry_restart.'//index//' '
         open(unit=987,file=trim(starfile),form='FORMATTED')
         open(unit=988,file=trim(restart_star),form='FORMATTED')
      CASE("wiggle_restart")
         starfile='starrysim.'//index//' '
         restart_star='starry_restart.'//index//' '
         open(unit=987,file=trim(starfile),form='FORMATTED')
         open(unit=988,file=trim(restart_star),form='FORMATTED')
      CASE DEFAULT
         print*,'Fixed star or indirect potential run.
     & No star ouput files'
      END SELECT

      CALL RITE(3,-1,  1,1,1,  1,1,1,  1,1,1)
!$OMP PARALLEL DEFAULT(SHARED)
      CALL CLEARUP
!$OMP END PARALLEL
      CALL RAD(ISOADI)

      DMAX=0.9*DEN
      CHGMAX=0.0001
      MIRP=2.d0*pi/ommax

      if (fluid_elements) then
        call Fluid_Setup()
      endif

      SELECT CASE (star_state)
            
      CASE("indirect")
         if (planet) THEN
            open(unit=22,file=planetfile,form='formatted')
            open(unit=24,file=planeten,form='formatted')
	    initial_r = sqrt(initial_x**2 + initial_y**2)
            r_hill = (r(jreq)*initial_r/Rdiskau)*
     &           ((mass_ratio/3.d0)**(1/3))
            smooth = (0.2d0*r_hill)**2
!$OMP PARALLEL DO SCHEDULE(STATIC) 
            do l=1,lmax
               do j=1,jmax2
                  do k=1,kmax2
                    phi(j,k,l) = phi(j,k,l)-starphi(j,k,l)-tinphi(j,k,l)
                  enddo
               enddo
            enddo
!$OMP END PARALLEL DO
            CALL planet_initialize(planetx,planety,planetvx,
     &           planetvy,planetax,planetay)
            planetr = sqrt(planetx**2+planety**2)
            planetphi = atan2(planety,planetx)
            if(planetphi .lt. 0.d0) planetphi = planetphi + 2.d0*pi

!$OMP PARALLEL DO SCHEDULE(STATIC) 
!$OMP&FIRSTPRIVATE(smooth,planetr,planetphi)
            do l=1,lmax
               do j=1,jmax2
                  do k=1,kmax2
                     phi(j,k,l) = phi(j,k,l)+((mass_planet*rhf(j)*
     &                    planetr*cos(((float(l)-0.5d0)*dtheta)-
     &                    planetphi))/planetr**3)-(mass_planet/(sqrt(
     &                    rhf(j)**2+planetr**2-2.d0*rhf(j)*planetr*
     &                    cos(((float(l)-0.5d0)*dtheta)-planetphi)+
     &                    zhf(k)**2+smooth)))+starphi(j,k,l)+
     &                    tinphi(j,k,l)
                  enddo
               enddo
            enddo
!$OMP END PARALLEL DO
            
!            write(22,57)itstep,time,planetx,planety,planetvx,planetvy,
!     &           planetax,planetay

         ENDIF
      END SELECT
 57   format(i8,1x,7(1pe20.12,1x))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------begin main loop-----------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(6,*) 'Beginning simulation'
      DO 90 ITSTEP=ITSTRT,ITSTOP
         INDX=ITSTEP-ITSTRT
         if(mod(ITSTEP,1).eq.0) write(6,10000) itstep,time,delt        !dkb
10000    format('step=',i8,'   time=',1pe26.13,'   delt=',1pe26.15)  !dkb

         IF(MOD(ITSTEP,IDIAG).EQ.0) then
            IHEAD=-1
            IPRINT=1
            CALL RITE(2,-1,ITSTOP,ITSTRT,ITSTOP,ISTOR,1,1,1,1,1)
         else 
            IPRINT=0
            IHEAD=0
         end if

C... Calculation of accreted mass and boundary outflow 
C... tmassini = initial disk mass 
C... tmass    = current disk mass
C... tmassadd = added mass by resetting bgrnd to rholmt
C... tmassacc = accreted mass
C... tmassout = mass flowing out of the top, bottom and outer boundary

         mdot=0.d0      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J,K,L)
!$OMP& SHARED(zof3n,delt,tmassout,dtheta,mdot)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:mdot)
         do j=jmin,jmax
            do l=1,lmax 
               mdot=mdot+T(j,kmax1,l)*dtheta*(R(j+1)**2-R(j)**2)
            enddo
         enddo
!$OMP END DO
         if(mdot.gt.0.0) then
!$OMP MASTER
           tmassout=tmassout+(mdot*DELT)
!$OMP END MASTER
         endif
!$OMP MASTER
         mdot=0.d0         
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:mdot)
         do k=2,kmax
            do l=1,lmax
               mdot=mdot+S(jmax1,k,l)*2.0*dtheta*R(jmax1)*(zof3n)
            enddo
         enddo
!$OMP END DO
         if(mdot.gt.0.0) then
!$OMP MASTER
           tmassout=tmassout+(mdot*DELT)
!$OMP END MASTER
         endif
!$OMP MASTER
         mdot=0.d0  
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:mdot)       
         do k=2,kmax
            do l=1,lmax
               mdot=mdot+S(jmin,k,l)*2.0*dtheta*R(jmin)*(zof3n)
            enddo  
         enddo
!$OMP END DO
!$OMP END PARALLEL
         tmassacc=tmassacc+(abs(mdot)*delt)


         write(17,20)itstep,time,tmass,tmassadd,tmassout,tmassacc,mdot
 20      format(i8,1x,6(1pe24.16,1x))

         CALL DELTA(MOD(ITSTEP,IDIAG).EQ.0)

         TIME=TIME+DELT

      if(fluid_elements) then 
         call fluidtrace(itstep,itstop)
      endif

!.............................................!
!.....START SECOND ORDER TIME INTEGRATION.....!
!.............................................!

C..(1) 1/2 SOURCE S, T, A, RHO, EPS.


         DELT=DELT*0.5

         CALL SOURCE

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,l)                           &
!$OMP&  SHARED(delt)
         CALL VELOCITY

         CALL VLIMIT

C..(2) STORE QUANTITIES.

!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
           DO K=1,KMAX2
             DO J=1,JMAX2
                SS(J,K,L)=S(J,K,L)
                TT(J,K,L)=T(J,K,L)
                AA(J,K,L)=A(J,K,L)
                RRHO(J,K,L)=RHO(J,K,L)
                EEPS(J,K,L)=EPS(J,K,L)
             ENDDO
           ENDDO
         ENDDO
C$OMP END DO
                      
C..(3) 1/2 FLUX S, T, A, RHO, EPS.

         CALL FLUX(S,T,A,RHO,EPS)

         CALL VELOCITY

         CALL VLIMIT

C..(4) 1 FLUX SS, TT, AA, RRHO, EEPS, 
C........BUT USE S, T, A, RHO, EEPS TO CALCULATE FLUXES.

!$OMP MASTER
         DELT=2.0*DELT
!$OMP END MASTER

         CALL FLUX(SS,TT,AA,RRHO,EEPS)

C..(5) UPDATE S, T, A, RHO, EPS.

!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
           DO K=1,KMAX2
             DO J=1,JMAX2
               S(J,K,L)=SS(J,K,L)
               T(J,K,L)=TT(J,K,L)
               A(J,K,L)=AA(J,K,L)
               RHO(J,K,L)=RRHO(J,K,L)
               EPS(J,K,L)=EEPS(J,K,L)
             ENDDO
           ENDDO
         ENDDO
!$OMP END DO

C..(6) 1/2 SOURCE S, T, A, RHO, EPS.

         CALL VELOCITY

         CALL VLIMIT

         CALL CLEARUP

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX
       do K = 1, KMAX2
        do J = 1, JMAX2
         phi(J,K,L) = 0.
        enddo
       enddo
      enddo
!$OMP END DO 
!$OMP END PARALLEL

      CALL RAD(ISOADI)

      REDGE=R(JMAX1)
      CALL BDYGEN(MAXTRM,ISYM,REDGE)

      CALL POT3(8,IPRINT)
      

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j,l)                           &
!$OMP&  SHARED(delt,tmassacc,zof3n,dtheta,tmass,mass_star)
        if (jmin.gt.2) then
!        recalculate the total mass, acm2000
!$OMP MASTER
          tmass=0.d0
          iphix=0.d0    ! reset for indirect potential
          iphiy=0.d0
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:tmass)
            do l=1,lmax
               do j=jmin,jmax
                  do k=2,kmax
                     tmass=tmass+rho(j,k,l)*0.5*dtheta*zof3n
     &                    *(r(j+1)**2-r(j)**2)
                  enddo
               enddo
            enddo
!$OMP END DO
!$OMP MASTER
            tmass=tmass*2.d0
!$OMP END MASTER

         SELECT CASE (star_state)

         CASE("fixed")
!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
               do j=1,jmax2
                  do k=1,kmax2
                     if(tmassacc.gt.0.0) then
!                    add potential of accreted mass 
                        tinphi(j,k,l) = -(tmassacc)
     &                       /(sqrt(RHF(j)*RHF(j)+zHF(k)*zHF(k)))
                     else
                        tinphi(j,k,l)=0.0
                     endif
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l)+ ! ELIMINATED diskphi. ACB.
     &                    tinphi(j,k,l)
                  enddo
               enddo
            enddo
!$OMP END DO

         CASE("indirect")
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:iphix,iphiy)
            do l=1,lmax
               do j=jmin,jmax
                  do k=2,kmax
                     iphix = iphix + (rho(j,k,l)*dtheta*zof3n*rhf(j)*
     &                    (r(j+1)**2-r(j)**2)*cos((float(l)-0.5d0)
     &                    *dtheta))/(sqrt(rhf(j)**2+zhf(k)**2))**3
                     iphiy = iphiy + (rho(j,k,l)*dtheta*zof3n*rhf(j)*
     &                    (r(j+1)**2-r(j)**2)*sin((float(l)-0.5d0)
     &                    *dtheta))/(sqrt(rhf(j)**2+zhf(k)**2))**3                     
                  enddo
               enddo
            enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
               do j=1,jmax2
                  do k=1,kmax2
                     if(tmassacc.gt.0.0) then
!                    add potential of accreted mass 
                        tinphi(j,k,l) = -(tmassacc)
     &                       /(sqrt(RHF(j)*RHF(j)+zHF(k)*zHF(k)))
                     else
                        tinphi(j,k,l)=0.0
                     endif
!!!!!!!!!!!!!!!! For the indirect potential !!!!!!!!!!!!!!!!!
                     ipot(j,k,l) = rhf(j)*cos((float(l)-0.5d0)*dtheta)
     &                    *iphix+rhf(j)*sin((float(l)-0.5d0)*dtheta)
     &                    *iphiy
                     phi(j,k,l)= phi(j,k,l)+ipot(j,k,l)      ! ELIMINATED diskphi. ACB.
                  enddo
               enddo
            enddo
!$OMP END DO

      IF(planet)THEN
!$OMP SINGLE                 
         CALL planet_integrate(planetx,planety,planetvx,planetvy
     &        ,planetax,planetay)
         planetr = sqrt(planetx**2+planety**2)
         planetphi = atan2(planety,planetx)
         if(planetphi .lt. 0.d0) planetphi = planetphi + 2.d0*pi

         IF(MOD(ITSTEP,IDIAG).eq.0)THEN
            disken = 0.d0
            CALL potential_energy(planetx,planety,disken)
         orbitenergy = 0.5d0*(((mass_star+tmassacc)*mass_planet)/
     &        (mass_star+mass_planet+tmassacc))*(planetvx**2+
     &        planetvy**2)-((mass_star+tmassacc)*(((mass_star+tmassacc)*
     &        mass_planet)/(mass_star+mass_planet+tmassacc))/(sqrt(
     &        planetx**2+planety**2)))-disken
         write(24,58)itstep,time,orbitenergy
 58      FORMAT(i8,1x,2(1pe20.12,1x))
         
         ENDIF
         r_hill = (r(jreq)*planetr/Rdiskau)*
     &        ((mass_ratio/3.d0)**(1/3))
         smooth = (0.2d0*r_hill)**2
         write(22,57)itstep,time,planetx,planety,planetvx,planetvy,
     &        planetax,planetay

!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC) 
!$OMP&FIRSTPRIVATE(smooth,planetr,planetphi) 
         do l=1,lmax
            do j=1,jmax2
               do k=1,kmax2
                  phi(j,k,l) = phi(j,k,l)+((mass_planet*rhf(j)*
     &                 planetr*
     &                 cos(((float(l)-0.5d0)*dtheta)-planetphi))/
     &                 planetr**3)-(mass_planet/(sqrt(rhf(j)**2+
     &                 planetr**2-2.d0*rhf(j)*planetr*
     &                 cos(((float(l)-0.5d0)*dtheta)-planetphi)+
     &                 zhf(k)**2+smooth)))+starphi(j,k,l)+
     &                 tinphi(j,k,l)
                  enddo
               enddo
            enddo
!$OMP END DO
      ELSE
!$OMP DO SCHEDULE(STATIC) 
         do l=1,lmax
            do j=1,jmax2
               do k=1,kmax2
                  phi(j,k,l) = phi(j,k,l)+starphi(j,k,l)+tinphi(j,k,l)
                  enddo
               enddo
            enddo
!$OMP END DO
      ENDIF
      
         CASE("wiggle")

            call wiggle(rho,rhf,zhf,vx,vy,fx,fy,delt,pi,rpstar,phi_star,
     &                  rof3n,zof3n,dtheta,JMAX,KMAX,LMAX,JMIN,.false.,
     &                   .false.)
            if(mod(INDX,10)==0)then
               write(987,'(9(1X,1pe15.8))') time,vx,vy,fx,fy,rpstar,
     &         phi_star,rpstar*cos(phi_star),rpstar*sin(phi_star)
            endif

!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
               do j=1,jmax2
                  do k=1,kmax2
                     starphi(j,k,l) = -(mass_star+tmassacc)
     &                       /(sqrt(rhf(J)*rhf(J)+rpstar*rpstar-2.
     & *rpstar*rhf(J)*cos((dble(L)-.5)*dtheta-phi_star)
     & +zhf(K)*zhf(K)))
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l) ! ELIMINATED diskphi. ACB.
                  enddo
               enddo
            enddo
!$OMP END DO

         CASE("wiggle_restart")

            call wiggle(rho,rhf,zhf,vx,vy,fx,fy,delt,pi,rpstar,phi_star,
     &                  rof3n,zof3n,dtheta,JMAX,KMAX,LMAX,JMIN,.false.,
     &                   .true.)
            if(mod(INDX,10)==0)then
               write(987,'(9(1X,1pe15.8))') time,vx,vy,fx,fy,rpstar,
     &         phi_star,rpstar*cos(phi_star),rpstar*sin(phi_star)
            endif

!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
               do j=1,jmax2
                  do k=1,kmax2
                     starphi(j,k,l) = -(mass_star+tmassacc)
     &                       /(sqrt(rhf(J)*rhf(J)+rpstar*rpstar-2.
     & *rpstar*rhf(J)*cos((dble(L)-.5)*dtheta-phi_star)
     & +zhf(K)*zhf(K)))
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l) ! ELIMINATED diskphi. ACB.
                  enddo
               enddo
            enddo
!$OMP END DO

         END SELECT

      endif

!$OMP MASTER
         DELT=DELT*0.5
!$OMP END MASTER

         CALL STATE
!$OMP END PARALLEL

         CALL SOURCE

C..(7) UPDATES DUE TO ENERGY EQUATION

!$OMP PARALLEL DEFAULT(SHARED)
         CALL STATE

         CALL VELOCITY

         CALL VLIMIT

         CALL CLEARUP

!$OMP END PARALLEL
         DELT=2.*DELT


C.............................................C
C.....COMPLETED ONE TIME STEP INTEGRATION.....C
C.............................................C




C...Written Output for result file...
c....Slices in j,k,l.................

         IF(MOD(ITSTEP,IDIAG).eq.0) then
            WRITE(3,10200) ITSTEP,TIME                                    !dkb
10200    FORMAT(//,' --------------------------------------  ',           !dkb
     &         'STEP=',i8,'   TIME= ',1pe15.8,                            !dkb
     &         '   -------------------------------------------------')    !dkb
            LHALF=LMAX/2
            CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  1,1,1)
c           CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  2,2,1)
c           CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  3,3,1)
c           CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  4,4,1)
c           CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  5,5,1)
c           CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  6,6,1)
c           CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  7,7,1)
c           CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  8,8,1)
            CALL RITE(1,-1,  2,2,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1,-1,  3,3,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1,-1,  6,6,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1,-1,  4,4,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1, 1,  8,8,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(1, 1,  16,16,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1, 1,  24,24,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(1, 1,  30,30,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(1, 1,  60,60,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(1, 1,  90,90,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(1, 1,  128,128,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(1, 1,  129,129,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1, 1,  180,180,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1, 1,  200,200,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1, 1,  300,300,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1, 1,  400,400,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1, 1,  500,500,1,  2,KMAX1,1,  1,1,1)
c           CALL RITE(1,-1,  5,5,1,  2,2,1,  1,LHALF,1)
            CALL RITE(1, 1,  2,2,1,  2,2,1,  1,LHALF,1)
c           CALL RITE(1, 1,  3,3,1,  2,2,1,  1,LHALF,1)
c           CALL RITE(1, 1,  4,4,1,  2,2,1,  1,LHALF,1)
c           CALL RITE(1,-1,  10,10,1,  2,2,1,  1,LHALF,1) 
C           CALL RITE(1,-1,  15,15,1,  2,2,1,  1,LHALF,1)
c           CALL RITE(1,-1,  20,20,1,  2,2,1,  1,LHALF,1)
c           CALL RITE(1,-1,  25,25,1,  2,2,1,  1,LHALF,1)
            CALL RITE(1,-1,  30,30,1,  2,2,1,  1,LHALF,1)
            CALL RITE(10,1,  1,1,1,  1,1,1,  1,1,1)
         END IF

         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENERATE OUTPUT FILES (RHO, TEMP), AND WRITE TO RECORD FILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (jmin.eq.2) then
            timc=int(10.d0*time/MIRP)
         else 
            timc=int(5.d0*time*omega(jmin+5,2,1)/pi)
         endif

         IF(MOD(itstep,IDIAG).EQ.0) THEN
            write(tim,'(i8.8)')itstep

            rhofile=trim(outpath)//'rho3d.'//tim
            OPEN(UNIT=14,FILE=rhofile,FORM='UNFORMATTED')
            write(14)rho
            write(14)time
            CLOSE(14)  

            tempfile=trim(outpath)//'temperat3d.'//tim
            OPEN(UNIT=23,FILE=tempfile,FORM='UNFORMATTED')
            write(23)tempK
            write(23)time
            CLOSE(23)  

C	    call qcalc(tim)

            write(21,100)itstep,time,totcool,totdflux,totheat,totirr,
     &                   etotfl,eflufftot

            jstart=jmin
            do j=jstart,jmax
               if (lambda(j,2,1).ne.0.0) then
                  tcool=eps(j,2,1)/(abs(lambda(j,2,1))*torp)
               else
                  tcool=0.0
               endif
               if (divflux(j,2,1).ne.0.0) then
                  tflux=eps(j,2,1)/(abs(divflux(j,2,1))*torp)
               else
                  tflux=0.0
               endif
               if (hgamma(j,2,1).ne.0.0) then
                  theat=eps(j,2,1)/(hgamma(j,2,1)*torp)
               else
                  theat=0.0
               endif
               if (igamma(j,2,1).ne.0.0) then
                  tirr=eps(j,2,1)/(igamma(j,2,1)*torp)
               else
                  tirr=0.0
               endif

               write(21,101)j,eps(j,2,1),lambda(j,2,1),divflux(j,2,1)
     &              ,hgamma(j,2,1),igamma(j,2,1),tcool,tflux,theat,tirr
     &              ,tau(j,2,1,1),TempK(j,2,1),TeffK(j,1),TphK(j,1)
            enddo

            write(21,105)itstep,time,totcool,totdflux,totheat,totirr

            do i=1,3
               if (i.eq.1) j=30
               if (i.eq.2) j=100
               if (i.eq.3) j=200
               do k=2,kmax
                  if (lambda(j,k,1).ne.0.0) then
                     tcool=eps(j,k,1)/(abs(lambda(j,k,1))*torp)
                  else
                     tcool=0.0
                  endif
                  if (divflux(j,k,1).ne.0.0) then
                     tflux=eps(j,k,1)/(abs(divflux(j,k,1))*torp)
                  else
                     tflux=0.0
                  endif
                  if (hgamma(j,k,1).ne.0.0) then
                     theat=eps(j,k,1)/(hgamma(j,k,1)*torp)
                  else
                     theat=0.0
                  endif
                  if (igamma(j,k,1).ne.0.0) then
                     tirr=eps(j,k,1)/(igamma(j,k,1)*torp)
                  else
                     tirr=0.0
                  endif
                  write(21,106)j,k,eps(j,k,1),lambda(j,k,1),divflux(j
     &                 ,k,1),hgamma(j,k,1),igamma(j,k,1),tcool,tflux
     &                 ,theat,tirr,tau(j,k,1,1),TempK(j,k,1)
               enddo
            enddo

         END IF

 100     format(///,'STEP',i9,' TIME',1pe24.16,' TOTAL COOLING',1pe24.16
     &        ,' TOTAL DIVFLUX',1pe24.16,' TOTAL HEATING',1pe24.16,
     &        ' TOTAL IRRADIATION',1pe24.16,' TOTAL FLUXLOSS',1pe24.16
     &        ,' TOTAL EFLUFF LOSS',1pe24.16,
     &          //,3X,'J',6X,'EPS',9X,'LAMBDA',6X,'DIVFLUX',6X
     &        ,'HGAMMA',7X,'IGAMMA',5X,'Tcool(orp)',1X,
     &        'Tdflx(orp)',1X,'Theat(orp)',1X,'Tirr(orp)',5X,'TAU',6X,
     &        'Tmid(K)',4X,'Teff(K)',5X,'Tph(K)'
C    &        ,2X,'SDen(g/cm2)'
     &        )


 105     format(/,'STEP',i9,' TIME',1pe24.16,' TOTAL COOLING',1pe24.16
     &        ,' TOTAL DIVFLUX',1pe24.16,' TOTAL HEATING',1pe24.16
     &        ,' TOTAL IRRADIATION',1pe24.16,//,3X,'J',4X,'K',6X,'EPS'
     &        ,9X,'LAMBDA',6X,'DIVFLUX',6X,'HGAMMA',7X,'IGAMMA',5X
     &        ,'Tcool(orp)',1X,'Tdflx(orp)',1X,'Theat(orp)',2X
     &        ,'Tirr(orp)',4X,'TAU',6X,'Temp(K)')
                  
 101     format(i4,1x,5(1pe12.5,1x),8(1pe10.2,1x))
 106     format(i4,1x,i4,1x,5(1pe12.5,1x),6(1pe10.2,1x))

C...Store Coefficient Information....................
C...Use to extract modes and rho in Equatorial Plane.
         
 1040    format(i3,1x,i3,2(1x,1pe13.4))
         

         IF(MOD(ITSTEP,IDIAG).EQ.0) THEN

c...Equatorial fourier mode information (M=1-LMAX/2)

            CALL RITE(0,IHEAD,ITSTEP,INDX,1,1,1,1,ISTRDN,MSTORE,MS)
            LMWR9 = LMAX/2      
            WRITE(9,103)ITSTEP,TIME,(((COEF(J,2,L,I),L=1,LMWR9 )
     &          ,I=1,2),J=2,JMAX)
 103        FORMAT('COEF AT STEP ',I8,' TIME',1PE13.4,/,(1P10E13.5))
            WRITE(9,124)(RHF(I),I=2,JMAX)
 124        FORMAT(' HALF RADII:',/,(1P10E13.5))
 102        CONTINUE

c...Vertical Fourier information (at moment for m=1-4)

c            mdnum = 4
c            do k=2,kmax
c               WRITE(18,203)k,TIME,(((COEF(J,k,L,I),L=1,mdnum ),I=1,2),
c     &              J=2,JMAX)
 203           FORMAT('COEF AT    K',I5,' TIME',1PE13.4,/,(1P10E13.5))
c               WRITE(18,224)(RHF(I),I=2,JMAX)
 224           FORMAT(' HALF RADII:',/,(1P10E13.5))
c            end do

c...Calculate total Fourier components m=1-8

            do j=2,jmax2
               massf(j)=0.d0
               massc(j)=0.d0
            end do
            
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j,l,m)                         &
!$OMP&  SHARED(zof3n,dtheta)
!$OMP DO SCHEDULE(STATIC)
            do m=1,8
               sumeven(m)=0.d0
               sumodd(m)=0.d0
               totcoef(m)=0.d0
               do i=1,10
                  totcoefcym(i,m)=0.d0
               end do
            end do
!$OMP END DO NOWAIT
!
!....Calculate SQRT(a**2 +b**2)
!....don't forget normalizations! 2x for other half of
!....disk. Normalization by A(0) hidden because A(0)=
!....total mass = 1.
!
!$OMP DO SCHEDULE(STATIC)
             do l=1,lmax
                theta(l)=dtheta*(l-1)
             end do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
             do j=2,jmax1
                volume(j)=0.5d0*dtheta*zof3n*(r(j+1)**2.-r(j)**2.)
             end do
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
             do m=1,8

               do j=2,jmax
               do k=2,kmax
               do l=1,lmax
                  sumeven(m)=sumeven(m)+rho(j,k,l)
     &                 *volume(j)*cos(m*theta(l))
                  sumodd(m)=sumodd(m)+rho(j,k,l)
     &                 *volume(j)*sin(m*theta(l))
               end do
               end do
               end do

               totcoef(m)=
     &              4.*SQRT(sumeven(m)**2. +sumodd(m)**2.)


            end do
!$OMP END DO nowait
!$OMP END PARALLEL

            write(13,104)time,(totcoef(m),m=1,8)

         END IF

C...Diagnostic Information for COM and M=1 power...
         X    = 0.d0
         Y    = 0.d0
         DISP = 0.d0

         CALL CENTMASS(X,Y,DISP,ALF)
         
         avgm = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)                            &
!$OMP&  SHARED(jreq) REDUCTION(+:avgm)
         do j=2,JREQ
            avgm=avgm + coef(j,2,1,1)
         enddo
!$OMP END PARALLEL DO
         avgm = avgm/(JREQ-1)
         write(19,1111) itstep,time,x,y,disp,alf,avgm
 1111    format(i8,6(1x,1pe13.4))

         IF(MOD(ITSTEP,IDIAG).eq.0) then
            WRITE(3,40)X,Y,DISP,ALF
 40         FORMAT(' CENTER OF MASS:   X=',1PE12.4,'   Y=',1PE12.4,           
     &           '   DISPLACEMENT=',1PE12.4,'   ANGLE=',0PF7.2)            
         ENDIF

 104     format(1pe13.4,8(1x,1pe13.4))

                      !t_rite = t_rite+dwalltime00(2) !TIMING

         if((MOD(itstep,ISTOR).eq.0))then

            if(ITSTEP .ne. ITSTOP) then
               close(9)
               close(19)
               close(13)
               close(15)
               close(17)
               close(21)
               close(3)
            
               write(index,'(i8.8)')itstep+ISTOR
               resfile =trim(outpath)//'rslts.'//index
               coeffile=trim(outpath)//'coefs.'//index
               modefile=trim(outpath)//'modes.'//index
               centfile=trim(outpath)//'c_o_m.'//index
               tcoeffile=trim(outpath)//'tcoef.'//index
               tcoeftenfile=trim(outpath)//'tctot.'//index
               massfile=trim(outpath)//'massflow.'//index
               epsfile=trim(outpath)//'coolheat.'//index
               planetfile=trim(outpath)//'planet_pos.'//index
               planeten=trim(outpath)//'planet_en.'//index

               SELECT CASE (star_state)
               CASE("wiggle")
                  close(987)
                  close(988)
                  starfile=trim(outpath)//'starrysim.'//index//' '
                  restart_star=trim(outpath)//'starry_restart.'
     &                         //index//' '
                  open(unit=987,file=trim(starfile),form='FORMATTED')
                  open(unit=988,file=trim(restart_star),
     &                 form='FORMATTED')
               CASE("wiggle_restart")
                  close(987)
                  close(988)
                  starfile=trim(outpath)//'starrysim.'//index//' '
                  restart_star=trim(outpath)//'starry_restart.'
     &                         //index//' '
                  open(unit=987,file=trim(starfile),form='FORMATTED')
                  open(unit=988,file=trim(restart_star),
     &                 form='FORMATTED')
               END SELECT

               OPEN(unit=9,file=coeffile)
               OPEN(unit=19,file=centfile)
               open(unit=13,file=tcoeffile)
               open(unit=15,file=tcoeftenfile)
               open(unit=3, file=resfile)
               open(unit=17,file=massfile,form='formatted')
               open(unit=21,file=epsfile,form='formatted')
               IF (planet) THEN
                  close(22)
                  close(24)
                  open(unit=22,file=planetfile,form='formatted')
                  open(unit=24,file=planeten,form='formatted')
               ENDIF
               
      
            endif
         ENDIF
         
         IF(ITSTEP.eq.ITSTOP)THEN

            write(index,'(i8.8)')itstop
            savedfile=trim(outpath)//'saved.'//index
            open(unit=8,file=savedfile,form='unformatted')
            WRITE(8) S
            WRITE(8) T
            WRITE(8) A
            WRITE(8) RHO
            WRITE(8) EPS
            WRITE(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
            if (jmin.gt.2) write(8) tmassini,tmass,tmassadd,tmassout
     &        ,tmassacc,totcool,totdflux,totheat,totirr,etotfl,eflufftot
            if (planet) write(8) planetx,planety,planetvx,planetvy
            close(8)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,l)
            do l=1,lmax
               do j=1,jmax2
                  do k=1,kmax2
                     tauross(j,k,l)=tau(j,k,l,1)
                  enddo
               enddo
            enddo
!$OMP END PARALLEL DO

            epsfull=trim(outpath)//'coolheat_full.'//index
            open(unit=29,file=epsfull,form='unformatted')
            write(29) divflux
            write(29) lambda
            write(29) hgamma
            write(29) igamma
            write(29) tauross
            write(29) TempK
            write(29) TeffK
            write(29) TphK
            write(29) time
            close(29)

            open(unit=29,file=trim(outpath)//"gamma1."//index,
     &           form='unformatted')
            write(29) gamma1
            write(29) time
            close(29)

         endif
         
 90   CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------end main loop-------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (fluid_elements) then 
        call fluid_writerestart(itstop)
      endif
      
      SELECT CASE (star_state)
      
      CASE("wiggle")

           write(988,'(8(1X,1pe15.8))') time,vx,vy,rpstar,
     &     phi_star,rpstar*cos(phi_star),rpstar*sin(phi_star),delt

      CASE("wiggle_restart")  

           write(988,'(8(1X,1pe15.8))') time,vx,vy,rpstar,
     &     phi_star,rpstar*cos(phi_star),rpstar*sin(phi_star),delt   

      END SELECT

      ITSTEP = ITSTEP - 1

 119  format(3(1x,1pe16.8))


C.....Finished with Current Run Here.....C
C.........Run Diagnostics Follow.........C
     

1113  format(i3,1x,i3,1x,1pe13.4,1x,1pe13.4)
1115  format(i2,1x,i2,1x,1pe16.8,1x,1pe16.8,1x,1pe16.8,1x,1pe16.8)

      CALL CENTMASS(X,Y,DISP,ALF)
      CALL RITE(2,-1,ITSTOP,ITSTRT,ITSTOP,ISTOR,1,1,1,1,1)


      close(17)
      close(21)
      close(22)
      close(987)
      close(988)

      STOP
      END





