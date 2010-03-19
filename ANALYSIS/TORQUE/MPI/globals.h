!-*-f90-*-
C***********************************************************************
C     
      real(KIND=8), dimension(JMAX2)    :: sigma_fc,
     &        cs_fc, omega_fc, tgrav, treyn,gt1,
     &        mass_fc,kepler_fc,rtv,oav,vrav, cs_sig

      real(KIND=8) :: vphiav,econv,rdisk,mstar,kconst

      real(KIND=8), dimension(JMAX2) :: masss,masse,mdot

      COMMON /alphadisk/sigma_fc,gt1,
     &        cs_fc, omega_fc, tgrav, treyn,
     &        mass_fc,kepler_fc,rtv,oav,vrav,
     &        cs_sig, masss,masse,mdot,
     &        vphiav,econv,rdisk,mstar,kconst
 

      REAL*8 KONST,NPRIME,xn,gamma,toverw
      COMMON /PTROPE/XN,GAMMA,KONST,NPRIME,TOVERW

      real*8 rcloud,constp,delt,bdytem,den,time,cormas,epscen
      COMMON /BLOK7/RCLOUD,CONSTP,DELT,BDYTEM,DEN,TIME,CORMAS,epscen 

      real*8 rholmt,epslmt,dumlmt
      COMMON /RELIMITS/RHOLMT,EPSLMT,DUMLMT

      integer jin,itype
      real*8 tmassini,tmassadd,tmassout,tmassacc
      real*8 starphi,diskphi,tinphi
      COMMON /GAP/starphi(jmax2,kmax2,lmax),
     &     tinphi(jmax2,kmax2,lmax),
     &     diskphi(jmax2,kmax2,lmax),
     &     tmassini,tmassadd,tmassout,tmassacc,
     &     jin,itype


      REAL*8 JN,s,t,a,u,w,omega
      COMMON /EOM/
     &     S(JMAX2,KMAX2,LMAX),
     &     T(JMAX2,KMAX2,LMAX),
     &     A(JMAX2,KMAX2,LMAX),
     &     U(JMAX2,KMAX2,LMAX),
     &     W(JMAX2,KMAX2,LMAX),
     &     JN(JMAX2,KMAX2,LMAX),
     &     OMEGA(JMAX2,KMAX2,LMAX)

      real*8 p,cv,eps
      real*8 phi,rho,tempc,rhosave
      COMMON /STATES/
     &     P(JMAX2,KMAX2,LMAX),
     &     CV(JMAX2,KMAX2,LMAX),
     &     EPS(JMAX2,KMAX2,LMAX),
     &     TEMPC(JMAX2,KMAX2,LMAX), ENON


      COMMON /POIS/
     &     PHI(POT3JMAX2,POT3KMAX2,LMAX),
     &     RHO(JMAX2,KMAX2,LMAX),
     &     RHOSAVE(JMAX2,KMAX2,LMAX)


c     COMMON /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
c     COMMON /TIMEST/INDX,ISOADI,ALLOW,DMAX,CHGMAX
c     COMMON /COEFS/COEF(POT3JMAX2,POT3KMAX2,LMAX2,2)
c     COMMON /LIMIT/SOUND
c     COMMON /ITS/ITSTRT,ITSTOP,ITSTEP

      real*8 Msyscgs,PKcgs,Tconv,Sconv,Dconv,Pconv,sigma
      COMMON /CONVERT/
     &     Msyscgs,PKcgs,Tconv,Sconv,
     &     Dconv,Pconv,sigma   

      real*8 r,z,rhf,zhf,rof3n,zof3n
      COMMON /GRID/jreq,kzpol,
     &     R(pot3jmax2),
     &     Z(POT3KMAX2),
     &     RHF(POT3JMAX2),
     &     ZHF(POT3KMAX2),
     &     ROF3N,ZOF3N

C$OMP THREADPRIVATE(/GRID/)
      









