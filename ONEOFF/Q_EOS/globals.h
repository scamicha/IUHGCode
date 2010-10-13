!-*-f90-*-
!***********************************************************************
!

      real*8 rcloud,constp,delt,bdytem,den,time,cormas,epscen
      COMMON /BLOK7/RCLOUD,CONSTP,DELT,BDYTEM,DEN,TIME,CORMAS,epscen 

      real*8 rholmt,epslmt,dumlmt
      COMMON /RELIMITS/RHOLMT,EPSLMT,DUMLMT

      integer jin,itype
      real*8 tmassini,tmassadd,tmassout,tmassacc
      real*8 starphi,diskphi,tinphi
      COMMON /GAP/starphi(jmax2,kmax2,lmax),                    &             
           tinphi(jmax2,kmax2,lmax),                            &
           diskphi(jmax2,kmax2,lmax),                           &
           tmassini,tmassadd,tmassout,tmassacc,                 &
           jin,itype


      REAL*8 s,t,a
      COMMON /EOM/S(JMAX2,KMAX2,LMAX),T(JMAX2,KMAX2,LMAX),      &
           A(JMAX2,KMAX2,LMAX)

      real*8 p,cv,eps
      real*8 phi,rho,tempc
      COMMON /STATES/P(JMAX2,KMAX2,LMAX),CV(JMAX2,KMAX2,LMAX),  &
           EPS(JMAX2,KMAX2,LMAX),TEMPC(JMAX2,KMAX2,LMAX)


      COMMON /POIS/PHI(POT3JMAX2,POT3KMAX2,LMAX),RHO(JMAX2,KMAX2,LMAX)

      real*8 Msyscgs,PKcgs,Tconv,Sconv,Dconv,Pconv,sigma,       &
           rhoconv,engconv,bkmpcode
      COMMON /CONVERT/Msyscgs,PKcgs,Tconv,Sconv,                &
           Dconv,Pconv,sigma,rhoconv,engconv,bkmpcode

      real*8 r,z,rhf,zhf,rof3n,zof3n
      COMMON /GRID/jreq,R(pot3jmax2),Z(POT3KMAX2),              &
           RHF(POT3JMAX2),ZHF(POT3KMAX2),ROF3N,ZOF3N

      real*8 :: temptable,engtable,gammatable,muc
      common /engtables/temptable(TTABLE),engtable(TTABLE),     &
              gammatable(TTABLE),muc

      real*8 etotfl, efl, eflufftot, gamma1
      COMMON /etally/ etotfl, efl, eflufftot,gamma1(JMAX2,KMAX2,LMAX)

      real*8 tempk,lambda,TeffK,TphK,divflux,hgamma,igamma,tauross
      COMMON /COOLING/tempk(JMAX2,KMAX2,LMAX),                  &
              lambda(JMAX2,KMAX2,LMAX),TeffK(JMAX2,LMAX),       &
              TphK(JMAX2,LMAX),divflux(JMAX2,KMAX2,LMAX),       &
              hgamma(JMAX2,KMAX2,LMAX),igamma(JMAX2,KMAX2,LMAX),&
              tauross(JMAX2,KMAX2,LMAX)

!$OMP THREADPRIVATE(/GRID/)
      









