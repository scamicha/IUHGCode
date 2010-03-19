      SUBROUTINE SOURCE
      IMPLICIT real*8 (a-h,o-z)      

      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'
      INCLUDE 'units.h'

      integer jstart
      REAL*8 VT(JMAX2,KMAX2,LMAX),dume,limiter,surflux
      COMMON /SOURCESHARED/VT,surflux

      real*8 epred,tpred,ppred

      DR=ROF3N
      DZ=ZOF3N

      HCOUNT = 0
      HHIT   = 0 !these are used to tally heating limiter
      
c  Should be ok to run the following loops over jmax1,kmax1 instead of
c  pot3jmax1,pot3kmax1 since we are only interested in sourcing the
c  regions where there is mass; the potential on the outer part of
c  the big grid does not really matter.

      limiter = den*phylim 

      if(jmin.gt.2) then
         jstart=jmin
      else   
         jstart=2
      endif 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k,LP,LM,l)                  &
!$OMP&  SHARED(dtheta,dz,delt,DR,jstart)
      DO L=1,LMAX
        LM=L-1
        LP=L+1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        DO K=2,KMAX1
 
C.....SOURCE S.
          DO J=jstart+1,JMAX1
            S(J,K,L)=S(J,K,L)+DELT*
     &          (R(J)*(RHO(J,K,L)+RHO(J-1,K,L))*
     &          ((OMEGA(J,K,L)+OMEGA(J-1,K,L))**2)/8.
     &          -(RHO(J,K,L)+RHO(J-1,K,L))*(PHI(J,K,L)-PHI(J-1,K,L))
     &          /(2.0*DR)-(P(J,K,L)-P(J-1,K,L))/DR)
          END DO

C......SOURCE T.
          DO J=jstart,JMAX1
            T(J,K,L)=T(J,K,L)-DELT*
     &         (0.5*(RHO(J,K,L)+RHO(J,K-1,L))*(PHI(J,K,L)-PHI(J,K-1,L))
     &         /DZ+(P(J,K,L)-P(J,K-1,L))/DZ)
          END DO
 
C......SOURCE A.
          DO J=jstart,JMAX1
            A(J,K,L)=A(J,K,L)-DELT*(RHO(J,K,L)*(PHI(J,K,LP)-
     &           PHI(J,K,LM))+P(J,K,LP)-P(J,K,LM))/(2.0*DTHETA)
          END DO

        END DO
      END DO
!$OMP END PARALLEL DO

      CALL RADTRANS

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(dume,ppred,tpred,epred,j,k,LP,LM,&
!$OMP& l)                                                               &
!$OMP&  SHARED(etotfl,efl,rholmt,epslmt,pconv,muc,rhoconv,HHIT,         &
!$OMP& HCOUNT,dtheta,dz,delt,DR,jstart,zof3n,rof3n) REDUCTION(+:surflux,&
!$OMP& totirr,totdflux,totheat,totcool)
      CALL VELOCITY
      CALL VLIMIT
      CALL AVISC
  
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            LM=L-1
            LP=L+1
            IF(L.EQ.1) LM=LMAX
            IF(L.EQ.LMAX) LP=1
            DO  K=2,KMAX1

               DO  J=jstart+1,JMAX1
                  S(J,K,L)=S(J,K,L)-DELT*(RHO(J,K,L)*QRR(J,K,L)
     &                 -RHO(J-1,K,L)*QRR(J-1,K,L))/(DR)
               END DO

               DO  J=jstart,JMAX1
                  T(J,K,L)=T(J,K,L)-DELT*(RHO(J,K,L)*QZZ(J,K,L)
     &                 -RHO(J,K-1,L)*QZZ(J,K-1,L))/(DZ)
               END DO
         
               DO  J=jstart,JMAX1
                  A(J,K,L)=A(J,K,L)-DELT*(RHO(J,K,LP)*QTT(J,K,LP)
     &                 -RHO(J,K,LM)*QTT(J,K,LM))/(2.d0*DTHETA) 
               END DO

               DO J=jstart,JMAX1
                  VT(J,K,L)=0.5d0*RHF(J)*(OMEGA(J,K,L)+OMEGA(J,K,LM))
               END DO

            END DO
         END DO
!$OMP END DO

C.......SOURCE EPS. (Includes artificial viscosity, L cooling, divflux,
c...... and irradiation.)

!$OMP DO SCHEDULE(STATIC) REDUCTION(+:HHIT,HCOUNT)
         DO L=1,LMAX
            LM=L-1
            LP=L+1
            IF(L.EQ.1) LM=LMAX
            IF(L.EQ.LMAX) LP=1
            DO K=2,KMAX
               DO J=jstart,JMAX

C    Artificial viscosity

                  Hgamma(J,K,L)=-RHO(j,k,l)*
     &                 (QRR(J,K,L)*(U(J+1,K,L)-U(J,K,L))/DR
     &                 +QZZ(J,K,L)*(W(J,K+1,L)-W(J,K,L))/DZ
     &                 +QTT(J,K,L)*(VT(J,K,LP)-VT(J,K,L))
     &                 /(DTHETA*RHF(J)))

C    Introduce a limiter just like in Igamma

                  if ((Hgamma(j,k,l).gt.0.0)) then

                        HCOUNT = HCOUNT + 1

                     if((eps(j,k,l)/Hgamma(j,k,l)).lt.(theatlmt
     &                 *delt)) then

                          Hgamma(j,k,l)=eps(j,k,l)/(theatlmt*delt)
                  
                          HHIT = HHIT + 1

                     endif

                  endif

C    Depending on the type of cooling...


                  if ((ictype.eq.1).or.(ictype.eq.2)) then

                     Divflux(j,k,l)=0.0

                  else if (ictype.eq.3.or.ictype.eq.5) then

                     Lambda(j,k,l)=0.0

                  else if (ictype.eq.4) then

                     if (tau(j,k,l,1).lt.(2.0/3.0)) then 
                        Divflux(j,k,l)=0.0
                     else
                        Lambda(j,k,l)=0.0
                    endif

                  endif

                  epred = eps(J,K,L) - p(J,K,L) * (
     &                 (r(J+1)*u(J+1,K,L)-r(J)*u(J,K,L))/(rhf(J)*rof3n)
     &               + (vt(J,K,LP)-vt(J,K,L))/(dtheta*rhf(J))
     &               + (w(J,K+1,L)-w(J,K,L))/zof3n)*delt
     &               +
     &                 ( Hgamma(j,k,l) + Igamma(j,k,l) - Lambda(j,k,l)
     &               -   Divflux(j,k,l) ) * delt

                  call TempFindSpec(epred,rho(J,K,L),tpred)

                  ppred = rho(J,K,L)*rhoconv
     &                  * bkmpcgs/(muc*pconv)*tpred
                  
                  eps(J,K,L) = eps(J,K,L) - 0.5d0*(p(J,K,L)+ppred) * (
     &                 (r(J+1)*u(J+1,K,L)-r(J)*u(J,K,L))/(rhf(J)*rof3n)
     &               + (vt(J,K,LP)-vt(J,K,L))/(dtheta*rhf(J))
     &               + (w(J,K+1,L)-w(J,K,L))/zof3n)*delt
     &               +  
     &                 ( Hgamma(j,k,l) + Igamma(j,k,l) - Lambda(j,k,l) 
     &               -   Divflux(j,k,l) ) * delt               

                  eps(J,K,L)=max(eps(J,K,L),epslmt)

                  if(rho(j,k,l).le.rholmt) eps(j,k,l)=epslmt

               enddo
            enddo
         enddo
!$OMP END DO         

! master switch added by aaron to tally limited heating cells.

!$OMP MASTER
       dume = dble(HHIT)/dble(HCOUNT)*100.

       if (dume > 1.)print *, " PERCENT TOTAL LIMITED HEATING CELLS ", 
     & dume

      etotfl = etotfl + efl 
      efl = 0.
!$OMP END MASTER


C....Calculate total heating and total cooling 

!$OMP DO SCHEDULE(STATIC)
      do l=1,lmax
         do k=2,kmax
           do j=jcool,jmax
               totcool=totcool+Lambda(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               totheat=totheat+Hgamma(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               totdflux=totdflux+Divflux(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               totirr=totirr+Igamma(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               if((tau(j,k,l,1).ge.(2.0/3.0)).and.(tau(j,k+1,l,1).lt.(2
     &              .0/3.0))) surflux=surflux+Radflux(j,k,l,2)*0.5
     &              *dtheta*(r(j+1)**2-r(j)**2)*delt
               
            enddo
         enddo
      enddo
!$OMP END DO NOWAIT


C....Set quantities around the z-axis.
      
      if (jmin.gt.2) then

!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO K=2,KMAX2
               do j=1,jmin1
                  S(j,K,L)  = 0.0
                  T(j,K,L)  = 0.0
                  A(j,K,L)  = 0.0
                  U(j,k,l)  = 0.0
                  W(j,k,l)  = 0.0
                  EPS(j,K,L)= epslmt
               enddo
            END DO
         END DO
!$OMP END DO NOWAIT

      else   
         
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            LP=1+MOD(L-1+LMAX/2,LMAX)
            DO K=2,KMAX2
               S(1,K,L)  = -S(3,K,LP)
               T(1,K,L)  = T(2,K,LP)
               A(1,K,L)  = A(2,K,LP)
               EPS(1,K,L)= EPS(2,K,LP)
            END DO
         END DO
!$OMP END DO NOWAIT

      endif


!$OMP DO SCHEDULE(STATIC) 
      DO L=1,LMAX

C.....Set S and U on the z-axis.
        DO K=2,KMAX
          U(2,K,L)=0.d0
          S(2,K,L)=0.d0
        END DO

C....Set quantities on the bottom of the grid.
        DO J=1,JMAX2
          S(J,1,L)  = S(J,2,L)
          T(J,1,L)  = -T(J,3,L)
          A(J,1,L)  = A(J,2,L)
          EPS(J,1,L)= EPS(J,2,L)
          W(J,2,L)  = 0.d0
          T(J,2,L)  = 0.d0
        END DO

      END DO
!$OMP END DO 

      CALL CLEANUP
!$OMP END PARALLEL

      RETURN
      END

C*********************************************************

      SUBROUTINE RADTRANS
      IMPLICIT real*8 (a-h,o-z)      
      
      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'
      INCLUDE 'units.h'
      
      integer jstart,ipta,kfit, kfitmax
      REAL*8 xmu,sum,Tfit,Teff,Teff4,Tbd4, Tenv4, tirrlmt,limiter
      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),absfrac
      REAL*8 Oross_ph, Oplck_env(jmax2,kmax2,lmax), Tfit_old
     &     ,Oplck_ph(jmax2,kmax2,lmax), Oirr, tau_face(kmax2,4)
      REAL*8 Tbdry4, tau_f(jmax2,3)
      real*8 dtau,dtau1
      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum, Oplck_env, 
     &      Oplck_ph
      
      DR=ROF3N
      DZ=ZOF3N
      
      limiter = den*phylim 

      if(jmin.gt.2) then
         jstart=jmin
      else   
         jstart=2
      endif 

C... New cooling (acm 2002)   

C    To convert from polytropic units of temperature to Kelvin
C    Use T(K)=T(pu)*Tconv
C    Opacities in cgs, so always use Oross or Oabs or Oplck times Sconv.

C    In this next loop, read opacities from table using the pollack 
C    function or the dalessio subroutine.
	
	Tenv4 = (TenvK/Tconv)**4	
c... Added by Kai: Envelope heating term

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Teff4,absfrac,dtau1,dtau,k,j,l,  &
!$OMP&Oirr)                                                             &
!$OMP&  SHARED(Dconv,Tconv,sigma,Sconv,dz,Pconv,                        &
!$OMP& limiter,jstart)
!$OMP DO SCHEDULE(STATIC)
      do l=1,lmax
         do j=jstart,jmax
            TeffK(j,l)=0.0
            TphK(j,l)=0.0
            surfcgs(j,l)=0.0
         enddo

         do k=1,kmax1
           do J = JSTART, JMAX 
               tau(j,k,l,1)=0.0
               tau(j,k,l,2)=0.0
               tau(j,k,l,3)=0.0
               tau(j,k,l,4)=0.0
               l_tau_z(J,K,L) = 0.d0
               Lambda(j,k,l)=0.0
               Hgamma(j,k,l)=0.0
               Igamma(j,k,l)=0.0
               Oross(j,k,l)=0.0
               Oplck(j,k,l)=0.0
               Oplck_env(j,k,l)=0.0
	       Oplck_ph(j,k,l)=0.0
               Oabs(j,k,l)=0.0
               Otot(j,k,l)=0.0
               Radflux(j,k,l,1)=0.0
               Radflux(j,k,l,3)=0.0
               Divflux(j,k,l)=0.0

               if (rho(j,k,l).ge.limiter) then
                  
                if (TempK(j,k,l).lt.Tbgrnd) TempK(j,k,l)=Tbgrnd
                  call dalessio(TenvK, (P(j,k,l)*Pconv), Oross(j,k,l),
     &                 Oplck_env(j,k,l), Oirr, Otot(j,k,l))
		   
                  call dalessio(TempK(j,k,l),(P(j,k,l)*Pconv),Oross(j,k
     &                 ,l),Oplck(j,k,l),Oabs(j,k,l),Otot(j,k,l))
                 Oross(j,k,l) = fk*Oross(j,k,l)
                 Oplck(j,k,l) = fk*Oplck(j,k,l)
                 Oplck_env(j,k,l) = fk*Oplck_env(j,k,l)                  
               endif
            enddo
          enddo

C     The following loop calculates vertical optical depths and the surface
c     density.  All taus are cell centered.
! the commented lines could be used to control maximum dtau change if such
! an occasion becomes necessary

          do k=kmax,2,-1
            do J = JSTART,JMAX
               if (tau(j,k+1,l,1).eq.0.0) then
                  if (rho(j,k,l).ge.limiter) then
                     dtau = oross(J,K,L)*rho(J,K,L)*sconv*dz*0.5d0
!                     if ( dtau > 5. ) dtau = 5.
                     tau(j,k,l,1)=dtau

                     dtau = oplck(J,K,L)*rho(J,K,L)*sconv*dz*0.5d0
!                     if ( dtau > 5. ) dtau = 5.
                     tau(j,k,l,2)=dtau

                     dtau = oplck_env(J,K,L)*rho(J,K,L)*sconv*dz*0.5d0
!                     if ( dtau > 5. ) dtau = 5.
                     tau(j,k,l,4)=dtau

                     surfcgs(j,l)=rho(j,k,l)*Sconv*dz/2.0
                  endif
               else      
                  dtau = oross(J,K+1,L)*rho(J,K+1,L)*sconv*dz*.5d0
                  dtau1= oross(J,K,L)*rho(J,K,L)*sconv*dz*0.5d0 
!                  if ( dtau > 5. ) dtau = 5.
!                  if ( dtau1 > 5. ) dtau1 = 5.
                  tau(j,k,l,1)=tau(J,K+1,L,1) + dtau + dtau1

                  dtau = oplck(J,K+1,L)*rho(J,K+1,L)*sconv*dz*.5d0
                  dtau1= oplck(J,K,L)*rho(J,K,L)*sconv*dz*0.5d0                      
!                  if ( dtau > 5. ) dtau = 5.
!                  if ( dtau1 > 5. ) dtau1 = 5.
                  tau(j,k,l,2)=tau(J,K+1,L,2) + dtau + dtau1

                  dtau = oplck_env(J,K+1,L)*rho(J,K+1,L)*sconv*dz*.5d0
                  dtau1= oplck_env(J,K,L)*rho(J,K,L)*sconv*dz*0.5d0                      
!                  if ( dtau > 5. ) dtau = 5.
!                  if ( dtau1 > 5. ) dtau1 = 5.
                  tau(j,k,l,4)=tau(J,K+1,L,4) + dtau + dtau1

                  surfcgs(j,l)=surfcgs(j,l) + ((rho(j,k+1,l)+rho(j,k,l))
     &                 *Sconv*dz/2.0)
               endif 
            enddo
         enddo

         if ( ICTYPE == 69 ) then
           do J = JSTART,JMAX
             l_tau_z(J,2,L) = tau(J,2,L,1) + 0.5d0*oross(J,2,L)*
     &         rho(J,2,L)*sconv*dz
           enddo

           do k=kmax,2,-1
             do J = JSTART,JMAX
               if (tau(j,k+1,l,1).eq.0.0) then

                  if (rho(j,k,l).ge.limiter) then
                     dtau = (l_tau_z(J,2,L)*oross(J,K,L)+oplck(J,K,L)
     &                    / l_tau_z(J,2,L))*rho(J,K,L)*sconv*dz*0.5d0 /
     &                    ( l_tau_z(J,2,L) + 1./l_tau_z(J,2,L))
!                     if ( dtau > 5. ) dtau = 5.
                     tau(j,k,l,1)=dtau

                  endif
 
               else      

                  dtau =   (l_tau_z(J,2,L)*oross(J,K+1,L)
     &                 + oplck(J,K+1,L)/l_tau_z(J,2,L)) / 
     &                 ( l_tau_z(J,2,L) + 1./l_tau_z(J,2,L)) *
     &                 rho(J,K+1,L) * sconv * dz * 0.5d0 

                  dtau1=   (l_tau_z(J,2,L)*oross(J,K,L)
     &                 + oplck(J,K,L)/l_tau_z(J,2,L)) /
     &                 ( l_tau_z(J,2,L) + 1./l_tau_z(J,2,L)) *
     &                 rho(J,K,L) * sconv * dz * 0.5d0

!                  if ( dtau > 5. ) dtau = 5.
!                  if ( dtau1 > 5. ) dtau1 = 5.
                  tau(j,k,l,1)=tau(J,K+1,L,1) + dtau + dtau1

               endif 
             enddo
          enddo
         endif

! I have taken out the following code because it is obsolete.  I am
! leaving it commented because it might need to be resurrected in
! a much more efficient format.

! BEGIN OBSOLETE CODE

C     The following loop calculates horizontal optical depths for the
C     irradiation.  Since the heating by irradiation affects one or two
c     cells, I am dividing the opacities by 'irop', defined in units.h.
c ... Annie        
!         do k=2,kmax
c            tau(jstart-1,k,l,3)=0.0
c            do j=jstart,jmax 
c               if (tau(j-1,k,l,3).eq.0.0) then
c                  if (rho(j,k,l).ge.limiter)
c     &               tau(j,k,l,3)=(Oabs(j,k,l)/irop)*rho(j,k,l)*Sconv
c     &                 *dr/2.0
c               else
c                  tau(j,k,l,3)=tau(j-1,k,l,3) + (((Oabs(j-1,k,l)/irop)
c     &                 *rho(j-1,k,l)+(Oabs(j,k,l)/irop)*rho(j,k,l))
c     &                 *Sconv*dr/2.0)
c               endif
c
c            enddo
c ... Change tau(j,k,l,3) to horizontal optical depth start from the 
c     outer edge of the grid: for Envelope heating (irop not used) 
c ... Kai
!            tau(jmax,k,l,3)=0.0
!          do K = 2, KMAX
!            do j=jmax,jstart, -1 
!              if (tau(j+1,k,l,3).eq.0.0) then             
!                 if (rho(j,k,l).ge.limiter)
!     &            tau(j,k,l,3)=Oplck_env(j,k,l)*rho(j,k,l)*Sconv*dr
!     &              /2.0
!              else
!                  tau(j,k,l,3)=tau(j+1,k,l,3) + (Oplck_env(j+1,k,l)*rho
!     &              (j+1,k,l)+Oplck_env(j,k,l)*rho(j,k,l))*Sconv*dz/2.0
!              endif
!            enddo
!         enddo

! END OBSOLETE CODE
   
      enddo    
!$OMP END DO

C      print*,'ICTYPE', ictype,irtype,jcool,jirr

C     call HEIGHTS ! obsolete: may need to be resurrected.

C     Irradiation according to D'Alessio et al. 1998.  Avoid heating 
C     anything under radial zone 28 because of artificial gap.  That 
c     material would be shadowed in a real disk. This is model specific!
C     See units.h for irop.  Absfrac is 1 - albedo = abs/(abs+scatt),
C     used as approximation since we do not treat scattering. 
 
                    
      IF (irtype.eq.1) THEN

      PRINT *, " IRTYPE==1: FUNCTIONALITY DISABLED "
      stop

!$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
            do j=jstart,jmax
               do k=2,kmax
                  if((rho(j,k,l).ge.limiter).and.(j.gt.jirr)) then

                     absfrac = Oabs(j,k,l)/Otot(j,k,l)
                     Igamma(j,k,l)=sigma*((Tstar/Tconv)**4)*
     &                    (Oabs(j,k,l)/irop)*Sconv*rho(j,k,l)*
     &                    ((Rstar*Rsuncgs/Dconv)**2)/
     &                    (rhf(j)**2+zhf(k)**2)*
     &                    exp(-tau(j,k,l,3))*absfrac
                  endif

C     Regulate the irradiation 
C     If tirr is lower than tirrlmt, make Igamma an upper limit.

                  if (Igamma(j,k,l).lt.1e-40) Igamma(j,k,l)=0.0
                  
                  if ((Igamma(j,k,l).gt.0.0).and.
     &                 ((eps(j,k,l)/Igamma(j,k,l)).lt.(theatlmt
     &                 *torp))) Igamma(j,k,l)=eps(j,k,l)/(theatlmt*torp)
                  
               enddo
            enddo
         enddo
!$OMP END DO NOWAIT

      ENDIF 

C     IF ICTYPE = 1 (Constant Tcool)

      IF (ICTYPE.EQ.1) THEN

!$OMP DO SCHEDULE(STATIC)  
      do l=1,lmax
         do j=jcool,jmax
            Teff4=0.0
            do k=2,kmax

               if (rho(j,k,l).ge.limiter) then
                  lambda(j,k,l)=eps(j,k,l)/(cct*torp)
               endif

               if(lambda(j,k,l).le.0.0) lambda(j,k,l)=0.0
               
               Teff4 = Teff4 + lambda(j,k,l)*dz/sigma

C     Store temperatures in kelvin.  These are pure estimates
C     since this cooling doesn't use opacities or temperatures or anything.

C     Weighted average for photospheric temperature (at tau=2/3)
               
               if ((tau(j,k,l,1).ge.(2.0/3.0)).and.(tau(j,k+1,l,1).lt
     &              .(2.0/3.0))) TphK(j,l)=10.0**((log10(TempK(j,k,l))+
     &              log10(TempK(j,k+1,l)))/2.0)
               
            enddo
            
C     Effective temperature is the sum of the cooling everywhere
            
            TeffK(j,l)=(Teff4**0.25) * Tconv
            
         enddo  
      enddo

!$OMP END DO NOWAIT

      ENDIF
!$OMP END PARALLEL
      
C     IF ICTYPE = 2 (Eddington atm)

      IF (ICTYPE.EQ.2) THEN

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Teff4,k,Tfit,sum,j,l)         &
!$OMP&  SHARED(sigma,gamma,dz,limiter,Tconv,jstart)
         do l=1,lmax
            do j=jstart,jmax
               sum=0.0
               Tfit = TempK(j,2,l)/Tconv

               do k=2,kmax
                  if (rho(j,k,l).ge.limiter) 
     &                 sum=sum + dz*rho(j,k,l)**gamma
               enddo
            
C     Define effective temperatures. 
            
               if (tau(j,2,l,1).lt.(2.0/3.0)) then
                  Teff4=Tfit**4 * tau(j,2,l,2)
               else
                  Teff4=Tfit**4 * (4.0/(3.0*(tau(j,2,l,1)+2
     &                 .0/3.0)))
C                  TphK(j,l)=(Teff4**0.25) * Tconv  
               endif   

C     Store temperatures in kelvin. 

               TeffK(j,l)=(Teff4**0.25) * Tconv

               if ((sum.gt.0.0).and.(rho(j,2,l).gt.0.0).and.(j.ge.jcool)
     &              ) lambda(j,2,l)=(sigma*Teff4)/(rho(j,2,l)*sum)
               
               do k=2,kmax
                  if (rho(j,k,l).ge.limiter) then
                     
                     lambda(j,k,l)=lambda(j,2,l)                        
                     
                     Lambda(j,k,l)=lambda(j,k,l)*rho(j,2,l)*
     &                    rho(j,k,l)**gamma
                     
                     if ((tau(j,k,l,1).ge.(2.0/3.0)).and.(tau(j,k+1,l,1)
     &                    .lt.(2.0/3.0))) TphK(j,l)=10.0
     &                    **((log10(TempK(j,k,l))+log10(TempK(j,k+1,l)))
     &                    /2.0)
                     
                  endif

C     Regulate lambda: if tcool is lower than tcoollmt, make lambda an
C     upper limit.
                  
                  if ((Lambda(j,k,l).ne.0.0).and.
     &                 (abs(eps(j,k,l)/Lambda(j,k,l)).lt.(tcoollmt
     &                 *torp))) Lambda(j,k,l)=eps(j,k,l)/(tcoollmt
     &                 *torp)

               enddo

            enddo
         enddo
!$OMP END PARALLEL DO
        print *," OLD LIMITER SET FOR THIS ROUTINE "

      ENDIF

C     IF ICTYPE = 3 (Diffusion approximation only)

      IF (ICTYPE.EQ.3) THEN

         do l=1,lmax
            lm=l-1
            if(l.eq.1) lm=lmax
            lp=l+1
            if(l.eq.lmax) lp=1
            do j=jstart,jmax
               Teff4=0.0
               do k=kmax,2,-1
                  
                  if (rho(j,k,l).ge.limiter) then
                     
                     if(tau(j,k,l,1).lt.(10.0)) then 
                        
                        Radflux(j,k,l,1) = 0.0
                        Radflux(j,k,l,2) = 0.0
                        Radflux(j,k,l,3) = 0.0
                     else                        
                        Radflux(j,k,l,1)=fluxlmdf(TempK(j,k,l),TempK(j-1
     &                       ,k,l),Oross(j,k,l),Oross(j-1,k,l),rho(j,k,l
     &                       ),rho(j-1,k,l),j,1)
                        Radflux(j,k,l,2)=fluxlmdf(TempK(j,k,l),TempK(j,k
     &                       -1,l),Oross(j,k,l),Oross(j,k-1,l),rho(j,k,l
     &                       ),rho(j,k-1,l),j,2)
                        Radflux(j,k,l,3)=fluxlmdf(TempK(j,k,lp),TempK(j
     &                       ,k,l),Oross(j,k,lp),Oross(j,k,l),rho(j,k,lp
     &                       ),rho(j,k,l),j,3)

C     Boundary conditions for the fluxes
                        
                        if (tau(j,k+1,l,1).lt.(10.0)) then
                          Radflux(j,k+1,l,2)=max(Radflux(j,k,l,2),0.0d0)
c                          Radflux(j,k+1,l,2)= fluxlmdf(TempK(j,k+1,l),
c     &                    TempK(j,k,l),Oross(j,k+1,l),Oross(j,k,l),rho(
c     &                    j,k+1,l),rho(j,k,l),j,2)

                          TeffK(j,l)=((Radflux(j,k+1,l,2)/sigma)**0.25)
     &                          * Tconv
                          TphK(j,l) = TeffK(j,l)
                          kfit = k
                        endif
                        if (tau(j-1,k,l,1).lt.(10.0))
     &                       Radflux(j,k,l,1)= -Radflux(j,kfit+1,l,2)
                        if (tau(j+1,k,l,1).lt.(10.0))
     &                       Radflux(j+1,k,l,1)=Radflux(j,kfit+1,l,2)
                        if (tau(j,k,lm,1).lt.(10.0))
     &                       Radflux(j,k,lm,3)=-Radflux(j,kfit+1,l,2)
                        if (tau(j,k,lp,1).lt.(10.0))
     &                       Radflux(j,k,l,3)=Radflux(j,kfit+1,l,2)
                                             
C     Calculate the divergence of fluxes.
                     
                        Divflux(j,k,l)=((r(j+1)*Radflux(j+1,k,l,1))-
     &                       (r(j)*Radflux(j,k,l,1)))/(rhf(j)*dr)  + 
     &                       (Radflux(j,k+1,l,2)-Radflux(j,k,l,2))/dz +
     &                       (Radflux(j,k,l,3)-Radflux(j,k,lm,3))
     &                       /(rhf(j)*dtheta)
                        
                     endif
                  endif
C     Regulate the divergence of the flux. 
C     If tdflux is lower than tcoollmt, make divflux an upper limit.
                  
                  if ((Divflux(j,k,l).ne.0.0).and.
     &                 (abs(eps(j,k,l)/Divflux(j,k,l)).lt.(tcoollmt
     &                 *torp))) Divflux(j,k,l)=eps(j,k,l)/(tcoollmt
     &                 *torp)*Divflux(j,k,l)/abs(Divflux(j,k,l))
c .... ???
                                                       
               enddo   
            enddo
         enddo
         
        print *," OLD LIMITER SET FOR THIS ROUTINE "

      ENDIF

C     IF ICTYPE = 4 (Diffusion approximation + shining atmosphere)

      IF (ICTYPE.EQ.4) THEN
         
         do l=1,lmax
            lm=l-1
            if(l.eq.1) lm=lmax
            lp=l+1
            if(l.eq.lmax) lp=1
            do j=jstart,jmax
               Tfit = TempK(j,2,l)/Tconv
               kfit=2
               Tbd4=0.0
               do k=kmax,2,-1
                  
                  if (rho(j,k,l).ge.limiter) then
                     
C    For optically thin, need lambdas and temperature from material
C    shining down on the disk
                     
                     if (tau(j,k,l,1).lt.(2.0/3.0)) then
                        
                       Lambda(j,k,l)=Oplck(j,k,l)*Sconv*rho(j,k,l
     &                       )*sigma*((TempK(j,k,l)/Tconv)**4)*4.0

                           Tbd4=Tbd4 + lambda(j,k,l)*dz/(2.0*sigma)

C     For optically thick, need fluxes
                        
                     else if (tau(j,k,l,1).ge.(2.0/3.0)) then

C     For the boundary, define Tfit, kfit.                        

                        if (tau(j,k+1,l,1).lt.(2.0/3.0)) then
                           Tfit=TempK(j,k,l)/Tconv
                           kfit=k
c                           lambda(j,k+1,l) = 0.d0
                        endif


                        Radflux(j,k,l,1)=fluxlmdf(TempK(j,k,l),TempK(j-1
     &                       ,k,l),Oross(j,k,l),Oross(j-1,k,l),rho(j,k,l
     &                       ),rho(j-1,k,l),j,1)
                        Radflux(j,k,l,2)=fluxlmdf(TempK(j,k,l),TempK(j,k
     &                       -1,l),Oross(j,k,l),Oross(j,k-1,l),rho(j,k,l
     &                       ),rho(j,k-1,l),j,2)
                        Radflux(j,k,l,3)=fluxlmdf(TempK(j,k,lp),TempK(j
     &                       ,k,l),Oross(j,k,lp),Oross(j,k,l),rho(j,k,lp
     &                       ),rho(j,k,l),j,3)
                     endif

                  endif
                  
               enddo

C     Define the temperatures and the flux at the boundary
C     between optical depths.  I call this the photospheric temperature,
C     even thought it might not be interpolated for tau= 2/3.  The
C     temperature seen by the observer adds the Tbd and the photospheric
C     temperature. I call this the effective temperature. 

               if (tau(j,2,l,1).lt.(2.0/3.0)) then
                  Teff4=2.0*Tbd4
                  TphK(j,l)=0.d0
               else
	          Tbd4=Tbd4 - lambda(j,kfit+1,l)*dz/(2.0*sigma)
	          lambda(j,kfit+1,l) = 0.d0                  
                  Radflux(j,kfit+1,l,2)=sigma*4.0*(Tfit**4-Tbd4
     &                -Tenv4*exp(-tau(j,kfit+1,l,4)))/(
     &                3.0*(tau(j,kfit,l,1)+2.0/3.0))
c .. still using tau @ ctr even though it really absorbs by whole cell
c ... Envelope heating term: downward so -
                  Tbdry4 = Radflux(j,kfit+1,l,2)/sigma
                  TphK(j,l)=(Tbdry4 +Tbd4+Tenv4*
     &                exp(-tau(j,kfit+1,l,4)))**0.25 * Tconv

c                  if (TphK(j,l).le.0.0) then
c                     print*, j, l, 'TphK <= 0!'
c                     TphK(j,l) = 0.d0
c                  endif
c ... This shouldn't ever happen: for debugging
                  
                  if (Radflux(j,kfit+1,l,2).lt.0.0) 
     &               TphK(j,l)=-TphK(j,l)
c ....    	     print*, 'F_bdry is downward.'       
               Tfit_old = TempK(j,kfit+1,l)
               TempK(j,kfit+1,l) = (3.0/4.0*Tbdry4*(tau(j,kfit+1,l,1)
     &            +2.0/3.0)+Tbd4+Tenv4*exp(-tau(j,kfit+1,l,4)))
     &              **0.25 * Tconv
c ... Eddington layer: reset to BC solution

                  Teff4=(TphK(j,l)/Tconv)**4*exp(-tau(j,kfit+1,l,2))
     &                    + Tbd4
               endif                

C     Store effective temperatures in kelvin
               
               TeffK(j,l)=(Teff4**0.25)* Tconv

               do k=2,kmax

C     Calculate heating of the atmosphere by the flux at the boundary

                 if (rho(j,k,l).ge.limiter) then

                    if (tau(j,k,l,1).lt.(2.0/3.0)) then
c ... if tau(j,2) < 2/3 then TphK(j,l) = 0 => Fb heating = 0
                              tau_face(k,2) = tau(j,k,l,2) + Oplck(j,k,
     &                            l)*Sconv*rho(j,k,l)*dz/2.0
                          tau_face(k,4) = max((tau(j,k,l,4) - Oplck_env
     &                         (j,k,l)*rho(j,k,l)*Sconv*dz/2.0), 0.0d0)
                          tau_f(j,3) = max((tau(j,k,l,3) - Oplck_env(j,
     &                        k,l)*rho(j,k,l)*Sconv*dr/2.0), 0.0d0)
                       if(tau(j,2,l,1).lt.(2.0/3.0).or.(k.gt.(kfit+1)))
     &                   then
     	                      lambda(j,k,l)=lambda(j,k,l)- rho(j,k,l)*
     &                            Sconv*sigma*(Oplck(j,k,l)*(TphK(j,l)
     &                            /Tconv)**4*exp(-(tau_face(kfit+2,2)     
     &                            -tau_face(k,2))) + Oplck_env(j,k,l)*
     &                            Tenv4*(exp(-tau_face(k,4))+exp(-tau_f
     &                            (j,3))))
c ... Start heating from cell kfit+2
c                            if (k.le.kfitmax) lambda(j,k,l)=lambda(j,k
c     &                          ,l)-rho(j,k,l)*Sconv*sigma*Oplck_env(
c     &                          j,k,l)*Tenv4*exp(-tau_f(j,3))      
c .. Envelope heating in the horizontal direction: up to kfitmax
                          endif                              
C ... Smooth lambda radially here: turned on for all k

                       if(tau(j,2,l,1).lt.(2.0/3.0).or.(k.gt.(kfit+1))) 
     &                     then
                             if(lambda(j-2,k,l).ne.0.0) then
                              if(lambda(j-1,k,l).ne.0.0) then
                               lambda(j,k,l)=
     &                         (lambda(j-2,k,l)+ lambda(j-1,k,l) +
     &                          lambda(j,k,l))/3.0
                              else
                               lambda(j,k,l)= (lambda(j-2,k,l)+
     &                           lambda(j,k,l))/2.0
                              endif
                            else if(lambda(j-1,k,l).ne.0.0) then
                             lambda(j,k,l)= (lambda(j-1,k,l)+
     &                          lambda(j,k,l))/2.0
                           endif
                        endif    

                     else
                        
C     Boundary conditions for the fluxes in r and theta


                        if (tau(j-1,k,l,1).lt.(2.0/3.0))
     &                       Radflux(j,k,l,1)=-Radflux(j,kfit+1,l,2)
                        if (tau(j+1,k,l,1).lt.(2.0/3.0))
     &                       Radflux(j+1,k,l,1)=Radflux(j,kfit+1,l,2)
                        if (tau(j,k,lm,1).lt.(2.0/3.0))
     &                       Radflux(j,k,lm,3)=-Radflux(j,kfit+1,l,2)
                        if (tau(j,k,lp,1).lt.(2.0/3.0))
     &                       Radflux(j,k,l,3)=Radflux(j,kfit+1,l,2)
c		    Radflux(j,k,l,1) = 0.0
c		    Radflux(j+1,k,l,1) = 0.0
c	            Radflux(j,k,l,3) = 0.0
c		    Radflux(j,k,lm,3) = 0.0 ...turn them off
                        
C     Calculate the divergence of fluxes.
                        
                        
                        Divflux(j,k,l)=((r(j+1)*Radflux(j+1,k,l,1))-
     &                       (r(j)*Radflux(j,k,l,1)) )/(rhf(j)*dr)  + 
     &                       (Radflux(j,k+1,l,2)-Radflux(j,k,l,2))/dz +
     &                       (Radflux(j,k,l,3)-Radflux(j,k,lm,3))
     &                       /(rhf(j)*dtheta)
                     endif
                     
                  endif 

               enddo
            enddo
         enddo


      LO=0
      CCOUNTER = 0
      HCOUNTER = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,l)                           &
!$OMP&  SHARED(delt,limiter) REDUCTION(+:ccounter,hcounter,LO)
!$OMP DO SCHEDULE(STATIC)
        do L = 1, LMAX
          do K = 2, KMAX1
            do J = JMIN, JMAX

             if (rho(J,K,L) >= limiter ) LO = LO + 1

              if (abs(divflux(J,K,L))
     &           .gt. eps(J,K,L)/(tcoollmt*delt))then
                 if (divflux(J,K,L) < 0.) then
                 divflux(J,K,L) = -eps(J,K,L)/(tcoollmt*delt)
                 if (rho(J,K,L) >= limiter ) HCOUNTER = HCOUNTER + 1
                 else
                   if(divflux(J,K,L) > eps(J,K,L)/(10*delt)) then
                   divflux(J,K,L) = eps(J,K,L)/(10*delt)
                   if (rho(J,K,L) >= limiter ) CCOUNTER = CCOUNTER + 1
                   endif
                 endif
              endif

           enddo
         enddo
       enddo
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
        do L = 1, LMAX
          do K = 2, KMAX1
            do J = JMIN, JMAX

             if (rho(J,K,L) >= limiter ) LO = LO + 1

              if (abs(lambda(J,K,L))
     &           .gt. eps(J,K,L)/(tcoollmt*delt))then
                 if (lambda(J,K,L) < 0.) then
                 lambda(J,K,L) = -eps(J,K,L)/(tcoollmt*delt)
                 if (rho(J,K,L) >= limiter ) HCOUNTER = HCOUNTER + 1
                 else
                   if(lambda(J,K,L) > eps(J,K,L)/(10*delt)) then
                   lambda(J,K,L) = eps(J,K,L)/(10*delt)
                   if (rho(J,K,L) >= limiter ) CCOUNTER = CCOUNTER + 1
                   endif
                 endif
              endif

           enddo
         enddo
       enddo
!$OMP END DO nowait
!$OMP END PARALLEL
        dumb=dble(hcounter)/dble(LO)*100.
        if (dumb>1.) then
          write(6,"(A,1pe12.3)")                                        &
     &    "->Shade Diag: percent relevant h-limited cells ",dumb
        endif
        dumb=dble(ccounter)/dble(LO)*100.
        if (dumb>1.) then
          write(6,"(A,1pe12.3)")                                        &
     &    "->Shade Diag: percent relevant c-limited cells ",dumb
        endif


      ENDIF

      if ( ictype == 69 ) call Hybrid()

      RETURN
      END

      FUNCTION POLLACK(TK)
      IMPLICIT real*8 (a-h,o-z)
      real*8 TK

C     Pollack et al. (1994) rosseland opacities

      if(TK.lt.80.0) then
         pollack=(TK**2)/3200.0
      else if(TK.lt.170.0) then
         pollack=-2.0 + 0.050*TK
      else if(TK.lt.180.0) then
         pollack=62.60 - 0.330*TK
      else if(TK.lt.270.0) then
         pollack=-1.0 + 0.0233*TK
      else if(TK.lt.300.0) then
         pollack=8.0 - 0.010*TK
      else if(TK.lt.425.0) then
         pollack=1.88 + 0.0104*TK
      else if(TK.lt.440.0) then
         pollack=128.13 - 0.2867*TK
      else if(TK.lt.670.0) then
         pollack=0.57 + 0.0033*TK
      else if(TK.lt.700.0) then
         pollack=19.50 - 0.0250*TK
      else if(TK.lt.1300.0) then
         pollack=-0.33 + 0.0033*TK
      else if(TK.lt.1350.0) then
         pollack=24.80 - 0.0160*TK
      else if(TK.lt.1450.0) then
         pollack=46.40 - 0.0320*TK
      else
         pollack=0.01
      endif

      return

      END


      FUNCTION XMMW(TK,Pcgs)

      implicit real*8 (a-h,o-z)
      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'

C     By courtesy of Paola D'Alessio (2002)  
c     limits in T&P: 10**0.5 to 1.e7 K, 1.e-12 to 1.e9 din/cm^2
c     tables in log T, log P, 
      real*8  o1(2),o2(2)
      real*8 TK, Pcgs, xmmw
      
      xlt=dlog10(TK)
      xlp=dlog10(Pcgs)
     
c     Out of bounds
      if (xlp.lt.Ptab(1)) xlp=Ptab(1)
      if (xlp.gt.Ptab(itable)) xlp=Ptab(itable)
      if (xlt.lt.Ttab(1)) xlt=Ttab(1)
      if (xlt.gt.Ttab(itable)) xlt=Ttab(itable)
 
      it1=dint((xlt-Ttab(1))/(Ttab(itable)-Ttab(1))*
     &     dfloat(itable-1)+1.)
      it2=it1+1
      if (it2.gt.itable) it2=it1
      
      ip1=dint((xlp-Ptab(1))/(Ptab(itable)-Ptab(1))*
     &     dfloat(itable-1)+1.)
      ip2=ip1+1
      if (ip2.gt.itable) ip2=ip1
      
      l=1

      do it=it1,it2
         i=1
         do ip=ip1,ip2
            o1(i)=XMMWtable(ip,it)
            i=i+1
C           P
         enddo 
         
         if (ip2.ne.ip1) then
            o2(l)=o1(1)+(o1(2)-o1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
         else
            o2(l)=o1(1)
         endif
         
         l=l+1
c	 T
      enddo 

      if (it2.ne.it1) then
         xmmw=o2(1)+(o2(2)-o2(1))*(xlt-Ttab(it1))/(Ttab(it2)-Ttab(it1))
      else
         xmmw=o2(1)
      endif

      RETURN
      END

      
      SUBROUTINE DALESSIO(TK,Pcgs,Ors,Opl,Oab,Oas)

      implicit real*8 (a-h,o-z)
      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'     

C     By courtesy of Paola D'Alessio (2002)  
c     limits in T&P: 10**0.5 to 1.e7 K, 1.e-12 to 1.e9 din/cm^2
c     tables in log T, log P, 
      real*8 o1(2),o2(2),p1(2),p2(2),g1(2),g2(2),s1(2),s2(2)
      real*8 TK,Pcgs,Ors,Opl,Oab,Oas

      Ors=0.0
      Opl=0.0
      Oab=0.0
      Oas=0.0
            
      xlt=dlog10(TK)
      xlp=dlog10(Pcgs)

c     Out of bounds
      
      if (xlt.gt.Ttab(itable)) xlt=Ttab(itable)
      if (xlp.lt.Ptab(1))xlp=Ptab(1)
      if (xlp.gt.Ptab(itable))xlp=Ptab(itable)
      if (xlt.lt.Ttab(1)) xlt=Ttab(1)
      
      it1=dint((xlt-Ttab(1))/(Ttab(itable)-Ttab(1))*
     &     dfloat(itable-1)+1.)
      it2=it1+1
      
      if (it2.gt.itable) it2=it1
      
      ip1=dint((xlp-Ptab(1))/(Ptab(itable)-Ptab(1))*
     &     dfloat(itable-1)+1.)
      ip2=ip1+1

      if (ip2.gt.itable) ip2=ip1
      

      l=1
      do it=it1,it2
         j=1
         do ip=ip1,ip2
            o1(j)=ROStable(ip,it)
            p1(j)=PLKtable(ip,it)
            g1(j)=IRRtable(ip,it)
            s1(j)=SCAtable(ip,it)
            j=j+1
c	    P
         enddo 
         
         if (ip2.ne.ip1) then
            o2(l)=o1(1)+(o1(2)-o1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
            p2(l)=p1(1)+(p1(2)-p1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
            g2(l)=g1(1)+(g1(2)-g1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
            s2(l)=s1(1)+(s1(2)-s1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
         else
            o2(l)=o1(1)
            p2(l)=p1(1)
            g2(l)=g1(1)
            s2(l)=s1(1)
         endif
         
         l=l+1
c	 T
      enddo 
      
      if (it2.ne.it1) then
         Ors=o2(1)+(o2(2)-o2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
         Opl=p2(1)+(p2(2)-p2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
         Oab=g2(1)+(g2(2)-g2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
         Oas=s2(1)+(s2(2)-s2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
      else
         Ors=o2(1)
         Opl=p2(1)
         Oab=g2(1)
         Oas=s2(1)
      endif
      
      Ors=10.**Ors
      Opl=10.**Opl
      Oab=10.**Oab
      Oas=10.**Oas

      if (xlt.ge.6.) then
         Ors=0.3429
      endif

      RETURN
      END
      
      FUNCTION FLUXLMDF(T0,T1,Ors,Ors1,Rho0,Rho1,jpos,dim)
      
      IMPLICIT real*8 (a-h,o-z)
      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'     
      INCLUDE 'units.h'

      real*8 T0,T1,Ors,Ors1,Rho0,Rho1,y,beta,ddim,dtau
      integer jpos,dim

      if (dim.eq.1) then
         ddim=rof3n
      else if (dim.eq.2) then
         ddim=zof3n
      else if (dim.eq.3) then
         ddim=rhf(jpos)*dtheta
      endif

      dtau = (ors+ors1)*sconv*(rho0+rho1)*ddim*.25d0
! if optical depth limiting is required, uncomment
!      if (dtau > 10.) dtau = 10.
      dtau = dtau*4

C     this y-beta thing is the flux limiter
      
      y=32.0*abs(T0-T1)/
     &     ( dtau*(T0+T1) )

      beta=(2.0 + y) / (6.0 + 3.0*y + y**2)
c      beta= 1.0/3.0  
c ... Turn the flux limiter off

      Fluxlmdf=-( 8.0*sigma*beta*(T0+T1)**3*(T0-T1) ) /
     &     ( dtau*(Tconv**4) )      
      
      RETURN
      END
c
      SUBROUTINE HEIGHTS ! THIS IS KEPT HERE FOR BONUS POINTS.
          ! IF YOU FIND IT, YOU GET THE ARCHAEOLOGY AWARD.
      
      IMPLICIT real*8 (a-h,o-z)
      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'     
      real*8 limiter
      limiter = den*phylim 
      do j=2,jmax
         do k=2,kmax
            if ((rho(j,k,1).ge.limiter).and.
     &           (rho(j,k+1,1).lt.limiter))print*,j,'klim ',k
            if ((tau(j,k,1,1).ge.(2.0/3.0)).and.
     &           (tau(j,k+1,1,1).lt.(2.0/3.0)))print*,j,' kross ',k
            if ((tau(j,k,1,2).ge.(2.0/3.0)).and.
     &           (tau(j,k+1,1,2).lt.(2.0/3.0)))print*,j,' kplck ',k
            if ((tau(j,k,1,3).ge.(2.0/3.0)).and.
     &           (tau(j,k+1,1,3).lt.(2.0/3.0)))print*,j,'kirr ',k
         enddo
      enddo
C      do k=2,kmax
C         do j=2,jmax
C            if ((rho(j,k,1).ge.limiter).and.
C     &           (rho(j-1,k,1).lt.limiter))print*,j,'jlim ',k 
C         enddo
C      enddo

      RETURN
      END      

      subroutine TempFind() ! ACB find the temperature based on the energy/particle 
                            ! temperature relation. The temperature is interpolated
                            ! from the temperature table temptable, which is 
                            ! calculated in setup based on X, Y, Z composition
                            ! and other assumptions.
         implicit none
         include 'hydroparam.h'
         include 'globals.h'
         include 'units.h'
        
         integer J,K,L,I 
         real*8 limiter
         real*8 eng,dummy

         limiter = den*phylim 

!$OMP DO SCHEDULE(STATIC) REDUCTION(+:eflufftot)
         do L = 1, LMAX
           do K = 1, KMAX2
             do J = 1, JMAX2

               eng = 0.d0
            
               if (rho(J,K,L) >= limiter) then

                 eng = eps(J,K,L)/rho(J,K,L)*engconv

                 if (eng < engtable(1)) then

                   tempk(J,K,L) = temptable(1)/engtable(1)*eng
                   eps(J,K,L) = engtable(1)/engconv*rho(J,K,L)

                   I = 2

                   gamma1(J,K,L) = gammatable(I-1) +
     & (gammatable(I)-gammatable(I-1))/(temptable(I)-temptable(I-1))
     & * (tempk(J,K,L)-temptable(I-1))

                 else

                   I = 2
                   do while ( engtable(I) < eng .and. I < TTABLE )
                    I = I+1
                   enddo

                   tempk(J,K,L) = temptable(I-1) +
     & (temptable(I)-temptable(I-1))/(engtable(I)-engtable(I-1))*
     & (eng - engtable(I-1))

                   gamma1(J,K,L) = gammatable(I-1) +
     & (gammatable(I)-gammatable(I-1))/(temptable(I)-temptable(I-1))
     & * (tempk(J,K,L)-temptable(I-1))


                 endif

               else

                 I = 2
                 do while ( temptable(I) < tbgrnd .and. I < TTABLE )
                  I = I+1
                 enddo

                 tempk(J,K,L) = tbgrnd

                 dummy   = engtable(I-1) + (engtable(I)-engtable(I-1))
     &                   / (temptable(I)-temptable(I-1))
     &                   * (tbgrnd-temptable(I-1))

                 dummy   = dummy*rho(J,K,L)*rhoconv/pconv

                 eflufftot = eflufftot + (eps(J,K,L)-dummy)*rhf(J)
     &                     * rof3n*dtheta*zof3n
                 eps(J,K,L) = dummy

                 gamma1(J,K,L) = gammatable(I-1) +
     & (gammatable(I)-gammatable(I-1))/(temptable(I)-temptable(I-1))
     & * (tempk(J,K,L)-temptable(I-1))


               endif

             enddo
           enddo
         enddo
!$OMP END DO

      return

      end subroutine TempFind

      subroutine TempFindSpec(eps,rho,temp)! like TempFind, but for a single cell. 

         implicit none

         include 'hydroparam.h'

         integer I
         real*8 eng,eps,rho,temp
         real*8 Msyscgs,PKcgs,Tconv,Sconv,Dconv,Pconv,sigma,rhoconv,
     &       engconv,bkmpcode
          COMMON /CONVERT/
     &     Msyscgs,PKcgs,Tconv,Sconv,
     &     Dconv,Pconv,sigma,rhoconv,engconv,bkmpcode
         real*8 :: temptable,engtable,gammatable,muc
         common /engtables/temptable(TTABLE),engtable(TTABLE),
     &           gammatable(TTABLE),muc



         eng = eps/rho*engconv

         if (eng < engtable(1)) then

            temp = temptable(1)/engtable(1)*eng
            eps = engtable(1)/engconv*rho

         else

            I = 2
            do while ( engtable(I) < eng .and. I < TTABLE )
              I = I+1
            enddo

            temp = temptable(I-1) +
     & (temptable(I)-temptable(I-1))/(engtable(I)-engtable(I-1))*
     & (eng - engtable(I-1))

         endif

      return

      end subroutine TempFindSpec
