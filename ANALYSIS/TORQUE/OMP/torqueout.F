      subroutine torqueout(AMCOUNT)
      implicit none
      include 'hydroparam.h'
      include 'globals.h'

      real(KIND=8), dimension(JMAX2,KMAX2,LMAX) :: gforce_phi, vr, vphi,
     &                                              treyn3d,
     &                                             delvr,delvphi 

      real(KIND=8), dimension(JMAX2,LMAX) :: sigma_2d
      real(KIND=8), dimension(JMAX2,KMAX2) :: treyn2d,vcol

      real(KIND=8) :: 
     &  kephi,dummy0,dummy1,dummy2,Aomega,rconv, vconv, limiter,tc,
     &  deriv,eta, dummy, dummy3, dummy4, dummy5, dummy6,timediff

      real(KIND=8), save :: times
      real*8 DTHETA, PI, GRAV, BGDEN
      COMMON /BLOK6/DTHETA,PI,GRAV,BGDEN
       
 
      integer :: AMCOUNT,J, K, L, I,LM,LP,JREQ,TINITK,KZPOL,
     &           ENON,LMM,LPP,O,H

      mstar  = 1.d0/(1.d0-tmassini)
      rdisk  = 40.d0 
      gamma  = 1.666666d0
      Aomega = sqrt(6.67d-8*0.5d0*1.99d33)
      
      kconst = rdisk*mstar**(0.333333333333333)/(r(JREQ-2)*7.93d-3)
      econv = 2.24d48*mstar**(7.d0/3.d0)/kconst
      rconv = rdisk/r(JREQ-2)*1.49d13
      tconv = 3.55d3/mstar*kconst**(1.5d0)
      vconv = 3.34d2*mstar**(2.d0/3.d0)/sqrt(kconst)
!      print *, tconv


      rholmt = 1.d-12*den
      limiter= 1.d-7*den

      IF(AMCOUNT.EQ.0)THEN
        print *, ' RHOLMT is ',rholmt, ' THIS IS HARDCODED INTO TORQUEOU
     &T'
      ENDIF
 
      pi = acos(-1.d0)
      dtheta = 2.d0*pi/dble(LMAX)

      oav      = 0.d0
      vrav     = 0.d0
      sigma_2d = 0.d0
      vcol = 0.d0

      do L = 1, LMAX ! first we will define the velocities
                     ! and the gravitational phi force.
                     ! the velocities are defined on the 
                     ! cell faces, as done in the hydro code.
        LP = L+1
        LM = L-1
        if (L .eq. 1) LM = LMAX
        if (L .eq. LMAX) LP = 1
        do J = 2, JMAX1
          do K = 2, KMAX1
            
             if ( rho(J,K,L) > limiter ) then
                gforce_phi(J,K,L) = (phi(J,K,LP)-phi(J,K,LM))/(dtheta
     &               * 2.d0 * rhf(J) )
             else
                gforce_phi(J,K,L) = 0.d0
             endif
           
             IF (AMCOUNT.EQ.0)THEN
                if (rho(J,K,L).gt.rholmt.and.rho(J+1,K,L).gt.rholmt)then

                   vr(J,K,L)      = 2.d0*s(J,K,L)/(rhosave(J,K,L)+
     &                  rhosave(J-1,K,L))
                 
                   vphi(J,K,L)    = (a(J-1,K,L)+a(J,K,L))/
     &                  ((rhosave(J-1,K,L)+rhosave(J,K,L))*r(J))

                else

                   vr(J,K,L)   = 0.d0
                   vphi(J,K,L) = 0.d0
                 
                endif
             ENDIF

          enddo
       enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create midplane terms. this includes the omega, sound speed, and 
! sigma.  cc is for cell center, fc is for face center.  the mass
! up to a face is mass_fc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(AMCOUNT.EQ.0)THEN
      
         do J = 3, JMAX 
            kepler_fc(J) = sqrt((1.d0-tmassini)/r(J)**3)
            dummy  = 0.d0
            dummy0 = 0.d0
            dummy1 = 0.d0
            dummy2 = 0.d0
            dummy3 = 0.d0
            dummy4 = 0.d0
            dummy6 = 0.d0
            cs_sig(J)   = 0.d0
            omega_fc(J) = 0.d0

            do L = 1, LMAX
               do K = 2, KMAX

                  dummy0 = dummy0 + (rhosave(J,K,L) + rhosave(J-1,K,L))
     &                 * zof3n

                  sigma_2d(J,L) =  sigma_2d(J,L) + 
     &                 (rhosave(J,K,L)+rhosave(J-1,K,L))*zof3n

                  dummy6 = dummy6+(gamma*(gamma-1.d0))*(eps(J,K,L)+ 
     &                 eps(J-1,K,L) )*zof3n

                  if ( 0.5d0*( rhosave(J,K,L) + rhosave(J-1,K,L) ) >= 
     &                 limiter ) then
                     dummy4 = dummy4 + 2.d0 * s(J,K,L) 

                     dummy3 = dummy3 + ( a(J,K,L) + a(J-1,K,L) ) /
     &                    ( r(J)**2 )
                  endif

               enddo

               dummy1 = dummy1 + sqrt((gamma*(gamma-1.d0))*(eps(J,2,L) + 
     &              eps(J-1,2,L) )/( rhosave(J,2,L) + 
     &              rhosave(J-1,2,L) ) )

               dummy2 = dummy2 + ( a(J,2,L) + a(J-1,2,L) ) / 
     &              ( r(J)**2 * ( rhosave(J,2,L) + rhosave(J-1,2,L)))

!               dummy6 = dummy6 + ((gamma*(gamma-1.d0))*(eps(J,2,L) + 
!     &              eps(J-1,2,L) )/( rhosave(J,2,L) + 
!     &              rhosave(J-1,2,L) ) )*sigma_2d(J,L)


            enddo
           
            sigma_fc(J) = dummy0 / dble(LMAX)   

            mass_fc(J)  = sigma_fc(J)*pi*(r(J)**2-r(J-1)**2)
     &           + mass_fc(J-1)

            cs_fc(J)    = dummy1 / dble(LMAX)
            cs_sig(J)   = cs_sig(J) + dummy6 / dble(LMAX)

            omega_fc(J) = omega_fc(J) + dummy2 / dble(LMAX)

            oav(J) = dummy3/( dble(LMAX)* sigma_fc(J) )*zof3n ! dummy
            vrav(J)= dummy4/( dble(LMAX)* sigma_fc(J) )*zof3n ! dummy
 
            masse(J) = omega_fc(J)/(3.d0*pi*cs_sig(J)) 

         enddo


         do J = 3, JMAX1
            do K = 2, KMAX1
               do L = 1, LMAX

                  if (0.5d0*(rhosave(J,K,L) + rhosave(J-1,K,L)) >= 
     &                 limiter) then

                     vcol(J,K) = vcol(J,K) + 0.5d0*(rho(J-1,K,L)+
     &                    rho(J,K,L))* (vphi(J,K,L)-oav(J)*r(J))
     &                    * (vr(J,K,L)-vrav(J) ) ! velocity corelation function
                                                 ! but average it over phi
                  endif
               enddo         

               treyn2d(J,K) = vcol(J,K)/dble(LMAX)

            enddo
         enddo

         do J = 2, JMAX1
            do K = 2, KMAX1     ! integrate vertically

               treyn(J) = treyn(J) + 2.d0*treyn2d(J,K)*zof3n

            enddo
         enddo
      
      ENDIF

      do J = 3, JMAX
         gt1(J) = 0.d0
         do L = 1, LMAX
            do K = 2, KMAX      ! back to gravity do the integral

               gt1(J) = gt1(J) + rhf(J-1) * rho(J-1,K,L) * 
     &              gforce_phi(J-1,K,L)
               
            enddo
         enddo
      enddo

      return
      end

