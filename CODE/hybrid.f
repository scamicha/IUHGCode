      subroutine Hybrid()
      implicit none
      include 'hydroparam.h'
      include 'units.h'
      include 'globals.h'

      real(KIND=8) ::  dummy , dumb 
      real(KIND=8) :: bcint, limiter,aconst,bconst,cconst, mu
      real(KIND=8) :: FluxLmDf,sum

      integer :: J,K,L,I,LM,LP

      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),absfrac,
     &      oplck_env(JMAX2,KMAX2,LMAX)

      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum,oplck_env

      mu = 1.d0/sqrt(3.d0)
      efl = 0.
      limiter = den*phylim 
      LO=0
      CCOUNTER = 0
      HCOUNTER = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(LM,dumb,cconst,bconst,aconst,I,  &
!$OMP& LP,J,K,L)                                                        &
!$OMP&  SHARED(ccounter,hcounter,LO,delt,dtheta,rof3n,limiter,zof3n,    &
!$OMP& sconv,sigma,pi,mu,tconv) REDUCTION(+:efl)
!$OMP DO SCHEDULE(STATIC)
      do L =  1, LMAX 
       do K = 1, KMAX2
         do J = 1, JMAX2
           tempk(J,K,L) = tempk(J,K,L)/tconv ! convert temps to code units
           tau(J,K,L,1) = tau(J,K,L,1)/mu ! convert to ray optical depths
           tau(J,K,L,4) = 0.
           temporary(J,K,L) = 0.d0
           divflux(J,K,L) = 0.d0
           intensity_z(J,K,L)=0.d0
           intensity_in_z(J,K,L)=0.d0
           int_temp(J,K,L) = 0.d0
           dtau_z(J,K,L) = 0.d0
         enddo
       enddo
       do J = 1, JMAX2
         init_int_in(J,L) = (tenvk/tconv)**4/pi*sigma
         KFITA(J,L) = 2
       enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC) ! set optical depths of positive plane in l_tau
      do L = 1, LMAX 
        do K = KMAX1, 2, -1
          do J = JMIN, JMAX

            ! interpolate between rosseland and planck means.

            if (l_tau_z(J,2,L) > 0.) then
              dtau_z(J,K,L) = rho(J,K,L)*sconv*zof3n *
     &          (l_tau_z(J,2,L)*oross(J,K,L) + oplck(J,K,L)/
     &           l_tau_z(J,2,L)) / (mu*(l_tau_z(J,2,L) +
     &           1.d0/l_tau_z(J,2,L)))
            else
              dtau_z(J,K,L) = 0.
            endif

! only set the following if you need to limit the opacties for some reason
!            if (dtau_z(J,K,L) > 10./mu) dtau_z(J,K,L)=10./mu
!            dtau_z(J,K,L) = rho(J,K,L)*oross(J,K,L)*sconv*zof3n/mu

            tau(J,K,L,4)   = tau(J,K+1,L,4) + 
     &              rho(J,K,L)*oplck_env(J,K,L)*sconv*zof3n/mu 
            l_tau_z(J,K,L) = l_tau_z(J,K+1,L) + dtau_z(J,K,L)

            sfunc(J,K,L) = sigma*tempk(J,K,L)**4/pi

          enddo
        enddo
        do J = JMIN, JMAX
           sfunc(J,1,L) = sfunc(J,2,L)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)! set optical depths of positive plane in l_tau
      do L = 1, LMAX
        do K = KMAX, 3, -1
          do J = JMIN, JMAX

            if ( rho(J,K,L) >= limiter 
     &          .and. rho(J,K+1,L) >= limiter 
     &          .and. rho(J,K-1,L) >= limiter ) then
       
              dsdt(J,K,L) =  2.d0*sigma*tempk(J,K,L)**3 * 
     &            (tempk(J,K-1,L)-tempk(J,K+1,L))/(pi*dtau_z(J,K,L) )


            else
              dsdt(J,K,L) = 0.d0
            endif

          enddo
        enddo

        do J = JMIN, JMAX
          if ( rho(J,2,L) >= limiter 
     &          .and. rho(J,3,L) >= limiter 
     &          .and. rho(J,4,L) >= limiter ) then

            dsdt(J,2,L) = sigma*tempk(J,2,L)**3 * 
     &        ( 8.d0*tempk(J,2,L)-7.d0*tempk(J,3,L)-tempk(J,4,L)) 
     &        / (3.d0*pi*dtau_z(J,2,L))


          else

            dsdt(J,2,L) = 0.d0

          endif

          dsdt(J,1,L)=-dsdt(J,2,L)

        enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        do K = KMAX, 2, -1
          do J = JMIN, JMAX 

            if ( rho(J,K,L) >= limiter 
     &         .and. rho(J,K+1,L) >= limiter 
     &         .and. rho(J,K-1,L) >= limiter 
     &         .and. dsdt(J,K,L) /= 0.) then

              ddsdtt(J,K,L) = ( dsdt(J,K-1,L)-dsdt(J,K+1,L) ) 
     &                      / (2.d0*dtau_z(J,K,L) )


! if the second derivative becomes important, kill it. This is a stability 
! consideration
 
          if(abs(ddsdtt(J,K,L)) >= 2*abs(dsdt(J,K,L)))then
               ddsdtt(J,K,L) = 0.
            endif

                     
            else
              ddsdtt(J,K,L) = 0.d0
            endif

          enddo
        enddo
  
        do J = JMIN, JMAX

          if ( rho(J,2,L) >= limiter 
     &          .and. rho(J,3,L) >= limiter 
     &          .and. rho(J,4,L) >= limiter 
     &         .and. dsdt(J,2,L) /= 0.) then

            ddsdtt(J,2,L) = 0.5d0*ddsdtt(J,2,L) 
     &        + (dsdt(J,2,L)-dsdt(J,3,L))/(2.d0*dtau_z(J,2,L))


! if the second derivative becomes important, kill it. This is a stability
! consideration

        if(abs(ddsdtt(J,2,L)) >= 2*abs(dsdt(J,2,L)))then
               ddsdtt(J,2,L) = 0.
            endif


          else
   
            ddsdtt(J,2,L) = 0.d0

          endif

        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX
          do J = JMIN, JMAX                                           
            if ( tau(J,K,L,1) >= 1. .and. rho(J,K,L)>=limiter ) then
            KFITA(J,L) = K

                radflux(J,K,L,1) = FluxLmDf(tempk(J,K,L)*tconv,
     &                  tempk(J-1,K,L)*tconv,
     &                  dtau_z(J,K,L)*mu/(zof3n*sconv*rho(J,K,L)), 
     &                  dtau_z(J-1,K,L)*mu/(zof3n*sconv*rho(J-1,K,L)),
     &                  rho(J,K,L),rho(J-1,K,L),J,1)

                radflux(J,K,L,3) = FluxLmDf(tempk(J,K,LP)*tconv,
     &                  tempk(J,K,L)*tconv,
     &                  dtau_z(J,K,LP)*mu/(zof3n*sconv*rho(J,K,LP)), 
     &                  dtau_z(J,K,L)*mu/(zof3n*sconv*rho(J,K,L)),
     &                  rho(J,K,LP),rho(J,K,L),J,3)

            endif
          enddo
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        do K = KMAX1, 2, -1
          do J = JMIN, JMAX  
   
              I = K+1

              aconst  = l_tau_z(J,I-1,L)**2
     &                - (2.d0*l_tau_z(J,I-1,L)-2.d0)
              aconst  = aconst - (l_tau_z(J,I,L)**2
     &                - (2.d0*l_tau_z(J,I,L)-2.d0)) 
     &                * exp (-dtau_z(J,I-1,L))

              aconst  = aconst * 0.5d0*ddsdtt(J,I-1,L)

              bconst  =  l_tau_z(J,I-1,L) - 1.d0 
     &                - ( l_tau_z(J,I,L) - 1.d0 ) 
     &                * exp (-dtau_z(J,I-1,L))

              bconst  = bconst * (dsdt(J,I-1,L) 
     &                - ddsdtt(J,I-1,L)*tau(J,I-1,L,1) )

              cconst  = ( sfunc(J,I-1,L) 
     &                - tau(J,I-1,L,1)*dsdt(J,I-1,L) 
     &                + 0.5d0*tau(J,I-1,L,1)**2*ddsdtt(J,I-1,L) )
     &                * ( 1.d0 - exp (-dtau_z(J,I-1,L))) 

              dumb    = cconst + bconst + aconst

              if ( dumb < 0.d0 ) dumb = sfunc(J,I-1,L) 
     &                  * ( 1.d0 - exp ( -dtau_z(J,I-1,L) ) )

            intensity_in_z(J,K,L) = dumb + intensity_in_z(J,I,L)
     &                            * exp(-dtau_z(J,K,L))

               
          enddo
        enddo
        do K = KMAX1,2,-1
          do J = JMIN, JMAX
            intensity_in_z(J,K,L) = intensity_in_z(J,K,L) 
     &              + init_int_in(J,L)*exp(-tau(J,K,L,4))

          enddo
        enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC) ! calculate intensity for each cell
      do L = 1, LMAX
        do J = JMIN, JMAX 
          intensity_z(J,2,L) = intensity_in_z(J,2,L)
        enddo
        do K = 3 , KMAX1
          do J = JMIN, JMAX

            I = K-1

              aconst  = l_tau_z(J,I+1,L)**2
     &                + (2.d0*l_tau_z(J,I+1,L)+2.d0)
              aconst  = aconst - (l_tau_z(J,I,L)**2
     &                + (2.d0*l_tau_z(J,I,L)+2.d0)) 
     &                * exp (-dtau_z(J,I,L))
              aconst  = aconst * 0.5d0*ddsdtt(J,I,L)

              bconst  =  l_tau_z(J,I+1,L) + 1.d0 
     &                - ( l_tau_z(J,I,L) + 1.d0 ) 
     &                * exp (-dtau_z(J,I,L))
              bconst  = bconst * (dsdt(J,I,L) 
     &                - ddsdtt(J,I,L)*tau(J,I,L,1) )

              cconst  = ( sfunc(J,I,L) - tau(J,I,L,1)*dsdt(J,I,L) 
     &                + 0.5d0*tau(J,I,L,1)**2*ddsdtt(J,I,L) )
     &                * ( 1.d0 - exp ( -dtau_z(J,I,L)))

              dumb    = aconst + bconst + cconst

              if ( dumb < 0.d0 ) dumb = sfunc(J,I,L) 
     &                   * ( 1.d0 - exp ( -dtau_z(J,I,L) ) )

            intensity_z(J,K,L) = intensity_z(J,I,L)
     &                         * exp(-dtau_z(J,I,L)) + dumb

           enddo
         enddo
         do K = 2, KMAX1
           do J = JMIN, JMAX
             int_temp(J,K,L) = 2.d0*pi*mu*(intensity_z(J,K,L) 
     &                       - intensity_in_z(J,K,L))
           enddo

        enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX
        do J = JMIN, JMAX
          teffk(J,L)= tconv * ( SQRT(SQRT( pi
     & * intensity_z(J,KMAX1,L) / sigma )) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX 
          do J = JMIN, JMAX
            if ( tau(J,K,L,1) >= 1.  ) then
            if (rho(J,K,L) .ge. limiter) then ! must be optically thick regime. 
                                               ! heat using divF
                       ! set BC for r and theta

                if ( tau(J+1,K,L,1) < 1. ) then

                     radflux(J+1,K,L,1)=abs(int_temp(J,KFITA(J,L)+1,L))
                     efl  = efl  +
     &                   radflux(J+1,K,L,1)*r(J+1)*dtheta*zof3n*delt                    
         
                endif

                if ( tau(J-1,K,L,1) < 1. ) then

                     radflux(J,K,L,1)=-abs(int_temp(J,KFITA(J,L)+1,L))
                     efl  = efl  +
     &                   abs(radflux(J,K,L,1))*r(J)*dtheta*zof3n*delt

                endif

                if ( tau(J,K,LM,1) < 1. ) then

                     radflux(J,K,LM,3)=-abs(int_temp(J,KFITA(J,L)+1,L))
                     efl  = efl  +
     &                   abs(radflux(J,K,LM,3))*rof3n*zof3n*delt
        
                endif

                if ( tau(J,K,LP,1) < 1. ) then

                     radflux(J,K,L,3)=abs(int_temp(J,KFITA(J,L)+1,L))
                     efl  = efl  +
     &                   radflux(J,K,L,3)*rof3n*zof3n*delt

                endif

            endif
            endif ! end cooling if structure
            enddo
          enddo
        enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX 
          do J = JMIN, JMAX
            if ( tau(J,K,L,1) >= 1.  ) then
            if (rho(J,K,L) .ge. limiter) then ! must be optically thick regime. 
                                               ! heat using divF
 
                temporary(J,K,L) = ((r(J+1)*radflux(J+1,K,L,1))
     &                        - (r(J)*radflux(J,K,L,1)))
     &                        / (rhf(J)*rof3n)        
     &                        + (radflux(J,K,L,3)-radflux(J,K,LM,3))
     &                        / (rhf(J)*dtheta)

            endif
            endif
          enddo
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)  ! calculate divflux
      do L = 1, LMAX 
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX
          do J = JMIN, JMAX
               
          if (rho(J,K,L) >= limiter ) then
             divflux(J,K,L) =                                          
     &         (int_temp(J,K+1,L)-int_temp(J,K,L))/zof3n      
     &       + temporary(J,K,L)
          endif

          enddo
        enddo
      enddo
!$OMP END DO NOWAIT

C$OMP DO SCHEDULE(STATIC) REDUCTION(+:CCOUNTER,HCOUNTER,LO) ! calculate divflux
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
!$OMP END DO
!$OMP MASTER
       dumb = dble(hcounter)/dble(LO)*100.
       if (dumb > 1.)  write(6,"(A,1pe12.3)")
     & "->Hybrid Diag: percent relevant h-limited cells ",
     & dumb
       dumb = dble(ccounter)/dble(LO)*100.
       if (dumb > 1.)  write(6,"(A,1pe12.3)")
     & "->Hybrid Diag: percent relevant c-limited cells ",
     & dumb
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC)
      do L =  1, LMAX 
       do K = 1, KMAX2
         do J = 1, JMAX2
           tempk(J,K,L) = tempk(J,K,L)*tconv ! convert temps to Kelvin
           tau(J,K,L,1) = tau(J,K,L,1)*mu
         enddo
       enddo
      enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      return

      end subroutine Hybrid
  
