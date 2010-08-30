      subroutine init()
      
      INCLUDE 'hydroparam_mpi.h'
      INCLUDE 'globals_mpi.h'

      INTEGER :: jreq

      Msyscgs=Mstar*Msuncgs*(1.0+(tmassini/(1.0-tmassini)))
      
      PKcgs = ( (Rdiskau*AUcgs/r(jreq))**(3.0-xn) * 
     &     Gcgs**xn * Msyscgs**(xn-1.0) )**(1.0/xn)
      
      sigma = (sigmacgs/(bkmpcgs**4)) * 
     &     ( Gcgs**(15.0/2.0 - 3.0*xn) * 
     &     Msyscgs**(5.0 - 2.0*xn) * PKcgs**(xn/2.0)
     &     )**(1.0/(3.0-xn)) 
      
      Tconv = (1.0/bkmpcgs) *( Gcgs**3 * Msyscgs**2 * PKcgs**(-xn)
     &     )**(1.0/(3.0-xn))
      
      Sconv = (Gcgs**(2.0*xn) * Msyscgs**(1.0+xn) * PKcgs**(-2
     &     .0*xn))**(1.0/(3.0-xn))
      
      Pconv = (Gcgs**(3.0+3.0*xn) * Msyscgs**(2.0+2.0*xn) * PKcgs**(-4
     &     .0*xn))**(1.0/(3.0-xn))

      Dconv = (Gcgs**(-xn) * Msyscgs**(1.0-xn) * PKcgs**(xn))**(1.0/(3.0
     &     -xn))

      rhoconv = (Gcgs**3*Msyscgs**2*PKcgs**(-3))
     &        **(1.d0/(3.d0*gamma-4.d0))

      engconv = pconv/rhoconv

      bkmpcode = 1.d0

      return

      end subroutine init

! Written by A. C. Boley (updated 12 June 2007).
! See Pathria for questions.  This uses E = NkT^2 d ln Z / d T to calculate internal energies.
! However, the zero point energies are subtracted out (relevant for orthohydrogen
! and for the vibrational states).  The 16.78 in the metals
! term is from cameron 1968. This takes into account rotational states for 
! hydrogen and 1 vibrational state. zp is the parahydrogen partition function
! and dzpdt and ddzpdtt are its derivatives, *e for equilibrium, and *o for ortho.
      subroutine Initengtable()
        implicit none
        include 'hydroparam_mpi.h'
        include 'globals_mpi.h'
        
        integer J,I,jreq
        
        real*8, parameter :: b = 85.4d0, vib = 5987.d0 ! see Draine et al. (1983)
        real*8 f1,f2,f3,f4,dummy0,dummy1

        real*8 zp(TTABLE),dzpdt(TTABLE)
        real*8 zo(TTABLE),dzodt(TTABLE)
        real*8 ze(TTABLE),dzedt(TTABLE)
        real*8 cvsm(TTABLE)


        do I = 1, TTABLE
          temptable(I) = dble(I-1)*5.d0 + tbgrnd ! for TTABLE = 400 takes to T = 1998
        enddo

        zp = 0.d0
        zo = 0.d0
        ze = 0.d0
        dzpdt = 0.d0
        dzodt = 0.d0
        dzedt = 0.d0

        do I = 1, TTABLE
         
          do J = 0, 24, 2 ! I found convergence with 24 terms.  You could use less if you want,
                          ! but why? If you are REALLY paranoid, you might increase the terms.

            f1 = dble(2*J+1)
            f2 = dble(J*(J+1))
            f3 = dble(2*(J+1)+1)
            f4 = dble((J+1)*(J+2))

            zp(I) = zp(I) +f1*exp(-f2*b/temptable(I))

            zo(I) = zo(I) +3.d0*f3*exp(-f4*b/temptable(I))

            dzpdt(I) = dzpdt(I) + f1*f2*b*exp(-f2*b/temptable(I))
     &               / temptable(I)**2 

            dzodt(I) = dzodt(I) + 3.d0*f3*f4*b
     &               * exp(-f4*b/temptable(I))
     &               / temptable(I)**2 

          enddo

            ze(I)  = zp(I) + zo(I)

            dzedt(I) = dzpdt(I) + dzodt(I) 

        enddo

        select case(H2STAT)

        case(0)

        do I = 1, TTABLE

          engtable(I) = bkmpcgs*0.5d0*temptable(I)*(1.5d0
     &                +
     &                  temptable(I)*dzpdt(I)/zp(I)
     &                +
     &                  vib/temptable(I)
     &                * (exp(-vib/temptable(I))
     &                / (1.d0-exp(-vib/temptable(I)))))


           engtable(I) = engtable(I)*xabun + (bkmpcgs*yabun*.25d0 +
     &                   bkmpcgs*zabun/16.78d0)*temptable(I)*1.5d0                 !from boss 1984 3.118d7 and 1.247d8
                                                              ! z based on mu_z = 16.78

        enddo

        case(1)

        do I = 1, TTABLE

          engtable(I) = bkmpcgs*0.5d0*temptable(I)*(1.5d0
     &                +
     &                 temptable(I)*(dzodt(I)/zo(I)
     &                - 2.d0*b/temptable(I)**2)
     &                +
     &                  vib/temptable(I)
     &                * (exp(-vib/temptable(I))
     &                / (1.d0-exp(-vib/temptable(I)))))


           engtable(I) = engtable(I)*xabun + (bkmpcgs*yabun*.25d0 +
     &                   bkmpcgs*zabun/16.78d0)*temptable(I)*1.5d0                 !from boss 1984 3.118d7 and 1.247d8
                                                              ! z based on mu_z = 16.78

        enddo

        case(2)

        do I = 1, TTABLE

          engtable(I) = bkmpcgs*0.5d0*temptable(I)*(1.5d0
     &                +
     &                  temptable(I)*dzpdt(I)/zp(I)
     &                * ac/(ac + bc)
     &                +
     &                  temptable(I)*(dzodt(I)/zo(I)
     &                - 2.d0*b/temptable(I)**2)
     &                * bc/(ac + bc)
     &                +
     &                  vib/temptable(I)
     &                * (exp(-vib/temptable(I))
     &                / (1.d0-exp(-vib/temptable(I)))))

           engtable(I) = engtable(I)*xabun + (bkmpcgs*yabun*.25d0 +
     &                   bkmpcgs*zabun/16.78d0)*temptable(I)*1.5d0                 !from boss 1984 3.118d7 and 1.247d8
                                                              ! z based on mu_z = 16.78

        enddo

        case(3)

        do I = 1, TTABLE

          engtable(I) = bkmpcgs*0.5d0*temptable(I)*(1.5d0
     &                +
     &                  temptable(I)*dzedt(I)/ze(I)
     &                + 
     &                  vib/temptable(I)
     &                * (exp(-vib/temptable(I))
     &                / (1.d0-exp(-vib/temptable(I)))))

           engtable(I) = engtable(I)*xabun + (bkmpcgs*yabun*.25d0 +
     &                   bkmpcgs*zabun/16.78d0)*temptable(I)*1.5d0                 !from boss 1984 3.118d7 and 1.247d8
                                                              ! z based on mu_z = 16.78

        enddo

        case(4)

        do I = 1, TTABLE

           engtable(I) = bkmpcgs*temptable(I)/(2.d0*(gamma-1.d0))

           engtable(I) = engtable(I)*xabun + (bkmpcgs*yabun*.25d0 +
     &                   bkmpcgs*zabun/16.78d0)*temptable(I)*1.5d0                 !from boss 1984 3.118d7 and 1.247d8
                                                              ! z based on mu_z = 16.78
 
        enddo

        end select


        muc = 1.d0 / (xabun*.5d0 + yabun*.25d0 + zabun/16.78d0)

        cvsm(1) = (engtable(1)-engtable(2))/(temptable(1)-temptable(2))
        cvsm(TTABLE) = (engtable(TTABLE)-engtable(TTABLE-1))
     &             / (temptable(TTABLE)-temptable(TTABLE-1))

        do I = 2, TTABLE-1
 
          cvsm(I) = (engtable(I-1)-engtable(I+1))
     &            / (temptable(I-1)-temptable(I+1))

        enddo

        do I = 1, TTABLE

          gammatable(I) = 1.d0 + bkmpcgs/(muc*cvsm(I))

        enddo

        open(unit=456,file="engtable.dat")
        do I = 1, TTABLE
         write(456, "(5(1X,1pe15.8))")temptable(I),
     &              1.+1./(cvsm(I)*muc/bkmpcgs),
     &              engtable(I)*muc/bkmpcgs,
     &              dzodt(I)/zo(I)*temptable(I)**2,
     &              gammatable(I)
        enddo
        close(456)


        print "(A,1pe10.3)","INITENGTABLE-> muc is ", muc
      
      return

      end subroutine Initengtable

      subroutine TempFind() ! ACB find the temperature based on the energy/particle 
                            ! temperature relation. The temperature is interpolated
                            ! from the temperature table temptable, which is 
                            ! calculated in setup based on X, Y, Z composition
                            ! and other assumptions.
         implicit none
         include 'hydroparam_mpi.h'
         include 'globals_mpi.h'
        
         integer J,K,L,I ,jreq
         real*8 limiter
         real*8 eng,dummy
         real*8 DTHETA, PI, GRAV, BGDEN
         COMMON /BLOK6/DTHETA,PI,GRAV,BGDEN


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

         include 'hydroparam_mpi.h'

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
