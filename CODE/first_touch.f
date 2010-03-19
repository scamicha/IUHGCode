      subroutine first_touch()
!
!     touch every field array  to organize data placement
!
      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'

      REAL*8 mass,mass1,
     &       VINV(JMAX2),
     &       VINVH(JMAX2),
     &       RD(JMAX2),
     &       RHD(JMAX2),
     &       VR(JMAX2,KMAX2,LMAX),
     &       VZ(JMAX2,KMAX2,LMAX),
     &       VT(JMAX2,KMAX2,LMAX),
     &       VLOP(JMAX2,KMAX2,LMAX),
     &       FR(JMAX2,KMAX2,LMAX),
     &       FZ(JMAX2,KMAX2,LMAX),
     &       FT(JMAX2,KMAX2,LMAX)

      COMMON  /FLUXPRIVATE/VINV,VINH,RD,RHD
      COMMON  /FLUXSHARED/VR,VZ,VT,VLOP,FR,FZ,FT
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO K=1,KMAX2
          DO J=1,JMAX2
            starphi(J,K,L)=0.0
            tinphi(J,K,L)=0.0        
            S(J,K,L)=0.0        
            T(J,K,L)=0.0        
            A(J,K,L)=0.0        
            U(J,K,L)=0.0        
            W(J,K,L)=0.0        
            JN(J,K,L)=0.0        
            OMEGA(J,K,L)=0.0        
            P(J,K,L)=0.0
            CV(J,K,L)=0.0
            EPS(J,K,L)=0.0
            RHO(J,K,L)=0.0
            PHI(J,K,L)=0.0
            QRR(J,K,L)=0.0
            QZZ(J,K,L)=0.0
            QTT(J,K,L)=0.0
            HGAMMA(J,K,L)=0.0
            LAMBDA(J,K,L)=0.0
            do i=1,4
              TAU(J,K,L,i)=0.0
            enddo
            TEMPK(J,K,L)=0.0
            DIVFLUX(J,K,L)=0.0
            do i=1,3
              RADFLUX(J,K,L,i)=0.0
            enddo
            temporary(J,K,L)=0.0
            dsdt(J,K,L)=0.0
            sfunc(J,K,L)=0.0
            l_tau_z(J,K,L)=0.0
            dtau_z(J,K,L)=0.0
            intensity_in_z(J,K,L)=0.0
            intensity_z(J,K,L)=0.0
            ddsdtt(J,K,L)=0.0
            int_temp(J,K,L)=0.0
            VR(J,K,L)=0.0
            VZ(J,K,L)=0.0
            VT(J,K,L)=0.0
            VLOP(J,K,L)=0.0
            FR(J,K,L)=0.0
            FZ(J,K,L)=0.0
            FT(J,K,L)=0.0
          ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO nowait

!$OMP DO SCHEDULE(STATIC)
      DO J=1,JMAX2
         VINV(J)=0.0
         VINVH(J)=0.0
         RD(J)=0.0
         RHD(J)=0.0         
      ENDDO
!$OMP ENDDO nowait      

!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO J=1,JMAX2
          TEFFK(J,L)=0.0
          TPHK(J,L)=0.0
          SURFCGS(J,L)=0.0
          init_int_in(J,L)=0.0
          KFITA(J,L)=0.0
        ENDDO
      ENDDO
!$OMP ENDDO nowait
!$OMP END PARALLEL
      rpstar = 0.
      phi_star=0.
      return
      END
