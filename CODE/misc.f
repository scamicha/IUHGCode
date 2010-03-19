C***********************************************************************
      SUBROUTINE VELOCITY
      IMPLICIT real*8 (a-h,o-z)

      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'
      real*8 rhox(kmax2)
      save rhox
      integer jstart,LP,K

C...FROM MOMENTA FIND VELOCITIES

      if (jmin.gt.2) then
         jstart=jmin
      else
         jstart=2
      endif

!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            U(J,K,L)=2.*S(J,K,L)/(RHO(J,K,L)+RHO(J-1,K,L))
            W(J,K,L)=2.*T(J,K,L)/(RHO(J,K,L)+RHO(J,K-1,L))
            JN(J,K,L)=A(J,K,L)/RHO(J,K,L)
            OMEGA(J,K,L) = JN(J,K,L)/(RHF(J)**2)
          ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO nowait

!$OMP DO SCHEDULE(STATIC)
      DO K=2,KMAX2
        RHOX(K)=0.0
      ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RHOX)
        do L=1,lmax,1
          do K=2,kmax2,1
          RHOX(K)=RHOX(K)+rho(2,K,L)
          ENDDO
        ENDDO
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
        do K=2,kmax2,1    
        RHOX(K)=RHOX(K)/lmax
        ENDDO
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
        do L=1,lmax,1
          do K=2,kmax2,1    
          u(2,K,L)=s(2,K,L)/RHOX(K)
          ENDDO
        ENDDO
!$OMP ENDDO nowait

C...SET VELOCITIES AROUND Z-AXIS.
      if (jmin.gt.2) then
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO K=2,KMAX1
               do j=1,jmin1
                  U(1,K,L)    = 0.0
                  W(1,K,L)    = 0.0
                  OMEGA(1,K,L)= 0.0
                  JN(1,K,L)   = 0.0
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO NOWAIT        
      else  
!$OMP DO SCHEDULE(STATIC) 
         DO L=1,LMAX
            LP=1+MOD(L-1+LMAX/2,LMAX)
            DO K=2,KMAX1
               IF(L.LE.LMAX/2) U(2,K,L)=0.d0 !dkb- what's this for?
               U(1,K,L)    = -U(3,K,LP)
               W(1,K,L)    = W(2,K,LP)
               OMEGA(1,K,L)= OMEGA(2,K,LP)
               JN(1,K,L)   = JN(2,K,LP)
            ENDDO
         ENDDO
!$OMP END DO
      endif

C...SET VELOCITES BELOW THE EQUATORIAL PLANE.

      if (jmin.gt.2) then
         jstart=jmin
      else
         jstart=1
      endif

!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO J=jstart,JMAX1
          U(J,1,L)    = U(J,2,L)
          W(J,2,L)    = 0.d0
          W(J,1,L)    = -W(J,3,L)
          OMEGA(J,1,L)= OMEGA(J,2,L)
          JN(J,1,L)   = JN(J,2,L)
        ENDDO
      ENDDO
!$OMP END DO


      RETURN
      END
 
 

C***********************************************************************
      SUBROUTINE QCALC(TIM)
      IMPLICIT real*8 (a-h,o-z) 

      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'

      REAL*8 mass,mass1

      COMMON /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
      COMMON /TIMEST/INDX,ISOADI,ALLOW,DMAX,CHGMAX
      COMMON /COEFS/COEF(pot3JMAX2,pot3KMAX2,LMAX2,2)
      COMMON /LIMIT/SOUND


      CHARACTER np*2,tw*2
      CHARACTER tim*6,rhofile*40,outfile1*40
      CHARACTER coeffile*40,index*6,modefile*40,centfile*40  

C The following arrays are local to this subroutine, and thus should not 
C need to be placed in common. But the OpenMP version of the code may fail
C unless they are.
      real*8 qphi(jmax2,lmax),
     &       rsig(jmax2,lmax),
     &       kaphi(jmax2,lmax),
     &       sspd(jmax2,lmax),
     &       dphidr(jmax2,lmax)

      integer zlim(jmax2,lmax),jstart
      
      character outfile*40
      COMMON /QCALCPRIVATE/QPHI,RSIG,KAPHI,sspd,DPHIDR,ZLIM
C$OMP THREADPRIVATE(/QCALCPRIVATE/)
      
      if (jmin.gt.2) then
         jstart=jmin
      else
         jstart=2
      endif   

      print*,'ENTER QCALC'

      cutoff=1.e-10
      print*,'cut',cutoff
      dr=rof3n
 
      print*,'ENTER first loop'
      print*,jstart,jmax,lmax,kmax,gamma
C$OMP DO
      do j=jstart,jmax
         do l=1,lmax
            zlim(j,l)=0
            do k=2,kmax
               if(rho(j,k,l).gt.cutoff) zlim(j,l)=k
            end do
         end do
      end do
C$OMP END DO      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Calculate Q = sspd kappa/(pi G Rsig)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*,'ENTER 2nd loop'
      print*,jmax2,lmax
C$OMP DO
      do j=13,jmax
         print*,j
         rsig(j,1)=0.0
         print*,rsig(j,1)
      enddo   
C$OMP END DO  
C$OMP DO
      do j=jstart,jmax2
         print*,'break j',j
         do l=1,lmax
            print*,'l',l
            rsig(j,l)=0.0
            print*,'rsig',rsig(j,l)
            sspd(j,l)=SQRT((p(j,2,l)/eps(j,2,l)+1)*p(j,2,l)/rho(j,2,l))
            print*,'sound',sspd(j,l)
            do k=2,kmax
               rsig(j,l)=rsig(j,l) + rho(j,k,l)*rof3n 
            end do
            rsig(j,l)=2.0*rsig(j,l)
         end do
      end do
C$OMP END DO      

      print*,'ENTER 3rd loop'

C$OMP DO
      do j=jstart+1,jmax1
         do l=1,lmax
            dphidr(j,l)=(phi(j+1,2,l)-phi(j-1,2,l))/(2.*dr)
         end do
      end do
C$OMP END DO

      print*,'ENTER 4th loop'

C$OMP DO
      do l=1,lmax
         dphidr(jmax2,l)=dphidr(jmax1,l)
      end do
C$OMP END DO

      print*,'ENTER 5th loop'

C$OMP DO
      do l=1,lmax
         do j=jstart+1,jmax1
            kaphi(j,l)=SQRT((
     &           (dphidr(j+1,l)-dphidr(j-1,l))/(2.d0*dr) 
     &           + (3.d0/rhf(j))*dphidr(j,l)))
            qphi(j,l)=sspd(j,l)*kaphi(j,l)/(pi*rsig(j,l))
         end do
      end do
C$OMP END DO

      outfile='qands.'//tim
      outfile1='kandc.'//tim
      
      open(unit=4,file=outfile)
      
       print*,'ENTER 6th loop'

C$OMP DO
      do j=jstart+1,jmax1
         do l=1,lmax
            write(4,66)j,l,qphi(j,l),rsig(j,l),kaphi(j,l),sspd(j,l)
         end do
      end do
C$OMP END DO

      close(4)
      
 66   format(1x,i3,1x,i3,1x,4(1P10E13.5,1x))
      return
      end
