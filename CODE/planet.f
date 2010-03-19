      SUBROUTINE planet_initialize(planet_x,planet_y,planet_vx,
     &     planet_vy,planet_ax,planet_ay)
      IMPLICIT NONE      

      include 'hydroparam.h'
      include 'globals.h'
      include 'units.h'

      REAL*8 planet_x,planet_y,planet_vx,planet_vy,v_c,planet_ax,
     &     planet_ay
      REAL*8 initial_x,initial_y,initial_vx,initial_vy
      COMMON /PLANETINIT/ initial_x,initial_y,initial_vx,initial_vy

!     Set planet mass

      mass_planet = mass_ratio*mass_star

!     Get planet positions

      planet_x  = initial_x
      planet_y  = initial_y
      planet_vx = initial_vx
      planet_vy = initial_vy
      
      CALL planet_force(planet_x,planet_y,planet_ax,planet_ay)


      RETURN
      END

      SUBROUTINE planet_force(planet_x,planet_y,planet_ax,planet_ay)
      IMPLICIT NONE      

      INCLUDE 'hydroparam.h'
      INCLUDE 'globals.h'
      INCLUDE 'units.h'      

      INTEGER planet_j,planet_l,J,K,L,LP,LM,acc_j,acc_l
      REAL*8 planet_x,planet_y,planet_r,planet_phi,planet_ax,planet_ay
      REAL*8 planet_ar,planet_aphi,dr,dz,dphi,r_factor,phi_factor
      REAL*8 ar(2,2,2), aphi(2,2,2)
      REAL*8 mid_ar(2,2),mid_aphi(2,2)
      
      dr = rof3n
      dz = zof3n
      dphi = 2.d0*pi/FLOAT(LMAX)

!     First get planet position on the grid
      planet_r = sqrt((planet_x)**2+(planet_y)**2)
      planet_phi = atan2(planet_y,planet_x)
      if(planet_phi.lt.0.d0) planet_phi = planet_phi + 2.d0*pi
      planet_j = INT(planet_r/dr)+2
      planet_l = INT(planet_phi/dphi)
      if(planet_l.eq.0) planet_l = LMAX
      IF((planet_j.ge.JMAX+1).or.(planet_j.le.JMIN)) STOP
      
!     Collect relevant r and phi forces
      DO L=planet_l,planet_l+1
         acc_l = L-planet_l+1
         IF(L.eq.1)then 
            LM = LMAX
         ELSE 
            LM = L-1
         ENDIF
         IF(L.eq.LMAX+1)then
            LP = 1
         ELSE
            LP = L
         ENDIF
         DO J=planet_j,planet_j+1
            acc_j = J-planet_j+1
            DO K=1,2
               ar(acc_j,K,acc_l)=-((PHI(J,K+1,LM)-PHI(J-1,K+1,LM))+
     &              (PHI(J,K+1,LP)-PHI(J-1,K+1,LP)))/(2.d0*dr)
               aphi(acc_j,K,acc_l)=-(((PHI(J-1,K+1,LP)-PHI(J-1,K+1,LM))/
     &              rhf(J-1))+((PHI(J,K+1,LP)-PHI(J,K+1,LM))/rhf(J)))/
     &              (2.d0*dphi)
            ENDDO
            mid_ar(acc_j,acc_l) = ((ar(acc_j,1,acc_l)*zhf(3)**2)
     &           -(ar(acc_j,2,acc_l)*zhf(2)**2))/(zhf(3)**2-
     &           zhf(2)**2)
            mid_aphi(acc_j,acc_l) = ((aphi(acc_j,1,acc_l)*
     &           zhf(3)**2)-(aphi(acc_j,2,acc_l)*zhf(2)**2))/(
     &           zhf(3)**2-zhf(2)**2)
            
         ENDDO
      ENDDO
      
!     Interpolate to planet position using midplane matrices

      r_factor = (planet_r-r(planet_j))/dr
      IF(planet_l.eq.LMAX) THEN
         phi_factor = (planet_phi)/dphi
      ELSE
         phi_factor = (planet_phi-(dphi*(planet_l)))/dphi
      ENDIF
      IF (r_factor.gt.1.d0) r_factor = 1.d0
      IF (phi_factor.gt.1.d0) phi_factor = 1.d0
      planet_ar = (1.d0-phi_factor)*(((1.d0-r_factor)*mid_ar(1,1))+
     &     (r_factor*mid_ar(2,1)))+(phi_factor)*(((1.d0-r_factor)*
     &     mid_ar(1,2))+(r_factor*mid_ar(2,2)))
      planet_aphi = (1.d0-phi_factor)*(((1.d0-r_factor)*mid_aphi(1,1))+
     &     (r_factor*mid_aphi(2,1)))+(phi_factor)*(((1.d0-r_factor)*
     &     mid_aphi(1,2))+(r_factor*mid_aphi(2,2)))      

!     Get x,y forces from r,phi forces, adding in accelration from star

!     TEST ORBIT

!      print*,'planet_x = ',planet_x,' planet_y = ',planet_y
!      print*,'planet_r = ',planet_r,' planet_phi = ',planet_phi
!      print*,planet_j, planet_l
!      print*,'ar = ',planet_ar,' aphi = ',planet_aphi
!      print*,'r_factor = ',r_factor,' phi_factor = ',phi_factor
!      print*,''
!      print*,mid_ar
!      print*,mid_aphi
!      print*, ''
!      if(abs(planet_aphi).gt.1.d-7)then
!      STOP
!      ENDIF

      planet_ax =-(((mass_star+tmassacc+mass_planet)*planet_r*
     &     cos(planet_phi))/planet_r**3)+((planet_ar*cos(planet_phi))-
     &     (planet_aphi*sin(planet_phi)))
      planet_ay = -(((mass_star+tmassacc+mass_planet)*planet_r*
     &     sin(planet_phi))/planet_r**3)+((planet_ar*sin(planet_phi))+
     &     (planet_aphi*cos(planet_phi)))
      
      RETURN
      END
      
      

      SUBROUTINE planet_integrate(planetx,planety,planetvx,planetvy,
     &     planetax,planetay)
      IMPLICIT NONE
      
      include 'hydroparam.h'
      include 'globals.h'
      include 'units.h'

      REAL*8 planetx,planety,planetvx,planetvy,planetax,planetay
      REAL*8 new_ax,new_ay

!     Calculate new x and y values

      planetx = planetx + planetvx*delt + 0.5d0*planetax*delt**2
      planety = planety + planetvy*delt + 0.5d0*planetay*delt**2

      CALL planet_force(planetx,planety,new_ax,new_ay)

      planetvx = planetvx + 0.5d0*(planetax + new_ax)*delt
      planetvy = planetvy + 0.5d0*(planetay + new_ay)*delt

      planetax = new_ax
      planetay = new_ay

      RETURN
      END

      SUBROUTINE potential_energy(planet_x,planet_y,potential)
      
      IMPLICIT NONE

      include 'hydroparam.h'
      include 'globals.h'
      include 'units.h'

      
      INTEGER planet_j,planet_l,J,K,L,LP,JP
      REAL*8 planet_r,planet_phi,dr,dz,dphi,r_factor,phi_factor
      REAL*8 planet_x,planet_y,potential,delz
      
      dr = rof3n
      dz = zof3n
      dphi = 2.d0*pi/FLOAT(LMAX)
      delz = (zhf(3)**2-zhf(2)**2)

!     First get planet position on the grid
      planet_r = sqrt((planet_x)**2+(planet_y)**2)
      planet_phi = atan2(planet_y,planet_x)
      if(planet_phi.lt.0.d0) planet_phi = planet_phi + 2.d0*pi
      J = INT(planet_r/dr)+2
      L = INT(planet_phi/dphi)
      if(L.eq.0) L = LMAX      
      

      r_factor = (planet_r-r(J))/dr
      IF(L.eq.LMAX) THEN
         phi_factor = (planet_phi)/dphi
      ELSE
         phi_factor = (planet_phi-(dphi*(L)))/dphi
      ENDIF
      IF (r_factor.gt.1.d0) r_factor = 1.d0
      IF (phi_factor.gt.1.d0) phi_factor = 1.d0
      IF (r_factor.gt.0.5d0) THEN
         JP = J+1
         r_factor = r_factor - 0.5d0
      ELSE
         JP = J-1
         r_factor = 0.5d0 - r_factor
      ENDIF
      IF (phi_factor.gt.0.5d0) THEN
         LP = L+1
         phi_factor= phi_factor - 0.5d0
      ELSE
         LP = L-1
         phi_factor = 0.5d0 - phi_factor
      ENDIF
      IF (LP.eq.LMAX+1) LP = 1
      IF (LP.eq.0) LP = LMAX

      
      potential = (1.d0-phi_factor)*((1.d0-r_factor)*(((PHI(J,2,L)*
     &     zhf(3)**2)-(PHI(J,3,L)*zhf(2)**2))/delz)+
     &     (r_factor*(((PHI(JP,2,L)*zhf(3)**2)-(PHI(JP,3,L)*zhf(2)**2))
     &     /delz)))+(phi_factor)*(((1.d0-r_factor)*
     &     (((PHI(J,2,LP)*zhf(3)**2)-(PHI(J,3,LP)*zhf(2)**2))/delz))+
     &     (r_factor*(((PHI(JP,2,LP)*zhf(3)**2)-(PHI(JP,3,LP)*zhf(2)**2)
     &     )/delz)))

      RETURN
      END
