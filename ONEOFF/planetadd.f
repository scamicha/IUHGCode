      program planetadd

C      implicit none.  This program reads the saved file to make
C      the same file with quadruple azimuthal resolution
c      Compile f77 -o interpol interpol.f

      Integer  jmax,jmax1,jmax2,lmax,kmax,kmax1,kmax2
      parameter (jmax=512, jmax1=jmax+1, jmax2=jmax+2)
      parameter (kmax=64, kmax1=kmax+1, kmax2=kmax+2)
      parameter (lmax=512)
      integer jreq, j, k, l, lp
      real*8 rof3n,zof3n, power,logl,loglp
      real*8 delt,den,time,elost,sound,omcen,tmassini,tmass
      real*8 tmassadd,tmassout,tmassacc,totcool,totheat,totdflux,totirr
      real*8 etotfl,eflufftot,planetx,planety,planetvx,planetvy
      real*8 v_c
      real*8 rdiskau,planetr,planetphi
      real*8 rho(jmax2,kmax2,lmax),s(jmax2,kmax2,lmax),t(jmax2,kmax2
     &     ,lmax),a(jmax2,kmax2,lmax),eps(jmax2,kmax2,lmax)
      real*8 r(JMAX2)
      real*8 diskmass,dphi,pi

      pi   = 3.1415926535897932384626433832795028841971d0
      dphi = 2.d0*pi/FLOAT(LMAX)
      planetr   = 25.d0
      planetphi = 0.d0
      rdiskau   = 40.d0

      open(unit=11,file='saved.01000000',form='unformatted')
      print*,'Opening file'
      read(11)s
      read(11)t
      read(11)a
      read(11)rho
      read(11)eps
      read(11)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMCEN
      read(11)tmassini,tmass,tmassadd,tmassout,tmassacc,totcool
     &     ,totdflux,totheat,totirr,etotfl,eflufftot
      close(11)

      DO J=1,JMAX2
         R(J) = (J-2)*ROF3N
      ENDDO

      diskmass = 0.d0
      J=11
      DO WHILE ((R(J)).lt.(r(JREQ)*planetr/rdiskau))
         DO l=1,lmax
            DO k=2,kmax
               diskmass = diskmass+rho(j,k,l)*0.5*dphi*zof3n*(r(j+1)**2
     &              -r(j)**2)
            ENDDO
         ENDDO
         J=J+1
         print*,J
      ENDDO

      print*,diskmass, tmass, (1.d0-tmassini)+diskmass
         
      planetx = (((JREQ-2)*ROF3N)*planetr*cos(planetphi))/rdiskau
      planety = (((JREQ-2)*ROF3N)*planetr*sin(planetphi))/rdiskau

!     Get circular velocities

      v_c = sqrt(((1.d0-tmassini)+diskmass)/(planetr*r(JREQ)/rdiskau))
      planetvx = -v_c*sin(planetphi)
      planetvy = v_c*cos(planetphi)
 



      open(unit=8,file='savedplanet.01000000',form='unformatted')
      WRITE(8) S
      WRITE(8) T
      WRITE(8) A
      WRITE(8) RHO
      WRITE(8) EPS
      Write(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMCEN
      write(8) tmassini,tmass,tmassadd,tmassout,tmassacc,totcool
     &     ,totdflux,totheat,totirr,etotfl,eflufftot
      write(8) planetx,planety,planetvx,planetvy
      close(8)


      end
