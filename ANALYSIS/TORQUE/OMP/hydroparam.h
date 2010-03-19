c     This file includes the parameter statements for the 3dhyd.f code.
c     These can be modified as necessary to run different sized grids;
c     however, parameters starting with 'pot3' must be powers of 2.
c     
c     -rpl, Sep, 1996

c Do not currently support different resolutions in the l-direction.
      integer  jmax,jmax1,jmax2
      integer  kmax,kmax1,kmax2
      integer  lmax,lmax2
      integer  pot3jmax,pot3jmax1,pot3jmax2
      integer  pot3kmax,pot3kmax1,pot3kmax2
      integer  jmin,jmin1,jmin2

      parameter (kmax= 64, kmax1=kmax+1, kmax2=kmax+2)
      parameter(pot3kmax=64, pot3kmax1=pot3kmax+1, pot3kmax2=pot3kmax+2)

      parameter (jmax=512, jmax1=jmax+1, jmax2=jmax+2)
      parameter(pot3jmax=512,pot3jmax1=pot3jmax+1,pot3jmax2=pot3jmax+2)

      parameter (lmax=128, lmax2=lmax/2)

c  Minimum radial grid point, for cutting out central star.
      parameter (jmin=11,jmin1=jmin-1,jmin2=jmin-2)

c  Other parameters.
      integer hj,hk,hj1,hj2,hk1,hk2,itable,lrlmax,lrjmax,lrkmax,
     &        lrjmax1,lrkmax1
      parameter (hj=256,hk=256,hj1=hj+1,hk1=hk+1,hj2=hj+2,hk2=hk+2)
      parameter (lrlmax=64,lrjmax=32,lrkmax=8,lrjmax1=lrjmax+1)
      parameter (itable=100)


