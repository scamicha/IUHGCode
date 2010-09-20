c-*-fortran-*- 
c     This file includes the parameter statements for the dcomp.f code.
c     These can be modified as necessary to run different sized grids;
c     however, parameters starting with 'pot3' must be powers of 2.
c     

C     Strings are set to 80 characters, if you need a longer path
C     resize them here.

      character datadir*80,outdir*80,outfile*80,amfile*80
      INTEGER ISTART,IEND,ISKIP
      integer  jmax,jmax1,jmax2
      integer  kmax,kmax1,kmax2
      integer  lmax,lmax2
      integer  pot3jmax,pot3jmax1,pot3jmax2
      integer  pot3kmax,pot3kmax1,pot3kmax2
      integer  jmin,jmin1,jmin2

C     datadir is the location of the saved files, 
C     outdir the directory where outfile will be written

      parameter (datadir = '../../WAN_DATA/')
      parameter (outdir  = './')
      parameter (outfile = 'decompbaselineMIXdat.MPI.')

C     The following parameters are for the EOS
!    Use H2STAT to select what type of mixture you want.
!    0 = pure para
!    1 = pure ortho -- this is not really physical for astro applications
!    2 = set mixture
!    3 = equilibrium
!    4 = single gamma law
!    Be sure to set ac and bc for mixture (ac: para component, bc: ortho component)
!    also pick a metallicity by adjusting xabun, yabun, and zabun.
!    please make sure you use the same values from the hydro run!!!

!     CHANGE THESE
      integer,parameter :: H2STAT=2
      real*8, parameter :: ac = 1.d0, bc = 3.d0
      real*8, parameter :: xabun = .73d0,yabun=.25d0,zabun=.02d0
      real*8, parameter :: gamma = 1.66666666666667d0,xn = 1.5d0
!     PROBABLY DON'T CHANGE THESE
      real*8,PARAMETER :: Mstar=1.0,Rstar=2.0,Rdiskau=40.0,Tstar=4000.0
      real*8,PARAMETER :: Msuncgs=1.989e33, Rsuncgs=6.96e10, Tbgrnd=3.
      real*8,PARAMETER :: sigmacgs=5.670e-5, Gcgs=6.672e-8, phylim=1.e-8 
      real*8,PARAMETER :: bkmpcgs=8.254d7,gridlim=1.e-12, AUcgs=1.496e13


C     ISTART is the first file IEND the last one, and ISKIP the interval
C     between saved files. Also make sure to set JMAX, KMAX, and LMAX 
C     correctly. pot3jmax and pot3kmax should be the same as JMAX and KMAX

      parameter (ISTART=550000,IEND=2050000,ISKIP=5000)
      parameter (kmax= 64, kmax1=kmax+1, kmax2=kmax+2)
      parameter(pot3kmax=kmax,pot3kmax1=pot3kmax+1,pot3kmax2=pot3kmax+2)

      parameter (jmax=512, jmax1=jmax+1, jmax2=jmax+2)
      parameter(pot3jmax=jmax,pot3jmax1=pot3jmax+1,pot3jmax2=pot3jmax+2)

      parameter (lmax=512, lmax2=lmax/2)

c  Minimum radial grid point, for cutting out central star.
      parameter (jmin=11,jmin1=jmin-1,jmin2=jmin-2)

c  Other parameters.
      integer hj,hk,hj1,hj2,hk1,hk2,itable,lrlmax,lrjmax,lrkmax,
     &        lrjmax1,lrkmax1
      parameter (hj=256,hk=256,hj1=hj+1,hk1=hk+1,hj2=hj+2,hk2=hk+2)
      parameter (lrlmax=64,lrjmax=32,lrkmax=8,lrjmax1=lrjmax+1)
      parameter (itable=100)
      INTEGER, PARAMETER :: TTABLE=400
