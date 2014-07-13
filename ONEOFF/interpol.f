      program inter

C      implicit none.  This program reads the saved file to make
C      the same file with quadruple azimuthal resolution
c      Compile f77 -o interpol interpol.f

      Integer  jmax,jmax1,jmax2,newjmax2
      integer  lmax,newlmax
      integer  kmax,kmax1,kmax2,newkmax2
      parameter (jmax=512, jmax1=jmax+1, jmax2=jmax+2,newjmax2=514)
      parameter (kmax=64, kmax1=kmax+1, kmax2=kmax+2,newkmax2=66)
      parameter (lmax=512,newlmax=2048)
      integer jreq, j, k, l, lp
      real*8 rof3n,zof3n, power,logl,loglp
      real*8 delt,den,time,elost,sound,omcen,tmassini,tmass
      real*8 tmassadd,tmassout,tmassacc,totcool,totheat,totdflux,totirr
      real*8 etotfl,eflufftot
      real*8 rho(jmax2,kmax2,lmax),s(jmax2,kmax2,lmax),t(jmax2,kmax2
     &     ,lmax),a(jmax2,kmax2,lmax),eps(jmax2,kmax2,lmax)
      real*8 newrho(newjmax2,newkmax2,newlmax),news(newjmax2,newkmax2
     &     ,newlmax),newt(newjmax2,newkmax2,newlmax),newa(newjmax2
     &     ,newkmax2,newlmax),neweps(newjmax2,newkmax2,newlmax)

      open(unit=11,file='saved.00550000',form='unformatted')
      print*,'Opening file'
      read(11)s
      read(11)t
      read(11)a
      read(11)rho
      read(11)eps
      read(11)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMCEN
      read(11)tmassini,tmass,tmassadd,tmassout,tmassacc,totcool
     &     ,totdflux,totheat,totirr,etotfl,eflufftot
      print*,'Check format',rho(208,3,4),rho(235,3,2),time
      close(11)

      
      do j=1,newjmax2
         do k=1,newkmax2
            do l=1,lmax
               lp=l+1
               if(l.eq.lmax) lp=1

C              For rho

               logl=0.0
               loglp=0.0
               if(log10(rho(j,k,l)).ne.0.0) logl=log10(rho(j,k,l))
               if(log10(rho(j,k,lp)).ne.0.0) loglp=log10(rho(j,k,lp))

               power=(7.0/8.0)*logl+(1.0/8.0)*loglp
               newrho(j,k,l*4-3)=10.0**power
               
               power=(5.0/8.0)*logl+(3.0/8.0)*loglp
               newrho(j,k,l*4-2)=10.0**power

               power=(3.0/8.0)*logl+(5.0/8.0)*loglp
               newrho(j,k,l*4-1)=10.0**power
 
               power=(1.0/8.0)*logl+(7.0/8.0)*loglp
               newrho(j,k,l*4)=10.0**power

               if(power.eq.0.0) then
                  newrho(j,k,l*4-3)=0.0
                  newrho(j,k,l*4-2)=0.0
                  newrho(j,k,l*4-1)=0.0
                  newrho(j,k,l*4)=0.0
                  print*,'rho',j,k,l
               endif
                  

C              For S

               news(j,k,l*4-3)=(7.0/8.0)*s(j,k,l)+(1.0/8.0)*s(j,k,lp)
               
               news(j,k,l*4-2)=(5.0/8.0)*s(j,k,l)+(3.0/8.0)*s(j,k,lp)
               
               news(j,k,l*4-1)=(3.0/8.0)*s(j,k,l)+(5.0/8.0)*s(j,k,lp)
               
               news(j,k,l*4)  =(1.0/8.0)*s(j,k,l)+(7.0/8.0)*s(j,k,lp)
               
               
C              For T
               
               newt(j,k,l*4-3)=(7.0/8.0)*t(j,k,l)+(1.0/8.0)*t(j,k,lp)
               
               newt(j,k,l*4-2)=(5.0/8.0)*t(j,k,l)+(3.0/8.0)*t(j,k,lp)
               
               newt(j,k,l*4-1)=(3.0/8.0)*t(j,k,l)+(5.0/8.0)*t(j,k,lp)
               
               newt(j,k,l*4)  =(1.0/8.0)*t(j,k,l)+(7.0/8.0)*t(j,k,lp)
               
               
C              For A
               
               newa(j,k,l*4-3)=(7.0/8.0)*a(j,k,l)+(1.0/8.0)*a(j,k,lp)
               
               newa(j,k,l*4-2)=(5.0/8.0)*a(j,k,l)+(3.0/8.0)*a(j,k,lp)
               
               newa(j,k,l*4-1)=(3.0/8.0)*a(j,k,l)+(5.0/8.0)*a(j,k,lp)
               
               newa(j,k,l*4)  =(1.0/8.0)*a(j,k,l)+(7.0/8.0)*a(j,k,lp)

               
C              For EPS
               
               logl=0.0
               loglp=0.0
               if(log10(eps(j,k,l)).ne.0.0) logl=log10(eps(j,k,l))
               if(log10(eps(j,k,lp)).ne.0.0) loglp=log10(eps(j,k,lp))

               power=(7.0/8.0)*logl+(1.0/8.0)*loglp
               neweps(j,k,l*4-3)=10.0**power
                  
               power=(5.0/8.0)*logl+(3.0/8.0)*loglp
               neweps(j,k,l*4-2)=10.0**power

               power=(3.0/8.0)*logl+(5.0/8.0)*loglp
               neweps(j,k,l*4-1)=10.0**power
 
               power=(1.0/8.0)*logl+(7.0/8.0)*loglp
               neweps(j,k,l*4)=10.0**power             

               if(power.eq.0.0) then
                  neweps(j,k,l*4-3)=0.0
                  neweps(j,k,l*4-2)=0.0
                  neweps(j,k,l*4-1)=0.0
                  neweps(j,k,l*4)=0.0
                  print*,'eps',j,k,l
               endif

            enddo
         enddo
      enddo

      do j=1,newjmax2
         print*,j,newrho(j,5,1),news(j,5,1),newt(j,5,1),newa(j,5,1)
     &        ,neweps(j,5,1)
      enddo

!      print*,'rho',rho(100,7,lmax),rho(100,7,1)
!      print*,'new',newrho(100,7,509),newrho(100,7,510),newrho(100,7,511)
!     &     ,newrho(100,7,512)
      

!      print*,'a',a(100,7,lmax),a(100,7,1)
!      print*,'new',newa(100,7,509),newa(100,7,510),newa(100,7,511)
!     &     ,newa(100,7,512)

!      print*,'rho',rho(5,7,lmax),rho(5,7,1)
!      print*,'new',newrho(5,7,509),newrho(5,7,510),newrho(5,7,511)
!     &     ,newrho(5,7,512)
      

!      print*,'a',a(5,7,lmax),a(5,7,1)
!      print*,'new',newa(5,7,509),newa(5,7,510),newa(5,7,511)
!     &     ,newa(5,7,512)


      open(unit=8,file='savedLMAX2048.00550000',form='unformatted')
      WRITE(8) newS
      WRITE(8) newT
      WRITE(8) newA
      WRITE(8) newRHO
      WRITE(8) newEPS
      Write(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMCEN
      write(8) tmassini,tmass,tmassadd,tmassout,tmassacc,totcool
     &     ,totdflux,totheat,totirr,etotfl,eflufftot
      close(8)


      end
