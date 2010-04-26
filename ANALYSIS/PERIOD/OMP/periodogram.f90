!-*-f90-*-
PROGRAM PERIODOGRAM
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)
  integer :: jmax, kmax, lmax, jreq, istart, iend, iskip, numargs
  integer :: jmax2, kmax2, lmax2,modes,numfiles,jreq,fileiter
  real(double) :: ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,OMMAX,tmassini
  real(double) :: dr,pi,torp,aujreq
  real(double),dimension(:,:,:),allocatable :: rho,s,t,a,eps,omega
  real(double),dimesion(:,:,:),allocatable :: amplitude,angle
  real(double),dimension(:),allocatable :: omegavg,kappa,massum,timearr
  character :: rhodir*80, savedir*80,rhofile*94,savedfile*94
  character :: jmaxin*8,kmaxin*8,lmaxin*8,istartin*8,iendin*8,iskipin*8,modein*8
  character :: filenum*8,aujreqin*10

  pi = ACOS(-1.d0)
  torp = 1605.63

  numargs = IARGC()

  if(numargs.ne. 11) then
     print*,"Incorrect number of arguments"
     STOP
  ENDIF
  
  call getarg(1,jmaxin)
  call getarg(2,kmaxin)
  call getarg(3,lmaxin)
  call getarg(4,modein)
  call getarg(5,istartin)
  call getarg(6,iendin)
  call getarg(7,iskipin)
  call getarg(8,aujreqin)
  call getarg(9,outfile)
  call getarg(10,rhodir)
  call getarg(11,savedir)
  read(jmaxin,*)jmax
  read(kmaxin,*)kmax
  read(lmaxin,*)lmax
  read(modein,*)modes
  read(istartin,*)istart
  read(iendin,*)iend
  read(iskipin,*)iskip
  read(aujreqin,*)aujreq

  numfiles = ((iend-istart)/iskip)+1
  jmax1 = jmax+1
  kmax1 = kmax+1
  jmax2 = jmax+2
  kmax2 = kmax+2
  lmax2 = lmax/2

  allocate(timearr(numfiles))
  allocate(rho(jmax2,kmax2,lmax))
  allocate(s(jmax2,kmax2,lmax))
  allocate(t(jmax2,kmax2,lmax))
  allocate(a(jmax2,kmax2,lmax))
  allocate(eps(jmax2,kmax2,lmax))
  allocate(omega(jmax2,kmax2,lmax))
  allocate(omegavg(jmax2))
  allocate(kappa(jmax2))
  allocate(massum(jmax2))
  allocate(amplitude(jmax2,modes,numfiles))
  allocate(angle(jmax2,modes,numfiles))

  write (filenum,'(I8.8)')iend
  savedfile=trim(savedir)//'saved.'//filenum
  OPEN(UNIT=8, FILE=trim(savedfile),FORM='UNFORMATTED',  &
       STATUS="OLD",ERR=911)      
  READ(8) S
  READ(8) T
  READ(8) A
  READ(8) RHO
  READ(8) EPS
  READ(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
  read(8) tmassini
  CLOSE(8)

!$OMP PARALLEL DEFAULT(SHARED)

!$OMP DO
  DO j = 1,jmax2
     r(j)   = (DBLE(J)-2.d0)*dr
     rhf(j) = r(j)+(0.5d0*dr)
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=1,jmax2
     DO k=1,kmax2
        DO l=1,LMAX
           omega(j,k,l)=a(j,k,l)/(rhf(j)**2)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO DEFAULT(SHARED) REDUCTION(+:omegavg,massum)
  DO J=2,jmax1
     omegavg(j) = 0.d0
     massum(j) = 0.d0
     DO k=1,kmax2
        DO l=1,lmax
           omegavg(j) = omegavg(j) + omega(j,k,l)
           massum(j) = massum(j) + rho(j,k,l)
        ENDDO
     ENDDO
  ENDDO
!$OMP END DO
     
!$OMP DO DEFAULT(SHARED)
  DO J=2,jmax1
     omegavg(j) = omegavg(j)/massum(j)
     kappa(j) = sqrt(ABS(rhf(j)*omegavg(j)*(omegavg(j+1)- &
          omegavg(j-1))/dr + 4.d0*omegavg(j)**2))*(torp/(2.d0*pi))
     omegavg(j) = omegavg(j)*(torp/(2.d0*pi))
!!!!!!!!!!!!!!!!!!!!! If not needed put change rhf to AU here.........
     
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  DO I = 1,NUMFILES
     fileiter = (I*stepskip)+stepbegin
     write (filenum,'(I8.8)')fileiter
     rhofile = trim(rhodir)//'rho3d.'//filenum
     OPEN(UNIT=8, FILE=trim(rhofile),FORM='UNFORMATTED',  &
          STATUS='OLD',ERR=911)
     READ(8) RHO
     READ(8) time
     CLOSE(8)
     timearr(I)=time
     CALL GlobalCoef(rho,I,jmax,lmax,modes,amplitude,angle)
  ENDDO

  
     
subroutine GlobalCoef(rho,I,JMAX,LMAX,MMAX,cn_amp,cn_ang)
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)

  real(double), dimension(:,:,:) :: rho
  real(double), allocatable, dimension(:,:)   :: an,bn
  real(double), allocatable, dimension(:)     :: rho0
  real(double), dimension(:,:,:) :: cn_amp,cn_ang
  real(double) :: angle, pi

  integer :: J, K, L, JMAX, LMAX
  integer :: N, I,MMAX

! begin I/O

  pi = acos(-1.d0)

! allocate variables
  allocate(  an(JMAX+2,0:LMAX/2)               )
  allocate(  bn(JMAX+2,0:LMAX/2)               )
  allocate(rho0(JMAX+2)                      )
       
!$OMP PARALLEL

!$OMP DO DEFAULT(SHARED) REDUCTION(+:rho0)
     do J = 2, JMAX+1
        rho0(J) = 0.d0
        do L = 1, LMAX
           
           rho0(J) = rho0(J) + rho(J,2,L)

        enddo
     enddo
!$OMP ENDDO

     rho0 = rho0/dble(LMAX)
     
!$OMP DO DEFAULT(SHARED) PRIVATE(angle) FIRSTPRIVATE(pi) REDUCTION(+:an,bn)
     do J = 2, JMAX+1
        an(J,0) = rho0*dble(LMAX)
        do N = 1, MMAX
           an(J,N) = 0.d0
           bn(J,N) = 0.d0
           cn_amp(J,N,I) = 0.d0
           do L = 1, LMAX
              angle = (dble(L))/(LMAX)*2.d0*pi
              an(J,N) = an(J,N) + rho(J,0,L)*cos(dble(N)*angle)
              bn(J,N) = bn(J,N) + rho(J,0,L)*sin(dble(N)*angle)            
           enddo

           an(J,N) = an(J,N)/(dble(LMAX))
           bn(J,N) = bn(J,N)/(dble(LMAX))
           
! compute the amplitude
           cn_amp(J,N,I) = 2.d0*sqrt(an(J,N)**2+bn(J,N)**2)/rho0(J)
! compute the phase angle
           cn_ang(J,N,I) = atan2(an(J,N),bn(J,N))+pi

        enddo
     enddo
!$OMP ENDDO

!$OMP END PARALLEL

  deallocate(rho0)
  deallocate(an)
  deallocate(bn)
  RETURN
END subroutine Globalcoef
