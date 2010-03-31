!-*-f90-*-
PROGRAM PERIODOGRAM
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)
  integer :: jmax, kmax, lmax, jreq, istart, iend, iskip, numargs
  integer :: jmax2, kmax2, lmax2,modemax,numfiles,jreq
  real(double) :: ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,OMMAX,tmassini
  real(double) :: dr
  real(double),dimension(:,:,:),allocatable :: rho,s,t,a,eps,omega
  real(double),dimension(:),allocatable :: omegavg,kappa
  character :: rhodir*80, savedir*80,rhofile*80,savedfile*80
  character :: jmaxin*8,kmaxin*8,lmaxin*8,istartin*8,iendin*8,iskipin*8,modein*8
  character :: filenum*8

  numargs = IARGC()

  if(numargs.ne. ) then
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
  call getarg(8,outfile)
  call getarg(9,rhodir)
  call getarg(10,savedir)
  read(jmaxin,*)jmax
  read(kmaxin,*)kmax
  read(lmaxin,*)lmax
  read(modein,*)modemax
  read(istartin,*)istart
  read(iendin,*)iend
  read(iskipin,*)iskip

  numfiles = ((iend-istart)/iskip)+1
  jmax2 = jmax+2
  kmax2 = kmax+2
  lmax2 = lmax/2

  allocate(rho(jmax2,kmax2,lmax))
  allocate(s(jmax2,kmax2,lmax))
  allocate(t(jmax2,kmax2,lmax))
  allocate(a(jmax2,kmax2,lmax))
  allocate(eps(jmax2,kmax2,lmax))
  allocate(omega(jmax2,kmax2,lmax))
  allocate(omegavg(jmax2))
  allocate(kappa(jmax2))

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

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dr)
  DO j = 1,jmax2
     r(j)   = (DBLE(J)-2.d0)*dr
     rhf(j) = r(j)+(0.5d0*dr)
  ENDDO
!$OMP END PARALLEL DO




  
subroutine Globalcoef(rho,JMAX,LMAX,cn_amp,cn_ang)
  implicit none
  integer, parameter :: double = selected_real_kind(15,300)

  real(double), dimension(:,:,:) :: rho
  real(double), allocatable, dimension(:,:)   :: an,bn
  real(double), allocatable, dimension(:)     :: rho0
  real(double), dimension(:,:) :: cn_amp,cn_ang
  real(double) :: angle, pi

  integer :: J, K, L, JMAX, LMAX
  integer :: N, I

! begin I/O

  pi = acos(-1.d0)

! allocate variables
  allocate(  an(JMAX+2,LMAX/2)               )
  allocate(  bn(JMAX+2,LMAX/2)               )
  allocate(rho0(JMAX+2)                      )
       
!$OMP PARALLEL

!$OMP DO DEFAULT(SHARED) REDUCTION(+:rho0)
     do J = 2, JMAX+1
        rho0(J) = 0.d0
        do L = 1, LMAX
           
           rho0(J) = rho0(J) + rho(J,0,L)

        enddo
     enddo
!$OMP ENDDO

     rho0 = rho0/dble(LMAX)
     
!$OMP DO DEFAULT(SHARED) PRIVATE(angle,pi) REDUCTION(+:an,bn)
     do J = 2, JMAX+1
        do N = 1, LMAX/2
           an(J,N) = 0.d0
           bn(J,N) = 0.d0
           cn_amp(J,N) = 0.d0
           do L = 1, LMAX
              angle = (dble(L))/(LMAX)*2.d0*pi
              an(J,N) = an(J,N) + rho(J,0,L)*cos(dble(N)*angle)
              bn(J,N) = bn(J,N) + rho(J,0,L)*sin(dble(N)*angle)            
           enddo

           an(J,N) = an(J,N)/(dble(LMAX))
           bn(J,N) = bn(J,N)/(dble(LMAX))
           
! compute the amplitude
           cn_amp(J,N) = 2.d0*sqrt(an(J,N)**2+bn(J,N)**2)/rho0(J)
! compute the phase angle
           cn_ang(J,N) = atan2(an(J,N),bn(J,N))+pi

        enddo
     enddo
!$OMP ENDDO

!$OMP END PARALLEL

  deallocate(rho0)
  deallocate(an)
  deallocate(bn)
  RETURN
END subroutine Globalcoef
