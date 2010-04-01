PRO ifof, start, finish

  JMAX = 512
  KMAX = 64
  LMAX = 128
  swap_endian = 1
  rhoprefix = './RHOTEMP/rho3d.'
  rhofile = ''
  torp  = 1605.63
  dr  = 0.2021807d0
  jreq = 250L
  aujreq = 40.d0
  mstar = 1.d0    ; stellar mass in solar masses
  mratio = 0.14   ; disk/star mass ratio

  mtot  = mstar + (mratio*mstar)
  JMAX2 = JMAX+2
  KMAX2 = KMAX+2
  RHOSTART  = DBLARR(JMAX2,KMAX2,LMAX)
  RHOEND    = DBLARR(JMAX2,KMAX2,LMAX)
  MCYLSTART = DBLARR(JMAX2)
  MCYLEND   = DBLARR(JMAX2)
  MDOT      = DBLARR(JMAX2)
  IFOFBDY   = DBLARR(JMAX2)
  r         = DBLARR(JMAX2)
  rhf       = DBLARR(JMAX2)
  STARTIME  = 0.d0
  ENDTIME   = 0.d0
  PI        =  3.1415926535897932384626433832795028841971  

  rhofile   = STRING(rhoprefix,start,FORMAT='(A,I6.6)') 

  FOR j = 0,jmax+1 DO BEGIN
     rhf(j) = (DOUBLE(j)-0.5d0)*dr
     r(j)   = rhf(j)-0.5d0*dr
  ENDFOR

  IF(swap_endian eq 1) THEN BEGIN
     OPENR,lun2,rhofile,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
  ENDIF ELSE BEGIN
     OPENR,lun2,rhofile,/F77_UNFORMATTED,/GET_LUN
  ENDELSE

  READU,lun2,RHOSTART
  READU,lun2,STARTIME
  FREE_LUN,lun2

  rhofile   = STRING(rhoprefix,finish,FORMAT='(A,I6.6)') 

  IF(swap_endian eq 1) THEN BEGIN
     OPENR,lun2,rhofile,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
  ENDIF ELSE BEGIN
     OPENR,lun2,rhofile,/F77_UNFORMATTED,/GET_LUN
  ENDELSE

  READU,lun2,RHOEND
  READU,lun2,ENDTIME
  FREE_LUN,lun2    

  FOR J=1,JMAX DO BEGIN
     FOR K=1,KMAX DO BEGIN
        FOR L=0,LMAX-1 DO BEGIN
           MCYLSTART(J) = MCYLSTART(J) + RHOSTART(J,K,L)
           MCYLEND(J)   = MCYLEND(J)   + RHOEND(J,K,L)
        ENDFOR
     ENDFOR
     MCYLSTART(J) = MCYLSTART(J)*dr*pi*(r(j+1)^2-r(j)^2)/LMAX
     MCYLEND(J) = MCYLEND(J)*dr*pi*(r(j+1)^2-r(j)^2)/LMAX
     IF (J GT 1) THEN BEGIN
        MCYLSTART(J) = MCYLSTART(J) + MCYLSTART(J-1)
        MCYLEND(J)   = MCYLEND(J)   + MCYLEND(J-1)
     ENDIF
     MDOT(J) = MCYLEND(J) - MCYLSTART(J)
  ENDFOR

  rhf = rhf*AUJREQ/rhf(JREQ)
  COUNT = 0L
  FOR J=2,JMAX DO BEGIN
     IF (((MDOT(J) GT 0.d0) AND (MDOT(J-1) LE 0.d0)) OR ((MDOT(J)  $
               LT 0.d0) AND (MDOT(J-1) GE 0.d0))) THEN BEGIN
        IFOFBDY(COUNT) = RHF(J)
        COUNT++
     ENDIF
  ENDFOR

  IFOFBDY = IFOFBDY(0:COUNT-1)

  MDOT = MDOT*mtot          ; mass now in solar masses
  TINTERVAL = ENDTIME-STARTIME
  ent  = AUJREQ/(r(JREQ)*7.93d-3*mtot^(-1.d0/3.d0))
  TINTERVAL = TINTERVAL*ent^(1.5d0)*mtot^(-1.d0)*3.55d3/3.15576d7
  STARTIME = STARTIME/torp
  ENDTIME = ENDTIME/torp
  MDOT = MDOT/TINTERVAL

  print, IFOFBDY
  print, STARTIME,ENDTIME
      
  titlestr  =  STRING('Mdot, time=',startime,'--',endtime,' ORP',FORMAT='(A11,F5.2,A2,F5.2,A4)')
  str         = STRING(IFOFBDY, FORMAT='(F6.2)')
  subtitlestr = STRJOIN(str,' ')

  DEVICE,DECOMPOSED=0
  LOADCT,42
      
  PLOT,RHF(1:JMAX),MDOT(1:JMAX),XRANGE=[0,60],XSTYLE=1,TITLE=titlestr, $
       SUBTITLE= subtitlestr,COLOR=0,BACKGROUND=255
  
  READ,OUTPUT,PROMPT='Press 1 to save this image. '
  IF (OUTPUT EQ 1) THEN BEGIN
     READ,OUTFILE,PROMPT='Please enter the file name. '
     image=tvrd(true=3)
     write_jpeg,OUTFILE,image,QUALITY=100,TRUE=3
  ENDIF

END
