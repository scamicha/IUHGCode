PRO PERIODPLOT
  JMAX   = 0L
  JREQ   = 242L
  AUJREQ = 0.d0
  MODES  = 0L
  MODENUM= 0L
  NFREQ  = 0L
  CHLEV  = 1L
  TSTART = 0.d0
  TEND   = 0.d0
  VALMAX = 0.d0
  VALMIN = 0.d0
  VALMID = 0.d0
  NUMLEV =  DBLARR(2)
  YRANGE =  DBLARR(2)
  XRANGE =  DBLARR(2)

  INFILE = 'period_indirect.dat'

  openr,lun1,INFILE,/GET_LUN,/F77_UNFORMATTED
  readu,lun1,JMAX,MODES,NFREQ,TSTART,TEND,AUJREQ
  values = DBLARR(JMAX,MODES,NFREQ)
  freqs  = DBLARR(JMAX,MODES,NFREQ)
  omega  = DBLARR(JMAX)
  kappa  = DBLARR(JMAX)
  rads   = DBLARR(JMAX)
  readu,lun1,values
  readu,lun1,freqs
  readu,lun1,rads
  readu,lun1,omega
  readu,lun1,kappa

  titlestr  =  STRING('Periodogram Contour m=',modenum,'time=',tstart,'--',tend,' ORP',FORMAT='(A22,I2,A8,F5.2,A2,F5.2,A4)')
  modestring = STRING('There are ', MODES, '  modes available. Please choose one to plot. ',FORMAT='(A,I2,A)')

  READ,MODENUM,PROMPT=modestring
  IF ((MODENUM LT 1) OR (MODENUM GT MODES)) THEN BEGIN 
     PRINT, "You chose a mode outside the data set"
     STOP
  ENDIF

  VALMAX = MAX(values,MIN=VALMIN)
  VALMID = MEDIAN(values,/EVEN)
  IF (VALMIN LE 1.e-5) THEN VALMIN = 1.e-5

  WHILE (CHLEV EQ 1) DO BEGIN

     PRINT,"The minimum value of the contour is: ", VALMIN
     PRINT,"The maximum value of the contour is: ", VALMAX
     PRINT,"The median value of the contour is: ",  VALMID

     READ,NUMLEV,                                                $
          PROMPT="Please enter contour limits (min,max) "

     LEVEL     =  INDGEN(256,/FLOAT)
     LEVEL     =  LEVEL/255.*(NUMLEV(1)-NUMLEV(0))+NUMLEV(0)
     LOADCT, 43            

     READ,XRANGE,PROMPT="Please enter Pattern range to plot: "
     READ,YRANGE,PROMPT="Please enter radius range to plot: "
     
     DEVICE, DECOMPOSED=0
     window,0,xsize=640,ysize=640

     CONTOUR,values(*,MODENUM,*),freqs(2,MODENUM,*),rads,           $
             YRANGE=YRANGE,LEVELS=LEVEL,YSTYLE=1,TITLE=titlestr,    $
             color=0,background=255,XTITLE='Pattern Speed (1/ORP)', $
             YTITLE='R(AU)',XRANGE=XRANGE,XSTYLE=1,/FILL
                      
     OPLOT,omega,rads,color=0
     OPLOT,omega+(kappa/DOUBLE(modenum)),rads,color=0
     OPLOT,omega-(kappa/DOUBLE(modenum)),rads,color=0

     READ,CHLEV,PROMPT='Press 1 to change contour levels '

  ENDWHILE
     
  READ,OUTPUT,PROMPT='Press 1 to save this image. '
  IF (OUTPUT EQ 1) THEN BEGIN
     READ,OUTFILE,PROMPT='Please enter the file name. '
     image=tvrd(true=3)
     write_jpeg,OUTFILE,image,QUALITY=100,TRUE=3
  ENDIF

END
