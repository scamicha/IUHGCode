; This program makes and plots surface density profiles and will fit a
; power law to the surface density if you want it to

 PRO sigmaplot, filenum, MIN=RMIN, MAX=RMAX, SWAP=swap_endian, PREFIX=prefix, OVER=overplot, $
	STARTF=startfit, ENDF=endfit, ZOOM=sigmazoom

 !P.COLOR = 0
 !P.BACKGROUND = 255
; !P.MULTI = [0,0,3,0,0]

 IF N_Elements(swap_endian) EQ 0 THEN swap_endian = 0L
 IF N_Elements(RMIN) EQ 0 THEN RMIN = 10.d0
 IF N_Elements(RMAX) EQ 0 THEN RMAX = 40.d0
 IF N_Elements(overplot) EQ 0 THEN overplot = 0L
 IF N_Elements(startfit) EQ 0 THEN startfit = 15.d0
 IF N_Elements(endfit) EQ 0 THEN endfit = 40.d0
 IF N_Elements(sigmazoom) EQ 0 THEN sigmazoom = 0.7477962d0
 IF N_Elements(prefix) EQ 0 THEN prefix = './saved'

 JMAX    = 512L
 KMAX    = 64L
 LMAX    = 512L
 JREQ    = 0L
 fitsize = 0L
 sigconv = 0d0
 time    = 0d0
 tconv   = 1605.63d0
 startjfit=0L
 endjfit  = 0L
 pi      = 3.1415926353892d0
 au      = 40.d0
 jreqs   = 0L
 dr      = 0.d0
 dz      = 0.d0

 s          = DBLARR(jmax+2,kmax+2,lmax)
 t          = DBLARR(jmax+2,kmax+2,lmax)
 a          = DBLARR(jmax+2,kmax+2,lmax)
 rho        = DBLARR(jmax+2,kmax+2,lmax)
 eps        = DBLARR(jmax+2,kmax+2,lmax)
 r          = DBLARR(jmax+2)
 rhf        = DBLARR(jmax+2)
 z          = DBLARR(kmax+2)
 zhf        = DBLARR(kmax+2)
 q          = DBLARR(jmax+2,lmax)
 qav        = DBLARR(jmax+2)
 csound     = DBLARR(jmax+2,lmax)
 kappa      = DBLARR(jmax+2,lmax)
 omega      = DBLARR(jmax+2,kmax+2,lmax)
 q0         = DBLARR(jmax+2,lmax)
 q0av       = DBLARR(jmax+2)
 sigma      = DBLARR(jmax+2,lmax)
 rplot      = DBLARR(jmax+2) 
 qminarr    = DBLARR(jmax+2)
 dphi       = 2.d0*pi/FLOAT(lmax)
 stringtime = ""
 linetype   = 0
 qmin       = 0L
 qsub       = 0
 den        = 0d0
 elost      = 0d0
 sound      = 0d0 
 delt       = 0d0
 ommax      = 0d0
 tmassini   = 0d0

 sigmavg    = DBLARR(jmax+2)
 dphi       = 2.d0*pi/FLOAT(lmax)
 height     = DBLARR(jmax+2,lmax)
 heightavg  = DBLARR(jmax+2)
 result   = DBLARR(2)
 


 titlestr = ""
 modelstr = ""
 fitstr   = ""
 slopestr = ""

 DEVICE,DECOMPOSED=0
 
 WINDOW,0

 getfile = STRING(prefix,".",filenum,FORMAT = '(A,A1,I8.8)') 

 IF(swap_endian eq 1) THEN BEGIN
    OPENR,lun1,getfile,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
 ENDIF ELSE BEGIN
    OPENR,lun1,getfile,/F77_UNFORMATTED,/GET_LUN
 ENDELSE    

 READU,lun1,s
 READU,lun1,t
 READU,lun1,a
 READU,lun1,rho
 READU,lun1,eps
 READU,lun1,dr,dz,delt,time,elost,den,sound,jreq,ommax
 READU,lun1,tmassini,tmass,tmassadd,tmassout,tmassacc

 FREE_LUN,lun1

 conv       = au/FLOAT(jreq-1)
 sigconv    = (conv/dz)^(-2.)*(1.d0/(1.d0-tmassini))*8.87d6

 print, dz, jreq
 print, "tmassini = ", tmassini


 FOR j=1, jmax+1 DO BEGIN
    rplot(j)=((FLOAT(j-1)+0.5d0))*au/(FLOAT(jreq-2))
    rhf(j) = (FLOAT(J)-0.5d0)*dr
 ENDFOR

 rconv = au/rhf(jreq-1)

 FOR j=1,jmax-1 DO BEGIN
    FOR l=0,lmax-1 DO BEGIN
       sigma(j,l)  = 0.d0
       height(j,l) = 0.d0
       FOR k = 1, kmax-1 DO BEGIN
          sigma(j,l)=sigma(j,l)+rho(j,k,l)*dz
       ENDFOR
       rhomid = rho(j,2,l)-((rhf(2)/(rhf(1)-rhf(2)))*(rho(j,1,l)-rho(j,2,l)))
       height(j,l) = sigma(j,l)/(2.d0*rhomid)
    ENDFOR
 ENDFOR
   
 FOR j=1,jmax-1 DO BEGIN
    sigmavg(j)   = 0.d0
    heightavg(j) = 0.d0
    FOR l=0,lmax-1 DO BEGIN
       sigmavg(j)=sigmavg(j)+sigma(j,l)
       heightavg(j)=heightavg(j)+height(j,l)
    ENDFOR
    sigmavg(j)=2.0*sigconv*sigmavg(j)/FLOAT(lmax)
    IF(overplot ne 1) THEN BEGIN
;       sigmavg(j)=alog10(sigmavg(j))
    ENDIF
 ENDFOR
 heightavg = heightavg*rconv/FLOAT(lmax)/rplot
 print, "sig r 5 = ", sigmavg(32)/sigconv
 print, "sig r 10 = ", sigmavg(62)/sigconv


 IF(overplot eq 1) THEN BEGIN
    FOR j=1,jmax-1 DO BEGIN
       sigmavg(j)=alog10(sigmavg(j))
       IF(rplot(j) gt startfit and startjfit eq 0) THEN BEGIN
          startjfit = FIX(j)
       ENDIF
       IF(rplot(j) gt endfit and endjfit eq 0) THEN BEGIN
          endjfit = FIX(j)
       ENDIF
    ENDFOR
         
    fitsize  = endjfit-startjfit
    rfit     = DBLARR(fitsize)
    sigmafit = DBLARR(fitsize)
    sigmaover= DBLARR(jmax+2)
    prob     = 0d0
    chisq    = 0d0
    
    FOR j=startjfit,endjfit-1 DO BEGIN
       rfit(j-startjfit)=alog10(rplot(j))
       sigmafit(j-startjfit)=sigmavg(j)
    ENDFOR

    result = LINFIT(rfit,sigmafit,CHISQ=chisq,PROB=prob,$
                    SIGMA=error)

    print, result, chisq,prob,error
         
    FOR j=startjfit,endjfit DO BEGIN
       sigmaover(j)=result(0)+result(1)*(alog10(rplot(j)))   
    ENDFOR
 ENDIF

 titlestr = string("T = ",time/tconv,FORMAT='(A,f5.2)')
 fitstr   = string("Sigma fit from ",startfit," to ",endfit," AU",   $
                   FORMAT= '(A,f4.1,A,f4.1,A)')
 slopestr = string("Slope = ",result(1),FORMAT='(A,f6.3)')
 minsig=MIN(sigmavg)
 maxsig=MAX(sigmavg)
 PLOT,rplot(*),sigmavg(*),XSTYLE=1,YSTYLE=1,     $
      XRANGE=[rmin,rmax],     $
      YRANGE=[0,1.1*maxsig],XTITLE=['R(AU)'],       $
      YTITLE=['log !I10!N !7R !3(g/cm!E2!N)'],      $
      YCHARSIZE=2,XCHARSIZE=2,YMARGIN=[7,2],color=1
 
IF(overplot eq 1) THEN BEGIN
    oplot,rplot(*),sigmaover(*),MIN_VALUE=1,LINESTYLE=2,color=1
 ENDIF
;     XYOUTS,(0.78*sigmazoom*rplot(jmax)),maxsig,titlestr,/DATA
 XYOUTS,(0.78*sigmazoom*rplot(jmax)),1.15*maxsig,titlestr,/DATA
 XYOUTS,(0.78*sigmazoom*rplot(jmax)),1.10*maxsig,fitstr,/DATA
 XYOUTS,(0.78*sigmazoom*rplot(jmax)),1.05*maxsig,slopestr,/DATA

 WINDOW,1
 
 print, heightavg

 PLOT,rplot(*),heightavg(*),XSTYLE=1,YSTYLE=1,XRANGE=[rmin,rmax],color=1,YRANGE=[0,0.1]

; image=tvrd()
; outfile  = STRING('sigma.',rend,".jpg",FORMAT='(A,I0,A4)')
; WRITE_JPEG,outfile,image,quality=100

 END
