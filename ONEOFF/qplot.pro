; hopefully this program will do q plots

 PRO qplot
 
 colorcount = 0L
 colorstep  = 0L
 datadir    = '../../INITCOND/RAND_PERT/QPLOT/'
 prefix     = 'saved'
 rendstart  = 10000L
 rendend    = 30000L
 rendskip   = 10000L
 numfiles   = (rendend-rendstart)/rendskip + 1
 SWAP_ENDIAN= 1

 time1      = 12.d0
 time2      = 19.d0
 
 jmax       = 256
 kmax       = 64
 lmax       = 128
 au         = 40.d0

 xmin       = 10.d0
 xmax       = 40.d0
 ymin       = 0.d0
 ymax       = 4.d0
 gamma      = 1.666666666666667d0
 dphi       = 0d0
 delt       = 0d0
 time       = 0d0
 elost      = 0d0
 den        = 0d0
 sound      = 0d0
 ommax      = 0d0
 tmassini   = 0d0
 tmass      = 0d0
 tmassadd   = 0d0
 tmassout   = 0d0
 tmassacc   = 0d0
 totcool    = 0d0
 totdflux   = 0d0
 totheat    = 0d0
 totirr     = 0d0
 rauconv    = 0d0
 massin     = 0d0
 mconv      = 0d0
 dr         = 0d0
 dz         = 0d0
 jmax2      = jmax+2
 kmax2      = kmax+2
 tconv      = 1605.63d0
 pi         = 3.14159265358979323846264338327950d0

 r          = DBLARR(jmax2)
 rhf        = DBLARR(jmax2)
 q          = DBLARR(jmax2,numfiles)
 csound     = DBLARR(jmax2,lmax)
 csoundavg  = DBLARR(jmax2)
 kappa      = DBLARR(jmax2)
 omega      = DBLARR(jmax2,kmax2,lmax)
 omegavg    = DBLARR(jmax2)
 q0         = DBLARR(jmax2,numfiles)
 sigma      = DBLARR(jmax2,lmax)
 sigmavg    = DBLARR(jmax2)
 rplot      = DBLARR(jmax2) 
 qminarr    = DBLARR(jmax2)
 dphi = 2.d0*pi/FLOAT(lmax)
 colorstep  = 256/(((rendend-rendstart)/rendskip)+1)
 colorstep  = FIX(colorstep)
 colorcount = 1
 xtype      = 0.85*(xmax-xmin)+xmin
 ytype      = 0.95*(ymax-ymin)+ymin
 stringtime = ""
 linetype   = 0
 qmin       = 0L
 qsub       = 0
 counter    = 0L
 jreq       = 0L
 timearr    = DBLARR(numfiles)


 !P.COLOR=0
 !P.BACKGROUND=255
 
 DEVICE,DECOMPOSED=0
 
 LOADCT,39

 FOR file = rendstart,rendend,rendskip DO BEGIN
 
     getfile = STRING(datadir,prefix,".",file,FORMAT='(A,A,A1,I6.6)')
     
;     IF (FILE_TEST(getfile) eq 0) then continue
     
     if (file eq 20000) then BEGIN
         SWAP_ENDIAN = 0
     ENDIF ELSE BEGIN
         SWAP_ENDIAN = 1
     ENDELSE

     if (file eq 10000) then begin
         KMAX = 32
         KMAX2 = 34
         print, "small"
     ENDIF ELSE BEGIN
         KMAX = 64
         KMAX2 = 66
         print, "big"
     ENDELSE

     s          = DBLARR(jmax2,kmax2,lmax)
     t          = DBLARR(jmax2,kmax2,lmax)
     a          = DBLARR(jmax2,kmax2,lmax)
     rho        = DBLARR(jmax2,kmax2,lmax)
     eps        = DBLARR(jmax2,kmax2,lmax)

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

     FREE_LUN, lun1

     print, "Opened file: ",getfile

;     IF ((time/tconv lt time1) or (time/tconv gt time2)) then continue

     

     IF (counter eq 0) then begin
        FOR j = 0,jmax+1 DO BEGIN
           r(j) = FLOAT(j-1)*dr
           rhf(j) = r(j)+0.5d0*dr
        ENDFOR
        rplot(*)=rhf(*)*au/(FLOAT(jreq-1)*dr)
     ENDIF
     
     counter = counter +1
     print, "Counter = ",counter
     
     FOR j=1,jmax+1 DO BEGIN
        csoundavg(j) = 0.d0
        omegavg(j) = 0.d0
        sigmavg(j) = 0.d0
        FOR l=0,lmax-1 DO BEGIN
        sigma(j,l)=0.d0
             FOR k=1,kmax+1 DO BEGIN
                 omega(j,k,l)=a(j,k,l)/rho(j,k,l)/rhf(j)^2
                 sigma(j,l)=sigma(j,l)+rho(j,k,l)*dz
             ENDFOR
             csound(j,l)=(gamma*(gamma-1.d0)*eps(j,2,l)/rho(j,2,l))^.5
             sigma(j,l)=2.d0*sigma(j,l)
             csoundavg(j)=csoundavg(j)+csound(j,l)
             omegavg(j)=omegavg(j)+omega(j,2,l)
             sigmavg(j)=sigmavg(j)+sigma(j,l)
         ENDFOR
     ENDFOR

     csoundavg(*)=csoundavg(*)/FLOAT(LMAX)
     omegavg(*)=omegavg(*)/FLOAT(LMAX)
     sigmavg(*)=sigmavg(*)/FLOAT(LMAX)



     FOR j=1,jmax DO BEGIN
        kappa(j)=(ABS(rhf(j)*omegavg(j)*(omegavg(j+1)- $
                 omegavg(j-1))/dr + 4.d0*omegavg(j)^2))^0.5
;        print,'J = ',J,' Thing 1 = ',rhf(j)*omegavg(j)*(omegavg(j+1)- $
;                 omegavg(j-1))/dr,' Thing 2 = ',4.d0*omegavg(j)^2
        q(j)=kappa(j)*csoundavg(j)/(pi*sigmavg(j))
        q0(j)=omegavg(j)*csoundavg(j)/(pi*sigmavg(j))
     ENDFOR



     time=time/tconv
     stringtime = STRING("T = ",time,FORMAT='(A4,F6.3)')

     FOR j = 25,235 DO BEGIN
         qminarr(j-25) = q(j)
     ENDFOR

     FOR j=210,jmax+1 DO BEGIN
         qminarr(j) = 5.d0
     ENDFOR
 
     qmin = min(qminarr,qsub,/NAN)

;;     print, qmin, qsub,rplot(qsub+25) 
              
     IF(file eq rendstart) THEN BEGIN
         PLOT,rplot(*),q0(*),MAX_VALUE=100,MIN_VALUE=-1,XSTYLE=1, $
           COLOR=colorcount,YSTYLE=1,XRANGE=[xmin,xmax],YRANGE=[ymin,ymax],$
           LINESTYLE=linetype
;         OPLOT,rplot(*),q(*),MAX_VALUE=100,MIN_VALUE=-1,COLOR=colorcount, $
;           LINESTYLE=2
         XYOUTS,xtype,ytype,stringtime,COLOR=colorcount,/DATA
     ENDIF ELSE BEGIN
     
         OPLOT,rplot(*),q0(*),MAX_VALUE=100,MIN_VALUE=-1,COLOR=colorcount,$
         LINESTYLE=linetype
 ;        OPLOT,rplot(*),q(*),MAX_VALUE=100,MIN_VALUE=-1,COLOR=colorcount, $
 ;          LINESTYLE=2
          XYOUTS,xtype,ytype,stringtime,COLOR=colorcount,/DATA
     ENDELSE

     colorcount=colorcount+colorstep     
     ytype = ytype-0.05*(ymax-ymin)
;     linetype++

 ENDFOR

;  q  = q(*,0:counter-1)
;  q0 = q0(*,0:counter-1)

;  qavg = TOTAL(q,2,/DOUBLE)
;  q0avg = TOTAL(q0,2,/DOUBLE)
;  qavg = qavg/counter
;  q0avg = q0avg/counter
;  index = where((rplot gt xmin) and (rplot lt xmax))
;  print, rplot
;  print, index
;  qmean = MEAN(qavg(index),/DOUBLE)

;   print, qmean

 RETURN
 END

