; hopefully this program will do q plots

 PRO tcomega
 
 SWAP_ENDIAN= 1

 jmax       = 256
 kmax       = 32
 lmax       = 128
 au         = 40.d0

 xmin       = 20.d0
 xmax       = 35.d0
 ymin       = 0.d0
 ymax       = 20.d0
 gamma      = 1.666666666666667d0
 tcool      = 0.6d0
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

 s          = DBLARR(jmax2,kmax2,lmax)
 t          = DBLARR(jmax2,kmax2,lmax)
 a          = DBLARR(jmax2,kmax2,lmax)
 rho        = DBLARR(jmax2,kmax2,lmax)
 eps        = DBLARR(jmax2,kmax2,lmax)
 r          = DBLARR(jmax2)
 rhf        = DBLARR(jmax2)
 omega      = DBLARR(jmax2,lmax)
 omegavg    = DBLARR(jmax2)
 rplot      = DBLARR(jmax2) 
 dphi = 2.d0*pi/FLOAT(lmax)
 stringtime = ""
 getfile    = ""
 teststring = "./SAVED/saved.000007"
 linetype   = 0
 jreq       = 0L
 first_time = 1L

 tcool = tcool*tconv

 !P.COLOR=0
 !P.BACKGROUND=255
 
 DEVICE,DECOMPOSED=0
 
 LOADCT,39
 ROUND2: print,"going on"
 IF(first_time eq 1) then begin
    getfile = STRING("./SAVED/saved.000007",FORMAT='(A)')
 ENDIF ELSE BEGIN
    getfile = STRING("../../G1.4TC2/SAVED/saved.000001",FORMAT='(A)')
    SWAP_ENDIAN = 0
    tcool = 1.d0*tconv
 ENDELSE
     
;; IF (FILE_TEST(getfile) eq 0) then continue
 
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
 
 FOR j = 0,jmax+1 DO BEGIN
    r(j) = FLOAT(j-1)*dr
    rhf(j) = r(j)+0.5d0*dr
 ENDFOR
 rplot(*)=rhf(*)*au/(FLOAT(jreq-1)*dr)
  
     
 FOR j=1,jmax+1 DO BEGIN
    omegavg(j) = 0.d0
    FOR l=0,lmax-1 DO BEGIN
       omega(j,l)=a(j,1,l)/rho(j,1,l)/rhf(j)^2
       omegavg(j)=omegavg(j)+omega(j,l)
    ENDFOR
 ENDFOR
 omegavg(*)=omegavg(*)/FLOAT(LMAX)
 tcomega = omegavg*tcool


 IF(first_time eq 1) THEN BEGIN
    PLOT,rplot(*),tcomega(*),XSTYLE=1, $
         COLOR=1,YSTYLE=1,XRANGE=[xmin,xmax],YRANGE=[ymin,ymax],$
         LINESTYLE=linetype
    first_time = 0L
    GOTO,ROUND2
 ENDIF
 OPLOT,rplot(*),tcomega(*),color=225

 
     
;;          OPLOT,rplot(*),q0(*),MAX_VALUE=100,MIN_VALUE=-1,COLOR=colorcount,$
;;          LINESTYLE=linetype
;;          OPLOT,rplot(*),q(*),MAX_VALUE=100,MIN_VALUE=-1,COLOR=colorcount, $
;;            LINESTYLE=2
;; ;         XYOUTS,xtype,ytype,stringtime,COLOR=colorcount,/DATA
;;      ENDELSE

     ;; colorcount=colorcount+colorstep     
     ;; ytype = ytype-0.05*(ymax-ymin)
     ;; linetype++

 RETURN
 END

