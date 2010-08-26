PRO torquemode

DEVICE,DECOMPOSED=0
 
 JMAX2    = 0L
 LMAX2    = 0L
 NUMFILES = 0L
 COUNTER  = 0L
 ISTART   = 0L
 IEND     = 0L
 ISKIP    = 0L
 INFILE   = ''
 gammag   = 0.d0
 gammang  = 0.d0
 time1    = 0.d0
 time2    = 0.d0
 dummy    = 0.

 INFILE = 'indirectdecompP0.5dat.300000'

 IF(swap_endian eq 1) THEN BEGIN
    OPENR,lun1,INFILE,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
 ENDIF ELSE BEGIN
    OPENR,lun1,INFILE,/F77_UNFORMATTED,/GET_LUN
 ENDELSE

 readu,lun1,JMAX2,LMAX2,NUMFILES,COUNTER
 readu,lun1,ISTART,IEND,ISKIP
 timearr   = DBLARR(NUMFILES)
 tgrav     = DBLARR(JMAX2,LMAX2,NUMFILES)
 treyn     = DBLARR(JMAX2,NUMFILES)
 alphagrav = DBLARR(JMAX2,LMAX2,NUMFILES)
 alphareyn = DBLARR(JMAX2,NUMFILES)
 alphasum  = DBLARR(JMAX2,LMAX2,NUMFILES)

 readu,lun1,timearr
 readu,lun1,tgrav
 readu,lun1,treyn
 readu,lun1,alphagrav
 readu,lun1,alphareyn
 readu,lun1,alphasum
 CLOSE,lun1
 FREE_LUN,lun1

 sub = WHERE((timearr ge time1) and (timearr le time2),subcount)

 timearr   = timearr(sub)
 tgrav     = tgrav(*,*,sub
 treyn     = treyn(*,sub)
 alphagrav = alphagrav(*,*,sub) 
 alphareyn = alphareyn(*,sub) 
 alphasum  = alphasum(*,sub)  

 print, "Of ",NUMFILES,"files ",COUNTER," were used."
 print, "There were ",subcount," files in the time interval ",time1," to ",time2
 print, "LMAX2 = ",LMAX2,"JMAX2 = ",JMAX2
 print, "Istart, Iend, Iskip: ",ISTART,IEND,ISKIP
 tgravavg = DBLARR(JMAX2,LMAX2)
 tgravavg = TOTAL(tgrav,3,/DOUBLE)
 tgravavg  = tgravavg/DOUBLE(COUNTER)

 tgravavg  = tgravavg/1.d39
 tgrav1    = tgravavg(*,1)
 tgrav12   = tgravavg(*,1)+tgravavg(*,2)
 tgrav123  = tgravavg(*,1)+tgravavg(*,2)+tgravavg(*,3)
 tgrav1234 = tgravavg(*,1)+tgravavg(*,2)+tgravavg(*,3)+tgravavg(*,4)

; alphagrav    = ALOG10(alphagrav+1.d-12)
; alphagravavg = ALOG10(ABS(alphagravavg)+1.d-12)

 r = dblarr(jmax2)
 rhf = dblarr(jmax2)
 x1=3.d0
 x2=80.d0
 dz=double(au/jreq)
 
 r(0)=-dz
 
 FOR j=1,jmax2-1 DO BEGIN
    r(j)=r(j-1)+(dz)
ENDFOR


 PLOT,r(20:jmax2-3),tgravavg(20:jmax2-3,0),XSTYLE=1,XRANGE=[10,55], $
     YTITLE = 'x 10!E39!N erg', XTITLE ="AU",YCHARSIZE=1.4,  $
     YSTYLE=1,YRANGE=[0,7],color=1,background=255

 OPLOT,r(20:jmax2-3),tgrav1(20:jmax2-3),LINESTYLE=1,color=1
 OPLOT,r(20:jmax2-3),tgrav12(20:jmax2-3),LINESTYLE=2,color=1
 OPLOT,r(20:jmax2-3),tgrav123(20:jmax2-3),LINESTYLE=3,color=1
 OPLOT,r(20:jmax2-3),tgrav1234(20:jmax2-3),LINESTYLE=4,color=1

 SET_PLOT,"PS"
 DEVICE, FILENAME="torque_indirect_nobdy.eps", /COLOR,ENCAPSULATED=1

 PLOT,r(20:jmax2-3),tgravavg(20:jmax2-3,0),XSTYLE=1,XRANGE=[10,55], $
     YTITLE = 'x 10!E39!N erg', XTITLE ="AU",YCHARSIZE=1.4,  $
     YSTYLE=1,YRANGE=[0,7],color=1,background=255

 OPLOT,r(20:jmax2-3),tgrav1(20:jmax2-3),LINESTYLE=1
 OPLOT,r(20:jmax2-3),tgrav12(20:jmax2-3),LINESTYLE=2
 OPLOT,r(20:jmax2-3),tgrav123(20:jmax2-3),LINESTYLE=3
 OPLOT,r(20:jmax2-3),tgrav1234(20:jmax2-3),LINESTYLE=4


 END
 



 munch = " "

titlestr = ""
time1 = time1/tconv
time2 = time2/tconv

titlestr = string("Alpha between ",time1," and ",time2,FORMAT='(A,f5.2,A,f5.2)')

r = dblarr(jmax2)
rhf = dblarr(jmax2)
x1=3.d0
x2=80.d0
dz=double(au/jreq)

r(0)=-dz

for j=1,jmax2-1 DO BEGIN
r(j)=r(j-1)+(dz)
rhf(j)=r(j)+0.5d0*(dz)
endfor

gammag = 3.d0-(2.d0/gamma)
gammang= (3.d0*gamma-1)/(gamma+1.d0)

 FOR i = 0, jmax2-1 DO BEGIN
;  READF,1,munch,FORMAT='(A256)'
;  READS,munch,dummy,dummy,FORMAT='(F16,F16)'
  tcomega(i) = 2.d0*tconv*omegafc(i)
 ENDFOR

alphapg(*) = 4./(9.*gammag*(gammag-1)*tcomega(*))
alphapng(*)= 4.d0/(9.d0*gammang*(gammang-1)*tcomega(*))

grav(*) = grav(*)/1.d39
reyn(*) = reyn(*)/1.d39
mdott(*)=mdott(*)/1.d39
mdot(*) = mdot(*)*3.156d7*1.d6

alphar(*) = ALOG10(ABS(alphar(*)+1.d-12))
alphas(*) = ALOG10(ABS(alphas(*)+1.d-12))
alpham(*) = ALOG10(ABS(alpham(*)+1.d-12))
alphapg(*)= ALOG10(ABS(alphapg(*)+1.d-12))
alphapng(*)=ALOG10(ABS(alphapng(*)+1.d-12))

tevolve = 0.d0
massev  = 0.d0
FOR j = 1, jmax2-1 DO BEGIN

 IF ( r(j) GE 20. ) AND ( r(j) LE 26. ) THEN BEGIN

  tevolve = tevolve + ABS(mdot(J)*(masscy2(J)-masscy2(J-1))) 
  massev  = massev  + masscy2(J)-masscy2(J-1)

 ENDIF
ENDFOR

PRINT," evolution time estimate, ",massev^2/tevolve,(masscy2(240)-masscy2(120))/1.d-7,massev/1.d-7

sum(*)=grav(*) - reyn(*) 
alpha06(*) = -1.22

LOADCT,ctable

WINDOW,0
PLOT,r(20:jmax2-3),grav(20:jmax2-3),XSTYLE=1,YRANGE=[-1.5,5],XRANGE=[x1,50],$
  YCHARSIZE=1.4,THICK=3, YTITLE = "erg (X10!E39!N)", XTITLE="AU",YSTYLE=1
OPLOT,r(20:jmax2-3),sum(20:jmax2-3)
OPLOT,r(20:jmax2-3),mdott(20:jmax2-3),color=200
OPLOT,r(20:jmax2-3),alpha06(20:jmax2-3),LINESTYLE=3
;PLOT,r(20:jmax2-3),mdot(20:jmax2-3),/NOERASE,/XSTYLE,YSTYLE=4,XRANGE=[x1,x2],LINESTYLE=2



WINDOW,1

PLOT,r(20:(jmax2-3)),alphas(20:(jmax2-3)),XSTYLE=1, $
     YTITLE = "log(abs[!7a!3])", XTITLE ="AU",YCHARSIZE=1.4,  $
     YRANGE=[-4,1],XRANGE=[5,45]
OPLOT,r(20:(jmax2-3)),alphag(20:(jmax2-3)),THICK=3
OPLOT,r(20:(jmax2-3)),alpham(20:(jmax2-3)),color=200
OPLOT,r(20:(jmax2-3)),alphapg(20:(jmax2-3)),THICK=2, LINESTYLE=2
OPLOT,r(20:(jmax2-3)),alphapng(20:(jmax2-3)),THICK=2, LINESTYLE=2
OPLOT,r(20:jmax2-3),alpha06(20:jmax2-3),LINESTYLE=3


WINDOW,2
PLOT,r(20:(jmax2-3)),masscy1(20:(jmax2-3))

WINDOW,3
PLOT,r(20:(jmax2-3)),mdot(20:(jmax2-3))



SET_PLOT, "PS"

DEVICE, FILENAME="torque_paperiii.eps", /COLOR,ENCAPSULATED=1

PLOT,r(20:jmax2-3),grav(20:jmax2-3),/XSTYLE,YRANGE=[-1.5,1.5],XRANGE=[x1,x2],YCHARSIZE=1.4, $
     THICK=3, YTITLE = "erg (X10!E39!N)", XTITLE="AU"
OPLOT,r(20:jmax2-3),sum(20:jmax2-3)
OPLOT,r(20:jmax2-3),mdott(20:jmax2-3),color=200
PLOT,r(20:jmax2-3),mdot(20:jmax2-3),/NOERASE,/XSTYLE,YSTYLE=4,XRANGE=[x1,x2],LINESTYLE=2
PLOTS,x1,0,/DATA
PLOTS,x2,0,/DATA,/CONTINUE,LINESTYLE=3


DEVICE, FILENAME="alphatcool2burst.eps", ENCAPSULATED=1

PLOT,r(20:(jmax2-3)),alphag(20:(jmax2-3)),/XSTYLE, $
     YTITLE = "log(abs[!7a!3])", XTITLE ="AU",YCHARSIZE=1.4,YRANGE=[-4,1]
OPLOT,r(20:(jmax2-3)),alphapg(20:(jmax2-3)),THICK=2, LINESTYLE=2
;OPLOT,r(20:(jmax2-3)),alphag(20:(jmax2-3)),THICK=3
;OPLOT,r(20:(jmax2-3)),alpham(20:(jmax2-3)),color=200
;OPLOT,r(20:(jmax2-3)),alphap(20:(jmax2-3)),THICK=2, LINESTYLE=2
;PLOTS,r(20),-1.6,/DATA
;PLOTS,r(jmax2-3),-1.6,/DATA,/CONTINUE,LINESTYLE=3
;PLOTS,r(20),-1.9,/DATA
;PLOTS,r(jmax2-3),-1.9,/DATA,/CONTINUE,LINESTYLE=3

DEVICE, FILENAME="alphaburst.eps",ENCAPSULATED=1

PLOT,r(20:(jmax2-3)),alphas(20:(jmax2-3)),XSTYLE=1, $
     YTITLE = "log(abs[!7a!3])", XTITLE ="AU",YCHARSIZE=1.4,  $
     YRANGE=[-3,-0.5],XRANGE=[10,55],TITLE=titlestr
OPLOT,r(20:(jmax2-3)),alphag(20:(jmax2-3)),THICK=3
OPLOT,r(20:(jmax2-3)),alpham(20:(jmax2-3)),color=200
OPLOT,r(20:(jmax2-3)),alphapg(20:(jmax2-3)),THICK=2, LINESTYLE=2
OPLOT,r(20:(jmax2-3)),alphapng(20:(jmax2-3)),THICK=2, LINESTYLE=2
OPLOT,r(20:jmax2-3),alpha06(20:jmax2-3),LINESTYLE=3


return
END
