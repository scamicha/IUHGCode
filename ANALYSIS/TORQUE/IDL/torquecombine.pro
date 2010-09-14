PRO torquecombine

 DEVICE,DECOMPOSED=0
 !P.BACKGROUND=255
 !P.COLOR=1

 timestart     = 12.d0
 timestop      = 19.d0
 gamma         = 1.66666666d0

 x1=5.d0
 x2=45.d0

 INFILE1   = 'torqueLMAX64dat.500000'
 INFILE2   = 'torquedecompP1dat.MPI.385007'
 INFILE3   = 'decompP1L256dat.MPI.650000'
 INFILE4   = 'decompP1L512dat.MPI.980000'

 tconv     = 1605.63
 JMAX1     = 0L
 LMAX1     = 0L
 FILES1    = 0L
 COUNT1    = 0L
 START1    = 0L
 STOP1     = 0L
 SKIP1     = 0L

 JMAX2     = 0L
 LMAX2     = 0L
 FILES2    = 0L
 COUNT2    = 0L
 START2    = 0L
 STOP2     = 0L
 SKIP2     = 0L

 JMAX3     = 0L
 LMAX3     = 0L
 FILES3    = 0L
 COUNT3    = 0L
 START3    = 0L
 STOP3     = 0L
 SKIP3     = 0L

 JMAX4     = 0L
 LMAX4     = 0L
 FILES4    = 0L
 COUNT4    = 0L
 START4    = 0L
 STOP4     = 0L
 SKIP4     = 0L

 openr,lun1,INFILE1,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
 openr,lun2,INFILE2,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
 openr,lun3,INFILE3,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN 
 openr,lun4,INFILE4,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN

 readu,lun1,JMAX1,LMAX1,FILES1,COUNT1
 readu,lun1,START1,STOP1,SKIP1

 readu,lun2,JMAX2,LMAX2,FILES2,COUNT2
 readu,lun2,START2,STOP2,SKIP2

 readu,lun3,JMAX3,LMAX3,FILES3,COUNT3
 readu,lun3,START3,STOP3,SKIP3

 readu,lun4,JMAX4,LMAX4,FILES4,COUNT4
 readu,lun4,START4,STOP4,SKIP4
 
 time1     = dblarr(count1)
 grav1     = dblarr(jmax1,lmax1,count1)
 reyn1     = dblarr(jmax1,count1)
 alphag1   = dblarr(jmax1,lmax1,count1)
 alphar1   = dblarr(jmax1,count1)
 alphas1   = dblarr(jmax1,lmax1,count1)
 alphapg1  = dblarr(jmax1)
 tcomega1  = dblarr(jmax1)
 alphapng1 = dblarr(jmax1)
 omega1    = DBLARR(jmax1,count1)
 alpha06   = DBLARR(jmax1)

 time2     = DBLARR(count2)
 grav2     = dblarr(jmax2,lmax2,count2)
 reyn2     = dblarr(jmax2,count2)
 alphag2   = dblarr(jmax2,lmax2,count2)
 alphar2   = dblarr(jmax2,count2)
 alphas2   = dblarr(jmax2,lmax2,count2)
 alphapg2  = dblarr(jmax2)
 tcomega2  = dblarr(jmax2)
 alphapng2 = dblarr(jmax2)
 omega2    = DBLARR(jmax2,count2)
 alpha062  = DBLARR(jmax2)

 time3     = dblarr(count3)
 grav3     = dblarr(jmax3,lmax3,count3)
 reyn3     = dblarr(jmax3,count3)
 alphag3   = dblarr(jmax3,lmax3,count3)
 alphar3   = dblarr(jmax3,count3)
 alphas3   = dblarr(jmax3,lmax3,count3)
 alphapg3  = dblarr(jmax3)
 tcomega3  = dblarr(jmax3)
 alphapng3 = dblarr(jmax3)
 omega3    = DBLARR(jmax3,count3)
 alpha063  = DBLARR(jmax3)

 time4     = dblarr(count4)
 grav4     = dblarr(jmax4,lmax4,count4)
 reyn4     = dblarr(jmax4,count4)
 alphag4   = dblarr(jmax4,lmax4,count4)
 alphar4   = dblarr(jmax4,count4)
 alphas4   = dblarr(jmax4,lmax4,count4)
 alphapg4  = dblarr(jmax4)
 tcomega4  = dblarr(jmax4)
 alphapng4 = dblarr(jmax4)
 omega4    = DBLARR(jmax4,count4)
 alpha064  = DBLARR(jmax4)

 gammag   = 0.d0
 gammang  = 0.d0

 readu,lun1,time1
 readu,lun1,grav1
 readu,lun1,reyn1
 readu,lun1,alphag1
 readu,lun1,alphar1
 readu,lun1,alphas1
 readu,lun1,omega1

 CLOSE,lun1

 readu,lun2,time2
 readu,lun2,grav2
 readu,lun2,reyn2
 readu,lun2,alphag2
 readu,lun2,alphar2
 readu,lun2,alphas2
 readu,lun2,omega2

 CLOSE,lun2

 readu,lun3,time3
 readu,lun3,grav3
 readu,lun3,reyn3
 readu,lun3,alphag3
 readu,lun3,alphar3
 readu,lun3,alphas3
 readu,lun3,omega3

 CLOSE,lun3

 readu,lun4,time4
 readu,lun4,grav4
 readu,lun4,reyn4
 readu,lun4,alphag4
 readu,lun4,alphar4
 readu,lun4,alphas4
 readu,lun4,omega4

 CLOSE,lun4

 print, alphag3(117,0,*)

 dummy = WHERE((time1 gt timestart) and (time1 lt timestop))
 timesub1   = time1(dummy)
 gravsub1   = grav1(*,*,dummy)
 reynsub1   = reyn1(*,dummy)
 alphagsub1 = alphag1(*,*,dummy)
 alpharsub1 = alphar1(*,dummy)
 alphassub1 = alphas1(*,*,dummy)
 omegasub1  = omega1(*,dummy)

 dummy = WHERE((time2 gt timestart) and (time2 lt timestop))
 timesub2   = time2(dummy)
 gravsub2   = grav2(*,*,dummy)
 reynsub2   = reyn2(*,dummy)
 alphagsub2 = alphag2(*,*,dummy)
 alpharsub2 = alphar2(*,dummy)
 alphassub2 = alphas2(*,*,dummy)
 omegasub2  = omega2(*,dummy)

 dummy = WHERE((time3 gt timestart) and (time3 lt timestop))
 timesub3   = time3(dummy)
 gravsub3   = grav3(*,*,dummy)
 reynsub3   = reyn3(*,dummy)
 alphagsub3 = alphag3(*,*,dummy)
 alpharsub3 = alphar3(*,dummy)
 alphassub3 = alphas3(*,*,dummy)
 omegasub3  = omega3(*,dummy)

 dummy = WHERE((time4 gt timestart) and (time4 lt timestop))
 timesub4   = time4(dummy)
 gravsub4   = grav4(*,*,dummy)
 reynsub4   = reyn4(*,dummy)
 alphagsub4 = alphag4(*,*,dummy)
 alpharsub4 = alphar4(*,dummy)
 alphassub4 = alphas4(*,*,dummy)
 omegasub4  = omega4(*,dummy)
 
 dummy = N_ELEMENTS(timesub1)
 grav1 = TOTAL(gravsub1,3,/DOUBLE)
 grav1 = grav1/dummy
 reyn1 = TOTAL(reynsub1,2,/DOUBLE)
 reyn1 = reyn1/dummy
 alphag1 = TOTAL(alphagsub1,3,/DOUBLE)
 alphag1 = alphag1/dummy
 alphar1 = TOTAL(alpharsub1,2,/DOUBLE)
 alphar1 = alphar1/dummy
 alphas1 = TOTAL(alphassub1,3,/DOUBLE)
 alphas1 = alphas1/dummy
 omega1 = TOTAL(omegasub1,2,/DOUBLE)
 omega1 = omega1/dummy

 dummy = N_ELEMENTS(timesub2)
 grav2 = TOTAL(gravsub2,3,/DOUBLE)
 grav2 = grav2/dummy
 reyn2 = TOTAL(reynsub2,2,/DOUBLE)
 reyn2 = reyn2/dummy
 alphag2 = TOTAL(alphagsub2,3,/DOUBLE)
 alphag2 = alphag2/dummy
 alphar2 = TOTAL(alpharsub2,2,/DOUBLE)
 alphar2 = alphar2/dummy
 alphas2 = TOTAL(alphassub2,3,/DOUBLE)
 alphas2 = alphas2/dummy
 omega2 = TOTAL(omegasub2,2,/DOUBLE)
 omega2 = omega2/dummy

 dummy = N_ELEMENTS(timesub3)
 grav3 = TOTAL(gravsub3,3,/DOUBLE)
 grav3 = grav3/dummy
 reyn3 = TOTAL(reynsub3,2,/DOUBLE)
 reyn3 = reyn3/dummy
 alphag3 = TOTAL(alphagsub3,3,/DOUBLE)
 alphag3 = alphag3/dummy
 alphar3 = TOTAL(alpharsub3,2,/DOUBLE)
 alphar3 = alphar3/dummy
 alphas3 = TOTAL(alphassub3,3,/DOUBLE)
 alphas3 = alphas3/dummy
 omega3 = TOTAL(omegasub3,2,/DOUBLE)
 omega3 = omega3/dummy

 dummy = N_ELEMENTS(timesub4)
 grav4 = TOTAL(gravsub4,3,/DOUBLE)
 grav4 = grav4/dummy
 reyn4 = TOTAL(reynsub4,2,/DOUBLE)
 reyn4 = reyn4/dummy
 alphag4 = TOTAL(alphagsub4,3,/DOUBLE)
 alphag4 = alphag4/dummy
 alphar4 = TOTAL(alpharsub4,2,/DOUBLE)
 alphar4 = alphar4/dummy
 alphas4 = TOTAL(alphassub4,3,/DOUBLE)
 alphas4 = alphas4/dummy
 omega4 = TOTAL(omegasub4,2,/DOUBLE)
 omega4 = omega4/dummy

 r = dblarr(jmax1)
 rhf = dblarr(jmax1)

 dz=40.d0/250.d0

 r(0)=-dz

 gammag = 3.d0-(2.d0/gamma)
 gammang= (3.d0*gamma-1)/(gamma+1.d0)

 for j=1,jmax2-1 DO BEGIN
    r(j)=r(j-1)+(dz)
    rhf(j)=r(j)+0.5d0*(dz)
 endfor

 FOR i = 0, jmax2-1 DO BEGIN
;  READF,1,munch,FORMAT='(A256)'
;  READS,munch,dummy,dummy,FORMAT='(F16,F16)'
    tcomega1(i) = 2.d0*tconv*omega1(i)
    tcomega2(i) = 2.d0*tconv*omega2(i)
    tcomega3(i) = 2.d0*tconv*omega3(i)  
    tcomega4(i) = 2.d0*tconv*omega4(i)    
 ENDFOR

 alphapg1(*) = 4./(9.*gammag*(gammag-1)*tcomega1(*))
 alphapng1(*)= 4.d0/(9.d0*gammang*(gammang-1)*tcomega1(*))
 alphapg2(*) = 4./(9.*gammag*(gammag-1)*tcomega2(*))
 alphapng2(*)= 4.d0/(9.d0*gammang*(gammang-1)*tcomega2(*))
 alphapg3(*) = 4./(9.*gammag*(gammag-1)*tcomega3(*))
 alphapng3(*)= 4.d0/(9.d0*gammang*(gammang-1)*tcomega3(*))
 alphapg4(*) = 4./(9.*gammag*(gammag-1)*tcomega4(*))
 alphapng4(*)= 4.d0/(9.d0*gammang*(gammang-1)*tcomega4(*))
 
 dummy = WHERE((r gt x1) and (r lt x2))

 print, dummy
 print, r(dummy)
 print, alphag3(dummy,0)

 alphag1avg = MEAN(alphag1(dummy,0))
 alphag2avg = MEAN(alphag2(dummy,0))
 alphag3avg = MEAN(alphag3(dummy,0))
 alphag4avg = MEAN(alphag4(dummy,0))
 alphapg1avg= MEAN(alphapg1(dummy,0))
 alphapng1avg=MEAN(alphapng1(dummy,0))
 
 alphag1(*) = ALOG10(ABS(alphag1(*)+1.d-12))
 alphag2(*) = ALOG10(ABS(alphag2(*)+1.d-12))
 alphag3(*) = ALOG10(ABS(alphag3(*)+1.d-12))
 alphag4(*) = ALOG10(ABS(alphag4(*)+1.d-12))
 alphapg1(*)= ALOG10(ABS(alphapg1(*)+1.d-12))
 alphapng1(*)=ALOG10(ABS(alphapng1(*)+1.d-12))
 alphapg2(*)= ALOG10(ABS(alphapg2(*)+1.d-12))
 alphapng2(*)=ALOG10(ABS(alphapng2(*)+1.d-12))
 alphapg3(*)= ALOG10(ABS(alphapg3(*)+1.d-12))
 alphapng3(*)=ALOG10(ABS(alphapng3(*)+1.d-12))
 alphapg4(*)= ALOG10(ABS(alphapg4(*)+1.d-12))
 alphapng4(*)=ALOG10(ABS(alphapng4(*)+1.d-12))
 
print,"alpha  64 = ",alphag1avg
print,"alpha 128 = ",alphag2avg
print,"alpha 256 = ",alphag3avg
print,"alpha 512 = ",alphag4avg
print,"alpha gamng = ",alphapng1avg
print,"alpha gamg  = ",alphapg1avg

alpha06(*) = -1.22


LOADCT,39

WINDOW,0

PLOT,r(20:(jmax2-3)),alphag1(20:(jmax2-3),0),XSTYLE=1, $
     YTITLE = "log(abs[!7a!3])", XTITLE ="AU",YCHARSIZE=1.4,  $
    XRANGE=[5,45],THICK=3; YRANGE=[-4,1],XRANGE=[5,45]
OPLOT,r(20:(jmax2-3)),alphag2(20:(jmax2-3),0),THICK=2
OPLOT,r(20:(jmax2-3)),alphag3(20:(jmax2-3),0),THICK=1
OPLOT,r(20:(jmax2-3)),alphapg2(20:(jmax2-3)),THICK=1,linestyle=2
OPLOT,r(20:(jmax2-3)),alphapng2(20:(jmax2-3)),THICK=1,linestyle=2
OPLOT,r(20:(jmax2-3)),alpha06(20:(jmax2-3)),LINESTYLE=3
OPLOT,r(20:(jmax2-3)),alphag4(20:(jmax2-3),0),THICK=7
;OPLOT,r(20:(jmax2-3)),alpham(20:(jmax2-3)),color=200
;OPLOT,r(20:(jmax2-3)),alphapg(20:(jmax2-3)),THICK=2, LINESTYLE=2
;OPLOT,r(20:(jmax2-3)),alphapng(20:(jmax2-3)),THICK=2, LINESTYLE=2
;OPLOT,r(20:jmax2-3),alpha06(20:jmax2-3),LINESTYLE=3


;; SET_PLOT, "PS"

;; DEVICE, FILENAME="torquecombine.eps", /COLOR,ENCAPSULATED=1,/BOLD

;; PLOT,r(20:(jmax2-3)),alphag1(20:(jmax2-3)),XSTYLE=1, $
;;      YTITLE = "!5log(abs[!7a!5])", XTITLE ="!5AU",YCHARSIZE=1.4,  $
;;     XRANGE=[10,40],YRANGE=[-3,-0.5],THICK=6,COLOR=1, $; YRANGE=[-4,1],XRANGE=[5,45]
;;   TITLE="!5Effective !7a!5 averaged from 12-19 ORPs",XTHICK=4,YTHICK=4
;; OPLOT,r(20:(jmax2-3)),alphag2(20:(jmax2-3)),THICK=6,COLOR=75
;; OPLOT,r(20:(jmax2-3)),alphag3(20:(jmax2-3)),THICK=6,COLOR=150
;; OPLOT,r(20:(jmax2-3)),alphag4(20:(jmax2-3)),THICK=6,COLOR=225
;; OPLOT,r(20:(jmax2-3)),alphapg1(20:(jmax2-3)),THICK=1,linestyle=2
;; OPLOT,r(20:(jmax2-3)),alphapng1(20:(jmax2-3)),THICK=1,linestyle=2
;; OPLOT,r(20:(jmax2-3)),alpha06(20:(jmax2-3)),LINESTYLE=3
;; XYOUTS,30,-0.6,"!5LMAX=64  Avg !7a!5: 0.059",color=1
;; XYOUTS,30,-0.685,"!5LMAX=128 Avg !7a!5: 0.040",color=75
;; XYOUTS,30,-0.765,"!5LMAX=256 Avg !7a!5: 0.024",color=150
;; XYOUTS,30,-0.845,"!5LMAX=512 Avg !7a!5: 0.024",color=225
;; ;OPLOT,r(20:(jmax2-3)),alpham(20:(jmax2-3)),color=200
;; ;OPLOT,r(20:(jmax2-3)),alphapg(20:(jmax2-3)),THICK=2, LINESTYLE=2
;; ;OPLOT,r(20:(jmax2-3)),alphapng(20:(jmax2-3)),THICK=2, LINESTYLE=2
;; ;OPLOT,r(20:jmax2-3),alpha06(20:jmax2-3),LINESTYLE=3

return
END
