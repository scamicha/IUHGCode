PRO amdatcomb

samprate = 1L

mmax1   = 0L
count1  = 0L
mmax2   = 0L
count2  = 0L
mmax3   = 0L
count3  = 0L
mmax4   = 0L
count4  = 0L
amtotal1= 0.d0
amtotal2= 0.d0
amtotal3= 0.d0
amtotal4= 0.d0
amrel1  = 0.d0
amrel2  = 0.d0
amrel3  = 0.d0
amrel4  = 0.d0
dummy   = 0.d0
time1   = 12.d0
time2   = 18.d0

INFILE1 = ''
INFILE2 = ''
INFILE3 = ''
INFILE4 = ''

INFILE1 = STRING("P1am2AUL64.dat")
INFILE2 = STRING("P1am2AUL128.dat")
INFILE3 = STRING("P1am2AUL256.dat")
INFILE4 = STRING("P1am2AUL512.dat")
openr,lun1,INFILE1,/GET_LUN,/F77_UNFORMATTED,/SWAP_ENDIAN
readu,lun1,mmax1,count1
avgam1 = DBLARR(mmax1,count1)
avgammid1 = DBLARR(mmax1,count1)
timearr1 = DBLARR(count1)
readu,lun1,avgam1
readu,lun1,avgammid1
readu,lun1,timearr1
free_lun,lun1

openr,lun2,INFILE2,/GET_LUN,/F77_UNFORMATTED,/SWAP_ENDIAN
readu,lun2,mmax2,count2
avgam2 = DBLARR(mmax2,count2)
avgammid2 = DBLARR(mmax2,count2)
timearr2 = DBLARR(count2)
readu,lun2,avgam2
readu,lun2,avgammid2
readu,lun2,timearr2
free_lun,lun2

openr,lun3,INFILE3,/GET_LUN,/F77_UNFORMATTED,/SWAP_ENDIAN
readu,lun3,mmax3,count3
avgam3 = DBLARR(mmax3,count3)
avgammid3 = DBLARR(mmax3,count3)
timearr3 = DBLARR(count3)
readu,lun3,avgam3
readu,lun3,avgammid3
readu,lun3,timearr3
free_lun,lun3

openr,lun4,INFILE4,/GET_LUN,/F77_UNFORMATTED,/SWAP_ENDIAN
readu,lun4,mmax4,count4
avgam4 = DBLARR(mmax4,count4)
avgammid4 = DBLARR(mmax4,count4)
timearr4 = DBLARR(count4)
readu,lun4,avgam4
readu,lun4,avgammid4
readu,lun4,timearr4
free_lun,lun4

timeindex1 = WHERE((timearr1 gt time1)and(timearr1 lt time2))

amtot1 = DBLARR(mmax1)
amdev1 = DBLARR(mmax1)
amerrplus1 = DBLARR(mmax1)
amerrminus1 = DBLARR(mmax1)

FOR i = 1,mmax1 DO BEGIN
    amtot1(i-1) = MEAN(avgam1(i-1,timeindex1),/DOUBLE) 
    amdev1(i-1) = STDDEV(avgam1(i-1,timeindex1),/DOUBLE)
    amerrplus1(i-1) = amtot1(i-1)+amdev1(i-1)
    amerrminus1(i-1) = amtot1(i-1)-amdev1(i-1)
ENDFOR

 amtotal1 = TOTAL(amtot1(1:mmax1-1),/DOUBLE)
 dummy    = TOTAL(amtot1(1:6),/DOUBLE)
 amrev1   = dummy/amtotal1

 amtot1 = ALOG10(amtot1)
 amerrplus1 = ALOG10(amerrplus1)
 amerrminus1 = ALOG10(amerrminus1)

timeindex2 = WHERE((timearr2 gt time1)and(timearr2 lt time2))

amtot2 = DBLARR(mmax2)
amdev2 = DBLARR(mmax2)
amerrplus2 = DBLARR(mmax2)
amerrminus2 = DBLARR(mmax2)

FOR i = 1,mmax2 DO BEGIN
    amtot2(i-1) = MEAN(avgam2(i-1,timeindex2),/DOUBLE) 
    amdev2(i-1) = STDDEV(avgam2(i-1,timeindex2),/DOUBLE)
    amerrplus2(i-1) = amtot2(i-1)+amdev2(i-1)
    amerrminus2(i-1) = amtot2(i-1)-amdev2(i-1)
ENDFOR

 amtotal2 = TOTAL(amtot2(1:mmax2-1),/DOUBLE)
 dummy    = TOTAL(amtot2(1:6),/DOUBLE)
 amrev2   = dummy/amtotal2

 amtot2 = ALOG10(amtot2)
 amerrplus2 = ALOG10(amerrplus2)
 amerrminus2 = ALOG10(amerrminus2)

timeindex3 = WHERE((timearr3 gt time1)and(timearr3 lt time2))

amtot3 = DBLARR(mmax3)
amdev3 = DBLARR(mmax3)
amerrplus3 = DBLARR(mmax3)
amerrminus3 = DBLARR(mmax3)

FOR i = 1,mmax3 DO BEGIN
    amtot3(i-1) = MEAN(avgam3(i-1,timeindex3),/DOUBLE) 
    amdev3(i-1) = STDDEV(avgam3(i-1,timeindex3),/DOUBLE)
    amerrplus3(i-1) = amtot3(i-1)+amdev3(i-1)
    amerrminus3(i-1) = amtot3(i-1)-amdev3(i-1)
ENDFOR

 amtotal3 = TOTAL(amtot3(1:mmax3-1),/DOUBLE)
 dummy    = TOTAL(amtot3(1:6),/DOUBLE)
 amrev3   = dummy/amtotal3

 amtot3 = ALOG10(amtot3)
 amerrplus3 = ALOG10(amerrplus3)
 amerrminus3 = ALOG10(amerrminus3)

timeindex4 = WHERE((timearr4 gt time1)and(timearr4 lt time2))

amtot4 = DBLARR(mmax4)
amdev4 = DBLARR(mmax4)
amerrplus4 = DBLARR(mmax4)
amerrminus4 = DBLARR(mmax4)

FOR i = 1,mmax4 DO BEGIN
    amtot4(i-1) = MEAN(avgam4(i-1,timeindex4),/DOUBLE) 
    amdev4(i-1) = STDDEV(avgam4(i-1,timeindex4),/DOUBLE)
    amerrplus4(i-1) = amtot4(i-1)+amdev4(i-1)
    amerrminus4(i-1) = amtot4(i-1)-amdev4(i-1)
ENDFOR

 amtotal4 = TOTAL(amtot4(1:mmax4-1),/DOUBLE)
 dummy    = TOTAL(amtot4(1:6),/DOUBLE)
 amrev4   = dummy/amtotal4

 amtot4 = ALOG10(amtot4)
 amerrplus4 = ALOG10(amerrplus4)
 amerrminus4 = ALOG10(amerrminus4)

XAXIS = DINDGEN(mmax4)
XAXIS = XAXIS+1
XAXIS = ALOG10(XAXIS)
 
print,"Total Am LMAX = 64: ",amtotal1
print,"Total Am LMAX = 128: ",amtotal2
print,"Total Am LMAX = 256: ",amtotal3
print,"Total Am LMAX = 512: ",amtotal4
print,"Relative Am LMAX = 64: ",amrev1
print,"Relative Am LMAX = 128: ",amrev2
print,"Relative Am LMAX = 256: ",amrev3
print,"Relative Am LMAX = 512: ",amrev4

DEVICE,DECOMPOSED=0
LOADCT,39
SET_PLOT, 'PS'
device, filename='amspec_res.eps', /COLOR,ENCAPSULATED=1
PLOT, xaxis(0:mmax1-1),amtot1(0:mmax1-1),XSTYLE=1,YSTYLE=1,  $
  XRANGE=[0.25,2.1],YRANGE=[-4.0,0.0],background=255,color=1,  $
  XTITLE="!5log m",YTITLE="!5 log <A!Im!N>",XTHICK=2,YTHICK=2, $
  XCHARSIZE=1.2,YCHARSIZE=1.6, psym=6,symsize=0.75
ERRPLOT, xaxis,amerrminus1,amerrplus1,color=1
OPLOT, xaxis(0:mmax2-1),amtot2(0:mmax2-1),psym=6,symsize=0.75,color=75
OPLOT, xaxis(0:mmax3-1),amtot3(0:mmax3-1),psym=6,symsize=0.75,color=150
OPLOT, xaxis(0:mmax4-1),amtot4(0:mmax4-1),psym=6,symsize=0.75,color=225
XYOUTS,1.7,-0.2,"!5LMAX=64",color=1,CHARSIZE=1.3
XYOUTS,1.7,-0.35,"!5LMAX=128",color=75,CHARSIZE=1.3
XYOUTS,1.7,-0.5,"!5LMAX=256",color=150,CHARSIZE=1.3
XYOUTS,1.7,-0.65,"!5LMAX=512",color=225,CHARSIZE=1.3
ERRPLOT, xaxis,amerrminus2,amerrplus2,color=75
ERRPLOT, xaxis,amerrminus3,amerrplus3,color=150
ERRPLOT, xaxis,amerrminus4,amerrplus4,color=225
DEVICE,/CLOSE
END



image1 = tvrd(true=3)

WRITE_JPEG,'amcomb.jpg',image1,QUALITY=100,true=3
print, ' am averaged over ', timearr1(0), ' to ',timearr1(count1-1),' using ',  $
  count1,' files.',FORMAT='(A,F5.2,A,F5.2,A,I4,A)'
END
