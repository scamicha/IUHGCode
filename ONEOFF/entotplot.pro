PRO entotplot
  
  LOADCT, 43
  DEVICE, DECOMPOSED=0

  INFILE1  = 'Qplot_baseline.00545000'
  INFILE2  = 'Qplot_baseline.02050000'

  openr,lun1,INFILE1,/F77_UNFORMATTED,/GET_LUN
  openr,lun2,INFILE2,/F77_UNFORMATTED,/GET_LUN

  JMAX1  = 0L
  FILES1 = 0L
  COUNTER1 = 0L
  START1 = 0L
  STOP1  = 0L
  SKIP1  = 0L 

  JMAX2  = 0L
  FILES2 = 0L
  COUNTER2 = 0L
  START2 = 0L
  STOP2  = 0L
  SKIP2  = 0L 

  readu,lun1,JMAX1,FILES1,COUNTER1
  readu,lun1,START1,STOP1,SKIP1

  timearr1 = DBLARR(FILES1)
  rplot1   = DBLARR(JMAX1)
  qomega1  = DBLARR(JMAX1,FILES1)
  qkappa1  = DBLARR(JMAX1,FILES1)
  toteng1  = DBLARR(FILES1)
  columncool1 = DBLARR(JMAX1,FILES1)

  readu,lun1, timearr1
  readu,lun1, toteng1
  readu,lun1, rplot1
  readu,lun1, qomega1
  readu,lun1, qkappa1
  readu,lun1, columncool1


  readu,lun2,JMAX2,FILES2,COUNTER2
  readu,lun2,START2,STOP2,SKIP2

  timearr2 = DBLARR(FILES2)
  rplot2   = DBLARR(JMAX2)
  qomega2  = DBLARR(JMAX2,FILES2)
  qkappa2  = DBLARR(JMAX2,FILES2)
  toteng2  = DBLARR(FILES2)
  columncool2 = DBLARR(JMAX2,FILES2)

  readu,lun2, timearr2
  readu,lun2, toteng2
  readu,lun2, rplot2
  readu,lun2, qomega2
  readu,lun2, qkappa2
  readu,lun2, columncool2

  alltime = DBLARR(COUNTER1+COUNTER2)
  alleng  = DBLARR(COUNTER1+COUNTER2)

  print, FILES1,COUNTER1
  print, FILES2,COUNTER2

  print,toteng1(0),toteng1(FILES1-2)
  print,toteng2(0),toteng2(FILES2-2)

  FOR I=0L,COUNTER1-1 DO BEGIN
      alltime(I) = timearr1(I)
      alleng(I)  = toteng1(I)
  ENDFOR

  FOR I=0L,COUNTER2-1 DO BEGIN
      alltime(I+COUNTER1) = timearr2(I)
      alleng(I+COUNTER1) = toteng2(I)
  ENDFOR

  alleng = alleng/alleng(0)

  SET_PLOT,'PS'
  LOADCT,43
;  DEVICE,DECOMPOSED=0
;  WINDOW,0

  DEVICE,filename='baseline_toteng.eps',/COLOR,/ENCAPSULATED

  PLOT,alltime,alleng,color=1,backgroun=255,XRANGE=[0,21.5],YTITLE='Internal Energy',  $
            XTITLE='Time (ORPs)',XSTYLE=1,YSTYLE=1,YRANGE=[0.25,1.05],THICK=3

  DEVICE,/CLOSE

END
