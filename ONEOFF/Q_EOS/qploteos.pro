PRO qploteos
  
  LOADCT, 43
  DEVICE, DECOMPOSED=0

  timestart = 10.d0
  timeend   = 20.d0
  INFILE   = ''

  openr,lun1,INFILE,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN

  readu,lun1,JMAX,FILES,COUNT
  readu,lun1,START,STOP,SKIP

  timearr = DBLARR(FILES)
  rplot   = DBLARR(JMAX)
  qomega  = DBLARR(JMAX,FILES)
  qkappa  = DBLARR(JMAX,FILES)
  toteng  = DBLARR(FILES)
  columncool = DBLARR(JMAX,FILES)

  readu,lun1, timearr
  readu,lun1, toteng
  readu,lun1, rplot
  readu,lun1, qomega
  readu,lun1, qkappa
  readu,lun1, columncool

  timesub = WHERE(timearr gt 13.d0)
  columnsub = columncool(*,timesub)

  WINDOW,0

  PLOT,rplot(*),qomega(*,3),MAX_VALUE=100,MIN_VALUE=-1,XSTYLE=1,color=1,$
       background=255,YSTYLE=1
;       COLOR=colorcount,YSTYLE=1,XRANGE=[xmin,xmax],YRANGE=[ymin,ymax],$
;           LINESTYLE=linetype
  OPLOT,rplot(*),qkappa(*,3),MAX_VALUE=100,MIN_VALUE=-1,color=1;,COLOR=colorcount, $
;           LINESTYLE=2

  WINDOW,1

  PLOT,rplot(*),columnsub(*,0),XSTYLE=1,color=1,background=255,YSTYLE=1

  WINDOW,2

  PLOT,timesub(*),toteng(*),XSTYLE=1,color=1,background=255,YSTYLE=1
 
END
