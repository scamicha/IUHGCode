PRO qploteos
  
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

  readu,lun1, timearr
  readu,lun1, rplot
  readu,lun1, qomega
  readu,lun1, qkappa

  PLOT,rplot(*),qomega(*,3),MAX_VALUE=100,MIN_VALUE=-1,XSTYLE=1
;       COLOR=colorcount,YSTYLE=1,XRANGE=[xmin,xmax],YRANGE=[ymin,ymax],$
;           LINESTYLE=linetype
  OPLOT,rplot(*),qkappa(*,3),MAX_VALUE=100,MIN_VALUE=-1;,COLOR=colorcount, $
;           LINESTYLE=2
 
END
