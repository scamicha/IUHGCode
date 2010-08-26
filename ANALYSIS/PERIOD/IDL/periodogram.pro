PRO PERIODOGRAM


;!!!!!!! User defined parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  JMAX         =  512L
  KMAX         =  64L
  LMAX         =  128L
  JREQ         =  250L
  AUJREQ       =  40.d0
  STEPBEGIN    =  290000L
  STEPEND      =  445000L
  STEPSKIP     =  200L
  JTOP         =  500L
  MINT         =  0.d0
  MAXT         =  0.d0
  MINDT        =  0.01d0
  RMIN         =  3L
  RMAX         =  JTOP
  TORP         =  1605.63d0
  DR           =  0.2021807d0
  NPRIME       =  9.d0
  YMIN         =  0.4d0
  YMAX         =  0.8d0
  XMIN         =  0.d0
  XMAX         =  2.5d0
  ZMIN         =  -4.d0
  ZMAX         =  0.d0
  ANGLE        =  [30.d0,30.d0]
  COEFFILE     =  "coefs.tot"
  CONTOURLABEL =  1L
  SWAP_ENDIAN  =  1L
  PI           =  3.1415926535897932384626433832795028841971


;!!!!!! Following parameters only needed if rhoswitch = true!!!!!!!!!!!!

  RHOSWITCH = 0
  RHOPREFIX    = '../RHOTEMP/rho3d.'
  rhofile  = ''


;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;!!!!!! Following parameters if crswitch = true, this plots CR,ILR and OLR

  CRSWITCH = 1
  SAVEDPREFIX = '../SAVED/saved.'
  SAVEDFILE = ''

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
;!!!!!!!!!!!! If IFOFSWITCH = true, plot the inflow/outflow boundary as
;determined by the first and last files considered.

  IFOFSWITCH = 1

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  JMAX2     =  JMAX+2
  KMAX2     =  KMAX+2
  RNUM      =  0L
  MODE      =  0L
  AMPANG    =  0L
  COUNT     =  0L
  CONTIN    =  1L
  OUTPUT    =  0L
  IK        =  0L
  FILENUM   =  0L 
  ANALTYPE  =  0L
  BDYSWITCH =  0L
  RINGSWITCH=  0L
  NUMRINGS  =  0L
  TIME      =  0.d0
  DUMMY     =  0.d0
  DT        =  0.d0
  DTOT      =  0.d0
  TSTART    =  0.d0
  TEND      =  0.d0
  CONMAX    =  0.d0
  CONMIN    =  0.d0
  CONMID    =  0.d0
  RPLOT     =  0.d0
  ZRANGE    =  0.d0
  DUMMY2    =  0L
  MAXMODE   =  8L
  CHANGELEV =  1L
  JSIZE     =  0L
  ZOOMIN    =  0L
  JPLOT     =  0L
  SUMALL    =  0L
  LIGHTDIR  =  [0,0,1]
  NUMLEV    =  DBLARR(2)
  YRANGE    =  DBLARR(2)
  XRANGE    =  DBLARR(2)
  RHF       =  DBLARR(JMAX+2)
  R         =  DBLARR(JMAX+2)
  IF(CRSWITCH EQ 1) THEN BEGIN
      S         =  DBLARR(JMAX+2,KMAX+2,LMAX)
      T         =  DBLARR(JMAX+2,KMAX+2,LMAX)
      A         =  DBLARR(JMAX+2,KMAX+2,LMAX)
      RHO       =  DBLARR(JMAX+2,KMAX+2,LMAX)
      EPS       =  DBLARR(JMAX+2,KMAX+2,LMAX)
      OMEGA     =  DBLARR(JMAX+2,KMAX+2,LMAX)
      OMEGAVG   =  DBLARR(JMAX+2)
      KAPPA     =  DBLARR(JMAX+2)
      MASSUM    =  0.d0
  ENDIF
  IF(IFOFSWITCH EQ 1) THEN BEGIN
      RHOSTART  = DBLARR(JMAX+2,KMAX+2,LMAX)
      RHOEND    = DBLARR(JMAX+2,KMAX+2,LMAX)
      MCYLSTART = DBLARR(JMAX+2)
      MCYLEND   = DBLARR(JMAX+2)
      MDOT      = DBLARR(JMAX+2)
      IFOFBDY   = DBLARR(JMAX+2)
      STARTIME  = 0.d0
      ENDTIME   = 0.d0
  ENDIF
  CN        =  DBLARR(JMAX+2,(LMAX/2)+1,2)
  STORAGE   =  DBLARR((JMAX-1)*(LMAX))
  TINTERV   =  DBLARR(2)
  ZMINMAX   =  DBLARR(2)
  DUMMY1    =  ''
  TITLESTR  =  ''
  SUBTITLESTR= ''
  STR       =  ''
  OUTFILE   =  ''
  GROWTHFILE=  ''
  DEVICE,DECOMPOSED=0
  LOADCT,42

  RNUM          =  JMAX/10 +1
  ARRAYSIZE     =  (STEPEND-STEPBEGIN)/STEPSKIP+1
  IF(ARRAYSIZE EQ 0) THEN ARRAYSIZE = 1
  COEF          =  DBLARR(ARRAYSIZE,JMAX,LMAX/2,2)
  TIMEARR       =  DBLARR(ARRAYSIZE)
  TIMESUB       =  DBLARR(ARRAYSIZE)
  B             =  DBLARR(JTOP,ARRAYSIZE,MAXMODE,2)

  FOR j = 0,jmax+1 DO BEGIN
      rhf(j) = (DOUBLE(j)-0.5d0)*dr
      r(j)   = rhf(j)-0.5d0*dr
  ENDFOR


  IF (RHOSWITCH NE 1) THEN OPENR,lun1,COEFFILE,/GET_LUN
  FOR I = 0L,ARRAYSIZE-1L DO BEGIN
      IF (RHOSWITCH EQ 1) THEN BEGIN

          filenum    = (I*STEPSKIP)+STEPBEGIN
          IF (filenum GE 1000000) THEN filenum = filenum-999999
          rhofile   = STRING(rhoprefix,filenum,FORMAT='(A,I6.6)')
          dummy      = call_external('global_coefs.so','global_coef',JMAX,  $
                                KMAX,LMAX,TORP,DR,TIME,CN,rhofile)
          FOR L=0,LMAX/2-1 DO BEGIN
              FOR M=0,1 DO BEGIN
                  FOR J=1,JMAX-1 DO BEGIN
                      COEF(I,J,L,M) = CN(J,L+1,M)
                  ENDFOR
              ENDFOR
          ENDFOR
          TIMEARR(I) = TIME

      ENDIF ELSE BEGIN
          READF,lun1,TIME,format='(25X,1e11.4)'
          print, "TIME = ",TIME
          READF,lun1,STORAGE
          print, "GOT STORAGE"
          READF,lun1,DUMMY1
          print,"READ DUMMY1"
          FOR J=1,RNUM DO BEGIN
              READF,lun1,DUMMY
          ENDFOR
          print, "READ DUMMY"
          N = 0L
          FOR J=1,JMAX-1 DO BEGIN
              FOR M=0,1L DO BEGIN
                  FOR L=0,LMAX/2-1 DO BEGIN
                      COEF(I,J,L,M) = STORAGE(N)
                      N++
                  ENDFOR
              ENDFOR
          ENDFOR
          TIMEARR(I) = TIME
      ENDELSE 

  ENDFOR
  
  IF(CRSWITCH EQ 1) THEN BEGIN

      savedfile   = STRING(savedprefix,stepend,FORMAT='(A,I6.6)')       
      
      IF(swap_endian eq 1) THEN BEGIN
          OPENR,lun2,savedfile,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
      ENDIF ELSE BEGIN
          OPENR,lun2,savedfile,/F77_UNFORMATTED,/GET_LUN
      ENDELSE

      READU,lun2,S
      READU,lun2,T
      READU,lun2,A
      READU,lun2,RHO
      READU,lun2,EPS
      FREE_LUN,lun2
      
      FOR j= 0,jmax+1 DO BEGIN
          FOR k = 0,kmax+1 DO BEGIN
              FOR l = 0,lmax-1 DO BEGIN
                  omega(j,k,l)=a(j,k,l)/(rhf(j)^2)
              ENDFOR
          ENDFOR
      ENDFOR

      FOR j = 1,jmax DO BEGIN
          omegavg(j) = 0.d0
          massum = 0.d0
          FOR k = 0,kmax+1 DO BEGIN
              FOR l = 0,lmax-1 DO BEGIN
                  omegavg(j) = omegavg(j) + omega(j,k,l)
                  massum = massum+rho(j,k,l)
              ENDFOR
          ENDFOR
          omegavg(j) = omegavg(j)/massum
      ENDFOR
      FOR j = 1,jmax DO BEGIN

          kappa(j) = abs(rhf(j)*omegavg(j)*(omegavg(j+1)- $
                      omegavg(j-1))/dr + 4.d0*omegavg(j)^2)
          
      ENDFOR

      kappa = kappa^0.5
      kappa= kappa*torp/(2.d0*PI)
      omegavg = omegavg*torp/(2.d0*PI)
      rhf = rhf*AUJREQ/rhf(JREQ)

;      print, omegavg+(kappa/2.d0)
;      print, omegavg-(kappa/2.d0)

  ENDIF 

  IF (IFOFSWITCH EQ 1) THEN BEGIN

      rhofile   = STRING(rhoprefix,stepbegin,FORMAT='(A,I6.6)') 

      IF(swap_endian eq 1) THEN BEGIN
          OPENR,lun2,rhofile,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
      ENDIF ELSE BEGIN
          OPENR,lun2,rhofile,/F77_UNFORMATTED,/GET_LUN
      ENDELSE

      READU,lun2,RHOSTART
      READU,lun2,STARTIME
      FREE_LUN,lun2

      rhofile   = STRING(rhoprefix,stepend,FORMAT='(A,I6.6)') 

      IF(swap_endian eq 1) THEN BEGIN
          OPENR,lun2,rhofile,/F77_UNFORMATTED,/SWAP_ENDIAN,/GET_LUN
      ENDIF ELSE BEGIN
          OPENR,lun2,rhofile,/F77_UNFORMATTED,/GET_LUN
      ENDELSE

      READU,lun2,RHOEND
      READU,lun2,ENDTIME
      FREE_LUN,lun2    

      FOR J=1,JMAX DO BEGIN
          FOR K=1,KMAX DO BEGIN
              FOR L=0,LMAX-1 DO BEGIN
                  MCYLSTART(J) = MCYLSTART(J) + RHOSTART(J,K,L)
                  MCYLEND(J)   = MCYLEND(J)   + RHOEND(J,K,L)
              ENDFOR
          ENDFOR
          MCYLSTART(J) = MCYLSTART(J)*dr*pi*(r(j+1)^2-r(j)^2)/LMAX
          MCYLEND(J) = MCYLEND(J)*dr*pi*(r(j+1)^2-r(j)^2)/LMAX
          IF (J GT 1) THEN BEGIN
              MCYLSTART(J) = MCYLSTART(J) + MCYLSTART(J-1)
              MCYLEND(J)   = MCYLEND(J)   + MCYLEND(J-1)
          ENDIF
          MDOT(J) = MCYLEND(J) - MCYLSTART(J)
      ENDFOR

      COUNT = 0L
      FOR J=2,JMAX DO BEGIN
          IF (((MDOT(J) GT 0.d0) AND (MDOT(J-1) LE 0.d0)) OR ((MDOT(J)  $
               LT 0.d0) AND (MDOT(J-1) GE 0.d0))) THEN BEGIN
              IFOFBDY(COUNT) = RHF(J)
              COUNT++
          ENDIF
      ENDFOR

      IFOFBDY = IFOFBDY(0:COUNT-1)
      STARTIME = STARTIME/torp
      ENDTIME = ENDTIME/torp

      print, IFOFBDY
      print, STARTIME,ENDTIME
      
      titlestr  =  STRING('Mdot, time=',startime,'--',endtime,' ORP',FORMAT='(A11,F5.2,A2,F5.2,A4)')
      str         = STRING(IFOFBDY, FORMAT='(F6.2)')
      subtitlestr = STRJOIN(str,' ')
      
      WINDOW,2
      PLOT,RHF(1:JMAX),MDOT(1:JMAX),XRANGE=[0,60],XSTYLE=1,TITLE=titlestr, $
        SUBTITLE= subtitlestr,COLOR=0,BACKGROUND=255

      READ,OUTPUT,PROMPT='Press 1 to save this image. '
      IF (OUTPUT EQ 1) THEN BEGIN
          READ,OUTFILE,PROMPT='Please enter the file name. '
          image=tvrd(true=3)
          write_jpeg,OUTFILE,image,QUALITY=100,TRUE=3
      ENDIF

      READ,OUTPUT,PROMPT='Press 1 to save this data. '
      IF (OUTPUT EQ 1) THEN BEGIN
          READ,OUTFILE,PROMPT='Please enter the file name. '
          openw,lun1,OUTFILE,/GET_LUN,/F77_UNFORMATTED
          writeu,STARTIME,ENDTIME,COUNT
          writeu,MDOT
          writeu,IFOFBDY
          free_lun,lun1
      ENDIF          
  ENDIF    
      
      
  WHILE (CONTIN EQ 1) DO BEGIN
      
      READ,MINT,PROMPT="Enter minimum time in ORPs for period search: "
      READ,MAXT,PROMPT="Enter maximum time in ORPs for period search: "
      MAXT          =  MAXT*TORP
      MINT          =  MINT*TORP

      COUNT     = 0
      CHANGELEV = 1

      FOR I=0L,ARRAYSIZE-1L DO BEGIN
          IF (I NE 0) THEN DT = TIMEARR(I)-TIMEARR(I-1)
          IF (I EQ 0) THEN DT = TIMEARR(I)
          IF ((TIMEARR(I) GT MINT) AND (TIMEARR(I) LT MAXT) AND $
              (DT GT MINDT)) THEN BEGIN
              COUNT++
              FOR J=0,JTOP-1 DO BEGIN
                  FOR M=0,MAXMODE-1 DO BEGIN
                      FOR N=0,1 DO BEGIN
                          K = J + 1
                          B(J,COUNT-1,M,N) = COEF(I,K,M,N)
                      ENDFOR
                  ENDFOR
              ENDFOR
              TIMESUB(COUNT-1) = TIMEARR(I)       
          ENDIF
      ENDFOR

;      print,"ARRAYSIZE = ",ARRAYSIZE
      print,"COUNT = ", COUNT

      INTERACT      =  0L
      WIREFRAME     =  0L
      TSTART        =  TIMESUB(0)
      TEND          =  TIMESUB(COUNT-1)
      CX            =  DBLARR(2*COUNT,JTOP)
      MODESUM       =  DBLARR(2*COUNT)
      NPX           =  DBLARR(2*COUNT)
      DTOT          = (TIMESUB(COUNT-1) - TIMESUB(0))/COUNT
      BSUB          =  DBLARR(JTOP,COUNT)
      TAMP          =  DBLARR(COUNT)

      PRINT,'Type 1 for surface of amplitudes:'
      READ,ANALTYPE,PROMPT='Type 2 for periodogram:  '
      IF ((ANALTYPE EQ 1)) THEN AMPANG = 0L
      IF ((ANALTYPE EQ 2) ) THEN AMPANG = 1L
      IF ((ANALTYPE NE 1) AND (ANALTYPE NE 2)) THEN BEGIN
          PRINT,"Invalid choice, try again"
          PRINT,'Type 1 for surface of amplitudes:'
          READ,ANALTYPE,PROMPT='Type 2 for periodogram:  '
      ENDIF

      READ,MODE,PROMPT='Which mode would you like? '      

      LIGHT         =  1L
      YAXIS   = DBLARR(JTOP)
      BSUB(*,0:COUNT-1) = B(*,0:COUNT-1,MODE-1,AMPANG)

      IF (ANALTYPE EQ 1) THEN BEGIN
          WHILE (LIGHT EQ 1) DO BEGIN
              READ,ymin,PROMPT='ENTER MIN J: '
              READ,ymax,PROMPT='ENTER MAX J: '
              JSIZE         =  ymax-ymin+1
              BAMP          =  DBLARR(COUNT,JSIZE)
              YAMP          =  DBLARR(JSIZE)
              SHADES        =  INTARR(COUNT,JSIZE)

              READ,ZMINMAX,PROMPT='Please enter Amplitude range (zmin,zmax): '
              ZMIN = ZMINMAX(0)
              ZMAX = ZMINMAX(1)
              ZRANGE = ZMAX-ZMIN
              ZINC   = ZRANGE/255.

              lev      = INDGEN(256,/FLOAT)
              lev      = -lev + 255.
              lev      = lev/255.*(ZMINMAX(0)-ZMINMAX(1))+ZMINMAX(1)
              
              LOADCT, 42

              FOR J=0,JSIZE-1 DO BEGIN
                  FOR I=0,COUNT-1 DO BEGIN
                      IF (BSUB(J,I) LE 1.0e-10) THEN BEGIN
                          BAMP(I,J) = -10.0
                      ENDIF ELSE BEGIN
                          BAMP(I,J) = ALOG10(BSUB(J,I))
                          SHADES(I,J) = (BAMP(I,J) - ZMIN)/ZINC
                      ENDELSE
                      TAMP(I)   = TIMESUB(I)/TORP
                  ENDFOR
                  YAMP(J) = (J+2)*DOUBLE(AUJREQ)/DOUBLE(JREQ)
              ENDFOR
              ymin = ymin*DOUBLE(AUJREQ)/DOUBLE(JREQ)
              ymax = ymax*DOUBLE(AUJREQ)/DOUBLE(JREQ)
              temp     = INDGEN(350,/FLOAT)/350.*(ZMINMAX(0)-ZMINMAX(1))   $
                +ZMINMAX(1)
              colorbar = FLTARR(20,350)

              FOR i = 0, 19 DO colorbar(i,*) = temp(*)
              temp     = INDGEN(20)


              READ,ANGLE,PROMPT='Please enter rotation angles AX,AZ (30,30 is standard): '
              window,0,xsize=640,ysize=640
              DEVICE,DECOMPOSED=0

              SHADE_SURF,BAMP,TAMP,YAMP,XCHARSIZE=1.75,                     $
                YCHARSIZE=1.75,ZCHARSIZE=2.0,YRANGE=[ymin,ymax],            $
                ZRANGE=[zmin,zmax],AX=ANGLE(0),AZ=ANGLE(1),                 $
                XSTYLE=9,YSTYLE=9,ZSTYLE=9,SHADES=SHADES,                   $
                POSITION=[35,35,570 ,630],/DEVICE

              CONTOUR,colorbar,temp,colorbar(0,*),XSTYLE=5,YSTYLE=1,$
                LEVELS=lev,/FILL,/DEVICE,                          $
                POSITION = [591,40,611,600],                                $
                /NOERASE

              READ,LIGHT,PROMPT='Type 1 to reposition, any other to continue: '
          ENDWHILE
          READ,ZOOMIN,PROMPT='Care to compute growth rates? (1 for yes): '
          WHILE (ZOOMIN EQ 1) DO BEGIN
              READ,TINTERV,PROMPT='Please enter the time interval in ORPS for fitting (tstart,tend): '
              READ,GROWTHFILE,PROMPT='Enter filename to save growth rates: '
              OPENW,lun1,GROWTHFILE,/GET_LUN
              PRINTF,lun1,"Radius","Time Interval","Slope","Intercept","R-Value","RMS",   $
                FORMAT='(4X,A,2X,A,4X,A,2X,A,2X,a)'
              I = 0L
              STARTINDX = 0L
              ENDINDX   = 0L
              IF((TINTERV(0) LT TSTART/TORP) OR (TINTERV(1) GT TEND/TORP))   $
                THEN BEGIN
                  PRINT,"Sorry I didn't read in those values"
                  BREAK
              ENDIF
              WHILE (TAMP(I) LT TINTERV(0)) DO BEGIN
                  I++
              ENDWHILE
              STARTINDX = I
              WHILE (TINTERV(1) GT TAMP(I)) DO BEGIN
                  I++ 
              ENDWHILE
              ENDINDX = I
              IF (ENDINDX LE STARTINDX) THEN BEGIN
                  print, 'Invalid time interval'
                  BREAK
              ENDIF
              PRINT,STARTINDX,ENDINDX
              LENGTH   = ENDINDX-STARTINDX+1
              AFIT     = DBLARR(JTOP,LENGTH)
              TFIT     = DBLARR(LENGTH)
              X2       = DBLARR(LENGTH)
              Y2       = DBLARR(JTOP,LENGTH)
              XY       = DBLARR(JTOP,LENGTH)
              Sy       = DBLARR(JTOP)
              Sxy      = DBLARR(JTOP)              
              Sy2      = DBLARR(JTOP)              
              FITLIN   = 0.d0
              Sx       = 0.d0
              Sx2      = 0.d0
              slope    = 0.d0
              rcoef    = 0.d0
              rms      = 0.d0
              inter    = 0.d0
              minslope = -1.0e10
              maxslope = 1.0e10
              K = 0L 
              FOR I=STARTINDX, ENDINDX DO BEGIN
                  FOR J=0,JTOP-1 DO BEGIN
                      AFIT(J,K) = ALOG10(BSUB(J,I))
                      TFIT(K) = TAMP(I)
                  ENDFOR
                  K++
              ENDFOR


              X2    = TFIT^2   
              Y2    = AFIT^2
              XY    = REBIN(REFORM(TFIT,1,LENGTH),JTOP,LENGTH)*AFIT
              Sy    = TOTAL(AFIT,2)
              Sy2   = TOTAL(Y2,2)
              Sxy   = TOTAL(XY,2)
              Sx    = TOTAL(TFIT)
              Sx2   = TOTAL(X2)
              
              FOR J=0L,JTOP-1 DO BEGIN
                  slope = (LENGTH*Sxy(J)-Sx*Sy(J))/(LENGTH*Sx2-Sx^2)

                  IF (slope GT maxslope) THEN maxslope = slope
                  IF (slope LT minslope) THEN minslope = slope
                  inter = (Sy(J)-slope*Sx)/LENGTH
                  rcoef = (LENGTH*Sxy(J)-Sx*Sy(J))/(SQRT((LENGTH*Sx2-Sx^2)*  $
                            (LENGTH*Sy2(J)-Sy(J)^2)))
                  FOR I=0L,LENGTH-1 DO BEGIN
                      FITLIN = inter + slope*TFIT(I)
                      rms = rms + (AFIT(J,I)-FITLIN)^2
                  ENDFOR
                  rms = SQRT(rms/(LENGTH-1))
                  RSPOT = (J+2)*DOUBLE(AUJREQ)/DOUBLE(JREQ)
             
                  PRINT, RSPOT, TFIT(0),TFIT(LENGTH-1),slope,inter,     $
                    rcoef,rms,FORMAT='(2X,F5.2,2X,F5.2,"--",F5.2,2X,F8.4,2X,F9.4,2X,F9.4,2X,F8.4)'
                  PRINTF,lun1, RSPOT, TFIT(0),TFIT(LENGTH-1),slope,     $
                    inter,rcoef,rms,                                         $
                    FORMAT='(2X,F5.2,2X,F5.2,"--",F5.2,2X,F8.4,2X,F9.4,2X,F9.4,2X,F8.4)'
              ENDFOR
              READ,RPLOT,PROMPT='Enter particular radius to plot 0 to continue: '

              IF (RPLOT NE 0.0) THEN BEGIN
                  WHILE (RPLOT NE 0.0) DO BEGIN
                      J = 0L
                      WHILE (((J+2)*DOUBLE(AUJREQ)/DOUBLE(JREQ)) LT RPLOT) DO J++
                      JPLOT = J
                      WINDOW,1
                      PLOT,TFIT,AFIT(JPLOT,*)
                      READ,RPLOT,PROMPT='Enter particular radius to plot 0 to continue: '                      
                  ENDWHILE
              ENDIF
              READ,ZOOMIN,PROMPT='Type 1 for another time interval: '

          ENDWHILE
              
      ENDIF 
          

      IF (ANALTYPE EQ 2) THEN BEGIN

          TAMP(0:COUNT-1) = TIMESUB(0:COUNT-1)

          dummy2 = call_external('global_coefs.so','surfprog',COUNT,JTOP,     $
                                 MODE,TORP,AMPANG+1,RMIN,RMAX,IK,BSUB,        $
                                 TAMP,CX,NPX)
          XAXIS   = DBLARR(2*COUNT)
 
          FOR J=0,JTOP-1 DO BEGIN
              YAXIS(J) = (J+2)*DOUBLE(AUJREQ)/DOUBLE(JREQ)
          ENDFOR
      
          XAXIS   = NPX
          titlestr  =  STRING('Periodogram Contour m=',mode,'time=',tstart/torp,'--',tend/torp,' ORP',FORMAT='(A22,I2,A8,F5.2,A2,F5.2,A4)')
          CONMAX  = MAX(CX,MIN=CONMIN)
          CONMID  = MEDIAN(CX,/EVEN)
          IF (CONMIN LE 1.0e-5) THEN CONMIN = 1.0e-5
              
;              REDVEC   = [0,255,51,0,204,51,0,255,255,255]
;              GREENVEC = [0,255,153,51,0,0,255,255,102,51]
;              BLUEVEC  = [0,255,255,255,255,51,0,0,0,0]
 ;             TVLCT, REDVEC,GREENVEC,BLUEVEC

          WHILE (CHANGELEV EQ 1) DO BEGIN

              PRINT,"The minimum value of the contour is: ", CONMIN
              PRINT,"The maximum value of the contour is: ", CONMAX
              PRINT,"The median value of the contour is: ",  CONMID

              READ,NUMLEV,                                              $
                PROMPT="Please enter contour limits (min,max) "

              LEVEL     =  INDGEN(256,/FLOAT)
              LEVEL     =  LEVEL/255.*(NUMLEV(1)-NUMLEV(0))+NUMLEV(0)
              LOADCT, 42
              

              READ,XRANGE,PROMPT="Please enter Pattern range to plot: "
              READ,YRANGE,PROMPT="Please enter radius range to plot: "
                  
              DEVICE, DECOMPOSED=0
              window,0,xsize=640,ysize=640
          
              CONTOUR,CX,XAXIS,YAXIS,YRANGE=YRANGE,LEVELS=LEVEL,    $
                YSTYLE=1,TITLE=titlestr,color=0,background=255,     $
                XTITLE='Pattern Speed (1/ORP)',YTITLE='R(AU)',      $
                XRANGE=XRANGE,XSTYLE=1,/FILL
                      
              IF(CRSWITCH EQ 1) THEN BEGIN
                  OPLOT,omegavg,rhf,color=0
                  OPLOT,omegavg+(kappa/DOUBLE(mode)),rhf,color=0
                  OPLOT,omegavg-(kappa/DOUBLE(mode)),rhf,color=0
              ENDIF

              READ,CHANGELEV,PROMPT='Press 1 to change contour levels '
              
          ENDWHILE

          READ,OUTPUT,PROMPT='Press 1 to save this image. '
          IF (OUTPUT EQ 1) THEN BEGIN
              READ,OUTFILE,PROMPT='Please enter the file name. '
              image=tvrd(true=3)
              write_jpeg,OUTFILE,image,QUALITY=100,TRUE=3
          ENDIF
          
          READ,OUTPUT,PROMPT='Press 1 to save this data. '
          IF (OUTPUT EQ 1) THEN BEGIN
              READ,OUTFILE,PROMPT='Please enter the file name. '
              openw,lun1,OUTFILE,/GET_LUN,/F77_UNFORMATTED
              writeu,2*COUNT,JTOP
              writeu,CX
              writeu,XAXIS,YAXIS
              writeu,LEVEL
              IF(CRSWITCH EQ 1) THEN BEGIN
                  writeu,omegavg,kappa,rhf
              ENDIF
              free_lun,lun1
          ENDIF
      ENDIF

      READ,SUMALL,PROMPT= "Care to sum up radii? Press 1, any other to continue. "
      IF (SUMALL EQ 1) THEN BEGIN
          FOR I=0,2*COUNT-1 DO BEGIN
              MODESUM(I) = 0.d0
              FOR J=0,JTOP-1 DO BEGIN
                  IF(FINITE(CX(I,J))) THEN                        $
                    MODESUM(I) = MODESUM(I)+CX(I,J)
              ENDFOR
          ENDFOR

          WINDOW,1
          PLOT,XAXIS,MODESUM,BACKGROUND=255,COLOR=0,XRANGE=XRANGE
          
          READ,OUTPUT,PROMPT='Press 1 to save this image. '
          IF (OUTPUT EQ 1) THEN BEGIN
              READ,OUTFILE,PROMPT='Please enter the file name. '
              image=tvrd(true=3)
              write_jpeg,OUTFILE,image,QUALITY=100,TRUE=3
          ENDIF

      ENDIF
          
      READ,CONTIN,                                                      $
        PROMPT='To plot another mode press 1, any other will exit '
      
  ENDWHILE
      
END
      
