PRO COOLRATE

JMAX      = 512L
KMAX      = 64L
LMAX      = 512L
AUJREQ    = 40.d0
STEPBEGIN = 1000000L
STEPEND   = 1005000L
STEPSKIP  = 5000L
TORP      = 1605.63d0
SWAP_ENDIAN = 1L


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

JREQ      = 0L


FOR step=STEPBEGIN,STEPEND,STEPSKIP DO BEGIN

    savedfile = string(savedprefix,step,FORMAT='(A,I8,8)')
    
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
      READU,lun2,dr,dz,delt,time,elost,den,sound,jreq,ommax
      READU,lun2,tmassini,tmass,tmassadd,tmassout,tmassacc
      FREE_LUN,lun2
