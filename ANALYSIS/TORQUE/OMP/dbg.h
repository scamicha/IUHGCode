
      integer, parameter   ::  DBG_MAIN     =  1
      integer, parameter   ::  DBG_SOURCE   =  2
      integer, parameter   ::  DBG_FLUX     =  3
      integer, parameter   ::  DBG_VELOCITY =  4
      integer, parameter   ::  DBG_VLIMIT   =  5
      integer, parameter   ::  DBG_BDYGEN   =  6
      integer, parameter   ::  DBG_POT3     =  7
      integer, parameter   ::  DBG_ZAXPHI   =  8
      integer, parameter   ::  DBG_DELTA    =  9

      integer, parameter   ::  MAXDBGL = 50

      integer    dbgstep
      integer    dbgl
      integer    dbgf

      common /dbgblk/
     &   dbgstep,            ! time step
     &   dbgl,               ! call stack level
     &   dbgf(0:MAXDBGL-1)   ! subroutine call stack

c$omp threadprivate(/dbgblk/)
