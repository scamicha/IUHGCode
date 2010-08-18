#include <stdio.h>
#include "idl_export.h"

void global_coef(int argc, void *argv[])
{
  extern void Globalcoef(); //Fortran routine

#define MAX_OUT_LEN 12

  IDL_LONG *JMAX, *KMAX, *LMAX, *JMIN;
  double *TORP, *DR, *CN, *TIME;
  int filelen = MAX_OUT_LEN;
  IDL_STRING *prefix;
  

  JMAX   = (IDL_LONG *)     argv[0];
  KMAX   = (IDL_LONG *)     argv[1];
  LMAX   = (IDL_LONG *)     argv[2];
  TORP   = (double *)       argv[3];
  DR     = (double *)       argv[4];
  TIME   = (double *)       argv[5];
  CN     = (double *)       argv[6];
  prefix = (IDL_STRING *)   argv[7];

  filelen = prefix->slen;

  globalcoef_(IDL_STRING_STR(prefix),TORP,JMAX,KMAX,LMAX,DR,TIME,CN,filelen);
 
  return;
}

void surfprog(int argc, void *argv[])
{
  extern void surfprgm();  //Fortran routine

  IDL_LONG *COUNT, *MODE, *AMPANG, *RMIN, *RMAX, *JTOP, *IK;
  double *B, *TIMESUB, *TORP, *CX, *NPX;

  //  if(argc != 6) return;

  COUNT     = (IDL_LONG *)     argv[0];
  JTOP      = (IDL_LONG *)     argv[1];
  MODE      = (IDL_LONG *)     argv[2];
  TORP      = (double *)       argv[3];
  AMPANG    = (IDL_LONG *)     argv[4];
  RMIN      = (IDL_LONG *)     argv[5];
  RMAX      = (IDL_LONG *)     argv[6];
  IK        = (IDL_LONG *)     argv[7];
  B         = (double *)       argv[8];
  TIMESUB   = (double *)       argv[9];  
  CX        = (double *)       argv[10];
  NPX       = (double *)       argv[11];


  printf("%s","returning");
  surfprgm_(COUNT,JTOP,MODE,TORP,AMPANG,RMIN,RMAX,IK,B,TIMESUB,CX,NPX);

  return;
}
	    
