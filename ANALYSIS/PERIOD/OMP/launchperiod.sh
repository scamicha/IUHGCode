#!/bin/tcsh
#sets some environment variables then launches the job
#the KMP_STACKSIZE variable is specific to the Intel
#compiler. If you're using another compiler set 
#OMP_STACKSIZE instead. Also, you should have a 
#unlimited stack set in your .bashrc
# compile main code with ifort -o period -fast -openmp periodogram.f90

setenv OMP_NUM_THREADS 8
setenv KMP_STACKSIZE 10M

set jmax=512
set kmax=64
set lmax=128
set mode=8
set start=060000
set finish=320000
set skip=200
set tstart=12.57
set tend=19.5
set aujreq=40.0
set outfile=P.5am15AU.dat
set rhodir=../RHOTEMP/
set savedir=../SAVED/

./period $jmax $kmax $lmax $mode $start $finish $skip $tstart $tend $aujreq $outfile $rhodir $savedir > period.out
