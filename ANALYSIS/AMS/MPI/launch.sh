#!/bin/tcsh
#sets some environment variables then launches the job
#the KMP_STACKSIZE variable is specific to the Intel
#compiler. If you're using another compiler set 
#OMP_STACKSIZE instead. Also, you should have a 
#unlimited stack set in your .bashrc

setenv OMP_NUM_THREADS 8
setenv KMP_STACKSIZE 1G

set jmax=512
set jstart=90
set kmax=64
set lmax=128
set start=060000
set finish=320000
set skip=200
set outfile=P.5am15AU.dat
set indir=../RHOTEMP/

./amspec_mpi $jmax $kmax $lmax $start $finish $skip $outfile $jstart $indir > amspec.out
