#!/bin/tcsh
set jmax=512
set jstart=90
set kmax=64
set lmax=128
set start=060000
set finish=320000
set skip=200
set outfile=P.5am15AU.dat
if (! -x amspec) then
    ifort -O3 -pad -align all -nocheck -shared-intel -openmp -r8 -mtune=core2 -mcmodel=medium -convert big_endian -o amspec amspec.f
fi
./amspec $jmax $kmax $lmax $start $finish $skip $outfile $jstart
