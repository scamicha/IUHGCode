#!/bin/tcsh
#sets some environment variables then launches the job
#the KMP_STACKSIZE variable is specific to the Intel
#compiler. If you're using another compiler set 
#OMP_STACKSIZE instead. Also, you should have a 
#unlimited stack set in your .bashrc

setenv OMP_NUM_THREADS 8
setenv KMP_STACKSIZE 1G
./decompose > P0.5MPI.out
