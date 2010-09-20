#!/bin/tcsh
#sets some environment variables then launches the job
#the KMP_STACKSIZE variable is specific to the Intel
#compiler. If you're using another compiler set 
#OMP_STACKSIZE instead. Also, you should have a 
#unlimited stack set in your .bashrc

setenv OMP_NUM_THREADS 4
setenv KMP_STACKSIZE 500M
./decompose_mpi > indirectMPI.out
