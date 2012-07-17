#!/bin/bash

#module add libs/cuda/4.0.17
#module add mpi/gcc/openmpi/1.4.4

start=0
end=1  
step=1
current=$start



numproc=$1
 

while [ "$current" -le $end ]
do


        mpirun -np $numproc bin/smaug -o gather $current
        current=`expr $current + $step`

done
