#!/bin/sh
#
#PBS -N teste
#PBS -l walltime=30:00
#PBS -l nodes=1:r662:ppn=48
#PBS -q mei

module load gcc/5.3.0
gcc main.c -o matrix

for size in 48 144 600 2048 
do
    ./matrix $size >> out.csv
    echo "\n" >> out.csv
done

