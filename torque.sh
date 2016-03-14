#!/bin/bash
#PBS -N 24Procesors
#PBS -l nodes=1:ppn=24
#PBS -l walltime=720:00:00
#PBS -q bigmem
#PBS -m bae
#PBS -M jvergar2@gmail.com

cd $PBS_O_WORKDIR
gmsh Mesher_Damian/Entrada/untitled2.geo -2
cd Entrada
./malla.out
cd ..
./damian.out
