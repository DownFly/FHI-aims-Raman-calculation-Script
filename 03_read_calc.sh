#!/bin/sh 
#SBATCH -p thcp1
#SBATCH -N 1 
#SBATCH -n 64

module purge
module load mpich/mpi-n-gcc9.3.0
UCX_TLS=tcp


#set -x
n_atom_start=1
n_atom_end=`grep "atom" geometry.in | wc -l`

mkdir -p read_all_atoms

for i in `seq $n_atom_start  $n_atom_end`
do
cp calc_atom_$i/*/i*.out  ./read_all_atoms/
cp control.in  ./read_all_atoms/
cp geometry.in ./read_all_atoms/ 
done


#--------read--------
cd read_all_atoms/
/thfs1/home/shanghh/gaoyingxiang/05_aims/01_single_node/bin/02_read_aims_harmonic_raman.pl 0.0025 0.0025  polar  $n_atom_start  $n_atom_end   > output
cd ../

echo 'finished :) '
