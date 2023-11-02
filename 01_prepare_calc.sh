#!/bin/sh 
#SBATCH -p thcp1
#SBATCH -N 1 
#SBATCH -n 64

module purge
module load mpich/mpi-n-gcc9.3.0
UCX_TLS=tcp


echo
echo "starting"
echo 

#n_atoms={{n_atoms}} 
n_atoms=`grep atom geometry.in |wc -l` 
n_atoms_per_group=1

#-------i is the atom index--------
for i in `seq 1  $n_atoms`

do

n_atom_start=`echo  $n_atoms_per_group'*('$i'-1)' |bc -l`
n_atom_end=`echo  $n_atoms_per_group'*('$i'-1)+'$n_atoms_per_group'-1'  |bc -l`

echo $n_atom_start ,  $n_atom_end
 
mkdir calc_atom_$i
cd calc_atom_$i 
cp ../control.in  ./
cp ../geometry.in ./

/thfs1/home/shanghh/gaoyingxiang/05_aims/01_single_node/bin/01_calc_aims_harmonic_raman.pl 0.0025 0.0025  polar  $n_atom_start  $n_atom_end   > output

cd ../

done
