#!/bin/sh 

#set -x
echo
echo "starting"
echo 

n_atoms={{n_atoms}}
n_atoms_per_group=1

for i in `seq 1  $n_atoms`
do

cd calc_atom_$i

for i_dir in $(ls -l |grep ^d|awk '{print $9}')
do
cd $i_dir
test_finish=`grep 'Have a nice day'  $i_dir.out |awk '{print $1}'`

if [ "$test_finish" != "Have" ]; then
   echo '************************************************************************* '
   echo 'this work need to be resubmit'
   echo "atom$i $i_dir"

   job_num=`squeue | wc -l`                                                                          
          
   while [ $job_num -gt 1200 ]
   do
       job_num=`squeue | wc -l`
   done
cat > sub.sh <<EOF 
#!/bin/bash
#SBATCH -p thcp1
#SBATCH -N 1 
#SBATCH -n 64
#SBATCH -J 20+
#SBATCH -t 12:00:00
module load mpich/mpi-x-gcc9.3.0
UCX_TLS=tcp

srun  /thfs1/home/shanghh/gaoyingxiang/05_aims/01_single_node/bin/aims.191127.scalapack.mpi.x > $i_dir.out

EOF

rm -f control.in
cp ~/gaoyingxiang/07_work/generate/more_than_20/pre_control.in pre_control.in
sh ~/gaoyingxiang/07_work/generate/more_than_20/control.sh
dos2unix control.in
python3 ~/gaoyingxiang/07_work/generate/more_than_20/prop1_for_20+.py

chmod u+x sub.sh
sbatch sub.sh
   
echo ' '
else 
   echo "atom$i $i_dir finished"
   echo '====================================='
   echo ''
fi

cd ../
done 

cd ../
done
