
set -x
cd output

for i_dir in $(ls -l |grep ^d|awk '{print $9}')
do
cd $i_dir
test_finish=`grep 'Have a nice day'  output |awk '{print $1}'`

if [ "$test_finish" != "Have" ]; then
   echo '************************************************************************* '
   echo 'this work need to be resubmit'
   echo " $i_dir"
   pwd

echo ' '
else 
   echo " $i_dir finished"
   echo '====================================='
   echo ''
fi

cd ../
done



