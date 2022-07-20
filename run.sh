#!/usr/bin/bash

echo "Please, write number of jobs:"
read num_jobs
echo "Please, write number of divisions:"
read divisions

for i in `seq $num_jobs`; do 
  nice ./exec $i $num_jobs $divisions &
  sleep 1 && echo $i &
done

wait
./matrices $num_jobs $divisions &
wait
./save-sizes &

wait

echo Done

