#!/bin/bash

data_sets=(af_shell2.mtx cage14.mtx darcy003.mtx GL7d19.mtx kim2.mtx lp_osa_60_linear.mtx sls.mtx stormG2_1000.mtx wikipedia-20051105.mtx)
path=/root/spmv/data/test_matrices
#data_sets=(494_bus.mtx)
#path=../data
make clean
make csr-omp

if [ -f result.log ]
then
    rm result.log
fi
touch result.log
for thread_num in `seq 1 12`
do
    for i in `seq 0 $((${#data_sets[@]}-1))`
    do
        data_set=${data_sets[$i]}
        echo "thread num: $thread_num, test data set: $data_set" >> result.log
        ./csr_omp $path/$data_set $thread_num >> result.log
    done
done
