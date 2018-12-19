#!/bin/bash
if [ -f compile.log ]
then
    rm compile.log
fi
touch compile.log
for i in `seq 1 100`
do
    ./csc_omp ../data/494_bus.mtx >> compile.log
done
