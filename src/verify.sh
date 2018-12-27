#!/bin/bash
#
# Verify the correctness of basic implement
#

CSR=csr
CSC=csc
CSR_OMP=csr_omp
CSC_OMP=csc_omp
num_threads=12

build() {
    #rm $CSR $CSR_OMP $CSC $CSC_OMP
    make clean
    make csr VER=TRUE
    make clean
    make csr-omp VER=TRUE
    make clean
    make csc VER=TRUE
    make clean
    make csc-omp VER=TRUE
    make clean
}
verify() {
    if [ -f verify.log ]
    then
        rm verify.log
    fi
    touch verify.log
    ./csr ../data/494_bus.mtx $num_threads >> verify.log
    ./csr_omp ../data/494_bus.mtx $num_threads >> verify.log
    ./csc ../data/494_bus.mtx $num_threads >> verify.log
    ./csc_omp ../data/494_bus.mtx $num_threads >> verify.log
    expect_pass=4
    pass_count=`grep -rin PASS verify.log | wc -l`

    echo "===--------------- statistics info ---------------==="
    echo "Expect passed: $expect_pass"
    echo "Actual passed: $pass_count"
    if [ $pass_count != $expect_pass ]
    then
        diff=`expr $expect_pass - $pass_count`
        echo "failed: $diff"
    else
        echo "PASS!"
    fi
    #rm $CSC $CSR $CSR_OMP $CSC_OMP verify.log
}
build
verify
