#!/bin/bash
#
# Verify the correctness of basic implement
#

CSR=csr
CSC=csc
BCSR=bcsr
CSR_OMP=csr_omp
CSC_OMP=csc_omp
num_threads=12

LOG=verify.log
DATA_SET=../data/494_bus.mtx
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
    make bcsr VER=TRUE
    make clean
}
verify() {
    if [ -f $LOG ]
    then
        rm $LOG
    fi
    touch $LOG
    ./csr $DATA_SET $num_threads >>$LOG
    ./csr_omp $DATA_SET $num_threads >>$LOG
    ./csc $DATA_SET $num_threads >>$LOG
    ./csc_omp $DATA_SET $num_threads >>$LOG
    expect_pass=4
    pass_count=`grep -rin PASS $LOG | wc -l`

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
    rm $CSC $CSR $CSR_OMP $CSC_OMP $LOG
}
BCSR_LOG=bcsr_verify.log
bcsr_verify() {
    if [ -f $BCSR_LOG ]
    then
        rm $BCSR_LOG
    fi
    touch $BCSR_LOG
    for i in `seq 1 12`
    do
        for j in `seq 1 12`
        do
            ./bcsr $DATA_SET $i $j >>$BCSR_LOG
        done
    done
    expect_pass=144
    pass_count=`grep -rin PASS $BCSR_LOG | wc -l`

    echo "===--------------- bcsr statistics info ---------------==="
    echo "Expect passed: $expect_pass"
    echo "Actual passed: $pass_count"
    if [ $pass_count != $expect_pass ]
    then
        diff=`expr $expect_pass - $pass_count`
        echo "failed: $diff"
    else
        echo "PASS!"
    fi
    rm $BCSR $BCSR_LOG
}
build
verify
bcsr_verify
