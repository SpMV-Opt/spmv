#!/bin/bash
# -----------------------------------------------------------
# Description: Register Blocking CSR kernel auto generation.
# -----------------------------------------------------------

#set -x

KERNEL_SRC=kernel.cpp

file_gen() {
    if [ -f ${KERNEL_SRC} ]
    then
        rm -f ${KERNEL_SRC}
        touch ${KERNEL_SRC}
    fi
    echo "/*">>${KERNEL_SRC}
    echo " *   Register Blocking CSR implement.">>${KERNEL_SRC}
    echo " *">>${KERNEL_SRC}
    echo " */">>${KERNEL_SRC}
    echo >>${KERNEL_SRC}
    echo >>${KERNEL_SRC}
}

emit_func_header() {
    R=$1
    C=$2
    echo "void bcsr_${R}x${C}(const int &bm, const int *b_row_start, const int *b_col_idx," >> ${KERNEL_SRC}
    echo "const double *b_values, const double *x, double *y) {" >> ${KERNEL_SRC}
}

emit_func_body() {
    local R=${1}
    local C=${2}
    #echo ${R} ${C}
    echo "int i, j;" >>${KERNEL_SRC}
    echo "double" >>${KERNEL_SRC}

    local m=$(( ${R} - 1 ))
    # emit d0, d1, ..., d(r-1), x0, x1, ..., x(c-1)
    for i in `seq 0 $m`
    do
        echo "d${i}, ">>${KERNEL_SRC}
    done
    local t=$(( ${C} - 1 ))
    for i in `seq 0 $t`
    do
        if [ ${i} -eq ${t} ]
        then
            echo "x${i};">>${KERNEL_SRC}
        else
            echo "x${i}, ">>${KERNEL_SRC}
        fi
    done
    # emit outer loop
    echo "for (i = 0; i < bm; ++i) {">>${KERNEL_SRC}
    # init d0, d1, ..., d(r-1)
    for i in `seq 0 $m`
    do
        echo "d${i} = y[${R} * i + ${i}];">>${KERNEL_SRC}
    done
    # emit inner loop
    echo "for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += ${R} * ${C}) {">>${KERNEL_SRC}
    # init x0, x1, ..., x(c-1)
    for i in `seq 0 $t`
    do
        echo "x${i} = x[${C} * b_col_idx[j] + ${i}];">>${KERNEL_SRC}
    done
    # reduce d0, d1, ..., d(r-1)
    for i in `seq 0 $t`
    do
        for j in `seq 0 $m`
        do
            index=$(( $(( ${j} * ${C} )) + ${i} ))
            echo "d${j} += b_values[${index}] * x${i};">>${KERNEL_SRC}
        done
    done
    # write back d0, d1, ..., d(r-1)
    for i in `seq 0 $m`
    do
        echo "y[${R} * i + ${i}] = d${i};">>${KERNEL_SRC}
    done
    # end inner loop
    echo "}">>${KERNEL_SRC}
    # end outer loop
    echo "}">>${KERNEL_SRC}
}

emit_func_tail() {
    echo "}" >> ${KERNEL_SRC}
}

echo "Begin to generate kernel..."
file_gen
<<!
for _i in `seq 1 12`
do
    for _j in `seq 1 12`
    do
        echo ${_i} ${_j}
        emit_func_header ${_i} ${_j}
        emit_func_body ${_i} ${_j}
        emit_func_tail
        echo >>${KERNEL_SRC}
    done
done
!

emit_func_header 1 1
emit_func_body   1 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 2
emit_func_body   1 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 3
emit_func_body   1 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 4
emit_func_body   1 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 5
emit_func_body   1 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 6
emit_func_body   1 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 7
emit_func_body   1 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 8
emit_func_body   1 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 9
emit_func_body   1 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 10
emit_func_body   1 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 11
emit_func_body   1 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 1 12
emit_func_body   1 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 2 1
emit_func_body   2 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 2
emit_func_body   2 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 3
emit_func_body   2 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 4
emit_func_body   2 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 5
emit_func_body   2 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 6
emit_func_body   2 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 7
emit_func_body   2 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 8
emit_func_body   2 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 9
emit_func_body   2 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 10
emit_func_body   2 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 11
emit_func_body   2 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 2 12
emit_func_body   2 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 3 1
emit_func_body   3 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 2
emit_func_body   3 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 3
emit_func_body   3 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 4
emit_func_body   3 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 5
emit_func_body   3 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 6
emit_func_body   3 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 7
emit_func_body   3 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 8
emit_func_body   3 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 9
emit_func_body   3 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 10
emit_func_body   3 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 11
emit_func_body   3 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 3 12
emit_func_body   3 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 4 1
emit_func_body   4 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 2
emit_func_body   4 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 3
emit_func_body   4 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 4
emit_func_body   4 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 5
emit_func_body   4 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 6
emit_func_body   4 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 7
emit_func_body   4 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 8
emit_func_body   4 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 9
emit_func_body   4 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 10
emit_func_body   4 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 11
emit_func_body   4 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 4 12
emit_func_body   4 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 5 1
emit_func_body   5 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 2
emit_func_body   5 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 3
emit_func_body   5 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 4
emit_func_body   5 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 5
emit_func_body   5 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 6
emit_func_body   5 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 7
emit_func_body   5 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 8
emit_func_body   5 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 9
emit_func_body   5 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 10
emit_func_body   5 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 11
emit_func_body   5 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 5 12
emit_func_body   5 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 6 1
emit_func_body   6 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 2
emit_func_body   6 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 3
emit_func_body   6 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 4
emit_func_body   6 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 5
emit_func_body   6 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 6
emit_func_body   6 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 7
emit_func_body   6 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 8
emit_func_body   6 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 9
emit_func_body   6 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 10
emit_func_body   6 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 11
emit_func_body   6 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 6 12
emit_func_body   6 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 7 1
emit_func_body   7 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 2
emit_func_body   7 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 3
emit_func_body   7 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 4
emit_func_body   7 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 5
emit_func_body   7 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 6
emit_func_body   7 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 7
emit_func_body   7 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 8
emit_func_body   7 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 9
emit_func_body   7 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 10
emit_func_body   7 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 11
emit_func_body   7 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 7 12
emit_func_body   7 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 8 1
emit_func_body   8 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 2
emit_func_body   8 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 3
emit_func_body   8 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 4
emit_func_body   8 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 5
emit_func_body   8 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 6
emit_func_body   8 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 7
emit_func_body   8 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 8
emit_func_body   8 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 9
emit_func_body   8 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 10
emit_func_body   8 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 11
emit_func_body   8 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 8 12
emit_func_body   8 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 9 1
emit_func_body   9 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 2
emit_func_body   9 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 3
emit_func_body   9 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 4
emit_func_body   9 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 5
emit_func_body   9 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 6
emit_func_body   9 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 7
emit_func_body   9 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 8
emit_func_body   9 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 9
emit_func_body   9 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 10
emit_func_body   9 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 11
emit_func_body   9 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 9 12
emit_func_body   9 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 10 1
emit_func_body   10 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 2
emit_func_body   10 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 3
emit_func_body   10 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 4
emit_func_body   10 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 5
emit_func_body   10 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 6
emit_func_body   10 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 7
emit_func_body   10 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 8
emit_func_body   10 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 9
emit_func_body   10 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 10
emit_func_body   10 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 11
emit_func_body   10 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 10 12
emit_func_body   10 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 11 1
emit_func_body   11 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 2
emit_func_body   11 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 3
emit_func_body   11 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 4
emit_func_body   11 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 5
emit_func_body   11 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 6
emit_func_body   11 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 7
emit_func_body   11 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 8
emit_func_body   11 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 9
emit_func_body   11 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 10
emit_func_body   11 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 11
emit_func_body   11 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 11 12
emit_func_body   11 12
emit_func_tail
echo >>${KERNEL_SRC}

emit_func_header 12 1
emit_func_body   12 1
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 2
emit_func_body   12 2
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 3
emit_func_body   12 3
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 4
emit_func_body   12 4
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 5
emit_func_body   12 5
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 6
emit_func_body   12 6
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 7
emit_func_body   12 7
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 8
emit_func_body   12 8
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 9
emit_func_body   12 9
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 10
emit_func_body   12 10
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 11
emit_func_body   12 11
emit_func_tail
echo >>${KERNEL_SRC}
emit_func_header 12 12
emit_func_body   12 12
emit_func_tail
echo >>${KERNEL_SRC}

clang-format -i $KERNEL_SRC
echo "End kernel generating..."
