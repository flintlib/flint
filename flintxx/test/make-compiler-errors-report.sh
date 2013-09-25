#!/bin/bash

CC=$1
TESTS=$2
FLAGS="-Wall -pedantic"
TMPF="/tmp/make-report.tmp"

# make sure we are in the test directory
cd $(dirname $0)

if [ "$CC" = "" ] ; then
    CC=g++
fi

if [ "$TESTS" = "" ] ; then
    TESTS="TEST_FMPZXX_INIT_WRONG TEST_FMPZXX_INIT_2 TEST_FMPZXX_ASSIGN_WRONG TEST_FMPZXX_CONVERT_WRONG TEST_FMPZXX_REF_INIT_WRONG_1 TEST_FMPZXX_REF_INIT_WRONG_2 TEST_FMPZXX_SRCREF_ASSIGN TEST_FMPZXX_ARITH_WRONG TEST_FMPZXX_ARITH_WRONG_DEEP TEST_FMPZXX_ARITHFUNC_WRONG_NARGS TEST_FMPZXX_ARITHFUNC_WRONG_TYPE TEST_FMPZXX_ARITHFUNC_WRONG_TYPE2 TEST_PADICXX_FORGET_EVAL"
fi

INCS=$(make -C ../../ print-INCS | grep 'INCS=' | sed 's/^INCS=//')

for test in $TESTS ; do
    echo $test
    $CC -E $INCS $FLAGS t-compiler-errors.cc -D$test -DEXTRACTING_SAMPLE \
        | grep -v '^#' | sed 's/^    //'
    $CC -c -o /dev/null $INCS $FLAGS t-compiler-errors.cc -D$test 2>$TMPF
    echo "Compiler error output is $(cat $TMPF | wc -l) line(s), $(cat $TMPF | wc -c) characters"
    echo "------START COMPILER ERROR OUTPUT-------"
    cat $TMPF
    echo "------END COMPILER ERROR OUTPUT---------"
    echo
    echo
done
