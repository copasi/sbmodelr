#!/bin/bash

# Caution that this test fails too often due to failure in the statistical tests
# it may fail at a rate of 1/15 times...

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py  -t Y -k 0.1 --pn X 0.1 norm --pn Y 0.1 norm  -o Wolf6x6Y.cps ../sources/Selkov-Wolf-Heinrich.cps 9 9 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# test -o option by checking if output file was created
if ! [ -f Wolf6x6Y.cps ]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py Wolf6x6Y.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# get values of X_i_j
grep -Po "X_\d+,\d+\s+reactions\s+(\d+\.\d+)\s" Wolf6x6Y.summary.txt | awk '{ print $3 }' > X.csv

# test if they look normal --  removed because it fails too often
#../normality_test.py X.csv
#if [ "$?" = 1 ] ; then
#  printf 'FAIL %s\n' "${test}"
#  let "fail = $fail + 4"
#fi

# test if the mean is 4.91
../ttest-mean.py X.csv 4.91
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test if the stdev is 0.1*mean
../chisqtest-sd.py X.csv 0.491
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# get values of Y_i_j
grep -Po "Y_\d+,\d+\s+reactions\s+(\d+\.\d+)\s" Wolf6x6Y.summary.txt | awk '{ print $3 }' > Y.csv

# test if they look normal -- removed because it fails too often
#../normality_test.py Y.csv
#if [ "$?" = 1 ] ; then
#  printf 'FAIL %s\n' "${test}"
#  let "fail = $fail + 32"
#fi

# test if the mean is 0.77
../ttest-mean.py Y.csv 0.77
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# test if the stdev is 0.1*mean
../chisqtest-sd.py Y.csv 0.077
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm Wolf6x6Y.summary.txt Wolf6x6Y.cps *.csv output
fi

exit $fail
