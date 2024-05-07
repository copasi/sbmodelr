#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
../../sbmodelr  -t Y -k 0.1 --pn X 0.1 norm --pn Y 0.1 norm  -o Wolf4x4Y.cps ../sources/Selkov-Wolf-Heinrich.cps 4 4 > output

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# test -o option by checking if output file was created
if ! [ -f Wolf4x4Y.cps ]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py Wolf4x4Y.cps

# get values of X_i_j
grep -Po "X_\d+,\d+\s+reactions\s+(\d+\.\d+)\s" Wolf4x4Y.summary.txt | awk '{ print $3 }' > X.csv
# test if they look normal
../shapiro-wilk.py X.csv
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# test if the mean is 4.91
../ttest-mean.py X.csv 4.91
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test if the stdev is 0.1
../chisqtest-sd.py X.csv 0.1
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# get values of Y_i_j
grep -Po "Y_\d+,\d+\s+reactions\s+(\d+\.\d+)\s" Wolf4x4Y.summary.txt | awk '{ print $3 }' > Y.csv
# test if they look normal
../shapiro-wilk.py Y.csv
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# test if the mean is 0.77
../ttest-mean.py Y.csv 0.77
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# test if the stdev is 0.1
../chisqtest-sd.py Y.csv 0.1
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm Wolf4x4Y.summary.txt Wolf4x4Y.cps *.csv output
fi

exit $fail