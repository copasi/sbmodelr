#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py  --pn Ka 0.2 norm ../sources/BindingKa.cps 100 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py BindingKa_100.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that 100 Ka parameters exist
n=$(grep -Pc "Ka_\d+\s+fixed" BindingKa_100.summary.txt)
if ((n != 100))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# get values of Ka
grep -Po "Ka_\d+\s+fixed\s+(\d+\.\d+)" BindingKa_100.summary.txt | awk '{ print $3 }' > Ka.csv
# test if they look normal
../normality_test.py Ka.csv
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# test if the mean is 1000
../ttest-mean.py Ka.csv 1000
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test if the stdev is 0.2
../chisqtest-sd.py Ka.csv 0.2
if [ "$?" = 1 ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# this would calculate mean and stdved using only grep and awk !
#grep -Po "Ka_\d+\s+fixed\s+(\d+\.\d+)" BindingKa_1000.summary.txt | awk '{ sum+=$3; sumsq +=$3^2 } END { print sum/NR, NR*(sqrt((sumsq-sum^2/NR)/NR))/sum }'

# check that the transport between two units exists
#if ! grep -Pq "t_c_1-2\s+c_1 = c_2\s+Mass action \(reversible\)" BindingKa_2.summary.txt; then
#  printf 'FAIL %s\n' "${test}"
#  fail=1
#fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm BindingKa_100.summary.txt Ka.csv output *.cps
fi

exit $fail
