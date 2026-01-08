#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py  --pn Ka 0.2 uni ../sources/BindingKa.cps 100 1> output 2> /dev/null

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

# test that the minimum is above 1000-1000*0.2 = 800
m=$(awk 'NR == 1 || $1 < min {line = $0; min = $1}END{print line}' Ka.csv)
#printf 'MIN %s\n' "$m"
if (( $(echo "$m < 800" |bc -l) )); then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# test that the maximum is below 1000+1000*0.2 = 1200
m=$(awk 'NR == 1 || $1 > max {line = $0; max = $1}END{print line}' Ka.csv)
#printf 'MAX %s\n' "$m"
if (( $(echo "$m > 1200" |bc -l) )); then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# this would calculate mean and stdved using only grep and awk !
#grep -Po "Ka_\d+\s+fixed\s+(\d+\.\d+)" BindingKa_1000.summary.txt | awk '{ sum+=$3; sumsq +=$3^2 } END { print sum/NR, NR*(sqrt((sumsq-sum^2/NR)/NR))/sum }'


if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm BindingKa_100.summary.txt Ka.csv output *.cps
fi

exit $fail
