#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t B -n ../sources/twins.gv --pn A 0.3 uni ../sources/BIOMD0000000539_url.xml 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# check that output file is SBML
s=$(head -n 3 BIOMD0000000539_url_2.xml | grep -q "sbml xmlns")
if ((s != 0))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py BIOMD0000000539_url_2.xml >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# are there 9 kinetic functions?
if ! grep -Pq "Kinetic Functions\: 9 =" BIOMD0000000539_url_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that there are exactly two rate laws for transcription and translation_1
n=$(grep -Pc "transcription and translation_1_([12])\s+rhof_\1\*gB_\1" BIOMD0000000539_url_2.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there is one transport reaction
if ! grep -Pq "t_B_1-2\s+B_1 = B_2\s+Function for t_B_1-2" BIOMD0000000539_url_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# get values of [A_i]_0
grep -Po "A_[12]\s+reactions\s+(\d+\.\d+)" BIOMD0000000539_url_2.summary.txt | awk '{ print $3 }' > A.csv

# test that the minimum is above 40-40*0.3 = 28
m=$(awk 'NR == 1 || $1 < min {line = $0; min = $1}END{print line}' A.csv)
if (( $(echo "$m < 28" |bc -l) )); then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# test that the maximum is below 40+40*0.3 = 52
m=$(awk 'NR == 1 || $1 > max {line = $0; max = $1}END{print line}' A.csv)
if (( $(echo "$m > 52" |bc -l) )); then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm BIOMD0000000539_url_2.summary.txt output A.csv *.xml
fi

exit $fail
