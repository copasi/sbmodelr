#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py  -t Y -k 0.1 -o Wolf3x3x3Y.cps ../sources/Selkov-Wolf-Heinrich.cps 3 3 3 1> output 2> /dev/null
fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# test -o option by checking if output file was created
if ! [ -f Wolf3x3x3Y.cps ]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py Wolf3x3x3Y.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# was the name of the model written correctly?
if ! grep -Pq "a 3D set of 27 replicas \(3x3x3\) of Simple model of glycolytic oscillations" Wolf3x3x3Y.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# are there 135 reactions?
if ! grep -Pq "Reactions: 135 =" Wolf3x3x3Y.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# are there 55 Global Quantities?
if ! grep -Pq "Global Quantities: 55 =" Wolf3x3x3Y.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check that there is no Cell_0,0,0
if grep -Pq "Cell_0,0,0" Wolf3x3x3Y.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# check that there is a Cell_3,3,3
if ! grep -Pq "Cell_3,3,3" Wolf3x3x3Y.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm Wolf3x3x3Y.summary.txt Wolf3x3x3Y.cps output
fi

exit $fail
