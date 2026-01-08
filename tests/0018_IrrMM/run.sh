#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -n ../sources/1to2.gv --Hill-transport c ../sources/BindingKa.cps 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py BindingKa_2.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that transport reaction exists
if ! grep -Pq "t_c_1-2\s+c_1 -\> c_2\s+Hill Cooperativity" BindingKa_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that the Henri-Michaelis-Menten rate law is used
if ! grep -Pq "^Hill Cooperativity\s+V\*\(substrate\/Shalve\)\^h\/\(1\+\(substrate\/Shalve\)\^h\)" BindingKa_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that the Km constant was created
if ! grep -Pq "^Km_c_transport\s+fixed" BindingKa_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the Vmax constant was created
if ! grep -Pq "^Vmax_c_transport\s+fixed" BindingKa_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check that the Hill constant was created
if ! grep -Pq "^h_c_transport\s+fixed" BindingKa_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm BindingKa_2.summary.txt output *.cps
fi

exit $fail
