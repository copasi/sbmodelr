#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -d Calcium-Cyt -c 0.1 ../sources/CalciumSpiking.cps 2 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py CalciumSpiking_2x2.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that the coupling constant has been defined with the correct value
if ! grep -Pq "^k_Calcium-Cyt_coupling\s+fixed\s+0\.1" CalciumSpiking_2x2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are four ODEs connected correctly with diffusive term
n=$(grep -Pc "Calcium-Cyt_([12]\,[12])\s+ode\s+.*\+ Values\[k_Calcium-Cyt_coupling\] \* \( \[Calcium-Cyt_[12]\,[12]\] - \[Calcium-Cyt_\1\] \) \+ Values\[k_Calcium-Cyt_coupling\] \* \( \[Calcium-Cyt_[12]\,[12]\] - \[Calcium-Cyt_\1\] \)" CalciumSpiking_2x2.summary.txt )
if ((n != 4))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm CalciumSpiking_2x2.summary.txt output *.cps
fi

exit $fail
