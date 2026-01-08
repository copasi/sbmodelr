#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -d x -n ../sources/1-2-3.gv ../sources/FitzHugh-Nagumo.cps 3 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py FitzHugh-Nagumo_3.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is exactly one event
n=$(grep -Pc "event\s+Time > 2 && false" FitzHugh-Nagumo_3.summary.txt)
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that three z ODEs are targets of the event
if ! grep -Pq "^event\s+Time.*target.*Values\[z_1\].*Values\[z_2\].*Values\[z_3\]" FitzHugh-Nagumo_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that the coupling constant has been defined with the default value
if ! grep -Pq "^k_x_coupling\s+fixed\s+1\.0" FitzHugh-Nagumo_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there are three ODEs connected correctly with diffusive term
n=$(grep -Pc "x_([123])\s+ode.*\+ Values\[k_x_coupling\] \* \( Values\[x_[123]\] - Values\[x_\1\]" FitzHugh-Nagumo_3.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm FitzHugh-Nagumo_3.summary.txt output *.cps
fi

exit $fail
