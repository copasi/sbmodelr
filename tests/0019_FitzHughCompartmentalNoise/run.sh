#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -d x -n ../sources/1-2-3.gv --cn 0.2 norm ../sources/FitzHugh-Nagumo.cps 3 1> output 2> /dev/null

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

# check that two coupling constants have been created
n=$(grep -Pc "^k_x_coupling_[12]-[23]\s+fixed" FitzHugh-Nagumo_3.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are three ODEs connected correctly with diffusive term
n=$(grep -Pc "x_([123])\s+ode.*\+ Values\[k_x_coupling_[12]-[23]\] \* \( Values\[x_[123]\] - Values\[x_\1\]" FitzHugh-Nagumo_3.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm FitzHugh-Nagumo_3.summary.txt output *.cps
fi

exit $fail
