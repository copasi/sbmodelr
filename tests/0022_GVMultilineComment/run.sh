#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t b -n ../sources/net4mlcomm.gv ../sources/BindingKa.cps 4 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py BindingKa_4.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are 4 transport reactions
n=$(grep -Pc "^t_b_[12]-[234]\s+b_[12] = b_[234]" BindingKa_4.summary.txt)
if ((n != 4))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check if a reaction 3 -- 4 exists
if grep -Pq "^t_b_3-4\s+b_3 = b_4" BindingKa_4.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm BindingKa_4.summary.txt output *.cps
fi

exit $fail
