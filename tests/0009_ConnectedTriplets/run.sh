#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t b -n ../sources/3_fully_connected.dot ../sources/BindingKa.cps 3 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py BindingKa_3.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are exactly three transport reactions
n=$(grep -Pc "t_b_([123])-([123])\s+b_\1 = b_\2" BindingKa_3.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that the expressions for the koff_i assignments are correct
n=$(grep -Pc "koff_([123])\s+assignment\s+[0-9\.]+\s+Values\[kon_\1\]\.InitialValue \/ Values\[Ka_\1\]\.InitialValue" BindingKa_3.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm BindingKa_3.summary.txt output *.cps
fi

exit $fail
