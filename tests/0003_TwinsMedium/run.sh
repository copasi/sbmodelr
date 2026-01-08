#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -n ../sources/twins.gv -t c --add-medium ../sources/BindingKa.cps 2 1> output 2> /dev/null

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

# check that all compartments were created
n=$(grep -Pc "^compartment_[12]\s+fixed" BindingKa_2.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that medium exists
if ! grep -Pq "medium\s+fixed" BindingKa_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that transport reactions to the medium exist
n=$(grep -Pc "t_c_[12]-medium\s+c_[12] = c_medium\s+Mass action \(reversible\)" BindingKa_2.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

#check that all units were created
n=$(grep -Pc "^[abc]_[12]\s+reactions" BindingKa_2.summary.txt)
if ((n != 6))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check that the transport between two units exists
if ! grep -Pq "t_c_1-2\s+c_1 = c_2\s+Mass action \(reversible\)" BindingKa_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm BindingKa_2.summary.txt output *.cps
fi

exit $fail
