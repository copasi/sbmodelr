#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t a ../sources/ATPlocs.cps 2 >output 2>&1

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output *.cps > /dev/null 2>&1
fi

exit $fail
