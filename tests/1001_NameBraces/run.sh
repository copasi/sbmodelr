#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# run sbmodelr
sbmodelr -t a ../sources/ATPlocs.cps 2 >output 2>&1

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

exit $fail
