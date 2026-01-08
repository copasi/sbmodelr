#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py ../sources/BindingKa.cps 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# validate the model with CopasiSE
../CopasiSE -c . --nologo --validate BindingKa_2.cps > cpsout 2>&1
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output tc.txt cpsout *.cps > /dev/null 2>&1
fi

exit $fail

