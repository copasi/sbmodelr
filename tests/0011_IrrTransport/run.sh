#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py  -t X -n ../sources/1to2.gv -o Wolf2Irr.cps ../sources/Selkov-Wolf-Heinrich.cps 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# test -o option by checking if output file was created
if ! [ -f Wolf2Irr.cps ]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py Wolf2Irr.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that the transport reaction is irreversible
if ! grep -Pq "t_X_1-2\s+X_1 -> X_2\s+Mass action \(irreversible\)" Wolf2Irr.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm Wolf2Irr.summary.txt Wolf2Irr.cps output
fi

exit $fail
