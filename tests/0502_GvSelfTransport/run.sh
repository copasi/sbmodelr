#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t X -n ../sources/self.gv -o ISC.cps ../sources/Selkov-Wolf-Heinrich.cps 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# test -o option by checking if output file was created
if ! [ -f ISC.cps ]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py ISC.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is a transport reaction 1->2
if ! grep -Pq "t_X_1-2\s+X_1 -> X_2" ISC.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there is no transport reaction 2->2
if grep -Pq "t_X_2-2\s+X_2 -> X_2" ISC.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that there is only one transport constant
n=$(grep -Pc "^k_.+_transport\s+fixed" ISC.summary.txt )
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm ISC.summary.txt output *.cps
fi

exit $fail
