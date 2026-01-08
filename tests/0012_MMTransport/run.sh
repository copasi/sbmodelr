#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py --Hill-transport C --transport-Km 0.15 --transport-Vmax 0.5 -o IRC_2x2.cps ../sources/IrreversibleReactionChain.cps 2 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# test -o option by checking if output file was created
if ! [ -f IRC_2x2.cps ]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py IRC_2x2.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is no compartment_0,0
if grep -Pq "compartment_0,0" IRC_2x2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that there is a compartment_2,2
if ! grep -Pq "compartment_2,2" IRC_2x2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there are 4 transport reactions
n=$(grep -Pc "t_C_[12]\,[12]-[12]\,[12]\s+C_[12]\,[12] = C_[12]\,[12]\s+Hill Transport" IRC_2x2.summary.txt)
if ((n != 4))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check the rate law for transport is defined
if ! grep -Pq "^Hill Transport\s+V \* \( S \^ h \- P \^ h \) \/ \( Km \^ h \+ S \^ h \+ P \^ h \)" IRC_2x2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# check that VMax is defined and has correct value
if ! grep -Pq "Vmax_C_transport\s+fixed\s+0.50" IRC_2x2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# check that Km is defined and has correct value
if ! grep -Pq "Km_C_transport\s+fixed\s+0.15" IRC_2x2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm IRC_2x2.summary.txt IRC_2x2.cps output
fi

exit $fail
