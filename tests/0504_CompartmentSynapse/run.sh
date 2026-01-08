#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# run sbmodelr
$PYTH ../../src/sbmodelr.py -s blob -n ../sources/1to2to3.gv ../sources/blob.cps 3 1> output 2> /dev/null

# check that the an error is issued
if ! grep -q "ERROR: blob is a compartment ODE, but compartments cannot have synaptic links" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output
fi

exit $fail
