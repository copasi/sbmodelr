#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# test model that has scan over parameter sets
$PYTH ../../src/sbmodelr.py -t M1 ../sources/multimodel.cps 2 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: a scan of parameter sets exists in the original model but was not included in the new model." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

# check that the other warning is also issued
if ! grep -q " Warning: in Parameter scan task the scanned or sampled items were converted to those of the first unit only." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# test model with a cross section
$PYTH ../../src/sbmodelr.py -t Y ../sources/crossect.cps 2 2 > output 2>&1

# check that the correct warning is issued for cross section
if ! grep -q " Warning: the cross section task was updated to use \[X_1,1\] as variable." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that the correct warning is issued for parameter scan
if ! grep -q " Warning: in Parameter scan task the scanned or sampled items were converted to those of the first unit only." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test model with a sensitivity analysis
$PYTH ../../src/sbmodelr.py -t B ../sources/sens.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: sensitivies task is now using items of unit " output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# test model with optimization
$PYTH ../../src/sbmodelr.py -t X ../sources/opt.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: in Optimization task the objective function and the search parameters were converted to those of the first unit only." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# test model with parameter estimation
$PYTH ../../src/sbmodelr.py -t E1 ../sources/pe_tcs.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: Parameter Estimation task settings were not copied to the new model." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# check that the correct warning is issued for time course sensitivities
if ! grep -q " Warning: Time Course Sensitivities task settings were not copied to the new model." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output *.cps > /dev/null 2>&1
fi

exit $fail
