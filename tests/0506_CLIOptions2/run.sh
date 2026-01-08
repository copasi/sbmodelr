#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# test network file with dim > 1
$PYTH ../../src/sbmodelr.py -g G --grn-a 1 -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: network file is only relevant for dimension 1 but you chose dimension 2 (2x2)" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

# test network file with more units than specified
$PYTH ../../src/sbmodelr.py -g G -n ../sources/CoopSW8_001.dot ../sources/GeneExpressionUnit.cps 3 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: network file lists nodes with numbers outside \[1,3\]" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# test grn with undirected graph
$PYTH ../../src/sbmodelr.py -g G -n ../sources/twins.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: regulatory connections require a directed graph" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# test loading invalid file
$PYTH ../../src/sbmodelr.py -t G void.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: void.cps failed to load." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test invalid network file
$PYTH ../../src/sbmodelr.py -g G -n void.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: file void.gv not found" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# test grn-a value outside bounds
$PYTH ../../src/sbmodelr.py -g G --grn-a "-3" -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: '--grn-a' value must be in the interval \[-1,1\]" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# since we are still checking bounds for --grn-a, we'll use the same error code
$PYTH ../../src/sbmodelr.py -g G --grn-a 2 -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: '--grn-a' value must be in the interval \[-1,1\]" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# check --grn-h value is within bounds
$PYTH ../../src/sbmodelr.py -g G --grn-h 0 -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "sbmodelr: error: argument --grn-h: value not in range 1-10" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# continue issuing the same error code for --grn-h out of bounds
$PYTH ../../src/sbmodelr.py -g G --grn-h 11 -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "sbmodelr: error: argument --grn-h: value not in range 1-10" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# continue issuing the same error code for --grn-h out of bounds
$PYTH ../../src/sbmodelr.py -g G --grn-h 0.5 -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "sbmodelr: error: argument --grn-h: value not an integer" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# test adding medium without any species transported
$PYTH ../../src/sbmodelr.py --add-medium ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: no medium unit created because no species are being transported or ODEs coupled" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output *.cps > /dev/null 2>&1
fi

exit $fail
