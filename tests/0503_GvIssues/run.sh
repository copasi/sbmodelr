#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# test an empty network file
$PYTH ../../src/sbmodelr.py -t G -n ../sources/empty.gv -o t.cps ../sources/GeneExpressionUnit.cps 2 1> output 2> /dev/null

# check that the correct warning is issued
if ! grep -q " Warning: ../sources/empty.gv did not contain any valid edges, no connections added" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

# test a network file without edges
$PYTH ../../src/sbmodelr.py -t G -n ../sources/no-edges.gv -o t.cps ../sources/GeneExpressionUnit.cps 2 1> output 2> /dev/null

# check that the correct warning is issued
if ! grep -q " Warning: ../sources/no-edges.gv did not contain any valid edges, no connections added" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# test a network file with repeated edges
$PYTH ../../src/sbmodelr.py -t G -n ../sources/duplicates.gv -o t.cps ../sources/GeneExpressionUnit.cps 4 1> output 2> /dev/null

# check that the correct warning is issued
if ! grep -q " Warning: duplicate entry for connection 1 to 3, ignored" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# create model summary
../model_report.py t.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is exactly one transport reaction 1->3 (the duplicate)
n=$(grep -Pc "^t_G_1-3\s+G_1 = G_3\s+Mass action" t.summary.txt )
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test a network file with non-numeric nodes
$PYTH ../../src/sbmodelr.py -t G -n ../sources/non-numeric-nodes.gv -o t.cps ../sources/GeneExpressionUnit.cps 3 1> output 2> /dev/null

# check that the correct warning is issued
if ! grep -q " Warning: ../sources/non-numeric-nodes.gv did not contain any valid edges, no connections added" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm t.summary.txt output *.cps
fi

exit $fail
