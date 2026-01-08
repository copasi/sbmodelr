#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# undirected, 2 nodes with labels, etc.
$PYTH ../../src/sbmodelr.py -t G -o t.cps -n ../sources/test08_2.gv ../sources/GeneExpressionUnit.cps 3 1> output 2> /dev/null

if ! grep -q "created new model t.cps with a set of 3 replicas of ../sources/GeneExpressionUnit.cps" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
elif grep -q "Warning" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
elif grep -q "Error" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

# create model summary
../model_report.py t.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is only 1 transport reaction
n=$(grep -Pc "^t_G_[123]-[123]\s+G_[123] = G_[123]" t.summary.txt)
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

# directed, 6 nodes with labels, etc.
$PYTH ../../src/sbmodelr.py -g G --ignore-compartments -o t22.cps -n ../sources/test22_6.gv ../sources/GeneExpressionUnit.cps 6 1> output 2> /dev/null

if ! grep -q "created new model t22.cps with a set of 6 replicas of ../sources/GeneExpressionUnit.cps" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
elif grep -q "Warning" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
elif grep -q "Error" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py t22.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are exactly 4 regulated synthesis reactions
n=$(grep -Pc "^synthesis G_(\d)\s+\-\>\s+G_\1;\s+G_\d" t22.summary.txt)
if ((n != 4))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 2 ]]; then
   let "fail = $fail + 2"
  fi
fi

# directed, 16 nodes with labels, landscape orientation etc.
$PYTH ../../src/sbmodelr.py -g G --ignore-compartments -o t16.cps -n ../sources/test16_12.gv ../sources/GeneExpressionUnit.cps 12 1> output 2> /dev/null

if ! grep -q "created new model t16.cps with a set of 12 replicas of ../sources/GeneExpressionUnit.cps" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
elif grep -q "Warning" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
elif grep -q "Error" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# create model summary
../model_report.py t16.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are exactly 11 regulated synthesis reactions
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1;\s+G_\d+" t16.summary.txt)
if ((n != 11))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 4 ]]; then
   let "fail = $fail + 4"
  fi
fi

# directed, 15 nodes with labels, landscape orientation etc.
$PYTH ../../src/sbmodelr.py -g G --ignore-compartments -o t15.cps -n ../sources/test15_11.gv ../sources/GeneExpressionUnit.cps 11 1> output 2> /dev/null

if ! grep -q "created new model t15.cps with a set of 11 replicas of ../sources/GeneExpressionUnit.cps" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
elif grep -q "Warning" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
elif grep -q "Error" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# create model summary
../model_report.py t15.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are exactly 11 regulated synthesis reactions
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1;\s+G_\d+" t15.summary.txt)
if ((n != 11))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 8 ]]; then
   let "fail = $fail + 8"
  fi
fi

# check that there are exactly 10 reactions with regulated_by_2 kinetics
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1.*regulated_by_2" t15.summary.txt)
if ((n != 10))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 8 ]]; then
   let "fail = $fail + 8"
  fi
fi

# check that there is exactly 1 reaction with regulated_by_5 kinetics
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1.*regulated_by_5" t15.summary.txt)
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 8 ]]; then
   let "fail = $fail + 8"
  fi
fi

# directed 8 nodes with ranks
$PYTH ../../src/sbmodelr.py -g G --ignore-compartments -o t10.cps -n ../sources/test10_8.gv ../sources/GeneExpressionUnit.cps 8 1> output 2> /dev/null

if ! grep -q "created new model t10.cps with a set of 8 replicas of ../sources/GeneExpressionUnit.cps" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
elif grep -q "Warning" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
elif grep -q "Error" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi


# create model summary
../model_report.py t10.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are exactly 6 regulated synthesis reactions
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1;\s+G_\d+" t10.summary.txt)
if ((n != 6))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 16 ]]; then
   let "fail = $fail + 16"
  fi
fi

# check that there are exactly 2 reactions with regulated_by_2 kinetics
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1.*regulated_by_2" t10.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 16 ]]; then
   let "fail = $fail + 16"
  fi
fi

# check that there are exactly 4 reactions with regulated_by_1 kinetics
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1.*regulated_by_1" t10.summary.txt)
if ((n != 4))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 16 ]]; then
   let "fail = $fail + 16"
  fi
fi

# directed 12 nodes with groups
$PYTH ../../src/sbmodelr.py -g G --ignore-compartments -o t11.cps -n ../sources/test11_12.gv ../sources/GeneExpressionUnit.cps 12 1> output 2> /dev/null

if ! grep -q "created new model t11.cps with a set of 12 replicas of ../sources/GeneExpressionUnit.cps" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
elif grep -q "Warning" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
elif grep -q "Error" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# create model summary
../model_report.py t11.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are exactly 9 regulated synthesis reactions
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1;\s+G_\d+" t11.summary.txt)
if ((n != 9))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 32 ]]; then
   let "fail = $fail + 32"
  fi
fi

# check that there are exactly 4 reactions with regulated_by_2 kinetics
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1.*regulated_by_2" t11.summary.txt)
if ((n != 4))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 32 ]]; then
   let "fail = $fail + 32"
  fi
fi

# check that there are exactly 5 reactions with regulated_by_1 kinetics
n=$(grep -Pc "^synthesis G_(\d+)\s+\-\>\s+G_\1.*regulated_by_1" t11.summary.txt)
if ((n != 5))  ; then
  printf 'FAIL %s\n' "${test}"
  if [[ $fail -lt 32 ]]; then
   let "fail = $fail + 32"
  fi
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm *.summary.txt output *.cps
fi

exit $fail
