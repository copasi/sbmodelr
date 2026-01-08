#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -g G --pn G 0.8 uni --grn-a -1 --grn-h 4 -n ../sources/3circle.gv -o repressilator.cps --ignore-compartments ../sources/GeneExpressionUnit.cps 3 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py repressilator.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is one rate law for regulation by 1
if ! grep -Pq "^regulated_by_1\s+V \* \( \( 1 \+ \( 1 \+ a1 \) \* M1 \^ h1 \) \/ \( 1 \+ M1 \^ h1 \) \)" repressilator.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are three species
n=$(grep -Pc "^G_[123]\s+reactions\s+[0-9]+\.[0-9]+" repressilator.summary.txt)
if ((n != 3)); then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that min and max of initial concentrations are different (meaning noise worked!)
n=$(grep -Po "^G_[123]\s+reactions\s+([0-9]+\.[0-9]+\s)" repressilator.summary.txt | awk '(NR==1){Min=$3;Max=$3};(NR>=2){if(Min>$3) Min=$3;if(Max<$3) Max=$3} END {if(Min<Max) {print 0} else {print 1}}')
if ((n != 0)) ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that no new compartments were created
n=$(grep -Pc "^cell_[123]\s+fixed" repressilator.summary.txt)
if ((n != 0))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# checks that there are three synthesis reactions
n=$(grep -Pc "synthesis (G_[123])\s+-> \1\;\s+G_[123]\s+regulated_by_1" repressilator.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm repressilator.summary.txt output *.cps
fi

exit $fail
