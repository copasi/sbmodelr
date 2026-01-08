#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# run sbmodelr
$PYTH ../../src/sbmodelr.py -s v --synapse-link-g -n ../sources/1to2to3.gv -o IzLinked.cps ../sources/IzhikevichBurstingNeuron.cps 3 1> output 2> /dev/null

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py IzLinked.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is a master synapse g
n=$(grep -Pc "^g_c_v_synapse\s+fixed" IzLinked.summary.txt)
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are two g_c...synapse that are assignments
n=$(grep -Pc "^g_c_v_[12],[23]_synapse\s+assignment" IzLinked.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# run the model with CopasiSE
../CopasiSE -c . --nologo IzLinked.cps > cpsout 2>&1
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the output was created and contains last time point
n=$(tail -n 3 I-tc.txt | grep -Pc "^300\s+")
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm IzLinked.summary.txt I-tc.txt cpsout output *.cps
fi

exit $fail
