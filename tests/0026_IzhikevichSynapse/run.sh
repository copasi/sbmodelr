#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# run sbmodelr
$PYTH ../../src/sbmodelr.py -s v -n ../sources/1to2to3.gv -o I-1-2-3.cps ../sources/IzhikevichBurstingNeuron.cps 3 1> output 2> /dev/null

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py I-1-2-3.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are three events that were variable-dependent
n=$(grep -Pc "event_0000001_([123])\s+Values\[v_\1\] gt Values\[Vthresh_\1\]" I-1-2-3.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there is only one event that is time-dependent
n=$(grep -Pc "Stimulus\s+Time gt10\s+" I-1-2-3.summary.txt)
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that there are three ODEs for v
n=$(grep -Pc "v_([123])\s+ode\s+" I-1-2-3.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there are two ODEs for br_i_j (the post-synaptic bound receptor)
n=$(grep -Pc "br_v_([12])\,([23])\s+ode.*\( 1 \/ Values\[tau_r_v_synapse\] - 1 \/ Values\[tau_d_v_synapse\] \) \* \( 1 - Values\[br_v_\1\,\2\] \) \/ \( 1 \+ exp \( Values\[V0_v_synapse\] - Values\[v_\1\] \) \) - Values\[br_v_\1\,\2\] \/ Values\[tau_d_v_synapse\]" I-1-2-3.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# run the model with CopasiSE
../CopasiSE -c . --nologo I-1-2-3.cps > cpsout 2>&1
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# check that the output was created and contains last time point
n=$(tail -n 3 I-tc.txt | grep -Pc "^300\s+")
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm I-1-2-3.summary.txt I-tc.txt cpsout output *.cps
fi

exit $fail
