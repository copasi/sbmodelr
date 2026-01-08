#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# run sbmodelr
$PYTH ../../src/sbmodelr.py -s v -n ../sources/1to2to3.gv --cn 0.1 uni -o IzNoise.cps ../sources/IzhikevichBurstingNeuron.cps 3 1> output 2> /dev/null

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py IzNoise.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are two ODEs for br_i_j (the post-synaptic bound receptor)
n=$(grep -Pc "br_v_([12])\,([23])\s+ode" IzNoise.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are two Vsyn parameters
n=$(grep -Pc "^Vsyn_v_synapse_[12]-[23]\s+fixed" IzNoise.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that there are two tau_r parameters
n=$(grep -Pc "^tau_r_v_synapse_[12]-[23]\s+fixed" IzNoise.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there are two tau_d parameters
n=$(grep -Pc "^tau_d_v_synapse_[12]-[23]\s+fixed" IzNoise.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# run the model with CopasiSE
../CopasiSE -c . --nologo IzNoise.cps > cpsout 2>&1
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check that the output was created and contains last time point
n=$(tail -n 3 I-tc.txt | grep -Pc "^300\s+")
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm IzNoise.summary.txt I-tc.txt cpsout output *.cps
fi

exit $fail
