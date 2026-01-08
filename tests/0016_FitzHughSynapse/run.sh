#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -s x -n ../sources/1to2.gv ../sources/FitzHugh-Nagumo.cps 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py FitzHugh-Nagumo_2.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is exactly one event
n=$(grep -Pc "event\s+Time > 2 && false" FitzHugh-Nagumo_2.summary.txt)
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are two ODEs for x
n=$(grep -Pc "x_([12])\s+ode.*Values\[c_\1\] \* \( Values\[y_\1\] \+ Values\[x_\1\] - Values\[x_\1\] \^ 3 \/ 3 \+ Values\[z_\1\] \)" FitzHugh-Nagumo_2.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# failures in the synaptic constants all generate error code 8
# check that the specific conductivity constant has been defined with the default value
if ! grep -Pq "^g_c_x_1\,2_synapse\s+fixed\s+1\.0" FitzHugh-Nagumo_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the tau_d constant has been defined with the default value
if ! grep -Pq "^tau_d_x_synapse\s+fixed\s+10\.0" FitzHugh-Nagumo_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the tau_r constant has been defined with the default value
if ! grep -Pq "^tau_r_x_synapse\s+fixed\s+0\.5" FitzHugh-Nagumo_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the V0 constant has been defined with the default value
if ! grep -Pq "^V0_x_synapse\s+fixed\s+-20\.0" FitzHugh-Nagumo_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the Vsyn constant has been defined with the default value
if ! grep -Pq "^Vsyn_x_synapse\s+fixed\s+20\.0" FitzHugh-Nagumo_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there is an ODE for br_1_2 (the post-synaptic bound receptor)
n=$(grep -Pc "br_x_1\,2\s+ode.*\( 1 \/ Values\[tau_r_x_synapse\] - 1 \/ Values\[tau_d_x_synapse\] \) \* \( 1 - Values\[br_x_1\,2\] \) \/ \( 1 \+ exp \( Values\[V0_x_synapse\] - Values\[x_1\] \) \) - Values\[br_x_1\,2\] \/ Values\[tau_d_x_synapse\]" FitzHugh-Nagumo_2.summary.txt)
if ((n != 1))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

#check that ODE for x_2 has the synaptic current term
if ! grep -Pq "^x_2\s+ode.*\+ Values\[g_c_x_1\,2_synapse\] \* Values\[br_x_1\,2] \* \( Values\[Vsyn_x_synapse\]" FitzHugh-Nagumo_2.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm FitzHugh-Nagumo_2.summary.txt output *.cps
fi

exit $fail
