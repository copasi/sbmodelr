#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -s x -n ../sources/1to2to3.gv ../sources/Hindmarsh-Rose.cps 3 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py Hindmarsh-Rose_3.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there are three ODEs for x
n=$(grep -Pc "x_([123])\s+ode.*\[y_\1\] - Values\[a_\1\] \* \[x_\1\] \^ 3 \+ Values\[b_\1\] \* \[x_\1\] \^ 2 \+ Values\[I_\1\] - \[z_\1\]" Hindmarsh-Rose_3.summary.txt)
if ((n != 3))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# failures in the synaptic constants all generate error code 8
# check that the specific conductivity constant has been defined with the default value
if ! grep -Pq "^g_c_x_1\,2_synapse\s+fixed\s+1\.0" Hindmarsh-Rose_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the tau_d constant has been defined with the default value
if ! grep -Pq "^tau_d_x_synapse\s+fixed\s+10\.0" Hindmarsh-Rose_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the tau_r constant has been defined with the default value
if ! grep -Pq "^tau_r_x_synapse\s+fixed\s+0\.5" Hindmarsh-Rose_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the V0 constant has been defined with the default value
if ! grep -Pq "^V0_x_synapse\s+fixed\s+-20\.0" Hindmarsh-Rose_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that the Vsyn constant has been defined with the default value
if ! grep -Pq "^Vsyn_x_synapse\s+fixed\s+20\.0" Hindmarsh-Rose_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there are two ODEs for br_1_2 and br_2_3 (the post-synaptic bound receptor)
n=$(grep -Pc "br_x_([12])\,([23])\s+ode.*\( 1 \/ Values\[tau_r_x_synapse\] - 1 \/ Values\[tau_d_x_synapse\] \) \* \( 1 - Values\[br_x_\1\,\2\] \) \/ \( 1 \+ exp \( Values\[V0_x_synapse\] - \[x_\1\] \) \) - Values\[br_x_\1\,\2\] \/ Values\[tau_d_x_synapse\]" Hindmarsh-Rose_3.summary.txt)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

#check that ODEs for x_2 and x_3 have the synaptic current term
if ! grep -Pq "^x_2\s+ode.*\+ Values\[g_c_x_1\,2_synapse\] \* Values\[br_x_1\,2] \* \( Values\[Vsyn_x_synapse\]" Hindmarsh-Rose_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if ! grep -Pq "^x_3\s+ode.*\+ Values\[g_c_x_2\,3_synapse\] \* Values\[br_x_2\,3] \* \( Values\[Vsyn_x_synapse\]" Hindmarsh-Rose_3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm Hindmarsh-Rose_3.summary.txt output *.cps
fi

exit $fail
