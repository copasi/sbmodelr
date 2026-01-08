#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -d blob -c 0.05 --pn blob 0.9 uni ../sources/blob.cps 3 3 3 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py blob_3x3x3.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that the coupling constant has been defined with the correct value
if ! grep -Pq "^k_blob_coupling\s+fixed\s+0.05" blob_3x3x3.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are 27 compartments
n=$(grep -Pc "^blob_[123],[123],[123]\s+ode" blob_3x3x3.summary.txt )
if ((n != 27))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that the ode for blob_2,2,2 (the center unit) has 6 diffusive terms
n=$(grep -P "^blob_2,2,2\s+ode" blob_3x3x3.summary.txt | grep -Po "\+ Values\[k_blob_coupling\] \* \("| wc -l)
if ((n != 6))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that min and max of [blob]_0 are different (checks that --pn worked)
n=$(grep -Po "^blob_[123],[123],[123]\s+ode\s+(.[0-9]+\.[0-9]+\s)" blob_3x3x3.summary.txt | awk '(NR==1){Min=$3;Max=$3};(NR>=2){if(Min>$3) Min=$3;if(Max<$3) Max=$3} END {if(Min == Max) {print 1} else {print 0} }')
if ((n != 0)) ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# validate the model with CopasiSE
../CopasiSE -c . --nologo --validate blob_3x3x3.cps > cpsout 2>&1
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm blob_3x3x3.* cpsout output
fi

exit $fail
