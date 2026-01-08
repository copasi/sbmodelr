#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py -g G --cn 2 uni --grn-a -1 --grn-h 8 --grn-V 2 -n ../sources/CoopSW8_001.dot -o CoopSW8_001.cps --ignore-compartments ../sources/GeneExpressionUnit.cps 100 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# create model summary
../model_report.py CoopSW8_001.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that there is one rate law for regulation by 2
if ! grep -Pq "^regulated_by_2\s+V \* \( \( 1 \+ \( 1 \+ a1 \) \* M1 \^ h1 \) \/ \( 1 \+ M1 \^ h1 \) \) \* \( \( 1 \+ \( 1 \+ a2 \) \* M2 \^ h2 \) \/ \( 1 \+ M2 \^ h2 \) \)" CoopSW8_001.summary.txt; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are 100 species
n=$(grep -Pc "^G_\d+\s+reactions\s+\d+\.\d+" CoopSW8_001.summary.txt)
if ((n != 100)); then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that min and max of a_synth_G_i-j are within bounds
n=$(grep -Po "^a_synth_G_\d+-\d+\s+fixed\s+(.[0-9]+\.[0-9]+\s)" CoopSW8_001.summary.txt | awk '(NR==1){Min=$3;Max=$3};(NR>=2){if(Min>$3) Min=$3;if(Max<$3) Max=$3} END {if(Min<-1.0 || Max>1.0) {print 1} else {print 0} }')
if ((n != 0)) ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that min and max of h_synth_G_i-j are within bounds
n=$(grep -Po "^h_synth_G_\d+-\d+\s+fixed\s+(.[0-9]+\.[0-9]+\s)" CoopSW8_001.summary.txt | awk '(NR==1){Min=$3;Max=$3};(NR>=2){if(Min>$3) Min=$3;if(Max<$3) Max=$3} END {if(Min<1.0 || Max>10.0) {print 1} else {print 0} }')
if ((n != 0)) ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check that V_synth_G_i-j is positive
n=$(grep -Po "^V_synth_G_\d+\s+fixed\s+(.[0-9]+\.[0-9]+\s)" CoopSW8_001.summary.txt | awk '(NR==1){Min=$3};(NR>=2){if(Min>$3) Min=$3} END {if(Min<0.0) {print 1} else {print 0} }')
if ((n != 0)) ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# are there 300 reactions ?
if ! grep -Pq "Reactions\: 300 =" CoopSW8_001.summary.txt ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# check that no new compartments were created
n=$(grep -Pc "^cell_[123]\s+fixed" CoopSW8_001.summary.txt)
if ((n != 0))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm CoopSW8_001.summary.txt output *.cps
fi

exit $fail
