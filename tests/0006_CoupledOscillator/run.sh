#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# run sbmodelr
$PYTH ../../src/sbmodelr.py  -t X -k 0.1 --add-medium -o Wolf2X.cps ../sources/Selkov-Wolf-Heinrich.cps 2 1> output 2> /dev/null

fail=0

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# test -o option by checking if output file was created
if ! [ -f Wolf2X.cps ]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# create model summary
../model_report.py Wolf2X.cps >/dev/null
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  exit -1
fi

# check that the transport rate constant is set properly
grep -Po "k_X_transport\s+fixed\s+0.1" Wolf2X.summary.txt >/dev/null
if [ "$?" != 0  ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# run the file with CopasiSE
../CopasiSE -c . --nologo Wolf2X.cps > /dev/null 2>&1
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that output file has time course header with all species
grep -Po "# Time\s+X_1\s+X_2\s+Y_1\s+Y_2\s+X_medium" tc.txt >/dev/null
if [ "$?" != 0  ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check that output file has last time point with 2 values equal (X) and two others equal (Y) plus one more
grep -Po "100\t([\d\.]+)\t\1\t([\d\.]+)\t\2\t[\d\.]+" tc.txt >/dev/null
if [ "$?" != 0  ] ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm Wolf2X.summary.txt tc.txt Wolf2X.cps output > /dev/null 2>&1
fi

exit $fail
