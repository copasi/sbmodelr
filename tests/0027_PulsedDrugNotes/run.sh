#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# run sbmodelr
$PYTH ../../src/sbmodelr.py  -t drug --pn Drug_dose 0.5 uni -n ../sources/twins.gv -o PD2.cps ../sources/PulsedDrug.cps 2 1> output 2> /dev/null

# compare output and target
difference=$(diff output target_stdout)
if [[ $difference ]]; then
  printf 'FAIL %s\n' "${test}"
  fail=1
fi

# poor man's way of checking existence of notes: grepping the cps file...

# check that there are two notes for the compartment
n=$(grep -c "a location within which every species is well mixed" PD2.cps)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# check that there are two notes for the S species
n=$(grep -c "fixed substrate" PD2.cps)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# check that there are two notes for the Drug_dose global quantity
n=$(grep -c "<p>the <i>amount</i> of drug added each time</p></body>" PD2.cps)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# check that there are two notes for reactions R1_i
n=$(grep -c "first reaction" PD2.cps)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# check that there are two notes for the events DrugOn_i
n=$(grep -c "<p>event that turns <i>on</i> drug addition</p></body>" PD2.cps)
if ((n != 2))  ; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# validate the model with CopasiSE
../CopasiSE -c . --nologo --validate PD2.cps > /dev/null 2>&1
if ! [[ $? = 0 ]]; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm  PD2.cps output > /dev/null 2>&1
fi

exit $fail
