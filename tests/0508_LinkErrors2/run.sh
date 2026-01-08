#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# test compartment ode with --ignore_compartments
$PYTH ../../src/sbmodelr.py -d vesicle -n ../sources/twins.gv --ignore-compartments ../sources/shrink.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: vesicle is a compartment but ignore_compartments is set, nothing done" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

# test diffusive link on compartment that is not an ode
$PYTH ../../src/sbmodelr.py -d ER -n ../sources/twins.gv ../sources/CalciumSpiking.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: ER is a compartment but it is not of type ODE" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# test synaptic connection on a compartment
$PYTH ../../src/sbmodelr.py -s vesicle -n ../sources/1to2.gv ../sources/shrink.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: vesicle is a compartment ODE, but compartments cannot have synaptic links" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# test diffusive connection between the same compartment
$PYTH ../../src/sbmodelr.py -d vesicle -n ../sources/self.gv ../sources/shrink.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: diffusive coupling onto the same unit not allowed, ignoring 2 -> 2" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test diffusive link on a non-existent entity
$PYTH ../../src/sbmodelr.py -d foo ../sources/shrink.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: foo is not a valid model entity" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# test regulatory link on 2D grid
$PYTH ../../src/sbmodelr.py -g G ../sources/GeneExpressionUnit.cps 2 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: regulatory synthesis reactions can only be created with dimension 1 and through a network (use option -n)" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
  let "rr = 1"
fi

# test regulatory link on 2D grid
$PYTH ../../src/sbmodelr.py -g G ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: regulatory synthesis reactions can only be created with dimension 1 and through a network (use option -n)" output; then
  printf 'FAIL %s\n' "${test}"
  if [ "$rr" != 1 ] ; then
	let "fail = $fail + 32"
  fi
fi

# test regulatory connection without --ignore_compartments
$PYTH ../../src/sbmodelr.py -g G -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: option --ignore-compartments is often desirable for building regulatory networks" output; then
  printf 'FAIL %s\n' "${test}"
  # if we saw this error already don't add error code
  let "fail = $fail + 64"
fi

# test regulatory connection with inexistent species
$PYTH ../../src/sbmodelr.py -g H --ignore-compartments -n ../sources/1to2.gv ../sources/GeneExpressionUnit.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: Species H does not exist in the model, no regulatory synthesis reactions added" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output *.cps > /dev/null 2>&1
fi

exit $fail
