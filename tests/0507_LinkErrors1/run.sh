#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

fail=0

# test non-existent species
$PYTH ../../src/sbmodelr.py -t J -n ../sources/twins.gv ../sources/PulsedDrug.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: Species J does not exist in the model" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = 1"
fi

# test transport on fixed species
$PYTH ../../src/sbmodelr.py -t S -n ../sources/twins.gv ../sources/PulsedDrug.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: S is a species that does not depend on reactions, no transport reactions can be added" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# test transport between same unit
$PYTH ../../src/sbmodelr.py -t A -n ../sources/self.gv ../sources/PulsedDrug.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q " Warning: transport on the same unit not allowed, ignoring 2 -> 2" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# test synaptic link in grid configuration
$PYTH ../../src/sbmodelr.py -s v ../sources/IzhikevichBurstingNeuron.cps 2 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: 2D or 3D grids cannot have synaptic connections" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

# test synaptic link in undirected network
$PYTH ../../src/sbmodelr.py -s v -n ../sources/twins.gv ../sources/IzhikevichBurstingNeuron.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q " Warning: network was defined as undirected, but synapses will be added as directed connections" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 16"
fi

# test synaptic link on global quantity that is not an ODE
$PYTH ../../src/sbmodelr.py -s i -n ../sources/1to2.gv ../sources/IzhikevichBurstingNeuron.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: i is a global variable that is not an ODE" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 32"
fi

# test synaptic link with noise and linked g
$PYTH ../../src/sbmodelr.py -s i --cn 0.1 uni --synapse-link-g -n ../sources/1to2.gv ../sources/IzhikevichBurstingNeuron.cps 2 > output 2>&1

# check that the correct error is issued
if ! grep -q "ERROR: --cn and --synapse-link-g options cannot be used together, chose only one!" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 64"
fi

# test diffusive connection on species that is not an ODE
$PYTH ../../src/sbmodelr.py -d A -n ../sources/twins.gv ../sources/PulsedDrug.cps 2 > output 2>&1

# check that the correct warning is issued
if ! grep -q "ERROR: A is a species but it is not of type ODE" output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 128"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output *.cps > /dev/null 2>&1
fi

exit $fail
