#!/bin/bash

# work out our folder name
test=${PWD##*/}          # to assign to a variable
test=${test:-/}          # to correct for the case where PWD=/

# create a completely random file called rdninput.cps
#< /dev/urandom tr -dc "q[:print:]" | head -c40000 > rndinput.cps
basenc --z85 -w 0 /dev/urandom | head -c 40000 > rndinput.cps

fail=0

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t B rndinput.cps 2 1> output 2> /dev/null

if ! grep -Pq "ERROR: rndinput.cps failed to load." output; then
  printf 'FAIL %s\n' "${test}"
  fail = 1
fi

# create a cps file that has a header but is followed by random input
printf "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" > rndinput2.cps
printf "<!-- generated with COPASI 4.42 (Build 284) (http://www.copasi.org) at 2024-04-26T02:00:17Z -->\n" >>rndinput2.cps
printf "<?oxygen RNGSchema=\"http://www.copasi.org/static/schema/CopasiML.rng\" type=\"xml\"?>\n" >>rndinput2.cps
printf "<COPASI xmlns=\"http://www.copasi.org/static/schema\" versionMajor=\"4\" versionMinor=\"42\" versionDevel=\"284\" copasiSourcesModified=\"0\">\n" >>rndinput2.cps
cat rndinput.cps >> rndinput2.cps

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t B rndinput2.cps 2 1> output 2> /dev/null

if ! grep -Pq "ERROR: rndinput2.cps failed to load." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 2"
fi

# now try as SBML
mv rndinput.cps rndinput.xml

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t B rndinput.xml 2 1> output 2> /dev/null

if ! grep -Pq "ERROR: rndinput.xml failed to load." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 4"
fi

# create a cps file that has a header but is followed by random input
printf "<?xml version='1.0' encoding='UTF-8' standalone='no'?>\n" > rndinput2.xml
printf "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\" xmlns:layout=\"http://www.sbml.org/sbml/level3/version1/layout/version1\" level=\"3\" metaid=\"e13bd125-750b-48d1-9687-1e31f49d9880\" xmlns:render=\"http://www.sbml.org/sbml/level3/version1/render/version1\" render:required=\"false\" version=\"1\" layout:required=\"false\">\n" >> rndinput2.xml
cat rndinput.xml >> rndinput2.xml

# run sbmodelr
$PYTH ../../src/sbmodelr.py -t B rndinput2.xml 2 1> output 2> /dev/null

if ! grep -Pq "ERROR: rndinput2.xml failed to load." output; then
  printf 'FAIL %s\n' "${test}"
  let "fail = $fail + 8"
fi

if [ "$fail" = 0 ] ; then
  printf 'PASS %s\n' "${test}"
  rm output *.cps *.xml
fi

exit $fail
