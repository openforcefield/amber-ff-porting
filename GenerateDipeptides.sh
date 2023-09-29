#!/usr/bin/env bash

set -x
set -euo pipefail

ALLRES="ALA ARG ASH ASN ASP CYS GLH GLN GLU GLY HID HIE HIP ILE LEU LYN LYS MET PHE PRO SER"
ALLRES="${ALLRES} THR TRP TYR VAL CYX CYM"

TRMRES="ALA ARG ASN ASP CYS GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER"
TRMRES="${TRMRES} THR TRP TYR VAL CYX CYM"

# Systems containing CYX-CYX bonds will need minimization in order to correct long S-S bonds
cat > min.in << EOF
Initial minimization
 &cntrl
  ntx    = 1,       irest  = 0,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 100,     ntwx   = 0,       ntwv   = 0,      ntwe   = 0,

  ntf    = 1,       ntb    = 1,
  ntc    = 1,

  cut    = 14.0,

  imin   = 1,
  maxcyc = 1000,
  ncyc   = 100,
 &end
 &ewald
  nfft1 = 32, nfft2 = 32, nfft3 = 32,
 &end
EOF

# Make main-chain tripeptides
mkdir -p MainChain
cd MainChain
for RESA in ${ALLRES} ; do
  if [ -d ${RESA} ] ; then
    continue
  fi
  echo "${RESA}"
  mkdir -p ${RESA}
  cd ${RESA}/
  echo "source leaprc.protein.ff14SB" > tleap.in
  if [ ${RESA} == "CYX" ] ; then
    echo "x = sequence { ACE ${RESA} NME ACE CYX NME }" >> tleap.in
    echo "bond x.2.SG x.5.SG" >> tleap.in
  else
    echo "x = sequence { ACE ${RESA} NME }" >> tleap.in
  fi
  echo "set x box { 48.0 48.0 48.0 }" >> tleap.in
  echo "saveAmberParm x ${RESA}.prmtop ${RESA}.inpcrd" >> tleap.in
  echo "quit" >> tleap.in
  tleap -f tleap.in > tleap.out
  sander -O -i ../../min.in -o min.out -p ${RESA}.prmtop -c ${RESA}.inpcrd -r mincrd
  mv mincrd ${RESA}.inpcrd
  ambpdb -p ${RESA}.prmtop < ${RESA}.inpcrd > ${RESA}.pdb
  echo "source leaprc.protein.ff14SB" > tleap2.in
  echo "x = loadPdb \"${RESA}.pdb\"" >> tleap2.in
  if [ ${RESA} == "CYX" ] ; then
    echo "bond x.2.SG x.5.SG" >> tleap2.in
  fi
  echo "set x box { 48.0 48.0 48.0 }" >> tleap2.in
  echo "saveMol2 x ${RESA}.mol2 1" >> tleap2.in
  echo "quit" >> tleap2.in
  tleap -f tleap2.in > tleap2.out
  antechamber -i ${RESA}.mol2 -fi mol2 -o ${RESA}.mol2 -fo mol2 -at sybyl -dr no > ac.out
  cd ../
done
cd ../

# Make N-terminal tripeptides
mkdir -p NTerminal
cd NTerminal
for RESA in ${TRMRES} ; do
  if [ -d ${RESA} ] ; then
    continue
  fi
  echo "${RESA}"
  mkdir -p ${RESA}
  cd ${RESA}/
  echo "source leaprc.protein.ff14SB" > tleap.in
  if [ ${RESA} == "CYX" ] ; then
    echo "x = sequence { N${RESA} NME ACE CYX NME }" >> tleap.in
    echo "bond x.1.SG x.4.SG" >> tleap.in
  else
    echo "x = sequence { N${RESA} NME }" >> tleap.in
  fi
  echo "set x box { 48.0 48.0 48.0 }" >> tleap.in
  echo "saveAmberParm x ${RESA}.prmtop ${RESA}.inpcrd" >> tleap.in
  echo "quit" >> tleap.in
  tleap -f tleap.in > tleap.out
  sander -O -i ../../min.in -o min.out -p ${RESA}.prmtop -c ${RESA}.inpcrd -r mincrd
  mv mincrd ${RESA}.inpcrd
  ambpdb -p ${RESA}.prmtop < ${RESA}.inpcrd > ${RESA}.pdb
  echo "source leaprc.protein.ff14SB" > tleap2.in
  echo "x = loadPdb \"${RESA}.pdb\"" >> tleap2.in
  if [ ${RESA} == "CYX" ] ; then
    echo "bond x.1.SG x.4.SG" >> tleap2.in
  fi
  echo "set x box { 48.0 48.0 48.0 }" >> tleap2.in
  echo "saveMol2 x ${RESA}.mol2 1" >> tleap2.in
  echo "quit" >> tleap2.in
  tleap -f tleap2.in > tleap2.out
  antechamber -i ${RESA}.mol2 -fi mol2 -o ${RESA}.mol2 -fo mol2 -at sybyl -dr no > ac.out
  cd ../
done
cd ../

# Make C-terminal tripeptides
mkdir -p CTerminal
cd CTerminal
for RESA in ${TRMRES} ; do
  if [ -d ${RESA} ] ; then
    continue
  fi
  echo "${RESA}"
  mkdir -p ${RESA}
  cd ${RESA}/
  echo "source leaprc.protein.ff14SB" > tleap.in
  if [ ${RESA} == "CYX" ] ; then
    echo "x = sequence { ACE C${RESA} ACE CYX NME }" >> tleap.in
    echo "bond x.2.SG x.4.SG" >> tleap.in
  else 
    echo "x = sequence { ACE C${RESA} }" >> tleap.in
  fi
  echo "set x box { 48.0 48.0 48.0 }" >> tleap.in
  echo "saveAmberParm x ${RESA}.prmtop ${RESA}.inpcrd" >> tleap.in
  echo "quit" >> tleap.in
  tleap -f tleap.in > tleap.out
  sander -O -i ../../min.in -o min.out -p ${RESA}.prmtop -c ${RESA}.inpcrd -r mincrd
  mv mincrd ${RESA}.inpcrd
  ambpdb -p ${RESA}.prmtop < ${RESA}.inpcrd > ${RESA}.pdb
  echo "source leaprc.protein.ff14SB" > tleap2.in
  echo "x = loadPdb \"${RESA}.pdb\"" >> tleap2.in
  if [ ${RESA} == "CYX" ] ; then
    echo "bond x.2.SG x.4.SG" >> tleap2.in
  fi
  echo "set x box { 48.0 48.0 48.0 }" >> tleap2.in
  echo "saveMol2 x ${RESA}.mol2 1" >> tleap2.in
  echo "quit" >> tleap2.in
  tleap -f tleap2.in > tleap2.out
  antechamber -i ${RESA}.mol2 -fi mol2 -o ${RESA}.mol2 -fo mol2 -at sybyl -dr no > ac.out
  cd ../
done
cd ../

