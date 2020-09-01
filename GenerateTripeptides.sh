#!/bin/bash

ALLRES="ALA ARG ASH ASN ASP CYS GLH GLN GLU GLY HID HIE HIP ILE LEU LYN LYS MET PHE PRO SER"
ALLRES="${ALLRES} THR TRP TYR VAL CYX"

TRMRES="ALA ARG ASN ASP CYS GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER"
TRMRES="${TRMRES} THR TRP TYR VAL CYX"

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
  for RESB in ${ALLRES} ; do
    if [ -d ${RESA}_${RESB} ] ; then
      continue
    elif [ ${RESA} == "CYX" ] && [ ${RESB} == "CYX" ] ; then
      continue
    fi
    echo "${RESA} ${RESB}"
    mkdir -p ${RESA}_${RESB}
    cd ${RESA}_${RESB}/
    echo "source leaprc.protein.ff14SB" > tleap.in
    if [ ${RESA} == "CYX" ] || [ ${RESB} == "CYX" ] ; then
      echo "x = sequence { ACE ${RESA} ${RESB} NME ACE CYX NME }" >> tleap.in
    else
      echo "x = sequence { ACE ${RESA} ${RESB} NME }" >> tleap.in
    fi
    if [ ${RESA} == "CYX" ] ; then
      echo "bond x.2.SG x.6.SG" >> tleap.in
    elif [ ${RESB} == "CYX" ] ; then
      echo "bond x.3.SG x.6.SG" >> tleap.in
    fi
    echo "set x box { 48.0 48.0 48.0 }" >> tleap.in
    echo "saveAmberParm x ${RESA}_${RESB}.prmtop ${RESA}_${RESB}.inpcrd" >> tleap.in
    echo "quit" >> tleap.in
    tleap -f tleap.in > tleap.out
    sander -O -i ../../min.in -o min.out -p ${RESA}_${RESB}.prmtop -c ${RESA}_${RESB}.inpcrd \
          -r mincrd
    mv mincrd ${RESA}_${RESB}.inpcrd
    ambpdb -p ${RESA}_${RESB}.prmtop < ${RESA}_${RESB}.inpcrd > ${RESA}_${RESB}.pdb
    echo "source leaprc.protein.ff14SB" > tleap2.in
    echo "x = loadPdb \"${RESA}_${RESB}.pdb\"" >> tleap2.in
    if [ ${RESA} == "CYX" ] ; then
      echo "bond x.2.SG x.6.SG" >> tleap2.in
    elif [ ${RESB} == "CYX" ] ; then
      echo "bond x.3.SG x.6.SG" >> tleap2.in
    fi
    echo "set x box { 48.0 48.0 48.0 }" >> tleap2.in
    echo "saveMol2 x ${RESA}_${RESB}.mol2 1" >> tleap2.in
    echo "quit" >> tleap2.in
    tleap -f tleap2.in > tleap2.out
    antechamber -i ${RESA}_${RESB}.mol2 -fi mol2 -o ${RESA}_${RESB}.mol2 -fo mol2 -at sybyl \
		-dr no > ac.out
    cd ../
  done
done
cd ../

# Make N-terminal tripeptides
mkdir -p NTerminal
cd NTerminal
for RESA in ${TRMRES} ; do
  for RESB in ${ALLRES} ; do
    if [ -d ${RESA}_${RESB} ] ; then
      continue
    elif [ ${RESA} == "CYX" ] && [ ${RESB} == "CYX" ] ; then
      continue
    fi
    echo "${RESA} ${RESB}"
    mkdir -p ${RESA}_${RESB}
    cd ${RESA}_${RESB}/
    echo "source leaprc.protein.ff14SB" > tleap.in
    if [ ${RESA} == "CYX" ] || [ ${RESB} == "CYX" ] ; then
      echo "x = sequence { N${RESA} ${RESB} NME ACE CYX NME }" >> tleap.in
    else
      echo "x = sequence { N${RESA} ${RESB} NME }" >> tleap.in
    fi
    if [ ${RESA} == "CYX" ] ; then
      echo "bond x.1.SG x.5.SG" >> tleap.in
    elif [ ${RESB} == "CYX" ] ; then
      echo "bond x.2.SG x.5.SG" >> tleap.in
    fi
    echo "set x box { 48.0 48.0 48.0 }" >> tleap.in
    echo "saveAmberParm x ${RESA}_${RESB}.prmtop ${RESA}_${RESB}.inpcrd" >> tleap.in
    echo "quit" >> tleap.in
    tleap -f tleap.in > tleap.out
    sander -O -i ../../min.in -o min.out -p ${RESA}_${RESB}.prmtop -c ${RESA}_${RESB}.inpcrd \
          -r mincrd
    mv mincrd ${RESA}_${RESB}.inpcrd
    ambpdb -p ${RESA}_${RESB}.prmtop < ${RESA}_${RESB}.inpcrd > ${RESA}_${RESB}.pdb
    echo "source leaprc.protein.ff14SB" > tleap2.in
    echo "x = loadPdb \"${RESA}_${RESB}.pdb\"" >> tleap2.in
    if [ ${RESA} == "CYX" ] ; then
      echo "bond x.1.SG x.5.SG" >> tleap2.in
    elif [ ${RESB} == "CYX" ] ; then
      echo "bond x.2.SG x.5.SG" >> tleap2.in
    fi
    echo "set x box { 48.0 48.0 48.0 }" >> tleap2.in
    echo "saveMol2 x ${RESA}_${RESB}.mol2 1" >> tleap2.in
    echo "quit" >> tleap2.in
    tleap -f tleap2.in > tleap2.out
    antechamber -i ${RESA}_${RESB}.mol2 -fi mol2 -o ${RESA}_${RESB}.mol2 -fo mol2 -at sybyl \
                -dr no > ac.out
    cd ../
  done
done
cd ../

# Make C-terminal tripeptides
mkdir -p CTerminal
cd CTerminal
for RESA in ${ALLRES} ; do
  for RESB in ${TRMRES} ; do
    if [ -d ${RESA}_${RESB} ] ; then
      continue
    elif [ ${RESA} == "CYX" ] && [ ${RESB} == "CYX" ] ; then
      continue
    fi
    echo "${RESA} ${RESB}"
    mkdir -p ${RESA}_${RESB}
    cd ${RESA}_${RESB}/
    echo "source leaprc.protein.ff14SB" > tleap.in
    if [ ${RESA} == "CYX" ] || [ ${RESB} == "CYX" ] ; then
      echo "x = sequence { ACE ${RESA} C${RESB} ACE CYX NME }" >> tleap.in
    else
      echo "x = sequence { ACE ${RESA} C${RESB} }" >> tleap.in
    fi
    if [ ${RESA} == "CYX" ] ; then
      echo "bond x.2.SG x.5.SG" >> tleap.in
    elif [ ${RESB} == "CYX" ] ; then
      echo "bond x.3.SG x.5.SG" >> tleap.in
    fi
    echo "set x box { 48.0 48.0 48.0 }" >> tleap.in
    echo "saveAmberParm x ${RESA}_${RESB}.prmtop ${RESA}_${RESB}.inpcrd" >> tleap.in
    echo "quit" >> tleap.in
    tleap -f tleap.in > tleap.out
    sander -O -i ../../min.in -o min.out -p ${RESA}_${RESB}.prmtop -c ${RESA}_${RESB}.inpcrd \
          -r mincrd
    mv mincrd ${RESA}_${RESB}.inpcrd
    ambpdb -p ${RESA}_${RESB}.prmtop < ${RESA}_${RESB}.inpcrd > ${RESA}_${RESB}.pdb
    echo "source leaprc.protein.ff14SB" > tleap2.in
    echo "x = loadPdb \"${RESA}_${RESB}.pdb\"" >> tleap2.in
    if [ ${RESA} == "CYX" ] ; then
      echo "bond x.2.SG x.5.SG" >> tleap2.in
    elif [ ${RESB} == "CYX" ] ; then
      echo "bond x.3.SG x.5.SG" >> tleap2.in
    fi
    echo "set x box { 48.0 48.0 48.0 }" >> tleap2.in
    echo "saveMol2 x ${RESA}_${RESB}.mol2 1" >> tleap2.in
    echo "quit" >> tleap2.in
    tleap -f tleap2.in > tleap2.out
    antechamber -i ${RESA}_${RESB}.mol2 -fi mol2 -o ${RESA}_${RESB}.mol2 -fo mol2 -at sybyl \
                -dr no > ac.out
    cd ../
  done
done
cd ../

