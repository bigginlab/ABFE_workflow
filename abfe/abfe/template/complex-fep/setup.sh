#!/bin/bash

wd=$PWD

prefix=../../../final-topology/
num_windows=$2

restraint_windows=12
vdw_windows=21
coul_windows=11

mkdir -p simulation

# restraints
let max_window=${restraint_windows}-1
for i in $(seq 0 ${max_window})
do
  cp -r template/restraints ./simulation/restraints.${i}
  sed -i "s#<setup>#${prefix}#g" ./simulation/restraints.${i}/*/*.bash
  sed -i "s#<state>#${i}#g" ./simulation/restraints.${i}/*/*.mdp
  cp restraint.sub ./simulation/
  cp rest.pre-prod.sub ./simulation/
done

#vdw
let max_window=${vdw_windows}-1
for i in $(seq 0 ${max_window})
do
  cp -r template/vdw ./simulation/vdw.${i}
  sed -i "s#<setup>#${prefix}#g" ./simulation/vdw.${i}/*/*.bash
  sed -i "s#<state>#${i}#g" ./simulation/vdw.${i}/*/*.mdp
  cp vdw.sub ./simulation/ 
  cp vdw.pre-prod.sub ./simulation/
done

#coul
let max_window=${coul_windows}-1
for i in $(seq 0 ${max_window})
do
  cp -r template/coul ./simulation/coul.${i}
  sed -i "s#<setup>#${prefix}#g" ./simulation/coul.${i}/*/*.bash
  sed -i "s#<state>#${i}#g" ./simulation/coul.${i}/*/*.mdp
  cp coul.sub ./simulation/
  cp coul.pre-prod.sub ./simulation/
done

# Create three replica

mv simulation simulation-run1
