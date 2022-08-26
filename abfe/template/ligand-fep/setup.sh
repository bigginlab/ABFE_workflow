#!/bin/bash

wd=$PWD

prefix=../../../final-topology/

cp -r ../6_gmx-prep/ligand.hmr-gmx/final final-topology

mkdir -p simulation

# coul
windows=11
let max_window=${windows}-1
for i in $(seq 0 ${max_window})
do
  cp -r template/coul ./simulation/coul.${i}
  sed -i "s#<setup>#${prefix}#g" ./simulation/coul.${i}/*/*.bash
  sed -i "s#<state>#${i}#g" ./simulation/coul.${i}/*/*.mdp
  cp coul.sub ./simulation/
done

# vdw
windows=21
let max_window=${windows}-1
for i in $(seq 0 ${max_window})
do
  cp -r template/vdw ./simulation/vdw.${i}
  sed -i "s#<setup>#${prefix}#g" ./simulation/vdw.${i}/*/*.bash
  sed -i "s#<state>#${i}#g" ./simulation/vdw.${i}/*/*.mdp
  cp vdw.sub ./simulation/
done

mv simulation simulation-run1
for i in {2..5}
do
    cp -r simulation-run1 simulation-run${i}
done
