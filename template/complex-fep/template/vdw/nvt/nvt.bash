#!/bin/bash

setup_dir=<setup>

gmx grompp -f nvt.mdp -c ../em/enmin.gro -r ../em/enmin.gro -p ${setup_dir}/complex.top \
	   -o nvt.tpr -maxwarn 2 2>&1 | tee nvt.grompp.log

gmx mdrun -v -deffnm nvt -ntmpi 1 -ntomp 6

wait
