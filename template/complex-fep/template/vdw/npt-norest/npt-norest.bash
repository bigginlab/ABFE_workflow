#!/bin/bash

setup_dir=<setup>

gmx grompp -f npt-norest.mdp -c ../npt/npt.gro -t ../npt/npt.cpt -p ${setup_dir}/complex.top \
	   -o npt-norest.tpr -maxwarn 2 2>&1 | tee npt.grompp.log

gmx mdrun -v -deffnm npt-norest -ntmpi 1 -ntomp 6

wait
