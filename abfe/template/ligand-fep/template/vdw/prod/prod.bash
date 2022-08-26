#!/bin/bash

setup_dir=<setup>

gmx grompp -f prod.mdp -c ../npt-norest/npt-norest.gro -t ../npt-norest/npt-norest.cpt -p ${setup_dir}/ligand.hmr.top \
	   -o prod.tpr -maxwarn 2 2>&1 | tee prod.grompp.log

gmx mdrun -v -deffnm prod -ntmpi 1 -ntomp 6 2>&1 | tee run.log

wait
