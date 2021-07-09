#!/bin/bash

setup_dir=<setup>

gmx grompp -f npt.mdp -c ../nvt/nvt.gro -t ../nvt/nvt.cpt -r ../nvt/nvt.gro -p ${setup_dir}/ligand.hmr.top \
	   -o npt.tpr -maxwarn 2 2>&1 | tee npt.grompp.log

gmx mdrun -v -deffnm npt -ntmpi 1 -ntomp 6 2>&1 | tee run.log

wait
