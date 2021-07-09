#!/bin/bash

setup_dir=<setup>

gmx grompp -f em.mdp -c ${setup_dir}/ligand.equil.hmr.gro -p ${setup_dir}/ligand.hmr.top -o enmin.tpr -maxwarn 2 2>&1 | tee emin.grompp.log

gmx mdrun -v -deffnm enmin -ntmpi 1 -ntomp 6 2>&1 | tee run.log

wait
