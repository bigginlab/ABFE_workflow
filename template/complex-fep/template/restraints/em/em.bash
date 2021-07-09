#!/bin/bash

setup_dir=<setup>

gmx grompp -f em.mdp -c ${setup_dir}/complex.gro -p ${setup_dir}/complex.top -o enmin.tpr -maxwarn 2 2>&1 | tee emin.grompp.log

gmx mdrun -v -deffnm enmin -ntmpi 1 -ntomp 6

wait
