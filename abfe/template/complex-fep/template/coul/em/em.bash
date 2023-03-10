#!/bin/bash

setup_dir=<setup>

gmx grompp -f em.mdp -c ${setup_dir}/complex.gro -p ${setup_dir}/complex.top -o emin.tpr -maxwarn 2 2>&1 | tee emin.grompp.log

gmx mdrun -v -deffnm emin -ntmpi 1 -ntomp 6

wait
