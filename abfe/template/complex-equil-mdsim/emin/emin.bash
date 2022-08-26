#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J Equil
#SBATCH -p small

gmx grompp -f emin.mdp -c ../../final-topology/mcwat.gro \
	   -p ../../final-topology/fragment_complex.hmr.top -o emin.tpr -maxwarn 2 2>&1 | tee emin.grompp.log

gmx mdrun -v -deffnm emin -ntmpi 1 -ntomp 6 2>&1 | tee run.log

wait

