#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J Equil
#SBATCH -p small

gmx grompp -f nvt_heat.mdp -c ../emin/emin.gro -r ../emin/emin.gro \
	   -p ../../final-topology/fragment_complex.hmr.top -o nvt_heat.tpr -maxwarn 2 2>&1 | tee nvt.grompp.log

gmx mdrun -v -deffnm nvt_heat -ntomp 6 2>&1 | tee run.log

wait

