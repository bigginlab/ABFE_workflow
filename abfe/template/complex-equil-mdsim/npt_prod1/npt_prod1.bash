#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J Equil
#SBATCH -p small

gmx grompp -f npt_prod.mdp -c ../npt_equil2/npt_equil2.gro -t ../npt_equil2/npt_equil2.cpt \
	   -p ../../final-topology/fragment_complex.hmr.top -o npt_prod1.tpr -maxwarn 2 2>&1 | tee npt.grompp.log

gmx mdrun -v -deffnm npt_prod1 -ntomp 6 2>&1 | tee run.log

wait

