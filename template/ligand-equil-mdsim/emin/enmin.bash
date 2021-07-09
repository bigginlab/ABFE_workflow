#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J Equil
#SBATCH -p small

module load use.own
module load apps/gromacs2019.1

gmx grompp -f enmin.mdp -c ../../final-topology/ligand.gaff2.cubic.150mM.hmr.gro \
	   -p ../../final-topology/ligand.hmr.top -o enmin.tpr -maxwarn 2 2>&1 | tee enmin.grompp.log

gmx mdrun -v -deffnm enmin -ntmpi 1 -ntomp 12 -pin on 2>&1 | tee run.log

wait

