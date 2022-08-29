#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=144:00:00
#SBATCH -J Equil
#SBATCH -p small

gmx grompp -f npt_equil1.mdp -c ../nvt_heat/nvt_heat.gro -t ../nvt_heat/nvt_heat.cpt \
           -r ../nvt_heat/nvt_heat.gro -p ../../final-topology/fragment_complex.hmr.top \
	   -o npt_equil1.tpr -maxwarn 2 2>&1 | tee npt.grompp.log

gmx mdrun -v -deffnm npt_equil1 -ntomp 6 2>&1 | tee run.log

wait

