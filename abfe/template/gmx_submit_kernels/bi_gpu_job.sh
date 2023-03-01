#!/bin/bash

#----------- Start of Edit -------------
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --job-name=Complex_smk

#SBATCH --cpus-per-task=8
#SBATCH --ntasks=8

### Edit Gromacs files to process:
OMP_NUM_THREADS="${1}" #dummy variable
STEPNAME="${2}"
TOP="${3}"
STRUCTURE="${4}"
GROMACS_TPR="./${STEPNAME}.tpr"
CONFOUT="./${STEPNAME}.gro"
#----------- End of Edit -------------


# define global environment variables
export PATH="/usr/local/cuda/bin:/apps/prod/COMPCHEM/compchem/linux/mvapich2-2.3.2_sl7_slurm/bin:/opt/slurm/prod/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/cuda/lib64:/apps/prod/COMPCHEM/compchem/linux/mvapich2-2.3.2_sl7_slurm/lib:/opt/slurm/prod/lib64:$LD_LIBRARY_PATH"
export CUDA_HOME="/usr/local/cuda"

# find number of threads for OpenMP
# find number of MPI tasks per node
echo "SLURM_TASKS_PER_NODE:" $SLURM_TASKS_PER_NODE
echo "SLURM_JOB_CPUS_PER_NODE:" $SLURM_JOB_CPUS_PER_NODE
TPN=`echo $SLURM_TASKS_PER_NODE | cut -f 1 -d \(`
# find number of CPU cores per node
PPN=`echo $SLURM_JOB_CPUS_PER_NODE | cut -f 1 -d \(`
THREADS=$(($PPN/$TPN))

OMP_NUM_THREADS=$THREADS
echo "OMP_NUM_THREADS: " $OMP_NUM_THREADS
echo ""
echo ""

# define gromacs environment
source /apps/prod/COMPCHEM/gromacs/gromacs-2021_sl7/bin/GMXRC.bash
source /apps/prod/COMPCHEM/gcc/gcc-7.4.0-source.bash

export MV2_ENABLE_AFFINITY=0
export OMP_SCHEDULE=dynamic
export KMP_BLOCKTIME=0
export MDRUN="mdrun_mpi -ntmpi $OMP_NUM_THREADS"

#Grompp this:
gmx grompp -f ./${STEPNAME}.mdp -c ${TOP} -r ${TOP} -p ${STRUCTURE} -o ./${GROMACS_TPR} -maxwarn 2

# run gromacs command
srun --mpi=pmi2 $MDRUN -s $GROMACS_TPR -c $CONFOUT -deffnm ./${STEPNAME} 

exit 0