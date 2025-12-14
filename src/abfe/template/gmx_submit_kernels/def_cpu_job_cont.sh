#!/bin/bash
# This is a base shell template file

module add GROMACS/2023.1-foss-2022a-CUDA-11.7.0

export OMP_NUM_THREADS="${1}"
STEPNAME="${2}"
TOPOLOGY="${3}"
STRUCTURE="${4}"
CPT="${5}"

#opts
GROMACS_TPR="${STEPNAME}.tpr"
CONFOUT="${STEPNAME}.gro"

#Grompp this:
gmx grompp -f ./${STEPNAME}.mdp -c ${STRUCTURE} -r ${STRUCTURE} -p ${TOPOLOGY} -t ${CPT} -o ${GROMACS_TPR} -maxwarn 3

# run gromacs command
gmx mdrun -ntomp ${OMP_NUM_THREADS} -s $GROMACS_TPR -c $CONFOUT -deffnm ${STEPNAME} -ntmpi 1

exit 0