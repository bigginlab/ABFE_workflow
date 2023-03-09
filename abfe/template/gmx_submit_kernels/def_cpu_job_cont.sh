#!/bin/bash
# This is a base shell template file

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
gmx mdrun -ntomp ${OMP_NUM_THREADS} -ntmpi 1 -s $GROMACS_TPR -c $CONFOUT -deffnm ${STEPNAME}

exit 0