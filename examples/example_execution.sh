#!/bin/bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

python ../ABFE_Calculator.py -p ./CyclophilinD_data/input/receptor.pdb \
                    -l ./CyclophilinD_data/input/ligands \
                    -o ./CyclophilinD_data/abfe  \
                    -gpu -hybrid -nc 8 -submit
