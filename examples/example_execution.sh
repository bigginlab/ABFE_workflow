#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

ABFE_Calculator -p ./CyclophilinD_data/input/receptor.pdb \
                    -l ./CyclophilinD_data/input/ligands \
                    -o ./CyclophilinD_data/abfe \
                    -nr 3 \
                    -submit
