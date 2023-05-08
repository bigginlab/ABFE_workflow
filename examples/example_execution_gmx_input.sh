#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

ABFE_GMX_CLI -d  ./CyclophilinD_data/gromacs_in/input \
                    -o ./CyclophilinD_data/gromacs_in/abfe \
                    -nr 3 \
                    -submit
