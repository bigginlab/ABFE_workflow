#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

cli-abfe-gmx -d  ./CyclophilinD_data/gromacs_in \
                    -o ./abfe_gmx \
                    -pn CyclophilinD_gmx \
                    -nr 3
