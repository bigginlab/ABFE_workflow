#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

cli-abfe-gmx -d  ./CyclophilinD_data/gromacs_in \
                    -o /scratch/riesbenj/ABFE/CyclophilinD_gmx \
                    -pn CyclophilinD_gmx \
                    -nr 3 \
                    -njr 8 \
                    -njl 40 
