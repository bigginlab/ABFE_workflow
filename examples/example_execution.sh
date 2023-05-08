#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

ABFE_CLI -p ./CyclophilinD_data/minimal_input/receptor_protein.pdb \
          -l ./CyclophilinD_data/minimal_input/ligands \
          -o ./CyclophilinD_data/minimal_input/abfe \
          -nr 3 \
          -submit
