#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

cli-abfe -p ./CyclophilinD_data/minimal_input/receptor_protein.pdb \
         -l ./CyclophilinD_data/minimal_input/ligands \
         -nr 3 -njr 8 -njl 40
