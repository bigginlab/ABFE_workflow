#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.

conda activate abfe

cli-abfe -p ./data/CyclophilinD_min/receptor.pdb \
         -l ./data/CyclophilinD_min/ligands \
         -nr 3 -njr 8 -njl 40
