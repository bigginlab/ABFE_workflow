#!/usr/bin/env bash
# This script is an example execution running the example CyclophilinD abfe simulations.


cli-abfe-gmx -d  ./data/HSP90_gmx \
             -o abfe_HSP90_out \
             -pn HSP90_gmx \
		         -njr 30 -nr 3
