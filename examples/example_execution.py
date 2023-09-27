#!/usr/bin/env python3

"""
This script is an example execution running the example CyclophilinD abfe simulations.
"""
import os
import glob
from abfe import calculate_abfe

root_path  = os.path.dirname(__file__)+"/CyclophilinD_data"

in_root_path = root_path+"/minimal_input"
ligand_sdfs = glob.glob(in_root_path+"/ligands/*sdf")
receptor_pdb = in_root_path+"/receptor.pdb"

out_folder = in_root_path+"/abfe_py"

calculate_abfe(protein_pdb_path=receptor_pdb, 
               ligand_sdf_paths=ligand_sdfs,
               out_root_folder_path=out_folder,
               num_replicas=3,
               submit=True
               )
