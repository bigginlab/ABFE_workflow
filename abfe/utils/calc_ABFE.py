#!/usr/bin/env python3

import os
from abfe.orchestration.build_and_run_ligands import calculate_all_ligands
   
if __name__ == "__main__":
    orig_dir = os.getcwd()
    
    # IO:
    out_root_path = "/path/to/data/MCL-1"
    in_root_path = "/path/to/data/input/MCL-1"

    input_ligand_paths = [in_root_path+"/"+dir for dir in os.listdir(in_root_path)]
    print("input ligand dirs: ", input_ligand_paths)
    print("output root dir: ", out_root_path)

    # Options:
    n_cores=1
    num_jobs = 10
    num_replicas=1
    submit=False

    # Do Fun! 
    if(not os.path.isdir(out_root_path)): os.mkdir(out_root_path)
    calculate_all_ligands(input_ligand_paths=input_ligand_paths, out_root_path=out_root_path, 
                           num_replicas=num_replicas, n_cores=n_cores, submit=submit, num_jobs=num_jobs)

    os.chdir(orig_dir)