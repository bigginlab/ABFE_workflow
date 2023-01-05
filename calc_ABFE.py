#!/usr/bin/env python3

import os
from abfe.orchestration.build_and_run_ligands import calculate_all_ligands
   
if __name__ == "__main__":
    orig_dir = os.getcwd()
    
    # IO:
    out_root_path = os.path.abspath("./data/system1")
    in_root_path = os.path.abspath("./data/input/system1")

    input_ligand_paths = [in_root_path+"/"+dir for dir in os.listdir(in_root_path) if(os.path.isdir(in_root_path+"/"+dir))]
    print("input ligand dirs: ", input_ligand_paths)
    print("output root dir: ", out_root_path)

    # Options:
    n_cores_per_job = 8
    num_jobs = 40 # spawn n jobs in parallel
    num_replicas=1 # how often shall the calculation be repeated?
    submit=False  # directly submit to the slurm queue

    
    cluster_config ={
        "partition": "cpu",
        "time": "48:00:00",
        "mem": "20GB",
    }
    
    # Do Fun!
    if(not os.path.isdir(out_root_path)): os.mkdir(out_root_path)
    calculate_all_ligands(input_ligand_paths=input_ligand_paths, out_root_path=out_root_path,  num_max_thread = n_cores_per_job,
                           num_replicas=num_replicas,  submit=submit, num_jobs=num_jobs, cluster_config=cluster_config)

    os.chdir(orig_dir)
