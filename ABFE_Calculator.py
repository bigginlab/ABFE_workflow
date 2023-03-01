#!/usr/bin/env python3

import glob
import argparse

from abfe import calculate_abfe
import logging
loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)

if __name__ == "__main__":
        # ARGPARSE
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--protein_pdb_path", help='Input protein pdb file path', required=True, type=str)
    parser.add_argument('-l', "--ligand_sdf_dir",  help='Input ligand(s) sdf file path', required=True, type=str) 
    parser.add_argument('-o', "--output_dir_path",  help='Output approach folder', required=True, type=str)
    parser.add_argument('-c', "--cofactor_sdf_path", help='Input cofactor(s) sdf file path', required=False, default=None, type=str)
    parser.add_argument('-nc', "--number_of_cpus_per_job", help='Number of cpus per job', required=False, default=None, type=int)
    parser.add_argument('-nj', "--number_of_parallel_jobs", help='Number of jobs in parallel', required=False, default=40, type=int)
    parser.add_argument('-nr', "--number_of_replicates", help='Number of replicates', required=False, default=3, type=int)
    parser.add_argument('-submit',help='Will automatically submit the ABFE calculations', required=False, action='store_true')
    
    #Todo:
    # dont creat / 
    args = parser.parse_args()
  
    sdf_paths = glob.glob(args.ligand_sdf_dir+"/*.sdf")
    cluster_config ={
        "partition": "cpu",
        "time": "60:00:00",
        "mem": "5000",
    }

    calculate_abfe(protein_pdb_path=args.protein_pdb_path, ligand_sdf_paths=sdf_paths, out_root_folder_path=args.output_dir_path, cofactor_sdf_path=args.cofactor_sdf_path, 
                n_cores_per_job = args.number_of_cpus_per_job, num_jobs = args.number_of_parallel_jobs, num_replicas=args.number_of_replicates,
                submit=args.submit, use_gpu=False, hybrid_job=False,
                cluster_config=cluster_config)
