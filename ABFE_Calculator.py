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
    parser.add_argument('-nr', "--number_of_replicates", help='Number of replicates', required=False, default=3, type=int)
    parser.add_argument('-njr', "--number_of_parallel_receptor_jobs", help='Number of jobs in parallel for receptor workflow', required=False, default=None,
                        type=int)
    parser.add_argument('-njl', "--number_of_parallel_ligand_jobs", help='Number of jobs in parallel for ligand workflow', required=False, default=40, type=int)
    parser.add_argument('-ncl', "--number_of_cpus_per_ligand_job", help='Number of cpus per ligand job', required=False, default=8, type=int)
    parser.add_argument('-submit',help='Will automatically submit the ABFE calculations', required=False, action='store_true')


    args = parser.parse_args()
  
    sdf_paths = glob.glob(args.ligand_sdf_dir+"/*.sdf")
    

    cluster_config ={
        "partition": "cpu",
        "time": "60:00:00",
        "mem": "5000",
    }
    
    calculate_abfe(protein_pdb_path=args.protein_pdb_path, ligand_sdf_paths=sdf_paths, out_root_folder_path=args.output_dir_path, cofactor_sdf_path=args.cofactor_sdf_path,
                   n_cores_per_job = args.number_of_cpus_per_ligand_job, num_jobs_per_ligand= args.number_of_parallel_ligand_jobs,
                   num_jobs_receptor_workflow=args.number_of_parallel_receptor_jobs,
                   num_replicas=args.number_of_replicates,
                   submit=args.submit,
                   cluster_config=cluster_config)
