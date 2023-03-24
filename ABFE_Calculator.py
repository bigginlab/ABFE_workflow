#!/usr/bin/env python3

import glob, os
import argparse
import yaml

from abfe import calculate_abfe
import logging
loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)

if __name__ == "__main__":
        # ARGPARSE
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--protein_pdb_path", help='Input protein pdb file path', required=True, type=str)
    parser.add_argument('-l', "--ligand_mol_dir",  help='Input ligand(s) mol file path', required=True, type=str) 
    parser.add_argument('-o', "--output_dir_path",  help='Output approach folder', required=True, type=str)
    parser.add_argument('-c', "--cofactor_mol_path", help='Input cofactor(s) mol file path', required=False, default=None, type=str)
    parser.add_argument('-m', "--membrane_pdb_path",
                        help="Input membrane pdb file path. WARNING: The CRYST1 information of this PDB will be used for solvating the system."\
                            "The protein-membrane system MUST be aligned to the Z-axis (needed for pressure couple scheme)."\
                            "CHARMM-GUI is a good option to get this file.",
                        required=False, default=None, type=str)
    parser.add_argument('-hmr', "--hmr_factor",
                        help='The Hydrogen Mass Repartition factor to use. 4 fs of integration time step will be used no matter what hmf_factor is provided. '\
                            'Values greater than 2 are advised, if not the system may be unstable.',
                        required=False, default=3.0, type=float)
    parser.add_argument('-nr', "--number_of_replicates", help='Number of replicates', required=False, default=3, type=int)
    parser.add_argument('-njr', "--number_of_parallel_receptor_jobs", help='Number of jobs in parallel for receptor workflow', required=False, default=None,
                        type=int)
    parser.add_argument('-njl', "--number_of_parallel_ligand_jobs", help='Number of jobs in parallel for ligand workflow', required=False, default=40, type=int)
    parser.add_argument('-ncl', "--number_of_cpus_per_ligand_job", help='Number of cpus per ligand job', required=False, default=8, type=int)
    parser.add_argument('-sc', '--slrum_config',
                        help="This is the configuration YAML file of your Slrum cluster. If nothing is provided:\n\n"\
                        "   partition = cpu\n"\
                        "   time=60:00:00\n"\
                        "   mem=5000"
                        , required=False, default = None, type=str)
    parser.add_argument('-submit',help='Will automatically submit the ABFE calculations', required=False, action='store_true')


    args = parser.parse_args()
  
    mol_paths = glob.glob(os.path.join(args.ligand_mol_dir,"*.mol"))
    
    if args.slrum_config:
        with open(args.slrum_config, 'r') as c:
            cluster_config =  yaml.safe_load(c)
    else:
        cluster_config ={
            "partition": "cpu",
            "time": "60:00:00",
            "mem": "5000",
        }
    
    calculate_abfe(protein_pdb_path = args.protein_pdb_path,
                   ligand_sdf_paths = mol_paths,
                   out_root_folder_path=args.output_dir_path,
                   cofactor_mol_path = args.cofactor_mol_path,
                   membrane_pdb_path = args.membrane_pdb_path,
                   hmr_factor = args.hmr_factor,
                   n_cores_per_job = args.number_of_cpus_per_ligand_job,
                   num_jobs_per_ligand = args.number_of_parallel_ligand_jobs,
                   num_jobs_receptor_workflow =args.number_of_parallel_receptor_jobs,
                   num_replicas = args.number_of_replicates,
                   submit = args.submit,
                   cluster_config=cluster_config,
                   )
