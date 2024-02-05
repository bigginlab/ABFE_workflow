#!/usr/bin/env python3

import glob
import argparse

from abfe import calculate_abfe
from abfe.template import default_slurm_config_path

import logging
loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)


def main():
    # ARGPARSE
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--protein_pdb_path", help='Input protein pdb file path', required=True, type=str)
    parser.add_argument('-l', "--ligand_sdf_dir", help='Input ligand(s) sdf file path', required=True, type=str)
    parser.add_argument('-o', "--output_dir_path", help='Output approach folder', required=True, type=str)
    parser.add_argument('-c', "--cofactor_sdf_path", help='Input cofactor(s) sdf file path', required=False, default=None, type=str)
    parser.add_argument('-pn', "--project_name", help='name prefix of jobs, etc.', required=False, type=str, default="ABFE")
    parser.add_argument('-nr', "--number_of_replicates", help='Number of replicates', required=False, default=3, type=int)
    parser.add_argument('-njr', "--number_of_parallel_receptor_jobs", help='Number of jobs in parallel for receptor workflow', required=False, default=None,
                        type=int)
    parser.add_argument('-njl', "--number_of_parallel_ligand_jobs", help='Number of jobs in parallel for ligand workflow', required=False, default=40, type=int)
    parser.add_argument('-ncl', "--number_of_cpus_per_ligand_job", help='Number of cpus per ligand job', required=False, default=8, type=int)
    parser.add_argument('-sff', "--small_mol_ff", help='Force Field used for small mols', required=False, default="gaff", type=str)
    parser.add_argument('-nosubmit', help='Will automatically submit the ABFE calculations', required=False, action='store_false')
    parser.add_argument('-nogpu', help='shall gpus be used for the submissions? WARNING: Currently Not working', required=False, action='store_true')
    parser.add_argument('-nohybrid', help='hybrid flag executes complex jobs on gpu and ligand jobs on cpu (requires gpu flag) WARNING: Currently Not working',
                        required=False,
                        action='store_true')

    args = parser.parse_args()

    sdf_paths = glob.glob(args.ligand_sdf_dir + "/*.sdf")

    if (args.nosubmit):
        if (args.nogpu):
            cluster_config = json.load(open(f"{default_slurm_config_path}", "r"))
            cluster_config["Snakemake_job"]["queue_job_options"]["cpus-per-task"] = int(args.number_of_parallel_receptor_jobs)
            cluster_config["Sub_job"]["queue_job_options"]["cpus-per-task"] = int(args.number_of_cpus_per_ligand_job)

        else:
            cluster_config = json.load(open(f"{default_slurm_config_path}", "r"))
            cluster_config["Snakemake_job"]["queue_job_options"]["cpus-per-task"] = int(args.number_of_parallel_receptor_jobs)
            cluster_config["Sub_job"]["queue_job_options"]["cpus-per-task"] = int(args.number_of_parallel_ligand_jobs)
            cluster_config["Sub_job"]["queue_job_options"]["partition"] = "gpu"

    else:
        cluster_config = None

    res = calculate_abfe(protein_pdb_path=args.protein_pdb_path, ligand_sdf_paths=sdf_paths, out_root_folder_path=args.output_dir_path,
                         cofactor_sdf_path=args.cofactor_sdf_path, approach_name=args.project_name,
                         n_cores_per_job=args.number_of_cpus_per_ligand_job, num_jobs_per_ligand=args.number_of_parallel_ligand_jobs,
                         num_jobs_receptor_workflow=args.number_of_parallel_receptor_jobs,
                         num_replicas=args.number_of_replicates,
                         submit=args.nosubmit, use_gpu=not bool(args.nogpu), hybrid_job=not bool(args.nohybrid), small_mol_ff=args.small_mol_ff,
                         cluster_config=cluster_config)

    return res


if __name__ == "__main__":
    main()
