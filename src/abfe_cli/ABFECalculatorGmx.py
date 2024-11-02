#!/usr/bin/env python3

import json
import argparse

from abfe.calculate_abfe_gmx import calculate_abfe_gmx
from abfe.template import default_slurm_config_path

import logging

loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)


def main():
    # ARGPARSE
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--gmx_files_root_dir",
        help="Input folder containing all the gmx files in the structure <ligand_name>/[solvent, "
        "complex]/[*.gro, *.top]",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_dir_path",
        help="Output approach folder",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-c",
        "--cofactor_sdf_path",
        help="Input cofactor(s) sdf file path",
        required=False,
        default=None,
        type=str,
    )
    parser.add_argument(
        "-pn",
        "--project_name",
        help="name prefix of jobs, etc.",
        required=False,
        type=str,
        default="ABFE",
    )
    parser.add_argument(
        "-nr",
        "--number_of_replicates",
        help="Number of replicates",
        required=False,
        default=3,
        type=int,
    )
    parser.add_argument(
        "-njr",
        "--number_of_parallel_receptor_jobs",
        help="Number of jobs in parallel for receptor workflow",
        required=False,
        default=None,
        type=int,
    )
    parser.add_argument(
        "-njl",
        "--number_of_parallel_ligand_jobs",
        help="Number of jobs in parallel for ligand workflow",
        required=False,
        default=40,
        type=int,
    )
    parser.add_argument(
        "-ncl",
        "--number_of_cpus_per_ligand_job",
        help="Number of cpus per ligand job",
        required=False,
        default=8,
        type=int,
    )
    parser.add_argument(
        "-nosubmit",
        help="Will automatically submit the ABFE calculations",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-nogpu",
        help="shall gpus be used for the submissions? WARNING: Currently Not working",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-nohybrid",
        help="hybrid flag executes complex jobs on gpu and ligand jobs on cpu (requires gpu flag) WARNING: Currently Not working",
        required=False,
        action="store_true",
    )

    args = parser.parse_args()

    if args.nosubmit:
        cluster_config = None
    else:
        if args.nogpu:
            cluster_config = json.load(open(f"{default_slurm_config_path}", "r"))
            cluster_config["Snakemake_job"]["queue_job_options"]["cpus-per-task"] = int(
                args.number_of_parallel_receptor_jobs
            )
            cluster_config["Sub_job"]["queue_job_options"]["cpus-per-task"] = int(
                args.number_of_parallel_ligand_jobs
            )

        else:
            cluster_config = json.load(open(f"{default_slurm_config_path}", "r"))
            cluster_config["Snakemake_job"]["queue_job_options"]["cpus-per-task"] = int(
                args.number_of_parallel_receptor_jobs
            )
            cluster_config["Sub_job"]["queue_job_options"]["cpus-per-task"] = int(
                args.number_of_parallel_ligand_jobs
            )
            cluster_config["Sub_job"]["queue_job_options"]["partition"] = "gpu"

    print(args.gmx_files_root_dir)

    res = calculate_abfe_gmx(
        input_dir=args.gmx_files_root_dir,
        out_root_folder_path=args.output_dir_path,
        approach_name=args.project_name,
        n_cores_per_job=args.number_of_cpus_per_ligand_job,
        num_jobs_per_ligand=args.number_of_parallel_ligand_jobs,
        num_jobs_receptor_workflow=args.number_of_parallel_receptor_jobs,
        num_replicas=args.number_of_replicates,
        submit=not bool(args.nosubmit),
        use_gpu=not bool(args.nogpu),
        hybrid_job=not bool(args.nohybrid),
        cluster_config=cluster_config,
    )

    return res


if __name__ == "__main__":
    main()
