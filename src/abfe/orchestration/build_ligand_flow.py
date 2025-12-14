import os
import shutil
from typing import List

import abfe
from abfe import template
from abfe.orchestration import generate_conf, generate_snake, generate_scheduler

from abfe.conf import std_conf


def build_input(
    input_ligand_path: str,
    input_protein_path: str,
    input_cofactor_path: str,
    out_ligand_path: str,
):
    out_ligand_input_path = out_ligand_path + "/input"
    out_orig_ligand_input_path = out_ligand_input_path + "/orig_in"

    ## Generate folders
    for dir_path in [out_ligand_input_path, out_orig_ligand_input_path]:
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

    input_ligand_path = shutil.copyfile(
        input_ligand_path,
        out_orig_ligand_input_path + "/" + os.path.basename(input_ligand_path),
    )
    input_protein_path = shutil.copyfile(
        input_protein_path,
        out_orig_ligand_input_path + "/" + os.path.basename(input_ligand_path),
    )
    input_cofactor_path = (
        shutil.copyfile(
            input_cofactor_path,
            out_orig_ligand_input_path + "/" + os.path.basename(input_ligand_path),
        )
        if (input_cofactor_path is not None)
        else None
    )

    return out_ligand_input_path, (
        input_ligand_path,
        input_protein_path,
        input_cofactor_path,
    )


def build_replicas_simulation_flow(
    out_ligand_path: str,
    input_ligand_path: str,
    ligand_name: str,
    num_max_thread: int = 1,
    approach_name="",
    num_replicas: int = 3,
    cluster_config={},
    submit: bool = False,
    num_jobs=1,
    use_gpu: bool = True,
    hybrid_job: bool = True,
):
    print("Use_GPU: ", use_gpu)
    code_path = os.path.abspath(os.path.dirname(abfe.__file__))
    outs = []
    ligand_rep_name = ""
    for num_replica in range(1, num_replicas + 1):
        ligand_rep_name = ligand_name + "_rep" + str(num_replica)
        out_replica_path = out_ligand_path + "/" + str(num_replica)

        if not os.path.isdir(out_replica_path):
            os.mkdir(out_replica_path)

        # set global files:
        global_snake_path = out_replica_path + "/Snakefile.smk"
        conf_path = out_replica_path + "/snake_conf.json"

        generate_snake.generate_snake_file(
            out_file_path=global_snake_path, conf_file_path=conf_path
        )

        # build scheduler class
        scheduler = generate_scheduler.scheduler(
            out_dir_path=out_replica_path, cluster_config=cluster_config
        )
        # In a use_gpu and hybrid_job setting, we use cpu for ligand and gpu for complex
        if use_gpu and hybrid_job:
            job_configs = [
                {
                    "subdir": "ligand",
                    "snake_file_name": "Snakefile.smk",
                    "conf_file_name": "snake_conf.json",
                    "gpu": False,
                    "snake_job": "fep_ana_get_dg_ligand",
                    "job_name_suffix": "job_ligand.sh",
                },
                {
                    "subdir": "complex",
                    "snake_file_name": "Snakefile.smk",
                    "conf_file_name": "snake_conf.json",
                    "gpu": True,
                    "snake_job": "fep_ana_get_dg_complex",
                    "job_name_suffix": "job_complex.sh",
                },
            ]
        else:
            job_configs = [
                {
                    "subdir": ".",
                    "snake_file_name": global_snake_path,  # Use global snake path directly
                    "conf_file_name": conf_path,  # Use global conf path directly
                    "gpu": use_gpu,
                    "snake_job": None,  # Default target
                    "job_name_suffix": "job.sh",
                }
            ]

        replica_job_paths = []

        for config in job_configs:
            # Determine paths
            if config["subdir"] == ".":
                current_process_path = out_replica_path
                snake_path = config["snake_file_name"]
                conf_file_path = config["conf_file_name"]
            else:
                current_process_path = os.path.join(out_replica_path, config["subdir"])
                if not os.path.exists(current_process_path):
                    os.mkdir(current_process_path)
                snake_path = os.path.join(
                    current_process_path, config["snake_file_name"]
                )
                conf_file_path = os.path.join(
                    current_process_path, config["conf_file_name"]
                )

                # For subdirs, we need to generate specific snakefiles.
                # Use global confiuration for template filling
                generate_snake.generate_snake_file(
                    out_file_path=snake_path, conf_file_path=conf_path
                )

            # Generate configuration for this strand
            std_cont_key = (
                "gmx_kernel_gpu_cont" if config["gpu"] else "gmx_kernel_cpu_cont"
            )
            std_run_key = "gmx_kernel_gpu" if config["gpu"] else "gmx_kernel_cpu"

            generate_conf.generate_ligand_conf(
                out_path=conf_file_path,
                run_path=out_replica_path,
                num_sim_threads=num_max_thread,
                input_data_path=input_ligand_path,
                num_replica=num_replica,
                code_path=code_path,
                gmx_cont_kernel_path=std_conf[std_cont_key],
                gmx_run_kernel_path=std_conf[std_run_key],
                gmx_add_flag=std_conf["gmx_add_flag"],
            )

            # Update scheduler output path
            scheduler.out_job_path = os.path.join(
                current_process_path, config["job_name_suffix"]
            )

            # Generate job file
            job_file_path = scheduler.generate_job_file(
                cluster=cluster_config is not None,
                cluster_config=cluster_config,
                cluster_conf_path=os.path.join(
                    current_process_path, "cluster_conf.json"
                ),
                out_prefix=ligand_rep_name,
                num_jobs=num_jobs,
                snake_file_path=snake_path,
                snake_job=config["snake_job"],
            )
            replica_job_paths.append(job_file_path)

        scheduler.out_job_path = replica_job_paths

        scheduler._final_job_path = job_file_path
        _ = scheduler.generate_scheduler_file(
            out_prefix=f"{approach_name}_{ligand_rep_name}"
        )

        if submit:
            out = scheduler.schedule_run()
            print("submitted " + str(input_ligand_path), out)
            outs.append(out)

    if submit:
        return outs
    else:
        return None


def build_ligand_flows(
    input_ligand_paths: List[str],
    input_protein_path: str,
    input_cofactor_path: str,
    out_root_path: str,
    num_replicas: int,
    cluster_config: dict,
    num_jobs: int,
    num_max_thread: int,
    use_gpu: bool = True,
    hybrid_job: bool = True,
):
    job_ids = []
    for input_ligand_path in input_ligand_paths:
        ligand_name = os.path.splitext(os.path.basename(input_ligand_path))[0]
        out_ligand_path = out_root_path + "/" + str(ligand_name)
        print("\t\tLigand: ", ligand_name)

        if not os.path.exists(out_ligand_path):
            os.mkdir((out_ligand_path))

        out_ligand_input_path, _ = build_input(
            input_ligand_path=input_ligand_path,
            input_protein_path=input_protein_path,
            input_cofactor_path=input_cofactor_path,
            out_ligand_path=out_ligand_path,
        )

        build_replicas_simulation_flow(
            out_ligand_path=out_ligand_path,
            input_ligand_path=out_ligand_input_path,
            ligand_name=ligand_name,
            num_max_thread=num_max_thread,
            num_replicas=num_replicas,
            cluster_config=cluster_config,
            submit=False,
            num_jobs=num_jobs,
            use_gpu=use_gpu,
            hybrid_job=hybrid_job,
        )
