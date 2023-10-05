import os
import shutil
from typing import List

import abfe
from abfe import template
from abfe.orchestration import generate_conf, generate_snake, generate_scheduler


def build_input(input_ligand_path: str,
                input_protein_path: str,
                input_cofactor_path: str,
                out_ligand_path: str, ):
    out_ligand_input_path = out_ligand_path + "/input"
    out_orig_ligand_input_path = out_ligand_input_path + "/orig_in"

    ## Generate folders
    for dir_path in [out_ligand_input_path, out_orig_ligand_input_path]:
        if (not os.path.isdir(dir_path)):
            os.mkdir(dir_path)

    input_ligand_path = shutil.copyfile(input_ligand_path,
                                        out_orig_ligand_input_path + "/" + os.path.basename(input_ligand_path))
    input_protein_path = shutil.copyfile(input_protein_path, out_orig_ligand_input_path + "/" + os.path.basename(input_ligand_path))
    input_cofactor_path = shutil.copyfile(input_cofactor_path, out_orig_ligand_input_path + "/" + os.path.basename(input_ligand_path)) if (
            input_cofactor_path is not None) else None

    return out_ligand_input_path, (input_ligand_path, input_protein_path, input_cofactor_path)


def build_replicas_simulation_flow(out_ligand_path: str, input_ligand_path: str, ligand_name: str, n_cores: int = 1,
                                   num_max_thread: int = 1,
                                   num_replicas: int = 3, cluster_config={}, submit: bool = False, num_jobs=1,
                                   use_gpu: bool = True, hybrid_job: bool = True):
    code_path = os.path.abspath(os.path.dirname(abfe.__file__))

    outs = []
    ligand_rep_name =  ""
    for num_replica in range(1, num_replicas + 1):
        ligand_rep_name = ligand_name + "_rep" + str(num_replica)
        out_replica_path = out_ligand_path + "/" + str(num_replica)

        if (not os.path.isdir(out_replica_path)):
            os.mkdir(out_replica_path)

        # set global files:
        snake_path = out_replica_path + "/Snakefile.smk.smk.smk"
        conf_path = out_replica_path + "/snake_conf.json"

        generate_snake.generate_snake_file(out_file_path=snake_path,
                                           conf_file_path=conf_path)

        # build scheduler class
        scheduler = generate_scheduler.scheduler(out_dir_path=out_replica_path, n_cores=n_cores)
        ##############################################################################
        # A bit hacky
        if (use_gpu and hybrid_job):

            # Prepare ligand strand
            out_out_ligand_path = out_replica_path + "/ligand"
            if (not os.path.exists(out_out_ligand_path)):
                os.mkdir(out_out_ligand_path)

            snake_path = out_out_ligand_path + "/Snakefile.smk.smk.smk"
            ligand_conf_path = out_out_ligand_path + "/snake_conf.json"
            generate_snake.generate_snake_file(out_file_path=snake_path,
                                               conf_file_path=conf_path)

            generate_conf.generate_ligand_conf(out_path=ligand_conf_path,
                                               run_path=out_replica_path,
                                               num_sim_threads=num_max_thread,
                                               input_data_path=input_ligand_path,
                                               num_replica=num_replica,
                                               code_path=code_path,
                                               )

            scheduler.out_job_path = out_out_ligand_path + "/job_ligand.sh"
            if(not cluster_config is None):
                cluster_config['partition'] = "cpu"
                cluster =True
            else:
                cluster =False
            job_ligand_file_path = scheduler.generate_job_file(cluster=cluster, cluster_config=cluster_config,
                                                               cluster_conf_path=out_out_ligand_path + "/cluster_conf.json",
                                                               out_prefix=ligand_rep_name, num_jobs=num_jobs,
                                                               snake_job="fep_ana_get_dg_ligand")

            # Prepare complex strand
            out_complex_path = out_replica_path + "/complex"
            if (not os.path.exists(out_replica_path + "/complex")):
                os.mkdir(out_complex_path)

            snake_path = out_complex_path + "/Snakefile.smk.smk.smk"
            complex_conf_path = out_complex_path + "/snake_conf.json"
            generate_snake.generate_snake_file(out_file_path=snake_path,
                                               conf_file_path=conf_path)

            generate_conf.generate_ligand_conf(out_path=complex_conf_path,
                                               run_path=out_replica_path,
                                               num_sim_threads=num_max_thread,
                                               input_data_path=input_ligand_path,
                                               num_replica=num_replica,
                                               code_path=code_path,
                                               gmx_cont_kernel_path=template.gmx_submit_kernels_path + "/def_gpu_job_cont.sh",
                                               gmx_run_kernel_path=template.gmx_submit_kernels_path + "/def_gpu_job.sh"
                                               )

            scheduler.out_job_path = out_complex_path + "/job_complex.sh"
            if(not cluster_config is None): cluster_config['partition'] = "gpu"
            job_complex_file_path = scheduler.generate_job_file(cluster=cluster, cluster_config=cluster_config,
                                                                cluster_conf_path=out_complex_path + "/cluster_conf.json",
                                                                out_prefix=ligand_rep_name, num_jobs=num_jobs,
                                                                snake_job="fep_ana_get_dg_complex")

            # Global script
            generate_conf.generate_ligand_conf(out_path=conf_path,
                                               run_path=out_replica_path,
                                               num_sim_threads=num_max_thread,
                                               input_data_path=input_ligand_path,
                                               num_replica=num_replica,
                                               code_path=code_path
                                               )

            scheduler.out_job_path = out_replica_path + "/job.sh"
            if(not cluster_config is None): cluster_config['partition'] = "cpu"
            job_file_path = scheduler.generate_job_file(cluster=cluster, cluster_config=cluster_config,
                                                        cluster_conf_path=out_replica_path + "/cluster_conf.json",
                                                        out_prefix=ligand_rep_name, num_jobs=num_jobs)

            # Final settings
            scheduler.out_job_path = [job_ligand_file_path, job_complex_file_path]
            if(not cluster_config is None): cluster_config['partition'] = "gpu"
            ##############################################################################
        else:
            if (use_gpu):
                generate_conf.generate_ligand_conf(out_path=conf_path,
                                                   run_path=out_replica_path,
                                                   num_sim_threads=num_max_thread,
                                                   input_data_path=input_ligand_path,
                                                   num_replica=num_replica,
                                                   code_path=code_path,
                                                   gmx_cont_kernel_path=template.gmx_submit_kernels_path + "/def_gpu_job_cont.sh",
                                                   gmx_run_kernel_path=template.gmx_submit_kernels_path + "/def_gpu_job.sh"
                                                   )
            else:
                generate_conf.generate_ligand_conf(out_path=conf_path,
                                                   run_path=out_replica_path,
                                                   num_sim_threads=num_max_thread,
                                                   input_data_path=input_ligand_path,
                                                   num_replica=num_replica,
                                                   code_path=code_path
                                                   )

            job_file_path = scheduler.generate_job_file(cluster=cluster_config is not None, cluster_config=cluster_config,
                                                        cluster_conf_path=out_replica_path + "/cluster_conf.json",
                                                        out_prefix=ligand_rep_name, num_jobs=num_jobs)
            scheduler.out_job_path = [job_file_path]

        scheduler._final_job_path = job_file_path
        _ = scheduler.generate_scheduler_file(out_prefix=ligand_rep_name, )

        if (submit):
            out = scheduler.schedule_run()
            print("submitted " + str(input_ligand_path), out)
            outs.append(out)

    if (submit):
        return outs
    else:
        return None


def calculate_all_ligands(input_ligand_paths: List[str],
                          out_root_path: str,
                          num_replicas: int, cluster_config: dict, submit: bool, num_jobs: int,
                          num_max_thread: int, use_gpu: bool = True, hybrid_job: bool = True):
    """
        Deappreciated!

    """
    job_ids = []
    for input_ligand_path in input_ligand_paths:
        ligand_name = os.path.splitext(os.path.basename(input_ligand_path))[0]
        out_ligand_path = out_root_path + "/" + str(ligand_name)

        if (not os.path.exists(out_ligand_path)):
            os.mkdir((out_ligand_path))

        job_id = build_replicas_simulation_flow(out_ligand_path=out_ligand_path, input_ligand_path=input_ligand_path,
                                                num_max_thread=num_max_thread,
                                                num_replicas=num_replicas, cluster_config=cluster_config, submit=submit,
                                                num_jobs=num_jobs,
                                                use_gpu=use_gpu, hybrid_job=hybrid_job)
        job_ids.append(job_id)

    return job_ids


def build_ligand_flows(input_ligand_paths: List[str],
                       input_protein_path: str,
                       input_cofactor_path: str,
                       out_root_path: str,
                       num_replicas: int, cluster_config: dict, num_jobs: int,
                       num_max_thread: int, use_gpu: bool = True, hybrid_job: bool = True):
    job_ids = []
    for input_ligand_path in input_ligand_paths:
        ligand_name = os.path.splitext(os.path.basename(input_ligand_path))[0]
        out_ligand_path = out_root_path + "/" + str(ligand_name)
        print("\t\tLigand: ", ligand_name)

        if (not os.path.exists(out_ligand_path)):
            os.mkdir((out_ligand_path))

        out_ligand_input_path, _ = build_input(input_ligand_path=input_ligand_path,
                                               input_protein_path=input_protein_path,
                                               input_cofactor_path=input_cofactor_path,
                                               out_ligand_path=out_ligand_path)

        build_replicas_simulation_flow(out_ligand_path=out_ligand_path,
                                       input_ligand_path=out_ligand_input_path,
                                       ligand_name=ligand_name,
                                       num_max_thread=num_max_thread,
                                       num_replicas=num_replicas, cluster_config=cluster_config, submit=False,
                                       num_jobs=num_jobs,
                                       use_gpu=use_gpu, hybrid_job=hybrid_job)
