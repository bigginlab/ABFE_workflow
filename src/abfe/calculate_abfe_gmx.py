import glob
import os
import shutil
from typing import List
import logging

from abfe.orchestration.build_approach_flow import build_approach_flow
from abfe.orchestration.build_ligand_flow import build_replicas_simulation_flow
from abfe.scripts import final_receptor_results

log = logging.getLogger(__file__)
log.setLevel(logging.INFO)
def calculate_abfe_gmx(input_dir:str, out_root_folder_path: str, approach_name: str = "",
                   n_cores_per_job: int = 8, num_jobs_receptor_workflow: int = None, num_jobs_per_ligand: int = 40, num_replicas: int = 3,
                   submit: bool = False, use_gpu: bool = True, hybrid_job: bool = True, cluster_config: dict = {}):
    orig_dir = os.getcwd()
    conf = {}

    # IO:
    ## Input standardization
    conf["input_path"] = os.path.abspath(input_dir)
    conf['build_system'] = False
    conf["out_approach_path"] = os.path.abspath(out_root_folder_path)
    conf["small_mol_ff"] = "custom"

    ## Generate output folders
    for dir_path in [conf["out_approach_path"]]:
        if (not os.path.isdir(dir_path)):
            os.mkdir(dir_path)

    # Prepare Input / Parametrize
    os.chdir(conf["out_approach_path"])

    #get Ligands:
    conf["num_jobs"] = num_jobs_receptor_workflow if (num_jobs_receptor_workflow is not None) else len(os.listdir(conf["input_path"])) * num_replicas * 2
    conf["num_replica"] = num_replicas

    ligand_dirs = [conf["input_path"]+"/"+d for d in os.listdir(conf["input_path"]) if(os.path.isdir(conf["input_path"]+"/"+d))]
    print("Found ligands: ")
    log.info("Found ligands:")
    conf["ligand_names"] = []
    for ligand_dir in ligand_dirs:
        ligand_name = os.path.basename(ligand_dir)
        print("\t", ligand_name)
        log.info("\t"+ligand_name)

        new_lig_dir = conf["out_approach_path"] + "/" + ligand_name

        if(not os.path.exists(new_lig_dir)):
            os.mkdir(new_lig_dir)

        if(not os.path.exists(new_lig_dir+"/input")):
            shutil.copytree(ligand_dir, new_lig_dir+"/input")

        conf["ligand_names"].append(ligand_name)
        build_replicas_simulation_flow(out_ligand_path=new_lig_dir,
                                       input_ligand_path=new_lig_dir+"/input",
                                       ligand_name=ligand_name,
                                       num_max_thread=n_cores_per_job,
                                       num_replicas=num_replicas, cluster_config=cluster_config, submit=False,
                                       num_jobs=num_jobs_per_ligand,
                                       use_gpu=use_gpu, hybrid_job=hybrid_job)

    log.info("\tstarting preparing ABFE-Approach file structur: ", out_root_folder_path)
    expected_out_paths = int(num_replicas) * len(conf["ligand_names"])
    result_paths = glob.glob(conf["out_approach_path"] + "/*/*/dG*tsv")

    if (len(result_paths) != expected_out_paths):
        print("\tBuild approach struct")

        job_approach_file_path = build_approach_flow(approach_name=approach_name,
                                                     num_jobs=conf["num_jobs"],
                                                     conf=conf, submit=submit,
                                                     cluster_config=cluster_config)
        print("Do")
        print("\tSubmit Job - ID: ", job_approach_file_path)

    # Final gathering
    print("\tAlready got results?: " + str(len(result_paths)))
    if (len(result_paths) > 0):
        print("Trying to gather ready results", out_root_folder_path)
        final_receptor_results.get_final_results(out_dir=out_root_folder_path, in_root_dir=out_root_folder_path)

    print()
    os.chdir(orig_dir)
