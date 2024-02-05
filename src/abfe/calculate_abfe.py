import glob
import os
from typing import List

from abfe.orchestration.build_approach_flow import build_approach_flow
from abfe.orchestration.build_ligand_flow import build_ligand_flows
from abfe.scripts import final_receptor_results


def calculate_abfe(protein_pdb_path: str, ligand_sdf_paths: List[str], out_root_folder_path: str,
                   approach_name: str = "", cofactor_sdf_path: str = None,
                   n_cores_per_job: int = 8, num_jobs_receptor_workflow: int = None, num_jobs_per_ligand: int = 40, num_replicas: int = 3, small_mol_ff="openff",
                   submit: bool = False, use_gpu: bool = True, hybrid_job: bool = True, cluster_config: dict = {}):
    orig_dir = os.getcwd()
    conf = {}

    # IO:
    ## Input standardization
    conf["input_protein_pdb_path"] = os.path.abspath(protein_pdb_path)
    conf["input_ligands_sdf_path"] = [os.path.abspath(ligand_sdf_path) for ligand_sdf_path in ligand_sdf_paths]

    if (cofactor_sdf_path is not None):
        conf["input_cofactor_sdf_path"] = os.path.abspath(cofactor_sdf_path)
    else:
        conf["input_cofactor_sdf_path"] = None

    conf["out_approach_path"] = os.path.abspath(out_root_folder_path)

    ## Generate output folders
    for dir_path in [conf["out_approach_path"]]:
        if (not os.path.isdir(dir_path)):
            os.mkdir(dir_path)

    # Prepare Input / Parametrize
    os.chdir(conf["out_approach_path"])

    conf["ligand_names"] = [os.path.splitext(os.path.basename(sdf))[0] for sdf in conf["input_ligands_sdf_path"]]
    conf["num_jobs"] = num_jobs_receptor_workflow if (num_jobs_receptor_workflow is not None) else len(conf["ligand_names"]) * num_replicas * 2
    conf["num_replica"] = num_replicas
    conf['build_system'] = True
    conf["small_mol_ff"] = small_mol_ff

    print("Prepare")
    print("\tstarting preparing ABFE-ligand file structur")
    build_ligand_flows(input_ligand_paths=conf["input_ligands_sdf_path"],
                       input_protein_path=conf["input_protein_pdb_path"],
                       input_cofactor_path=conf["input_cofactor_sdf_path"],
                       out_root_path=conf["out_approach_path"],
                       num_max_thread=n_cores_per_job,
                       num_replicas=num_replicas, num_jobs=num_jobs_per_ligand,
                       cluster_config=cluster_config,
                       use_gpu=use_gpu, hybrid_job=hybrid_job)

    print("\tstarting preparing ABFE-Approach file structure: ", out_root_folder_path)
    expected_out_paths = int(num_replicas) * len(conf["ligand_names"])
    result_paths = glob.glob(conf["out_approach_path"] + "/*/*/dG*tsv")

    job_approach_file_path= None
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
