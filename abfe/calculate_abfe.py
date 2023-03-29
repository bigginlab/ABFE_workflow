import glob
import os
from typing import List

from abfe.orchestration.build_approach_flow import build_approach_flow
from abfe.orchestration.build_ligand_flow import build_ligand_flows
from abfe.scripts import final_receptor_results
def calculate_abfe(
        protein_pdb_path: str,
        ligand_mol_paths: List[str],
        out_root_folder_path: str,
        approach_name: str = "",
        cofactor_mol_path: str = None,
        membrane_pdb_path: str = None,
        hmr_factor: float = 3.0,
        n_cores_per_job: int = 8,
        num_jobs_receptor_workflow: int = None,
        num_jobs_per_ligand: int = 40,
        num_replicas: int = 3,
        submit: bool = False,
        cluster_config: dict = {}):
    orig_dir = os.getcwd()
    conf = {}

    if hmr_factor < 2:
        raise ValueError(f'hmr_factor must be equal or higher than 2 (provided {hmr_factor}) to avoid instability during MD simulations. ABFE_workflow uses dt = 4 fs')
    # IO:
    ## Input standardization
    conf["input_protein_pdb_path"] = os.path.abspath(protein_pdb_path)
    conf["input_ligand_mol_paths"] = [os.path.abspath(ligand_mol_path) for ligand_mol_path in ligand_mol_paths]

    if not conf["input_ligand_mol_paths"]:
        raise ValueError(f'There were not provided any ligands or they are not accessible on: {ligand_mol_paths}')

    if cofactor_mol_path:
        conf["input_cofactor_mol_path"] = os.path.abspath(cofactor_mol_path)
    else:
        conf["input_cofactor_mol_path"] = None
    if membrane_pdb_path:
        conf["input_membrane_pdb_path"] = os.path.abspath(membrane_pdb_path)
    else:
        conf["input_membrane_pdb_path"] = None

    conf["hmr_factor"] = hmr_factor
    

    conf["out_approach_path"] = os.path.abspath(out_root_folder_path)

    ## Generate output folders
    for dir_path in [conf["out_approach_path"]]:
        if (not os.path.isdir(dir_path)):
            os.mkdir(dir_path)

    # Prepare Input / Parametrize
    os.chdir(conf["out_approach_path"])

    conf["ligand_names"] = [os.path.splitext(os.path.basename(mol))[0] for mol in conf["input_ligand_mol_paths"]]
    conf["num_jobs"] = num_jobs_receptor_workflow if (num_jobs_receptor_workflow is not None) else len(conf["ligand_names"]) * num_replicas * 2
    conf["num_replica"] = num_replicas

    print("Prepare")
    print("\tstarting preparing ABFE-ligand file structure")

    build_ligand_flows(input_ligand_paths=conf["input_ligand_mol_paths"],
                       input_protein_path=conf["input_protein_pdb_path"],
                       input_cofactor_path=conf["input_cofactor_mol_path"],
                       input_membrane_path=conf["input_membrane_pdb_path"],
                       hmr_factor = conf["hmr_factor"],
                       out_root_path=conf["out_approach_path"],
                       num_max_thread=n_cores_per_job,
                       num_replicas=num_replicas, num_jobs=num_jobs_per_ligand,
                       cluster_config=cluster_config)

    print("\tstarting preparing ABFE-Approach file structur: ", out_root_folder_path)
    expected_out_paths = int(num_replicas) * len(conf["ligand_names"])
    result_paths = glob.glob(conf["out_approach_path"] + "/*/*/dG*csv")

    if (len(result_paths) != expected_out_paths):
        print("\tBuild approach struct")
        # TODO is it right to set an empty dict the user input of cluster_config here
        cluster_config = {}
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
