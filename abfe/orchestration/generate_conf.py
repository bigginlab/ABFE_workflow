import json
from typing import List

import numpy as np

from abfe import template


def generate_approach_conf(out_path: str,
                           input_ligand_mol_paths: str,
                           input_protein_pdb_path: str,
                           input_cofactor_mol_path: str,
                           input_membrane_pdb_path:str,
                           out_approach_path: str,
                           ligand_basenames: List[str],
                           num_replica: int
                           ):
    ## Ugly implementation every defined variable is added to conf! :)
    conf_settings = {k: v for k, v in locals().items() if not k.startswith("__")}

    with open(out_path, "w") as out_IO:
        json.dump(conf_settings, out_IO, indent=4)


def generate_ligand_conf(out_path: str, code_path: str, run_path: str, input_data_path: str, # TODO, it is not used
                         num_replica: int = 1, num_sim_threads: int = 8,
                         n_vdw_windows_complex: int = 21, n_rest_windows_complex: int = 11, n_coul_windows_complex: int = 11,
                         n_vdw_windows_ligand: int = 11, n_coul_windows_ligand: int = 11,
                         gmx_run_kernel_path: str = template.gmx_submit_kernels_path + "/def_cpu_job.sh",
                         gmx_cont_kernel_path: str = template.gmx_submit_kernels_path + "/def_cpu_job_cont.sh"):
    ## Ugly implementation every defined variable is added to conf! :)

    ## get all the window ids
    lam_vdw_complex_range = list(np.round(np.linspace(0, 1, n_vdw_windows_complex), 2))  # TODO, it is not used
    lam_coul_complex_range = list(np.round(np.linspace(0, 1, n_rest_windows_complex), 2))  # TODO, it is not used
    lam_rest_complex_range = list(np.round(np.linspace(0, 1, n_coul_windows_complex), 2))  # TODO, it is not used

    lam_vdw_ligand_range = list(np.round(np.linspace(0, 1, n_vdw_windows_ligand), 2))  # TODO, it is not used
    lam_coul_ligand_range = list(np.round(np.linspace(0, 1, n_coul_windows_ligand), 2))  # TODO, it is not used

    vdw_complex_windows = [f'vdw.{i}' for i in range(n_vdw_windows_complex)]
    rest_complex_windows = [f'restraints.{i}' for i in range(n_rest_windows_complex)]
    coul_complex_windows = [f'coul.{i}' for i in range(n_coul_windows_complex)]

    vdw_ligand_windows = [f'vdw.{i}' for i in range(n_vdw_windows_ligand)]
    coul_ligand_windows = [f'coul.{i}' for i in range(n_coul_windows_ligand)]

    complex_windows = vdw_complex_windows + rest_complex_windows + coul_complex_windows  # TODO, it is not used
    ligand_windows = vdw_ligand_windows + coul_ligand_windows  # TODO, it is not used

    # Parallel:
    num_sim_threads = num_sim_threads

    run_num = num_replica  # TODO, it is not used
    run_path = run_path
    num_retries = 3 # TODO, it is not used

    gmx_run_kernel_path = gmx_run_kernel_path
    gmx_cont_kernel_path = gmx_cont_kernel_path

    conf_settings = {k: v for k, v in locals().items() if not k.startswith("__")}

    with open(out_path, "w") as out_IO:
        json.dump(conf_settings, out_IO, indent=4)
