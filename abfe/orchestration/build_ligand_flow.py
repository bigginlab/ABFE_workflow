import os
import tarfile
from typing import List

from abfe.scripts.preparation import system_builder as sb

import abfe
from abfe.orchestration import generate_conf, generate_snake, generate_scheduler


# TODO I am here!
# I have to check because one of this snakes make the files and I already did it, so
# that snake it is unnecessary
def build_replicas_simulation_flow(out_ligand_path: str, input_ligand_path: str, ligand_name: str, n_cores: int = 1,
                                   num_max_thread: int = 1,
                                   num_replicas: int = 3, cluster_config={}, submit: bool = False, num_jobs=1):
    code_path = os.path.abspath(os.path.dirname(abfe.__file__))

    outs = []
    for num_replica in range(1, num_replicas + 1):
        ligand_rep_name = ligand_name + "_rep" + str(num_replica)
        out_replica_path = out_ligand_path + "/" + str(num_replica)

        if (not os.path.isdir(out_replica_path)):
            os.mkdir(out_replica_path)

        # set global files:
        snake_path = out_replica_path + "/Snakefile"
        conf_path = out_replica_path + "/snake_conf.json"

        generate_snake.generate_snake_file(out_file_path=snake_path,
                                           conf_file_path=conf_path)

        # build scheduler class
        scheduler = generate_scheduler.scheduler(out_dir_path=out_replica_path, n_cores=n_cores)

        generate_conf.generate_ligand_conf(out_path=conf_path,
                                           run_path=out_replica_path,
                                           num_sim_threads=num_max_thread,
                                           input_data_path=input_ligand_path,
                                           num_replica=num_replica,
                                           code_path=code_path
                                           )

        job_file_path = scheduler.generate_job_file(cluster=True, cluster_config=cluster_config,
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



def build_ligand_flows(
        input_ligand_paths:List[str],
        input_protein_path:str,
        input_cofactor_path:str,
        input_membrane_path:str,
        hmr_factor:float,
        out_root_path:str,
        num_replicas:int,
        cluster_config:dict,
        num_jobs:int,
        num_max_thread:int):
    """Prepare all the directories and files

    Parameters
    ----------
    input_ligand_paths : List[str]
        List of path of the ligands
    input_protein_path : str
        Path for the PDB protein file
    input_cofactor_path : str
        Path for the mol cofactor file
    input_membrane_path : str
        Path for the PDB Membrane file. It must have the correct CRYST1 definition, it will be used for the solvation phase
    hmr_factor : float
        Hydrogen Mass Repartition factor
    out_root_path : str
        Out dir for all the systems
    num_replicas : int
        Number of Replicas
    cluster_config : dict
        Configuration of the HPC-Cluster
    num_jobs : int
        How many simulations will be launch in parallel
    num_max_thread : int
        ?
    """
    # Initialize the files builder
    builder = sb.MakeInputs(
        protein_pdb=input_protein_path,
        membrane_pdb=input_membrane_path,
        cofactor_mol=input_cofactor_path,
        hmr_factor = hmr_factor,
        keep_tmp_files_on=None
    )

    for input_ligand_path in input_ligand_paths:
        ligand_name = os.path.splitext(os.path.basename(input_ligand_path))[0]
        out_ligand_path = os.path.join(out_root_path, ligand_name)
        out_ligand_input_path = os.path.join(out_ligand_path, 'input')

        # Create topologies and input files
        builder(ligand_mol=input_ligand_path,out_dir=out_ligand_input_path)
        # Archive original files      
        with tarfile.open(os.path.join(out_ligand_input_path, 'orig_in'), "w:gz") as tar:
            tar.add(input_ligand_path, arcname=os.path.basename(input_ligand_path))
            tar.add(input_protein_path,arcname=os.path.basename(input_protein_path))
            if input_cofactor_path:
                tar.add(input_cofactor_path,arcname=os.path.basename(input_cofactor_path))
            if input_membrane_path:
                tar.add(input_membrane_path,arcname=os.path.basename(input_membrane_path))

        # Build the replicas
        build_replicas_simulation_flow(out_ligand_path=out_ligand_path,
                                       input_ligand_path=out_ligand_input_path,
                                       ligand_name=ligand_name,
                                       num_max_thread=num_max_thread,
                                       num_replicas=num_replicas, cluster_config=cluster_config, submit=False,
                                       num_jobs=num_jobs)

