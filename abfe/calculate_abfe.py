
import os
import subprocess
import multiprocessing

from typing import List

from abfe.scripts.preparation.generate_ABFE_systems import prepare_input_files
from abfe.orchestration.build_and_run_ligands import calculate_all_ligands
from abfe.scripts import final_receptor_results


def calculate_abfe(protein_pdb_path:str, ligand_sdf_paths: List[str], out_root_folder_path:str, cofactor_sdf_path:str=None,
             n_cores_per_job:int = 8, num_jobs:int = 40, num_replicas:int=1, 
             submit:bool=False, use_gpu:bool=False, hybrid_job:bool=False, cluster_config:dict={}, force_parametrize:bool=False):
    orig_dir = os.getcwd()

    # IO:
    ## Input standardization
    protein_pdb_path = os.path.abspath(protein_pdb_path)
    ligand_sdf_paths = [os.path.abspath(ligand_sdf_path) for ligand_sdf_path in ligand_sdf_paths]

    out_root_folder_path = os.path.abspath(out_root_folder_path)
    if(cofactor_sdf_path is not None): cofactor_sdf_path = os.path.abspath(cofactor_sdf_path)

    ## Generate output folders
    input_dir = out_root_folder_path+"/input"
    for dir_path in [out_root_folder_path, input_dir]:
        if(not os.path.isdir(dir_path)):
            os.mkdir(dir_path)

    # Prepare Input / Parametrize
    os.chdir(out_root_folder_path)    

    print("Parameterize ")        
    for sdf in ligand_sdf_paths:
        out_folder = input_dir+"/"+os.path.splitext(os.path.basename(sdf))[0]+""
        if(not os.path.isdir(out_folder)):
            os.mkdir(out_folder)
        prepare_input_files(protein_pdb=protein_pdb_path, ligand_sdf=sdf, cofactor_sdf=cofactor_sdf_path, out_dir=out_folder)
    
    ## Check prepared_input
    print("Gather Input")
    input_ligand_paths = [out_root_folder_path+"/input/"+dir+"" for dir in os.listdir(out_root_folder_path+"/input") if(os.path.isdir(out_root_folder_path+"/input/"+dir+""))]
    print("input ligand dirs: ", input_ligand_paths)
    print("output root dir: ", out_root_folder_path)


    # Do ABFE
    print("starting ABFE: submit-", str(submit))
    job_ids = calculate_all_ligands(input_ligand_paths=input_ligand_paths, out_root_path=out_root_folder_path,  num_max_thread = n_cores_per_job,
                           num_replicas=num_replicas,
                           submit=submit, num_jobs=num_jobs, 
                           cluster_config=cluster_config, 
                           use_gpu=use_gpu, hybrid_job=hybrid_job)
    print("\t started ABFEs with job_ids: ", job_ids)
    os.chdir(out_root_folder_path)
    
    
    # Final gathering
    if(len(job_ids)>0):
        cmd = "python "+final_receptor_results.__file__+" -i "+out_root_folder_path+" -o "+out_root_folder_path
        subprocess.getoutput("sbatch -c1 --dependecy=afterok:"+":".join(map(str, job_ids))+" "+cmd)
    else:
        print("Trying to gather results", out_root_folder_path)
        final_receptor_results.get_final_results(out_dir=out_root_folder_path, in_root_dir=out_root_folder_path)
    os.chdir(orig_dir)
