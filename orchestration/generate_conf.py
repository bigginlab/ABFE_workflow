import json 
import numpy as np

def generate_ligand_conf(out_path:str, code_path:str, run_path:str,  input_data_path:str,
                         num_replica:int=1, num_sim_threads:int=8,
                         n_vdw_windows_complex:int=21, n_rest_windows_complex:int=11, n_coul_windows_complex:int=11,
                         n_vdw_windows_ligand:int=11,  n_coul_windows_ligand:int=11,):
  
    ## get all the window ids
    lam_vdw_complex_range = list(np.round(np.linspace(0,1, n_vdw_windows_complex),2))
    lam_coul_complex_range = list(np.round(np.linspace(0,1, n_rest_windows_complex),2))
    lam_rest_complex_range = list(np.round(np.linspace(0,1, n_coul_windows_complex),2))
    
    lam_vdw_ligand_range = list(np.round(np.linspace(0,1, n_vdw_windows_ligand),2))
    lam_coul_ligand_range = list(np.round(np.linspace(0,1, n_coul_windows_ligand),2))
    
    
    vdw_complex_windows = [f'vdw.{i}' for i in range(n_vdw_windows_complex)]
    rest_complex_windows = [f'restraints.{i}' for i in range(n_rest_windows_complex)]
    coul_complex_windows = [f'coul.{i}' for i in range(n_coul_windows_complex)]
    
    vdw_ligand_windows = [f'vdw.{i}' for i in range(n_vdw_windows_ligand)]
    coul_ligand_windows = [f'coul.{i}' for i in range(n_coul_windows_ligand)]
    
    
    complex_windows = vdw_complex_windows + rest_complex_windows + coul_complex_windows
    ligand_windows = vdw_ligand_windows + coul_ligand_windows
    
    #Parallel:
    num_sim_threads= num_sim_threads
    

    run_num = num_replica
    run_path = run_path
    num_retries = 3

    conf_settings = { k:v for k,v in locals().items() if not k.startswith("__")}
    
    out_IO = open(out_path,"w")
    json.dump(conf_settings, out_IO, indent=4)