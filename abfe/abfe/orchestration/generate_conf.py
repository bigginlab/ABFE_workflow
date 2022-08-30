import json 
import numpy as np

def generate_ligand_conf(out_path:str, code_path:str, run_path:str,  input_data_path:str,
                         num_replica:int=1,
                         n_vdw_windows:int=21, n_rest_windows:int=12, n_coul_windows:int=11):
  
    ## get all the window ids
    lam_vdw_range = list(np.round(np.linspace(0,1, n_vdw_windows),2))
    lam_coul_range = list(np.round(np.linspace(0,1, n_coul_windows),2))
    lam_rest_range = list(np.round(np.linspace(0,1, n_rest_windows),2))
    
    vdw_windows = [f'vdw.{i}' for i in range(n_vdw_windows)]
    rest_windows = [f'restraints.{i}' for i in range(n_rest_windows)]
    coul_windows = [f'coul.{i}' for i in range(n_coul_windows)]
    complex_windows = vdw_windows + rest_windows + coul_windows
    ligand_windows = vdw_windows + coul_windows
    
    #Parallel:
    num_sim_threads= 8

    run_num = num_replica
    run_path = run_path

    conf_settings = { k:v for k,v in locals().items() if not k.startswith("__")}

    out_IO = open(out_path,"w")
    json.dump(conf_settings, out_IO, indent=4)