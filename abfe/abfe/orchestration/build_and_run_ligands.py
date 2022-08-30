
import os
from typing import List

from abfe.orchestration import generate_conf, generate_snake, generate_scheduler


def build_run_ligand(out_root_path:str, input_ligand_path:str, 
                     num_replicas:int=3, n_cores:int=10, submit:bool=False, num_jobs=1):
    prefix = "_".join(input_ligand_path.split("/")[-2:])
    code_path = os.path.abspath(os.path.dirname(__file__)+"../../..")   #TODO: shaky approach but works for now!
    ligand_path = out_root_path+"/"+str(os.path.basename(input_ligand_path))
    
    if(not os.path.isdir(ligand_path)):
        os.mkdir(ligand_path)

    for num_replica in range(1, num_replicas+1):
        out_path = ligand_path+"/"+str(num_replica)
        if(not os.path.isdir(out_path)):
            os.mkdir(out_path)

        #set files:    
        snake_path = out_path+"/Snakefile"
        conf_path = out_path+"/snake_conf.json"

        generate_snake.generate_snake_file(out_file_path=snake_path, 
            conf_file_path=conf_path, 
            code_file_path=code_path)
        
        generate_conf.generate_ligand_conf(out_path=conf_path,
                                        run_path=out_path,
                                        input_data_path=input_ligand_path,
                                        num_replica=num_replica,
                                        code_path=code_path
                                        )
        
        scheduler = generate_scheduler.scheduler(out_dir_path=out_path, n_cores=n_cores)
        job_file_path = scheduler.generate_job_file(cluster=True, cluster_conf_path=out_path+"/cluster_conf.json", out_prefix=prefix, num_jobs=num_jobs)
        sched_file_path = scheduler.generate_scheduler_file( out_prefix=prefix,)
        print(job_file_path, sched_file_path)
        
        if(submit):
            out = scheduler.schedule_run()
            print("submitted "+str(input_ligand_path), out)
            return out
        return None
    
def calculate_all_ligands(input_ligand_paths:List[str], out_root_path:str,
                         num_replicas:int, n_cores:int, submit:bool, num_jobs:int):
    job_ids = []
    for ligand_id in input_ligand_paths:
        job_id = build_run_ligand(out_root_path=out_root_path, input_ligand_path=ligand_id,
                                           num_replicas=num_replicas, n_cores=n_cores, submit=submit, num_jobs=num_jobs)
        job_ids.append(job_id)
        
    return job_ids