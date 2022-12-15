import os, stat, json
import subprocess

class scheduler():
    
    def __init__(self, out_dir_path:str, n_cores:int=1, time:str="96:00:00", partition="cpu") -> None:
        self.n_cores=n_cores
        self.out_dir_path = out_dir_path
        self.out_job_path = out_dir_path+"/job.sh"
        self.out_scheduler_path = out_dir_path+"/scheduler.sh"
        self.time=time
        self.partition=partition
    
    def generate_scheduler_file(self, out_prefix):           
        file_str = "\n".join([
            "#!/bin/env bash",
            "",
            #"conda activate abfe",
            "",
            "sbatch -p cpu -c "+str(self.n_cores)+" -J "+str(out_prefix)+"_scheduler --time="+self.time+" "+self.out_job_path
            ])
        
        file_io = open(self.out_scheduler_path, "w")
        file_io.write(file_str)
        file_io.close()
        os.chmod(self.out_scheduler_path, stat.S_IRWXU+stat.S_IRGRP+stat.S_IXGRP+stat.S_IROTH+stat.S_IXOTH)
        
        return self.out_scheduler_path
        
    def generate_job_file(self, out_prefix, cluster_conf_path:str, cluster_config:dict=None, cluster=False, num_jobs:int=1, latency_wait:int=60):
        if(cluster):
            root_dir = os.path.dirname(cluster_conf_path)
            slurm_logs = os.path.dirname(cluster_conf_path)+"/slurm_logs"
            if(not os.path.exists(slurm_logs)): os.mkdir(slurm_logs)
            
            if(cluster and cluster_config is None):
                cluster_config ={
                    "partition": "cpu",
                    "time": "48:00:00",
                    "mem": "20GB",
                }
            if(not all([x in cluster_config for x in ["partition", "time", "mem"]])):
                raise ValueError("missing keys in cluster_config! at least give: [\"partition\", \"time\", \"mem\"] ", cluster_config)
            
            self.n_cores=1
            cluster_config.update({
                    "cpus-per-task": '{threads}',
                    "cores-per-socket": '{threads}',
                    "chdir": root_dir,                
                    "job-name": "\\\""+str(out_prefix)+".{name}.{jobid}\\\"",
                    "output":  "\\\""+slurm_logs+"/"+str(out_prefix)+"_{name}_{jobid}.out\\\"",
                    "error":  "\\\""+slurm_logs+"/"+str(out_prefix)+"_{name}_{jobid}.err\\\""
                })
            
            json.dump(cluster_config, open(cluster_conf_path, "w"), indent="  ")
            cluster_options = " ".join(["--"+key+"="+str(val)+" " for key, val in cluster_config.items()])
            
            #TODO: change this here, such each job can access resource from cluster-config!
            file_str = "\n".join([
                "#!/bin/env bash",
                "snakemake --cluster \"sbatch "+cluster_options+"\" --cluster-config "+cluster_conf_path+" --jobs "+str(num_jobs)+" --latency-wait "+str(latency_wait)+" --rerun-incomplete "
            ])
            
        else:
            file_str = "\n".join([
                "#!/bin/env bash",
                "snakemake -c "+str(self.n_cores)+"  --latency-wait "+str(latency_wait)+" --rerun-incomplete "
            ])
        
        file_io = open(self.out_job_path, "w")
        file_io.write(file_str)
        file_io.close()
        os.chmod(self.out_job_path, stat.S_IRWXU+stat.S_IRGRP+stat.S_IXGRP+stat.S_IROTH+stat.S_IXOTH)

        return self.out_job_path

    def schedule_run(self)->int:
        orig_path = os.getcwd()
        os.chdir(self.out_dir_path)
        out = subprocess.getoutput(self.out_scheduler_path)
        job_id = int(out.split()[-1])
        os.chdir(orig_path)
        
        return job_id
        
    def submit_run(self, out_prefix="ABFE", cluster=True)->int:
        self.generate_job_file(out_prefix, cluster=cluster)
        self.generate_scheduler_file(out_prefix)
        out = self.schedule_run()
        return out