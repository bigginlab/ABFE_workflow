import json
import os
import stat
import copy
import subprocess
from abfe.template import default_slurm_config_path


class scheduler:

    def __init__(
        self,
        out_dir_path: str,
        n_cores: int = 1,
        time: str = "96:00:00",
        partition: str = "cpu",
        cluster_config: dict = {},
    ) -> None:
        self.n_cores = n_cores
        self.out_dir_path = out_dir_path
        self.out_job_path = out_dir_path + "/job.sh"
        self.out_scheduler_path = out_dir_path + "/scheduler.sh"
        self.time = time
        self.partition = partition
        self.def_cluster_config = json.load(open(f"{default_slurm_config_path}", "r"))

        self.cluster_config = cluster_config
        if cluster_config is not None:
            for big_key, small_keys in self.def_cluster_config.items():
                if big_key in self.cluster_config:
                    for def_key, def_val in small_keys.items():
                        if def_key not in self.cluster_config[big_key]:
                            self.cluster_config[big_key][def_key] = def_val
                else:
                    self.cluster_config[big_key] = small_keys

    def generate_scheduler_file(self, out_prefix):
        if isinstance(self.out_job_path, str):
            self.out_job_path = [self.out_job_path]

        file_str = [
            "#!/bin/env bash",
            "",
        ]
        # SnakeMake Scheduler File:
        cluster_config = copy.deepcopy(self.cluster_config["Snakemake_job"])

        for i, job_path in enumerate(self.out_job_path):
            basename = os.path.basename(job_path).replace(".sh", "")

            cluster_config["queue_job_options"]["job-name"] = (
                str(out_prefix) + "_" + str(basename) + "_scheduler "
            )
            cluster_options = " ".join(
                [
                    "--" + key + "=" + str(val) + " "
                    for key, val in cluster_config["queue_job_options"].items()
                ]
            )

            file_str.extend(
                [
                    "",
                    "cd " + os.path.dirname(job_path),
                    "job"
                    + str(i)
                    + "=$("
                    + cluster_config["queue_submission_cmd"]
                    + " "
                    + cluster_options
                    + job_path
                    + ")",
                    "jobID" + str(i) + "=$(echo $job" + str(i) + " | awk '{print $4}')",
                    'echo "${jobID' + str(i) + '}"',
                ]
            )

        if len(self.out_job_path) > 1:
            cluster_config = copy.deepcopy(self.cluster_config["Snakemake_job"])
            cluster_config["job-name"] = str(out_prefix) + "_final_ana_scheduler "
            cluster_options = " ".join(
                [
                    "--" + key + "=" + str(val) + " "
                    for key, val in cluster_config["queue_job_options"].items()
                ]
            )

            dependency_key = cluster_config["queue_dependency"]["key"]
            dependency_value = cluster_config["queue_dependency"]["value"]
            dependency_sep = cluster_config["queue_dependency"]["sep"]

            deps = (
                f"--{dependency_key}={dependency_value}{dependency_sep}"
                + dependency_sep.join(
                    ["${jobID" + str(i) + "}(" for i in range(len(self.out_job_path))]
                )
            )
            file_str.append("\n")
            file_str.append(
                "echo "
                + ":".join(
                    ["${jobID" + str(i) + "}" for i in range(len(self.out_job_path))]
                )
            )
            file_str.append(
                cluster_config["queue_submission_cmd"]
                + " "
                + cluster_options
                + " "
                + deps
                + " "
                + self._final_job_path
            )

        file_str = "\n".join(file_str)
        file_io = open(self.out_scheduler_path, "w")
        file_io.write(file_str)
        file_io.close()
        os.chmod(
            self.out_scheduler_path,
            stat.S_IRWXU + stat.S_IRGRP + stat.S_IXGRP + stat.S_IROTH + stat.S_IXOTH,
        )

        return self.out_scheduler_path

    def generate_job_file(
        self,
        out_prefix,
        cluster_conf_path: str = None,
        cluster_config: dict = None,
        cluster=False,
        num_jobs: int = 1,
        latency_wait: int = 1000,
        snake_file_path=None,
        snake_job="",
    ):

        if snake_file_path is not None:
            snake_job = " -s " + snake_file_path + " " + snake_job

        if (
            cluster
            and self.cluster_config is not None
            and cluster_conf_path is not None
        ):
            root_dir = os.path.dirname(cluster_conf_path)
            slurm_logs = os.path.dirname(cluster_conf_path) + "/slurm_logs"
            if not os.path.exists(slurm_logs):
                os.mkdir(slurm_logs)

            if out_prefix == "":
                name = str(out_prefix) + "{name}.{jobid}"
                log = slurm_logs + "/" + str(out_prefix) + "{name}_{jobid}"
            else:
                name = str(out_prefix) + ".{name}.{jobid}"
                log = slurm_logs + "/" + str(out_prefix) + "_{name}_{jobid}"

            cluster_config = copy.deepcopy(self.cluster_config["Sub_job"])
            cluster_config["queue_job_options"].update(
                {
                    "cpus-per-task": "{threads}",
                    "cores-per-socket": "{threads}",
                    "chdir": root_dir,
                    "job-name": '\\"' + name + '\\"',
                    "output": '\\"' + log + '.out\\"',
                    "error": '\\"' + log + '.err\\"',
                }
            )

            json.dump(
                self.cluster_config["Sub_job"]["queue_job_options"],
                open(cluster_conf_path, "w"),
                indent="  ",
            )
            cluster_options = (
                " ".join(
                    [
                        "--" + key + "=" + str(val) + " "
                        for key, val in cluster_config["queue_job_options"].items()
                    ]
                )
                + " --parsable"
            )

            status_script_path = self.cluster_config["Sub_job"]["queue_status_script"]

            # TODO: change this here, such each job can access resource from cluster-config!
            file_str = "\n".join(
                [
                    "#!/bin/env bash",
                    'snakemake --cluster "'
                    + cluster_config["queue_submission_cmd"]
                    + " "
                    + cluster_options
                    + '" '
                    "--cluster-config " + cluster_conf_path + " "
                    "--cluster-status " + status_script_path + " "
                    '--cluster-cancel "' + cluster_config["queue_abort_cmd"] + '" '
                    "--jobs "
                    + str(num_jobs)
                    + " --latency-wait "
                    + str(latency_wait)
                    + " "
                    "--rerun-incomplete "
                    + snake_job
                    + " 1>  "
                    + str(out_prefix)
                    + ".out 2>"
                    + str(out_prefix)
                    + ".err",
                ]
            )
        elif cluster:
            raise ValueError("give cluster conf! ")
        else:

            file_str = "\n".join(
                [
                    "#!/bin/env bash",
                    "snakemake -c "
                    + str(self.n_cores)
                    + " -j "
                    + str(num_jobs)
                    + " --latency-wait "
                    + str(latency_wait)
                    + " --rerun-incomplete "
                    + snake_job,
                ]
            )

        file_io = open(self.out_job_path, "w")
        file_io.write(file_str)
        file_io.close()
        os.chmod(
            self.out_job_path,
            stat.S_IRWXU + stat.S_IRGRP + stat.S_IXGRP + stat.S_IROTH + stat.S_IXOTH,
        )

        return self.out_job_path

    def schedule_run(self) -> int:
        orig_path = os.getcwd()
        os.chdir(self.out_dir_path)
        out = subprocess.getoutput(self.out_scheduler_path)

        job_id = int(out.strip())
        os.chdir(orig_path)

        return job_id

    def submit_run(self, out_prefix="ABFE", cluster=True) -> int:
        self.generate_job_file(out_prefix, cluster=cluster)
        self.generate_scheduler_file(out_prefix)
        out = self.schedule_run()
        return out
