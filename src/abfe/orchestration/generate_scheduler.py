import json
import os
import stat
import copy
import subprocess
from pathlib import Path

from abfe.template import default_slurm_config_path


class scheduler:
    def __init__(
        self,
        out_dir_path: str,
        cluster_config: dict = {},
    ) -> None:
        self.out_dir_path = out_dir_path
        self.out_dir = Path(out_dir_path)
        # We assume out_job_path can change or be a list later, but start as a string path default
        self.out_job_path = str(self.out_dir / "job.sh")
        self.out_scheduler_path = str(self.out_dir / "scheduler.sh")
        self.cluster_config = cluster_config
        # Internal tracking for a final job, used in complex workflows
        self._final_job_path = None

    def generate_scheduler_file(self, out_prefix):
        # 1. Normalize job paths to list
        job_paths = self.out_job_path
        if isinstance(job_paths, str):
            job_paths = [job_paths]

        # 2. Extract configuration
        # Base config for the scheduler wrapper job itself (often runs on login node or light queue)
        sched_config = copy.deepcopy(self.cluster_config.get("Snakemake_job", {}))
        submission_cmd = sched_config.get("queue_submission_cmd", "sbatch")

        # 3. Build Script Content
        script_content = ["#!/bin/env bash", ""]

        job_ids_variable_names = []

        for i, job_path in enumerate(job_paths):
            job_path_obj = Path(job_path)
            # e.g., "job" from "job.sh"
            basename = job_path_obj.stem.replace(".sh", "")

            # Config for this particular submission step
            # We assume we use the same queue options as the base scheduler config
            # but with a modified job name.
            current_job_config = copy.deepcopy(sched_config)
            current_job_options = current_job_config.get("queue_job_options", {})
            current_job_options["job-name"] = f"{out_prefix}_{basename}_scheduler"

            # Create option flags: --key=val
            options_str = " ".join(
                [f"--{k}={v}" for k, v in current_job_options.items()]
            )

            # Shell commands to submit and capture ID
            var_name = f"jobID{i}"
            job_ids_variable_names.append(var_name)

            script_content.append("")
            script_content.append(f"cd {os.path.dirname(job_path)}")
            # Submit and capture output, then parse ID (assuming 4th word is ID, standard Slurm)
            script_content.append(
                f"job{i}=$({submission_cmd} {options_str} {job_path})"
            )
            script_content.append(f"{var_name}=$(echo $job{i} | awk '{{print $4}}')")
            script_content.append(f'echo "${{{var_name}}}"')

        # 4. Handle Final Dependency Job (if applicable)
        # If there are multiple jobs (e.g. ligand + complex), we schedule a final job that depends on them
        if len(job_paths) > 1 and getattr(self, "_final_job_path", None):
            final_job_config = copy.deepcopy(sched_config)
            final_job_config["queue_job_options"]["job-name"] = (
                f"{out_prefix}_final_ana_scheduler"
            )

            final_opts_str = " ".join(
                [
                    f"--{k}={v}"
                    for k, v in final_job_config.get("queue_job_options", {}).items()
                ]
            )

            # Dependency logic
            dep_info = final_job_config.get("queue_dependency", {})
            dep_key = dep_info.get("key", "dependency")
            dep_val = dep_info.get("value", "afterok")
            dep_sep = dep_info.get("sep", ":")

            # Construct dependency string (e.g. --dependency=afterok:123:124)
            ids_str = dep_sep.join([f"${{{var}}}" for var in job_ids_variable_names])
            dep_flag = f"--{dep_key}={dep_val}{dep_sep}{ids_str}"

            script_content.append("\n")
            # Echo all IDs together
            all_ids_echo = ":".join([f"${{{var}}}" for var in job_ids_variable_names])
            script_content.append(f'echo "{all_ids_echo}"')

            # Submit final job
            script_content.append(
                f"{submission_cmd} {final_opts_str} {dep_flag} {self._final_job_path}"
            )

        # 5. Write to file
        file_content = "\n".join(script_content)
        with open(self.out_scheduler_path, "w") as f:
            f.write(file_content)

        # 6. Make executable
        st = os.stat(self.out_scheduler_path)
        os.chmod(self.out_scheduler_path, st.st_mode | stat.S_IEXEC)

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
        # 1. Update snake job string if file path provided
        if snake_file_path is not None:
            snake_job = f" -s {snake_file_path} {snake_job}"

        script_content = []

        # 2. Case: Cluster Submission (Snakemake submits to Slurm)
        if cluster and self.cluster_config and cluster_conf_path:
            # Prepare paths
            cluster_conf_path_obj = Path(cluster_conf_path)
            root_dir = cluster_conf_path_obj.parent
            slurm_logs = root_dir / "slurm_logs"
            if not slurm_logs.exists():
                slurm_logs.mkdir()

            # Prepare job naming
            if out_prefix:
                job_name = str(out_prefix)
                log_base = str(slurm_logs / out_prefix)
            else:
                job_name = "job"
                log_base = str(slurm_logs / "job")

            # Create cluster config for this specific run
            sub_job_config = copy.deepcopy(self.cluster_config.get("Sub_job", {}))
            queue_opts = sub_job_config.get("queue_job_options", {})

            # Update typical Slurm fields
            # Using escaped quotes pattern from original code
            queue_opts.update(
                {
                    "chdir": str(root_dir),
                    "job-name": f'\\"{job_name}\\"',
                    "output": f'\\"{log_base}.out\\"',
                    "error": f'\\"{log_base}.err\\"',
                }
            )

            # Write cluster.json
            with open(cluster_conf_path, "w") as f:
                json.dump(queue_opts, f, indent="  ")

            # Construct Snakemake's cluster submission command
            # This is the command snakemake will call for EACH job.
            submission_cmd = sub_job_config.get("queue_submission_cmd", "sbatch")
            abort_cmd = sub_job_config.get("queue_abort_cmd", "scancel")
            status_script = sub_job_config.get("queue_status_script")

            # Options for the cluster command
            cluster_opts_str = " ".join([f"--{k}={v}" for k, v in queue_opts.items()])
            cluster_opts_str += " --parsable"

            cmd_parts = [
                "snakemake",
                f'--cluster "{submission_cmd} {cluster_opts_str}"',
                f"--cluster-config {cluster_conf_path}",
            ]

            if status_script:
                cmd_parts.append(f"--cluster-status {status_script}")

            cmd_parts.extend(
                [
                    f'--cluster-cancel "{abort_cmd}"',
                    f"--jobs {num_jobs}",
                    f"--latency-wait {latency_wait}",
                    "--rerun-incomplete",
                    snake_job.strip(),
                    f"1> {out_prefix}.out",
                    f"2> {out_prefix}.err",
                ]
            )

            script_content = ["#!/bin/env bash", " ".join(cmd_parts)]

        # 3. Case: Local Execution
        elif not cluster:
            cmd = (
                f"snakemake -c -j {num_jobs} "
                f"--latency-wait {latency_wait} "
                "--rerun-incomplete "
                f"{snake_job}"
            )
            script_content = ["#!/bin/env bash", cmd]

        # 4. Case: Invalid State
        else:
            # If cluster=True but no config provided, previous code raised error
            raise ValueError(
                "Cluster requested but no configuration or path available!"
            )

        # 5. Write Job File
        with open(self.out_job_path, "w") as f:
            f.write("\n".join(script_content))

        # Make executable
        st = os.stat(self.out_job_path)
        os.chmod(self.out_job_path, st.st_mode | stat.S_IEXEC)

        return self.out_job_path

    def schedule_run(self) -> int:
        orig_path = os.getcwd()
        os.chdir(self.out_dir_path)

        # Run the generated scheduler script
        # Using subprocess.getoutput to capture stdout easily
        out = subprocess.getoutput(self.out_scheduler_path)

        # Parse job ID (expecting integer output)
        try:
            job_id = int(out.strip())
        except ValueError:
            print(f"Warning: Could not parse job ID from output: '{out}'")
            job_id = 0

        os.chdir(orig_path)
        return job_id

    def submit_run(self, out_prefix="ABFE", cluster=True) -> int:
        self.generate_job_file(out_prefix, cluster=cluster)
        self.generate_scheduler_file(out_prefix)
        out = self.schedule_run()
        return out
