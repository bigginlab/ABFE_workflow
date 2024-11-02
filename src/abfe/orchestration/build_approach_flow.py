import os
from abfe.orchestration import generate_conf, generate_snake, generate_scheduler


def build_approach_flow(approach_name: str, num_jobs: int, conf: dict, cluster_config={}, submit=False):
    out_path = conf["out_approach_path"]
    snake_path = out_path + "/Snakefile.smk"
    approach_conf_path = out_path + "/snake_conf.json"

    if("input_ligands_sdf_path" in conf):
        generate_conf.generate_approach_conf(out_path=approach_conf_path,
                                             out_approach_path=conf["out_approach_path"],
                                             input_ligands_sdf_path=conf["input_ligands_sdf_path"],
                                             input_protein_pdb_path=conf["input_protein_pdb_path"],
                                             input_cofactor_sdf_path=conf["input_cofactor_sdf_path"],
                                             ligand_names=conf["ligand_names"],
                                             num_replica=conf['num_replica'],
                                             python_bin=os.environ["CONDA_PREFIX"] + "/bin/python",
                                             build_system=conf["build_system"],
                                             small_mol_ff = conf["small_mol_ff"]
                                             )
    else:
        generate_conf.generate_approach_conf(out_path=approach_conf_path,
                                             out_approach_path=conf["out_approach_path"],
                                             input_ligands_sdf_path=None,
                                             input_protein_pdb_path=None,
                                             input_cofactor_sdf_path=None,
                                             ligand_names=conf["ligand_names"],
                                             num_replica=conf['num_replica'],
                                             python_bin=os.environ["CONDA_PREFIX"] + "/bin/python",
                                             build_system=conf["build_system"],
                                             small_mol_ff = conf["small_mol_ff"],
                                             )

    generate_snake.generate_approach_snake_file(out_file_path=snake_path,
                                                conf_file_path=approach_conf_path,
                                                gmx=not conf["build_system"])

    scheduler = generate_scheduler.scheduler(out_dir_path=out_path, n_cores=num_jobs, cluster_config=cluster_config)
    scheduler.generate_job_file(cluster=False,
                                out_prefix=approach_name, num_jobs=num_jobs,
                                snake_file_path=snake_path,
                                snake_job="", cluster_config=cluster_config)

    scheduler.generate_scheduler_file(out_prefix="ABFE_approach" + approach_name, )

    if (submit):
        out = scheduler.schedule_run()
    else:
        out = None
    return out
