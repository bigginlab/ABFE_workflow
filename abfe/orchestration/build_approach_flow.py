from abfe.orchestration import generate_conf, generate_snake, generate_scheduler


def build_approach_flow(approach_name: str, num_jobs: int, conf: dict, cluster_config:dict, submit=False):
    out_path = conf["out_approach_path"]
    snake_path = out_path + "/Snakefile"
    approach_conf_path = out_path + "/snake_conf.json"

    generate_conf.generate_approach_conf(out_path=approach_conf_path,
                                         out_approach_path=conf["out_approach_path"],
                                         input_ligand_mol_paths=conf["input_ligand_mol_paths"],
                                         input_protein_pdb_path=conf["input_protein_pdb_path"],
                                         input_cofactor_mol_path=conf["input_cofactor_mol_path"],
                                         input_membrane_pdb_path=conf["input_membrane_pdb_path"],
                                         hmr_factor=conf["hmr_factor"],
                                         ligand_names=conf["ligand_names"],
                                         num_replica=conf['num_replica']
                                         )

    # TODO fix pass same info time partition and then as for the cluster config, ask from the initialization step.
    generate_snake.generate_approach_snake_file(out_file_path=snake_path,
                                                conf_file_path=approach_conf_path)

    scheduler = generate_scheduler.scheduler(out_dir_path=out_path, n_cores=num_jobs, cluster_config=cluster_config,)
    scheduler.generate_job_file(cluster=False,
                                out_prefix=approach_name, num_jobs=num_jobs,
                                snake_job="")

    scheduler.generate_scheduler_file(out_prefix="ABFE_approach" + approach_name, )

    if (submit):
        out = scheduler.schedule_run()
    else:
        out = None
    return out
