from abfe import scripts

#Final Check Job
approach_path = config["out_approach_path"]
ligand_names = config['ligand_names']
num_replica =  config['num_replica']
replica_list = list(map(str, range(1,1 + num_replica)))

rule gather_receptor_results:
    input:
        approach_path=approach_path,
        prior_result_paths=expand(approach_path + "/{ligand_names}/{replica}/dG_results.tsv",ligand_names=ligand_names,replica=replica_list, allow_missing=True)
    params:
        script_dir = scripts.root_path
    output:
        out_dG_File=approach_path+"/abfe_results.tsv",
    shell:
        "python {params.script_dir}/final_receptor_results.py --in_root_dir {input.approach_path} --out_dir {input.approach_path}"