from abfe import scripts

#Final Check Job
approach_path = config["approach_path"]
ligand_number = config["num_ligands"]

rule gather_receptor_results:
    input:
        approach_path=approach_path,
        prior_result_paths=expand(approach_path+"/ligand-{ligand_number}/dG_results.csv", state=ligand_number)
    params:
        script_dir = scripts.root_path
    output:
        out_dG_File=run_path+"/abfe_final_result.csv",

    shell:
        "python {params.script_dir}/final_receptor_results.py --in_root_dir {input.approach_path} --out_dir {input.run_path}"