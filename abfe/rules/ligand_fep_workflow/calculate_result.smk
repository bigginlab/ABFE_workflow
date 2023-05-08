run_path = config["run_path"]


rule fep_gather_dGs:
    input:
        complex_var=run_path+"/complex/fep/ana/dg_results.tsv",
        ligand_var=run_path+"/ligand/fep/ana/dg_results.tsv"
    params:
        conf_path = run_path+"/snake_conf.json",
        script_dir = scripts.root_path
    output:
        out_file_path=run_path+"/dG_results.tsv",
    shell:
        "python {params.script_dir}/free_energy/calculate_ABFE_ligand_dG.py --in_lig_path {input.ligand_var} --in_comp_path {input.complex_var} --out_csv_path {output.out_file_path}"

