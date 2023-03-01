from abfe import scripts

#Final Check Job
run_path = config["run_path"]
input_path = config['input_data_path']

input_protein_pdb = config['input_data_path']
input_ligand_sdf = config['input_data_path']
input_cofactor_sdf = config['input_data_path']

#ligand_number

rule build_system:
    input:
        protein_pdb=input_protein_pdb,
        ligand_sdf=input_ligand_sdf,
        cofactor_sdf = input_cofactor_sdf,
        output_dir= run_path+'/input'
    params:
        script_dir = scripts.root_path
    output:
        out_ligand_gro=input_path+"/ligand/ligand.gro",
        out_ligand_top=input_path+"/ligand/ligand.top",
        out_complex_gro=input_path+"/complex/complex.gro",
        out_complex_top=input_path+"/complex/complex.top"
    shell:
        "python {params.script_dir}/preparation-generate_ABFE_systems.py --ligand_sdf_dir {input.ligand_sdf_dir} --protein_pdb_path {input.protein_pdb}  --cofactor_sdf_path {input.cofactor_sdf} --output_dir_path {input.out_file_path}"