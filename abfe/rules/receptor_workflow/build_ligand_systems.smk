from abfe import scripts

#Final Check Job
approach_path = config["out_approach_path"]
ligand_names = config['ligand_names']

input_protein_pdb = config['input_protein_pdb_path']
input_ligand_sdfs = config['input_ligands_sdf_path']
input_cofactor_sdf = config['input_cofactor_sdf_path']


rule build_ligand_systems:
    input:
        protein_pdb= input_protein_pdb,
        ligand_sdf= expand("{input_ligand_sdf}", input_ligand_sdf=input_ligand_sdfs),
        output_dir= expand(approach_path+"/{ligand_name}/input",  ligand_name=ligand_names)
    params:
        cofactor_sdf= str(input_cofactor_sdf),
        script_dir = scripts.root_path
    output:
        out_ligand_gro=expand(approach_path+"/{ligand_name}/input/ligand/ligand.gro", ligand_name=ligand_names),
        out_ligand_top=expand(approach_path+"/{ligand_name}/input/ligand/ligand.top", ligand_name=ligand_names),
        out_complex_gro=expand(approach_path+"/{ligand_name}/input/complex/complex.gro", ligand_name=ligand_names),
        out_complex_top=expand(approach_path+"/{ligand_name}/input/complex/complex.top", ligand_name=ligand_names)
    shell:
        """
            python {params.script_dir}/preparation/generate_ABFE_systems.py --ligand_sdf_dir {input.ligand_sdf} \
            --protein_pdb_path {input.protein_pdb} \
            --cofactor_sdf_path {params.cofactor_sdf} --output_dir_path {input.output_dir}
        """