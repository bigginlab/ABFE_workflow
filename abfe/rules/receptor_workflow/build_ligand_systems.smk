from abfe import scripts

#Final Check Job
approach_path = config["out_approach_path"]
ligand_names = config['ligand_names']

input_protein_pdb = config['input_protein_pdb_path']
input_ligand_sdfs = config['input_ligands_sdf_path']
input_cofactor_sdf = config['input_cofactor_sdf_path']


rule gather_files:
    input:
        ligand_sdfs=expand("{input_ligand_sdfs}", input_ligand_sdfs=input_ligand_sdfs)
    params:
        out_dir = approach_path+"/orig_input"
    output:
        out_files = expand(approach_path+"/orig_input/{ligand_name}.sdf", ligand_name=ligand_names)

    shell:
        """
            mkdir {params.out_dir} -p
            ligand_files="{input.ligand_sdfs}"
            echo "$ligand_files"

            for ligand_sdf in "$ligand_files"
            do
                cp $ligand_sdf {params.out_dir}
            done
        """

rule build_ligand_system:
    input:
        # TODO Now the class can work with all the ligands at hte same time, but they must be .mol
        protein_pdb= input_protein_pdb,
        ligand_sdf=approach_path+"/orig_input/{ligand_name}.sdf",
        output_dir= approach_path+"/{ligand_name}/input"
    params:
        cofactor_sdf= str(input_cofactor_sdf),
        script_dir = scripts.root_path
    output:
        out_ligand_gro=approach_path+"/{ligand_name}/input/ligand/ligand.gro",
        out_ligand_top=approach_path+"/{ligand_name}/input/ligand/ligand.top",
        out_complex_gro=approach_path+"/{ligand_name}/input/complex/complex.gro",
        out_complex_top=approach_path+"/{ligand_name}/input/complex/complex.top"
    shell:
        # TODO Generalize in case of COF of membrnaes, right now is not taked into account
        """
            python {params.script_dir}/preparation/system_builder.py --ligand_mol_dir {input.ligand_mol} \
            --protein_pdb_path {input.protein_pdb} \
            --cofactor_sdf_path {params.cofactor_sdf} --output_dir_path {input.output_dir}
        """