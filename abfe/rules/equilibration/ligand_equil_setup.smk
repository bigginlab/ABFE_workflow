from abfe import template

run_path = config["run_path"]
input_path = config['input_data_path']

rule equil_setup_ligand:
    input:
        input_path+"/ligand"
    params:
        sim_dir=run_path+"/ligand",
        template_dir=template.ligand_equil_template_path
    output:
        ligand_top=directory(run_path+"/ligand/topology"),
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/topology/ligand.gro"
    shell:
        '''
            cp -r {params.template_dir} {params.sim_dir}/equil-mdsim
            cp -r {input}/* {params.sim_dir}/topology
        '''