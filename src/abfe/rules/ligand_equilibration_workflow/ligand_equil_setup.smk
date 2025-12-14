from abfe import template

run_path = config["run_path"]
input_path = config['input_data_path']

rule equil_setup_ligand:
    input:
        ligand_input=input_path+"/ligand"
    params:
        sim_dir=run_path+"/ligand",
        template_dir=template.ligand_equil_template_path
    output:
        ligand_top=directory(run_path+"/ligand/topology"),
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/topology/ligand.gro",
        enmin_mdp=run_path+"/ligand/equil-mdsim/emin/emin.mdp"
    shell:
        r'''
            set -euo pipefail

            mkdir -p {params.sim_dir}/equil-mdsim/emin
            mkdir -p {params.sim_dir}/topology

            cp -r {params.template_dir}/. {params.sim_dir}/equil-mdsim
            cp -r {input.ligand_input}/. {params.sim_dir}/topology
        '''
