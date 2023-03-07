from abfe import template

run_path = config["run_path"]
input_path = config['input_data_path']

rule equil_setup_complex:
    input:
        input_path+"/complex"
    params:
        sim_dir=run_path+"/complex",
        complex_template=template.complex_equil_template_path
    output:
        complex_top=directory(run_path+"/complex/topology"),
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/topology/complex.gro"
    shell:
        '''
            cp -r {params.complex_template} {params.sim_dir}/equil-mdsim
            cp -r {input}/* {params.sim_dir}/topology
        '''