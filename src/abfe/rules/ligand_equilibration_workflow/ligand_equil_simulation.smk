
run_path = config["run_path"]
num_sim_threads = config['num_sim_threads']

gromacs_run_script=config['gmx_run_kernel_path']
gromacs_cont_script=config['gmx_cont_kernel_path']


rule equil_run_ligand_emin:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/topology/ligand.gro"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/ligand/equil-mdsim/emin",
        gmx_template=gromacs_run_script
    output:
        gro=run_path+"/ligand/equil-mdsim/emin/emin.gro"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_emin.sh   
            chmod +x ./job_emin.sh
            ./job_emin.sh {params.nthreads} emin {input.top} {input.gro}
        '''

rule equil_run_ligand_nvt_heat:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/emin/emin.gro"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/ligand/equil-mdsim/nvt_heat",
        gmx_template=gromacs_run_script
    output:
        gro=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_nvt_heat.sh  
            chmod +x ./job_nvt_heat.sh
            ./job_nvt_heat.sh {params.nthreads} nvt_heat {input.top} {input.gro}
        '''

rule equil_run_ligand_npt_eq1:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/ligand/equil-mdsim/npt_equil1",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_npt_eq1.sh   
            chmod +x ./job_npt_eq1.sh
            ./job_npt_eq1.sh {params.nthreads} npt_equil1 {input.top} {input.gro} {input.cpt}
        '''

rule equil_run_ligand_npt_eq2:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/ligand/equil-mdsim/npt_equil2",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.cpt",
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_npt_eq2.sh   
            chmod +x ./job_npt_eq2.sh
            ./job_npt_eq2.sh {params.nthreads} npt_equil2 {input.top} {input.gro} {input.cpt}
        '''
