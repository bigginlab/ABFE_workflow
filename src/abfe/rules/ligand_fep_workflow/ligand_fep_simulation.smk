run_path = config["run_path"]
simulation_dir = run_path+"/ligand/fep/simulation"
top_dir = run_path+"/ligand/fep/fep-topology"

num_sim_threads = config['num_sim_threads']
num_retries = config['num_retries']

gromacs_run_script=config['gmx_run_kernel_path']
gromacs_cont_script=config['gmx_cont_kernel_path']


rule fep_run_ligand_emin:
    input:
        gro=top_dir+"/equil.gro",
        top=top_dir+"/ligand.top",
        mdp=simulation_dir+"/{state}/emin/emin.mdp"
    params:
        nthreads=num_sim_threads,
        run_dir=simulation_dir+"/{state}/emin",
        gmx_template=gromacs_run_script
    output:
        gro=simulation_dir+"/{state}/emin/emin.gro"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_emin.sh   
            chmod +x ./job_emin.sh
            ./job_emin.sh {params.nthreads} emin {input.top} {input.gro}
        '''

rule fep_run_ligand_nvt_heat:
    input:
        top=top_dir+"/ligand.top",
        mdp=simulation_dir+"/{state}/nvt/nvt.mdp",
        gro=simulation_dir+"/{state}/emin/emin.gro"
    params:
        nthreads=num_sim_threads,
        run_dir=simulation_dir+"/{state}/nvt",
        gmx_template=gromacs_run_script
    output:
        gro=simulation_dir+"/{state}/nvt/nvt.gro",
        cpt=simulation_dir+"/{state}/nvt/nvt.cpt"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_nvt.sh   
            chmod +x ./job_nvt.sh
            ./job_nvt.sh {params.nthreads} nvt {input.top} {input.gro}
        '''

rule fep_run_ligand_npt_eq1:
    input:
        top=top_dir+"/ligand.top",
        mdp=simulation_dir+"/{state}/npt/npt.mdp",
        gro=simulation_dir+"/{state}/nvt/nvt.gro",
        cpt=simulation_dir+"/{state}/nvt/nvt.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=simulation_dir+"/{state}/npt",
        gmx_template=gromacs_cont_script
    output:
        gro=simulation_dir+"/{state}/npt/npt.gro",
        cpt=simulation_dir+"/{state}/npt/npt.cpt"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_npt.sh   
            chmod +x ./job_npt.sh
            ./job_npt.sh {params.nthreads} npt {input.top} {input.gro} {input.cpt}
        '''

rule fep_run_ligand_npt_eq2:
    input:
        top=top_dir+"/ligand.top",
        mdp=simulation_dir+"/{state}/npt-norest/npt-norest.mdp",
        gro=simulation_dir+"/{state}/npt/npt.gro",
        cpt=simulation_dir+"/{state}/npt/npt.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=simulation_dir+"/{state}/npt-norest",
        gmx_template=gromacs_cont_script
    output:
        gro=simulation_dir+"/{state}/npt-norest/npt-norest.gro",
        cpt=simulation_dir+"/{state}/npt-norest/npt-norest.cpt"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_npt_norest.sh   
            chmod +x ./job_npt_norest.sh
           ./job_npt_norest.sh {params.nthreads} npt-norest {input.top} {input.gro} {input.cpt}
        '''

rule fep_run_ligand_prod:
    input:
        top=top_dir+"/ligand.top",
        mdp=simulation_dir+"/{state}/prod/prod.mdp",
        gro=simulation_dir+"/{state}/npt-norest/npt-norest.gro",
        cpt=simulation_dir+"/{state}/npt-norest/npt-norest.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=simulation_dir+"/{state}/prod",
        gmx_template=gromacs_cont_script
    output:
        gro=simulation_dir+"/{state}/prod/prod.gro",
        xvg=simulation_dir+"/{state}/prod/prod.xvg"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_prod.sh   
            chmod +x ./job_prod.sh
            ./job_prod.sh {params.nthreads} prod {input.top} {input.gro} {input.cpt}
        '''
