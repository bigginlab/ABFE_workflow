run_path = config["run_path"]
num_retries = config['num_retries']
num_sim_threads = config['num_sim_threads']

gromacs_run_script=config['gmx_run_kernel_path']
gromacs_cont_script=config['gmx_cont_kernel_path']


rule fep_run_complex_emin:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        gro=run_path+"/complex/fep/fep-topology/complex.gro",
        mdp=run_path+"/complex/fep/simulation/{state}/emin/emin.mdp",
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/fep/simulation/{state}/emin",
        gmx_template=gromacs_run_script
    output:
        gro=run_path+"/complex/fep/simulation/{state}/emin/emin.gro"
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

rule fep_run_complex_nvt_heat:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/nvt/nvt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/emin/emin.gro"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/fep/simulation/{state}/nvt",
        gmx_template=gromacs_run_script
    output:
        gro=run_path+"/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/nvt/nvt.cpt"
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

rule run_fep_complex_npt_eq1:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/npt/npt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/fep/simulation/{state}/npt",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/complex/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt/npt.cpt"
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

rule fep_run_complex_npt_eq2:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt/npt.cpt"
    params:
        nthreads=num_sim_threads,        
        run_dir=run_path+"/complex/fep/simulation/{state}/npt-norest",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.cpt"
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

rule fep_run_complex_prod:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/prod/prod.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/fep/simulation/{state}/prod",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/complex/fep/simulation/{state}/prod/prod.gro",
        xvg=run_path+"/complex/fep/simulation/{state}/prod/prod.xvg"
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
