run_path = config["run_path"]
num_sim_threads = config['num_sim_threads']
num_retries = config['num_retries']

rule fep_run_ligand_emin:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/em/em.mdp"
    params:
        nthreads=num_sim_threads,
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/em/emin.gro"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.fep_dir}/em/em.mdp -c {params.top_dir}/equil.gro \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/em/emin.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/em/emin -ntomp {params.nthreads}
        '''

rule fep_run_ligand_nvt_heat:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/em/emin.gro"
    params:
        nthreads=num_sim_threads,
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.cpt"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.fep_dir}/nvt/nvt.mdp -c {input.gro} -r {input.gro} \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/nvt/nvt.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/nvt/nvt -ntomp {params.nthreads}
        '''

rule fep_run_ligand_npt_eq1:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/npt/npt.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        nthreads=num_sim_threads,
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt/npt.cpt"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.fep_dir}/npt/npt.mdp -c {input.gro} -t {input.cpt} \
                    -r {input.gro} -p {params.top_dir}/ligand.top -o {params.fep_dir}/npt/npt.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/npt/npt -ntomp {params.nthreads}
        '''

rule fep_run_ligand_npt_eq2:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt/npt.cpt"
    params:
        nthreads=num_sim_threads,
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.fep_dir}/npt-norest/npt-norest.mdp -c {input.gro} -t {input.cpt} \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/npt-norest/npt-norest.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/npt-norest/npt-norest -ntomp {params.nthreads}
        '''

rule fep_run_ligand_prod:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/prod/prod.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    params:
        nthreads=num_sim_threads,
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/prod/prod.gro",
        xvg=run_path+"/ligand/fep/simulation/{state}/prod/prod.xvg"
    threads: num_sim_threads
    retries: num_retries
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.fep_dir}/prod/prod.mdp -c {input.gro} -t {input.cpt} \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/prod/prod.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/prod/prod -ntomp {params.nthreads}
        '''
