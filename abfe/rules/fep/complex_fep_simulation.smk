from abfe import scripts
run_path = config["run_path"]
complex_windows = config["complex_windows"]
n_coul_windows = config['n_coul_windows']
n_rest_windows = config['n_rest_windows']
n_vdw_windows = config['n_vdw_windows']

rule fep_run_complex_emin:
    input:
        mdp=run_path+"/complex/fep/simulation/{state}/em/em.mdp",
    params:
        fep_dir=run_path+"/complex/fep/simulation/{state}",
        top_dir=run_path+"/complex/fep/fep-topology/"
    output:
        gro=run_path+"/complex/fep/simulation/{state}/em/emin.gro"
    threads: 8
    shell:
        '''
            gmx grompp -f {params.fep_dir}/em/em.mdp -c {params.top_dir}/complex.gro \
                    -p {params.top_dir}/complex.top -o {params.fep_dir}/em/emin.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/em/emin -ntmpi {threads}
        '''

rule fep_run_complex_nvt_heat:
    input:
        mdp=run_path+"/complex/fep/simulation/{state}/nvt/nvt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/em/emin.gro"
    params:
        fep_dir=run_path+"/complex/fep/simulation/{state}",
        top_dir=run_path+"/complex/fep/fep-topology/"
    output:
        gro=run_path+"/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/nvt/nvt.cpt"
    threads: 8
    shell:
        '''
            gmx grompp -f {params.fep_dir}/nvt/nvt.mdp -c {input.gro} -r {input.gro} \
                    -p {params.top_dir}/complex.top -o {params.fep_dir}/nvt/nvt.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/nvt/nvt -ntmpi {threads}
        '''

rule run_fep_complex_npt_eq1:
    input:
        mdp=run_path+"/complex/fep/simulation/{state}/npt/npt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        fep_dir=run_path+"/complex/fep/simulation/{state}",
        top_dir=run_path+"/complex/fep/fep-topology/"
    output:
        gro=run_path+"/complex/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt/npt.cpt"
    threads: 8
    shell:
        '''
            gmx grompp -f {params.fep_dir}/npt/npt.mdp -c {input.gro} -t {input.cpt} \
                    -r {input.gro} -p {params.top_dir}/complex.top -o {params.fep_dir}/npt/npt.tpr -maxwarn 3
            gmx mdrun -deffnm {params.fep_dir}/npt/npt -ntmpi {threads}
        '''

rule fep_run_complex_npt_eq2:
    input:
        mdp=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt/npt.cpt"
    params:
        fep_dir=run_path+"/complex/fep/simulation/{state}",
        top_dir=run_path+"/complex/fep/fep-topology/"
    output:
        gro=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    threads: 8
    shell:
        '''
            gmx grompp -f {params.fep_dir}/npt-norest/npt-norest.mdp -c {input.gro} -t {input.cpt} \
                    -p {params.top_dir}/complex.top -o {params.fep_dir}/npt-norest/npt-norest.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/npt-norest/npt-norest -ntmpi {threads}
        '''

rule fep_run_complex_prod:
    input:
        mdp=run_path+"/complex/fep/simulation/{state}/prod/prod.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    params:
        fep_dir=run_path+"/complex/fep/simulation/{state}",
        top_dir=run_path+"/complex/fep/fep-topology/"
    output:
        gro=run_path+"/complex/fep/simulation/{state}/prod/prod.gro",
        xvg=run_path+"/complex/fep/simulation/{state}/prod/prod.xvg"
    threads: 8
    shell:
        '''
            gmx grompp -f {params.fep_dir}/prod/prod.mdp -c {input.gro} -t {input.cpt} \
                    -p {params.top_dir}/complex.top -o {params.fep_dir}/prod/prod.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/prod/prod -ntmpi {threads}
        '''
