from abfe import scripts
run_path = config["run_path"]
complex_windows = config["complex_windows"]
n_coul_windows = config['n_coul_windows']
n_rest_windows = config['n_rest_windows']
n_vdw_windows = config['n_vdw_windows']

rule run_fep_complex_emin:
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

rule run_fep_complex_nvt_heat:
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

rule run_fep_complex_npt_eq2:
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

rule run_fep_complex_prod:
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

rule gather_complex_xvg:
    input:
        xvg_loc=expand(run_path+"/complex/fep/simulation/{state}/prod/prod.xvg",
                       state=complex_windows)
    params:
        sub_loc=run_path+"/complex/fep/simulation",
        rest_max_windows=n_rest_windows,
        vdw_max_windows=n_vdw_windows,
        coul_max_windows=n_coul_windows
    output:
        xvg_dir=directory(run_path+"/complex/fep/simulation/xvgs")
    shell:
        '''
            mkdir -p {params.sub_loc}/xvg/restraints-xvg
            mkdir -p {params.sub_loc}/xvg/vdw-xvg
            mkdir -p {params.sub_loc}/xvg/coul-xvg

            # restraints
            let max_window={params.rest_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sub_loc}/restraints.${{i}}/prod/prod.xvg {params.sub_loc}/xvg/restraints-xvg/dhdl.${{i}}.xvg
            done

            # vdw
            let max_window={params.vdw_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sub_loc}/vdw.${{i}}/prod/prod.xvg {params.sub_loc}/xvg/vdw-xvg/dhdl.${{i}}.xvg
            done

            # coul
            let max_window={params.coul_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sub_loc}/coul.${{i}}/prod/prod.xvg {params.sub_loc}/xvg/coul-xvg/dhdl.${{i}}.xvg
            done
        '''

rule get_dg_complex:
    input:
        xvg_dir=run_path+"/complex/fep/simulation/xvgs"
    params:
        out_dir=run_path+"/complex/",
        script_dir = scripts.root_path
    output:
        complex_var=run_path+"/complex/dg_results.csv"
    shell:
        "python {params.script_dir}/alchemlyb-analysis-complex.py --xvgpath {input.xvg_dir} --outpath {params.out_dir}"

