from abfe import scripts

run_path = config["run_path"]
ligand_windows = config["ligand_windows"]
num_sim_threads = config['num_sim_threads']
n_coul_windows = config['n_coul_windows']
n_vdw_windows = config['n_vdw_windows']

rule run_fep_ligand_emin:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/em/em.mdp"
    params:
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/em/emin.gro"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.fep_dir}/em/em.mdp -c {params.top_dir}/equil.gro \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/em/emin.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/em/emin -ntmpi {threads}
        '''

rule run_fep_ligand_nvt_heat:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/em/emin.gro"
    params:
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.cpt"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.fep_dir}/nvt/nvt.mdp -c {input.gro} -r {input.gro} \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/nvt/nvt.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/nvt/nvt -ntmpi {threads}
        '''

rule run_fep_ligand_npt_eq1:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/npt/npt.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt/npt.cpt"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.fep_dir}/npt/npt.mdp -c {input.gro} -t {input.cpt} \
                    -r {input.gro} -p {params.top_dir}/ligand.top -o {params.fep_dir}/npt/npt.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/npt/npt -ntmpi {threads}
        '''

rule run_fep_ligand_npt_eq2:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt/npt.cpt"
    params:
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.fep_dir}/npt-norest/npt-norest.mdp -c {input.gro} -t {input.cpt} \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/npt-norest/npt-norest.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/npt-norest/npt-norest -ntmpi {threads}
        '''

rule run_fep_ligand_prod:
    input:
        mdp=run_path+"/ligand/fep/simulation/{state}/prod/prod.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    params:
        fep_dir=run_path+"/ligand/fep/simulation/{state}",
        top_dir=run_path+"/ligand/fep/fep-topology/"
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/prod/prod.gro",
        xvg=run_path+"/ligand/fep/simulation/{state}/prod/prod.xvg"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.fep_dir}/prod/prod.mdp -c {input.gro} -t {input.cpt} \
                    -p {params.top_dir}/ligand.top -o {params.fep_dir}/prod/prod.tpr -maxwarn 2
            gmx mdrun -deffnm {params.fep_dir}/prod/prod -ntmpi {threads}
        '''

rule gather_ligand_xvg:
    input:
        xvg_loc=expand(run_path+"/ligand/fep/simulation/{state}/prod/prod.xvg",
                       state=ligand_windows)
    params:
        sub_loc=run_path+"/ligand/fep/simulation",
        vdw_max_windows=n_coul_windows,
        coul_max_windows=n_coul_windows
    output:
        xvg_dir=directory(run_path+"/ligand/fep/simulation/xvgs")
    shell:
        '''
            mkdir -p {params.sub_loc}/xvg/vdw-xvg
            mkdir -p {params.sub_loc}/xvg/coul-xvg

            # vdw
            for i in $(seq 0 1 {params.vdw_max_windows})
            do
                cp {params.sub_loc}/vdw.${{i}}/prod/prod.xvg {params.sub_loc}/xvg/vdw-xvg/dhdl.${{i}}.xvg
            done

            # coul
            for i in $(seq 0 1 {params.coul_max_windows})
            do
                cp {params.sub_loc}/coul.${{i}}/prod/prod.xvg {params.sub_loc}/xvg/coul-xvg/dhdl.${{i}}.xvg
            done
        '''

rule get_dg_ligand:
    input:
        xvg_dir=run_path+"/ligand/fep/simulation/xvgs"
    params:
        out_dir=run_path+"/ligand/",
        script_dir=scripts.root_path
    output:
        complex_var=run_path+"/ligand/dg_results.csv"
    shell:
        "python {params.script_dir}/alchemlyb-analysis-ligand.py --xvgpath {input.xvg_dir} --outpath {params.out_dir}"
        
