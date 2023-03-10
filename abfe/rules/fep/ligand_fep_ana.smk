from abfe import scripts

run_path = config["run_path"]
ligand_windows = config["ligand_windows"]
n_coul_windows = config['n_coul_windows_ligand']
n_vdw_windows = config['n_vdw_windows_ligand']


# Ana
rule fep_ana_gather_ligand_xvg:
    input:
        xvg_loc=expand(run_path+"/ligand/fep/simulation/{state}/prod/prod.xvg",
                       state=ligand_windows)
    params:
        sim_loc=run_path+"/ligand/fep/simulation",
        ana_loc=run_path+"/ligand/fep/ana",
        vdw_max_windows=n_vdw_windows,
        coul_max_windows=n_coul_windows
    output:
        xvg_dir=directory(run_path+"/ligand/fep/ana/xvgs")
    shell:
        '''
            mkdir -p {params.ana_loc}/xvgs/vdw-xvg
            mkdir -p {params.ana_loc}/xvgs/coul-xvg

            # vdw
            let max_window={params.vdw_max_windows}
            for i in $(seq 0 $((max_window-1)))
            do
                cp {params.sim_loc}/vdw.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/vdw-xvg/dhdl.${{i}}.xvg
            done

            # coul
            let max_window={params.coul_max_windows}
            for i in $(seq 0 $((max_window-1)))
            do
                cp {params.sim_loc}/coul.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/coul-xvg/dhdl.${{i}}.xvg
            done
        '''

rule fep_ana_get_dg_ligand:
    input:
        xvg_dir=run_path+"/ligand/fep/ana/xvgs"
    params:
        conf_path = run_path+"/snake_conf.json",
        out_dir=run_path+"/ligand/fep/ana",
        script_dir=scripts.root_path
    output:
        complex_var=run_path+"/ligand/fep/ana/dg_results.csv"
    shell:
        "python {params.script_dir}/alchemlyb-analysis-ligand.py --xvgpath {input.xvg_dir} --confpath {params.conf_path} --outpath {params.out_dir}"