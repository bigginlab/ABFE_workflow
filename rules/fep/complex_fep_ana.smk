from abfe import scripts
run_path = config["run_path"]
complex_windows = config["complex_windows"]

n_coul_windows = config['n_coul_windows_complex']
n_rest_windows = config['n_rest_windows_complex']
n_vdw_windows = config['n_vdw_windows_complex']


#Ana
rule fep_ana_gather_complex_xvg:
    input:
        xvg_loc=expand(run_path+"/complex/fep/simulation/{state}/prod/prod.xvg",
                       state=complex_windows)
    params:
        sim_loc=run_path+"/complex/fep/simulation",
        ana_loc=run_path+"/complex/fep/ana",
        rest_max_windows=n_rest_windows,
        vdw_max_windows=n_vdw_windows,
        coul_max_windows=n_coul_windows
    output:
        xvg_dir=directory(run_path+"/complex/fep/ana/xvgs")
    shell:
        '''
            mkdir -p {params.ana_loc}/xvgs/restraints-xvg
            mkdir -p {params.ana_loc}/xvgs/vdw-xvg
            mkdir -p {params.ana_loc}/xvgs/coul-xvg

            # restraints
            let max_window={params.rest_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sim_loc}/restraints.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/restraints-xvg/dhdl.${{i}}.xvg
            done

            # vdw
            let max_window={params.vdw_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sim_loc}/vdw.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/vdw-xvg/dhdl.${{i}}.xvg
            done

            # coul
            let max_window={params.coul_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sim_loc}/coul.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/coul-xvg/dhdl.${{i}}.xvg
            done
        '''

rule fep_ana_get_dg_complex:
    input:
        xvg_dir=run_path+"/complex/fep/ana/xvgs"
    params:
        conf_path = run_path+"/snake_conf.json",
        out_dir=run_path+"/complex/fep/ana",
        script_dir = scripts.root_path
    output:
        complex_var=run_path+"/complex/fep/ana/dg_results.csv"
    shell:
        "python {params.script_dir}/alchemlyb-analysis-complex.py --xvgpath {input.xvg_dir}  --confpath {params.conf_path} --outpath {params.out_dir}"

