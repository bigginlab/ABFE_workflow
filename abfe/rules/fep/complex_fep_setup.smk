from abfe import template


run_path = config["run_path"]
complex_windows = config["complex_windows"]
lam_rest_range = config['lam_rest_range']
n_rest_windows = len(lam_rest_range)
lam_coul_range = config['lam_coul_range']
n_coul_windows = len(lam_coul_range)
lam_vdw_range = config['lam_vdw_range']
n_vdw_windows = len(lam_vdw_range)

rule setup_complex_fep:
    input:
        complex_top=run_path+"/complex/topology",
        boresch_top=run_path+"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro"
    params:
        sim_dir=run_path+"/complex/fep",
        restraint_range=" ".join(map(str, lam_rest_range)),
        restraint_windows=n_rest_windows,
        vdw_range=" ".join(map(str, lam_vdw_range)),
        vdw_windows=n_vdw_windows,
        coul_range=" ".join(map(str, lam_coul_range)),
        coul_windows=n_coul_windows,
        template_dir=template.complex_fep_template_path,
    output:
        fep_em=expand(run_path+"/complex/fep/simulation/{state}/em/em.mdp", state=complex_windows),
        fep_npt=expand(run_path+"/complex/fep/simulation/{state}/npt/npt.mdp", state=complex_windows),
        fep_npt_norest=expand(run_path+"/complex/fep/simulation/{state}/npt-norest/npt-norest.mdp", state=complex_windows),
        fep_nvt=expand(run_path+"/complex/fep/simulation/{state}/nvt/nvt.mdp", state=complex_windows),
        fep_prod=expand(run_path+"/complex/fep/simulation/{state}/prod/prod.mdp", state=complex_windows),
    shell:
        '''
            mkdir -p {params.sim_dir}/template
            cp -r {params.template_dir}/template/* {params.sim_dir}/template
            cp -r {input.complex_top} {params.sim_dir}/fep-topology
            cp {input.boresch_gro} {params.sim_dir}/fep-topology/complex.gro
            cat {params.sim_dir}/fep-topology/complex.top {input.boresch_top} > {params.sim_dir}/fep-topology/complex_boresch.top

            # create simulation directory
            mkdir -p {params.sim_dir}/simulation

            let max_window={params.restraint_windows}
            for i in $(seq 0 ${{max_window}})
            do
                mkdir -p {params.sim_dir}/simulation/restraints.${{i}}
                cp -r {params.sim_dir}/template/restraints/* {params.sim_dir}/simulation/restraints.${{i}}
                sed -i "s/<state>/${{i}}/g" {params.sim_dir}/simulation/restraints.${{i}}/*/*.mdp
                sed -i "s/<lamRange>/{params.restraint_range}/g" {params.sim_dir}/simulation/restraints.${{i}}/*/*.mdp
            done
            
            let max_window={params.vdw_windows}
            for i in $(seq 0 ${{max_window}})
            do
                mkdir -p {params.sim_dir}/simulation/vdw.${{i}}
                cp -r {params.sim_dir}/template/vdw/* {params.sim_dir}/simulation/vdw.${{i}}
                sed -i "s/<state>/${{i}}/g" {params.sim_dir}/simulation/vdw.${{i}}/*/*.mdp
                sed -i "s/<lamRange>/{params.vdw_range}/g" {params.sim_dir}/simulation/vdw.${{i}}/*/*.mdp
            done

            let max_window={params.coul_windows}
            for i in $(seq 0 ${{max_window}})
            do
                mkdir -p {params.sim_dir}/simulation/coul.${{i}}
                cp -r {params.sim_dir}/template/coul/* {params.sim_dir}/simulation/coul.${{i}}
                sed -i "s/<state>/${{i}}/g" {params.sim_dir}/simulation/coul.${{i}}/*/*.mdp
                sed -i "s/<lamRange>/{params.coul_range}/g" {params.sim_dir}/simulation/coul.${{i}}/*/*.mdp
            done
        '''