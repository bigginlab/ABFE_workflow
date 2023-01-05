from abfe import template

run_path = config["run_path"]
ligand_windows = config["ligand_windows"]
lam_coul_range = config['lam_coul_ligand_range']
lam_vdw_range = config['lam_vdw_ligand_range']
n_coul_windows = len(lam_coul_range)
n_vdw_windows = len(lam_vdw_range)

rule fep_setup_ligand:
    input:
        ligand_top=run_path+"/ligand/topology",
        equil_gro=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.gro",
    params:
        sim_dir=run_path+"/ligand/fep",
        vdw_windows=n_vdw_windows,
        vdw_range=" ".join(map(str, lam_vdw_range)),
        coul_windows=n_coul_windows,
        coul_range=" ".join(map(str, lam_coul_range)),
        template_dir = template.ligand_fep_template_path,
    output:
        fep_em=expand(run_path+"/ligand/fep/simulation/{state}/emin/emin.mdp", state=ligand_windows),
        fep_npt=expand(run_path+"/ligand/fep/simulation/{state}/npt/npt.mdp", state=ligand_windows),
        fep_npt_norest=expand(run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.mdp", state=ligand_windows),
        fep_nvt=expand(run_path+"/ligand/fep/simulation/{state}/nvt/nvt.mdp", state=ligand_windows),
        fep_prod=expand(run_path+"/ligand/fep/simulation/{state}/prod/prod.mdp", state=ligand_windows),
        fep_top=run_path+"/ligand/fep/fep-topology/ligand.top",
        fep_gro=run_path+"/ligand/fep/fep-topology/equil.gro"
    shell:
        '''
            mkdir -p {params.sim_dir}/template
            cp -r {params.template_dir}/template/* {params.sim_dir}/template

            cp -r {input.ligand_top} {params.sim_dir}/fep-topology
            cp {input.equil_gro} {params.sim_dir}/fep-topology/equil.gro

            # create simulation directory
            mkdir -p {params.sim_dir}/simulation
           
            let max_window={params.vdw_windows}
            for i in $(seq 0 $((max_window-1)))
            do
                mkdir -p {params.sim_dir}/simulation/vdw.${{i}}
                cp -r {params.sim_dir}/template/vdw/* {params.sim_dir}/simulation/vdw.${{i}}
                sed -i "s/<state>/${{i}}/g" {params.sim_dir}/simulation/vdw.${{i}}/*/*.mdp
                sed -i "s/<lamRange>/{params.vdw_range}/g" {params.sim_dir}/simulation/vdw.${{i}}/*/*.mdp
            done

            let max_window={params.coul_windows}
            for i in $(seq 0 $((max_window-1)))
            do
                mkdir -p {params.sim_dir}/simulation/coul.${{i}}
                cp -r {params.sim_dir}/template/coul/* {params.sim_dir}/simulation/coul.${{i}}
                sed -i "s/<state>/${{i}}/g" {params.sim_dir}/simulation/coul.${{i}}/*/*.mdp
                sed -i "s/<lamRange>/{params.coul_range}/g" {params.sim_dir}/simulation/coul.${{i}}/*/*.mdp
            done

        '''