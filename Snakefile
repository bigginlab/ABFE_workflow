import os

# system definitions
## get all the ligand ids
ligand_ids = [dir[1] for dir in os.walk('data/samples')]

## get all the window ids
vdw_windows = [f'vdw.{i}' for i in range(4)]
rest_windows = [f'restraints.{i}' for i in range(2)]
coul_windows = [f'coul.{i}' for i in range(1)]
all_windows = vdw_windows + rest_windows + coul_windows
ligand_windows = vdw_windows + coul_windows


## set the number of replica
replicas = 1

rule clean:
    shell:
        '''
        rm -r data/results
        rm -r data/sims
        '''

rule all:
    input:
        [expand("data/results/{ligand}/{run}/complex/dg_results.csv",
                 ligand=ligand_ids[0], run=range(replicas)),
         expand("data/results/{ligand}/{run}/ligand/dg_results.csv",
                 ligand=ligand_ids[0], run=range(replicas))]

rule setup_complex:
    input:
        "data/samples/{ligand}/complex-topology"
    params:
        sim_dir="data/sims/{ligand}/{run}/complex/"
    output:
        top="data/sims/{ligand}/{run}/complex/topology/complex.top",
        gro="data/sims/{ligand}/{run}/complex/topology/complex.gro"
    shell:
        '''
        cp -r template/complex-equil-mdsim {params.sim_dir}/equil-mdsim
        cp -r {input} {params.sim_dir}/topology
        '''

rule setup_ligand:
    input:
        "data/samples/{ligand}/ligand-topology"
    params:
        sim_dir="data/sims/{ligand}/{run}/ligand/"
    output:
        top="data/sims/{ligand}/{run}/ligand/topology/ligand.top",
        gro="data/sims/{ligand}/{run}/ligand/topology/ligand.gro"
    shell:
        '''
        cp -r template/ligand-equil-mdsim {params.sim_dir}/equil-mdsim
        cp -r {input} {params.sim_dir}/topology
        '''

include: "rules/complex_equil.smk"
include: "rules/ligand_equil.smk"

rule setup_complex_fep:
    input:
        complex_top="data/samples/{ligand}/complex-topology",
        boresch_top="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_gro="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro"
    params:
        sim_dir="data/sims/{ligand}/{run}/complex/fep",
        restraint_windows=12,
        vdw_windows=21,
        coul_windows=11,
    output:
        fep_em=directory(expand("data/sims/{{ligand}}/{{run}}/complex/fep/simulation/{state}/em/em.mdp",
                                state=all_windows))
    shell:
        '''
        cp -r template/complex-fep {params.sim_dir}
        cp -r {input.complex_top} {params.sim_dir}/fep-topology
        cp {input.boresch_gro} {params.sim_dir}/fep-topology/complex.gro
        cat {params.sim_dir}/fep-topology/complex.top {input.boresch_top} > {params.sim_dir}/fep-topology/complex.top

        # create simulation directory
        mkdir -p {params.sim_dir}/simulation

        let max_window=${params.restraint_windows}-1
        for i in $(seq 0 ${{max_window}})
        do
            cp -r {params.sim_dir}/template/restraints {params.sim_dir}/simulation/restraints.${{i}}
            sed -i "s#<state>#${{i}}#g" {params.sim_dir}/simulation/restraints.${{i}}/*/*.mdp
        done

        let max_window=${params.vdw_windows}-1
        for i in $(seq 0 ${{max_window}})
        do
            cp -r {params.sim_dir}/template/vdw {params.sim_dir}/simulation/vdw.${{i}}
            sed -i "s#<state>#${{i}}#g" {params.sim_dir}/simulation/vdw.${{i}}/*/*.mdp
        done

        let max_window=${params.coul_windows}-1
        for i in $(seq 0 ${{max_window}})
        do
            cp -r {params.sim_dir}/template/coul {params.sim_dir}/simulation/coul.${{i}}
            sed -i "s#<state>#${{i}}#g" {params.sim_dir}/simulation/coul.${{i}}/*/*.mdp
        done
        '''

rule setup_ligand_fep:
    input:
        ligand_top="data/samples/{ligand}/ligand-topology",
        equil_gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil2/npt_equil2.gro"
    params:
        sim_dir="data/sims/{ligand}/{run}/ligand/fep",
        vdw_windows=21,
        coul_windows=11
    output:
        fep_em=directory(expand("data/sims/{{ligand}}/{{run}}/ligand/fep/simulation/{state}/em/em.mdp",
                                state=ligand_windows))
    shell:
        '''
        cp -r template/ligand-fep {params.sim_dir}
        cp -r {input.ligand_top} {params.sim_dir}/fep-topology
        cp {input.equil_gro} {params.sim_dir}/fep-topology/equil.gro

        # create simulation directory
        mkdir -p {params.sim_dir}/simulation

        let max_window=${params.vdw_windows}-1
        for i in $(seq 0 ${{max_window}})
        do
            cp -r {params.sim_dir}/template/vdw {params.sim_dir}/simulation/vdw.${{i}}
            sed -i "s#<state>#${{i}}#g" {params.sim_dir}/simulation/vdw.${{i}}/*/*.mdp
        done

        let max_window=${params.coul_windows}-1
        for i in $(seq 0 ${{max_window}})
        do
            cp -r {params.sim_dir}/template/coul {params.sim_dir}/simulation/coul.${{i}}
            sed -i "s#<state>#${{i}}#g" {params.sim_dir}/simulation/coul.${{i}}/*/*.mdp
        done
        '''

include: "rules/complex_fep.smk"
include: "rules/ligand_fep.smk"

rule results:
    input:
        complex_var="data/sims/{ligand}/{run}/complex/dg_results.csv",
        ligand_var="data/sims/{ligand}/{run}/ligand/dg_results.csv"
    output:
        complex_dg="data/results/{ligand}/{run}/complex/dg_results.csv",
        ligand_dg="data/results/{ligand}/{run}/ligand/dg_results.csv"
    shell:
        "cat {input.complex_var} > {output.complex_dg}"
        "cat {input.ligand_var} > {output.ligand_dg}"
