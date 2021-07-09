
rule run_fep_ligand_min:
    input:
        fep_em="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/em/em.mdp"
    params:
        fep_dir="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/ligand/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/em/enmin.gro"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/em/em.mdp -c {params.top_dir}/equil.gro \
                   -p {params.top_dir}/ligand.top -o {params.fep_dir}/em/enmin.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/em/enmin -ntmpi 1 -ntomp {threads}
        '''

rule run_fep_ligand_heat:
    input:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/em/enmin.gro"
    params:
        fep_dir="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/ligand/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/nvt/nvt.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/nvt/nvt.mdp -c {input.gro} -r {input.gro} \
                   -p {params.top_dir}/ligand.top -o {params.fep_dir}/nvt/nvt.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/nvt/nvt -ntomp {threads}
        '''

rule run_fep_ligand_equil1:
    input:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        fep_dir="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/ligand/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt/npt.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/npt/npt.mdp -c {input.gro} -t {input.cpt} \
                   -r {input.gro} -p {params.top_dir}/ligand.top -o {params.fep_dir}/npt/npt.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/npt/npt -ntomp {threads}
        '''

rule run_fep_ligand_equil2:
    input:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt/npt.cpt"
    params:
        fep_dir="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/ligand/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/npt-norest/npt-norest.mdp -c {input.gro} -t {input.cpt} \
                   -p {params.top_dir}/ligand.top -o {params.fep_dir}/npt-norest/npt-norest.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/npt-norest/npt-norest -ntomp {threads}
        '''

rule run_fep_ligand_prod:
    input:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    params:
        fep_dir="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/ligand/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/prod/prod.gro",
        xvg="data/sims/{ligand}/{run}/ligand/fep/simulation/{state}/prod/prod.xvg"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/prod/prod.mdp -c {input.gro} -t {input.cpt} \
                   -p {params.top_dir}/ligand.top -o {params.fep_dir}/prod/prod.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/prod/prod -ntomp {threads}
        '''

rule gather_ligand_xvg:
    input:
        xvg_loc=expand("data/sims/{{ligand}}/{{run}}/ligand/fep/simulation/{state}/prod/prod.xvg",
                       state=ligand_windows)
        # need to add gro input
    params:
        sub_loc="data/sims/{ligand}/{run}/ligand/fep/simulation",
        vdw_max_windows=3,
        coul_max_windows=0
    output:
        xvg_dir=directory("data/sims/{ligand}/{run}/ligand/fep/simulation/xvgs")
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
        xvg_dir="data/sims/{ligand}/{run}/ligand/fep/simulation/xvgs"
    params:
        out_dir="data/sims/{ligand}/{run}/ligand/"
    output:
        complex_var="data/sims/{ligand}/{run}/ligand/dg_results.csv"
    shell:
        "python scripts/alchemlyb-analysis-ligand.py --xvgpath {input.xvg_dir} --outpath {params.out_dir}"
