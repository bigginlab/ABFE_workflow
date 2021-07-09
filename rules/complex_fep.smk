
rule run_fep_complex_min:
    input:
        fep_em="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/em/em.mdp"
    params:
        fep_dir="data/sims/{ligand}/{run}/complex/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/complex/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/em/enmin.gro"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/em/em.mdp -c {params.top_dir}/complex.gro \
                   -p {params.top_dir}/complex.top -o {params.fep_dir}/em/enmin.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/em/enmin -ntmpi 1 -ntomp {threads}
        '''

rule run_fep_complex_heat:
    input:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/em/enmin.gro"
    params:
        fep_dir="data/sims/{ligand}/{run}/complex/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/complex/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/nvt/nvt.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/nvt/nvt.mdp -c {input.gro} -r {input.gro} \
                   -p {params.top_dir}/complex.top -o {params.fep_dir}/nvt/nvt.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/nvt/nvt -ntomp {threads}
        '''

rule run_fep_complex_equil1:
    input:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        fep_dir="data/sims/{ligand}/{run}/complex/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/complex/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt/npt.gro",
        cpt="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt/npt.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/npt/npt.mdp -c {input.gro} -t {input.cpt} \
                   -r {input.gro} -p {params.top_dir}/complex.top -o {params.fep_dir}/npt/npt.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/npt/npt -ntomp {threads}
        '''

rule run_fep_complex_equil2:
    input:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt/npt.gro",
        cpt="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt/npt.cpt"
    params:
        fep_dir="data/sims/{ligand}/{run}/complex/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/complex/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/npt-norest/npt-norest.mdp -c {input.gro} -t {input.cpt} \
                   -p {params.top_dir}/complex.top -o {params.fep_dir}/npt-norest/npt-norest.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/npt-norest/npt-norest -ntomp {threads}
        '''

rule run_fep_complex_prod:
    input:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    params:
        fep_dir="data/sims/{ligand}/{run}/complex/fep/simulation/{state}",
        top_dir="data/sims/{ligand}/{run}/complex/fep/fep-topology/"
    output:
        gro="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/prod/prod.gro",
        xvg="data/sims/{ligand}/{run}/complex/fep/simulation/{state}/prod/prod.xvg"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.fep_dir}/prod/prod.mdp -c {input.gro} -t {input.cpt} \
                   -p {params.top_dir}/complex.top -o {params.fep_dir}/prod/prod.tpr -maxwarn 2
        gmx mdrun -deffnm {params.fep_dir}/prod/prod -ntomp {threads}
        '''

rule gather_complex_xvg:
    input:
        xvg_loc=expand("data/sims/{{ligand}}/{{run}}/complex/fep/simulation/{state}/prod/prod.xvg",
                       state=all_windows)
    params:
        sub_loc="data/sims/{ligand}/{run}/complex/fep/simulation",
        rest_max_windows=1,
        vdw_max_windows=3,
        coul_max_windows=0
    output:
        xvg_dir=directory("data/sims/{ligand}/{run}/complex/fep/simulation/xvgs")
    shell:
        '''
        mkdir -p {params.sub_loc}/xvg/restraints-xvg
        mkdir -p {params.sub_loc}/xvg/vdw-xvg
        mkdir -p {params.sub_loc}/xvg/coul-xvg

        # restraints
        for i in $(seq 0 1 {params.rest_max_windows})
        do
            cp {params.sub_loc}/restraints.${{i}}/prod/prod.xvg {params.sub_loc}/xvg/restraints-xvg/dhdl.${{i}}.xvg
        done

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

rule get_dg_complex:
    input:
        xvg_dir="data/sims/{ligand}/{run}/complex/fep/simulation/xvgs"
    params:
        out_dir="data/sims/{ligand}/{run}/complex/"
    output:
        complex_var="data/sims/{ligand}/{run}/complex/dg_results.csv"
    shell:
        "python scripts/alchemlyb-analysis-complex.py --xvgpath {input.xvg_dir} --outpath {params.out_dir}"

