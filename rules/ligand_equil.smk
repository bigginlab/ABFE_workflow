
rule run_ligand_min:
    input:
        top="data/sims/{ligand}/{run}/ligand/topology/ligand.top",
        gro="data/sims/{ligand}/{run}/ligand/topology/ligand.gro"
    params:
        run_dir="data/sims/{ligand}/{run}/ligand/equil-mdsim/enmin"
    output:
        gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/enmin/enmin.gro"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/enmin.mdp -c {input.gro} \
                   -p {input.top} -o {params.run_dir}/enmin.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/enmin -ntmpi 1 -ntomp {threads}
        '''

rule run_ligand_heat:
    input:
        top="data/sims/{ligand}/{run}/ligand/topology/ligand.top",
        gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/enmin/enmin.gro"
    params:
        run_dir="data/sims/{ligand}/{run}/ligand/equil-mdsim/nvt_heat"
    output:
        gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt="data/sims/{ligand}/{run}/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/nvt_heat.mdp -c {input.gro} -r {input.gro} \
                   -p {input.top} -o {params.run_dir}/nvt_heat.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/nvt_heat -ntomp {threads}
        '''

rule run_ligand_equil1:
    input:
        top="data/sims/{ligand}/{run}/ligand/topology/ligand.top",
        gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt="data/sims/{ligand}/{run}/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt"
    params:
        run_dir="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil1"
    output:
        gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/npt_equil1.mdp -c {input.gro} -t {input.cpt} \
                   -r {input.gro} -p {input.top} -o {params.run_dir}/npt_equil1.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/npt_equil1 -ntomp {threads}
        '''

rule run_ligand_equil2:
    input:
        top="data/sims/{ligand}/{run}/ligand/topology/ligand.top",
        gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt"
    params:
        run_dir="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil2"
    output:
        gro="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt="data/sims/{ligand}/{run}/ligand/equil-mdsim/npt_equil2/npt_equil2.cpt",
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/npt_equil2.mdp -c {input.gro} -t {input.cpt} \
                   -p {input.top} -o {params.run_dir}/npt_equil2.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/npt_equil2 -ntomp {threads}
        '''
