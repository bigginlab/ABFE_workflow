
rule run_complex_min:
    input:
        top="data/sims/{ligand}/{run}/complex/topology/complex.top",
        gro="data/sims/{ligand}/{run}/complex/topology/complex.gro",
    params:
        run_dir="data/sims/{ligand}/{run}/complex/equil-mdsim/enmin"
    output:
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/enmin/enmin.gro"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/enmin.mdp -c {input.gro} \
                   -p {input.top} -o {params.run_dir}/enmin.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/enmin -ntmpi 1 -ntomp {threads}
        '''

rule run_complex_heat:
    input:
        top="data/sims/{ligand}/{run}/complex/topology/complex.top",
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/enmin/enmin.gro"
    params:
        run_dir="data/sims/{ligand}/{run}/complex/equil-mdsim/nvt_heat"
    output:
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt="data/sims/{ligand}/{run}/complex/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/nvt_heat.mdp -c {input.gro} -r {input.gro} \
                   -p {input.top} -o {params.run_dir}/nvt_heat.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/nvt_heat -ntomp {threads}
        '''

rule run_complex_equil1:
    input:
        top="data/sims/{ligand}/{run}/complex/topology/complex.top",
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt="data/sims/{ligand}/{run}/complex/equil-mdsim/nvt_heat/nvt_heat.cpt"
    params:
        run_dir="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil1"
    output:
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/npt_equil1.mdp -c {input.gro} -t {input.cpt} \
                   -r {input.gro} -p {input.top} -o {params.run_dir}/npt_equil1.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/npt_equil1 -ntomp {threads}
        '''

rule run_complex_equil2:
    input:
        top="data/sims/{ligand}/{run}/complex/topology/complex.top",
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil1/npt_equil1.cpt"
    params:
        run_dir="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil2"
    output:
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/npt_equil2.mdp -c {input.gro} -t {input.cpt} \
                   -p {input.top} -o {params.run_dir}/npt_equil2.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/npt_equil2 -ntomp {threads}
        '''

rule run_complex_prod:
    input:
        top="data/sims/{ligand}/{run}/complex/topology/complex.top",
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
    params:
        run_dir="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_prod1"
    output:
        tpr="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_prod1/npt_prod1.tpr",
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_prod1/npt_prod1.gro",
        xtc="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_prod1/npt_prod1.xtc",
    threads: 8
    shell:
        '''
        gmx grompp -f {params.run_dir}/npt_prod.mdp -c {input.gro} -t {input.cpt} \
                   -p {input.top} -o {params.run_dir}/npt_prod1.tpr -maxwarn 2
        gmx mdrun -deffnm {params.run_dir}/npt_prod1 -ntomp {threads}
        '''

rule run_trjconv:
    input:
        tpr="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_prod1/npt_prod1.tpr",
        xtc="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_prod1/npt_prod1.xtc"
    params:
        run_dir="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/"
    output:
        xtc="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/npt_prod1_center.xtc"
    shell:
        '''
        echo 0 | gmx trjconv -s {input.tpr} -f {input.xtc} -o {params.run_dir}/whole.xtc -pbc whole
        echo 0 | gmx trjconv -s {input.tpr} -f {params.run_dir}/whole.xtc -o {params.run_dir}/nojump.xtc -pbc nojump
        gmx trjconv -s {input.tpr} -f {params.run_dir}/nojump.xtc -o {output.xtc} -pbc mol -center -ur compact << EOF
        1
        0
        EOF
        rm {params.run_dir}/whole.xtc {params.run_dir}/nojump.xtc
        '''

rule run_boresch:
    input:
        tpr="data/sims/{ligand}/{run}/complex/equil-mdsim/npt_prod1/npt_prod1.tpr",
        xtc="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/npt_prod1_center.xtc"
    params:
        run_dir="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/"
    output:
        gro="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        dG_off="data/sims/{ligand}/{run}/complex/equil-mdsim/boreschcalc/dG_off.dat"
    shell:
        '''
        python scripts/boresch.py --top {input.tpr} --trj {input.xtc} --output {params.run_dir}
        '''
