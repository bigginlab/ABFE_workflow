from abfe import scripts

run_path = config["run_path"]
num_sim_threads = config['num_sim_threads']

rule equil_run_complex_emin:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/topology/complex.gro",
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/emin"
    output:
        gro=run_path+"/complex/equil-mdsim/emin/emin.gro"
    threads: num_sim_threads
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.run_dir}/emin.mdp -c {input.gro} \
                    -p {input.top} -o {params.run_dir}/emin.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/emin -ntomp {params.nthreads}
        '''

rule equil_run_complex_nvt_heat:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/emin/emin.gro"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/nvt_heat"
    output:
        gro=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: num_sim_threads
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.run_dir}/nvt_heat.mdp -c {input.gro} -r {input.gro} \
                    -p {input.top} -o {params.run_dir}/nvt_heat.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/nvt_heat -ntomp {params.nthreads}
        '''

rule equil_run_complex_npt_eq1:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_equil1"
    output:
        gro=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: num_sim_threads
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.run_dir}/npt_equil1.mdp -c {input.gro} -t {input.cpt} \
                    -r {input.gro} -p {input.top} -o {params.run_dir}/npt_equil1.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/npt_equil1 -ntomp {params.nthreads}
        '''

rule equil_run_complex_npt_eq2:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_equil2"
    output:
        gro=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
    threads: num_sim_threads
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.run_dir}/npt_equil2.mdp -c {input.gro} -t {input.cpt} \
                    -p {input.top} -o {params.run_dir}/npt_equil2.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/npt_equil2 -ntomp {params.nthreads}
        '''

rule equil_run_complex_prod:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_prod1"
    output:
        tpr=run_path+"/complex/equil-mdsim/npt_prod1/npt_prod1.tpr",
        gro=run_path+"/complex/equil-mdsim/npt_prod1/npt_prod1.gro",
        xtc=run_path+"/complex/equil-mdsim/npt_prod1/npt_prod1.xtc",
    threads: num_sim_threads
    shell:
        '''
            export OMP_NUM_THREADS={params.nthreads}
            gmx grompp -f {params.run_dir}/npt_prod.mdp -c {input.gro} -t {input.cpt} \
                    -p {input.top} -o {params.run_dir}/npt_prod1.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/npt_prod1 -ntomp {params.nthreads}
        '''

rule equil_run_complex_trjconv:
    input:
        tpr=run_path+"/complex/equil-mdsim/npt_prod1/npt_prod1.tpr",
        xtc=run_path+"/complex/equil-mdsim/npt_prod1/npt_prod1.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/"
    output:
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/npt_prod1_center.xtc"
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

rule equil_run_complex_get_boresch_restraints:
    input:
        tpr=run_path+"/complex/equil-mdsim/npt_prod1/npt_prod1.tpr",
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/npt_prod1_center.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/",
        code_path = scripts.root_path
    output:
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top=run_path+"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off=run_path+"/complex/equil-mdsim/boreschcalc/dG_off.dat"
    shell:
        '''
            python {params.code_path}/boresch.py --top {input.tpr} --trj {input.xtc} --outpath {params.run_dir}
        '''
