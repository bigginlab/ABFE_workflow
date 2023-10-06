from abfe import scripts

run_path = config["run_path"]
num_sim_threads = config['num_sim_threads']

gromacs_run_script=config['gmx_run_kernel_path']
gromacs_cont_script=config['gmx_cont_kernel_path']


rule equil_run_complex_emin:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/topology/complex.gro",
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/emin",
        gmx_template=gromacs_run_script
    output:
        gro=run_path+"/complex/equil-mdsim/emin/emin.gro"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_emin.sh   
            chmod +x ./job_emin.sh
            ./job_emin.sh {params.nthreads} emin {input.top} {input.gro} 
        '''

rule equil_run_complex_nvt_heat:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/emin/emin.gro"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/nvt_heat",
        gmx_template=gromacs_run_script
    output:
        gro=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_nvt_heat.sh   
            chmod +x ./job_nvt_heat.sh
            ./job_nvt_heat.sh {params.nthreads} nvt_heat  {input.top} {input.gro}
        '''

rule equil_run_complex_npt_eq1:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_equil1",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_npt_eq1.sh   
            chmod +x ./job_npt_eq1.sh
            ./job_npt_eq1.sh {params.nthreads} npt_equil1 {input.top} {input.gro} {input.cpt}
        '''

rule equil_run_complex_npt_eq2:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_equil2",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_npt_eq2.sh   
            chmod +x ./job_npt_eq2.sh
            ./job_npt_eq2.sh {params.nthreads} npt_equil2 {input.top} {input.gro} {input.cpt}
        '''

rule equil_run_complex_prod:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_prod",
        gmx_template=gromacs_cont_script
    output:
        tpr=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.tpr",
        gro=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.gro",
        xtc=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.xtc"
    threads: num_sim_threads
    shell:
        '''
            set -e
            cd {params.run_dir}
            cp {params.gmx_template} ./job_npt_prod.sh   
            chmod +x ./job_npt_prod.sh
            ./job_npt_prod.sh  {params.nthreads} npt_prod {input.top} {input.gro} {input.cpt}
        '''

rule equil_run_complex_trjconv:
    input:
        tpr=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/"
    output:
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/npt_prod_center.xtc"
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
        tpr=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/npt_prod_center.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/",
        code_path = scripts.root_path
    output:
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top=run_path+"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off=run_path+"/complex/equil-mdsim/boreschcalc/dG_off.dat"
    shell:
        '''
            python {params.code_path}/preparation/generate_boresch_restraints.py --top {input.tpr} --trj {input.xtc} --outpath {params.run_dir}
        '''
