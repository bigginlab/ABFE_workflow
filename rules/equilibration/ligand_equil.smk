run_path = config["run_path"]
num_sim_threads = config['num_sim_threads']

rule run_ligand_emin:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/topology/ligand.gro"
    params:
        run_dir=run_path+"/ligand/equil-mdsim/emin"
    output:
        gro=run_path+"/ligand/equil-mdsim/emin/emin.gro"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.run_dir}/emin.mdp -c {input.gro} \
                    -p {input.top} -o {params.run_dir}/emin.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/emin -ntmpi {threads}
        '''

rule run_ligand_nvt_heat:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/emin/emin.gro"
    params:
        run_dir=run_path+"/ligand/equil-mdsim/nvt_heat"
    output:
        gro=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.run_dir}/nvt_heat.mdp -c {input.gro} -r {input.gro} \
                    -p {input.top} -o {params.run_dir}/nvt_heat.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/nvt_heat -ntmpi {threads}
        '''

rule run_ligand_npt_eq1:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt"
    params:
        run_dir=run_path+"/ligand/equil-mdsim/npt_equil1"
    output:
        gro=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.run_dir}/npt_equil1.mdp -c {input.gro} -t {input.cpt} \
                    -r {input.gro} -p {input.top} -o {params.run_dir}/npt_equil1.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/npt_equil1 -ntmpi {threads}
        '''

rule run_ligand_npt_eq2:
    input:
        top=run_path+"/ligand/topology/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt"
    params:
        run_dir=run_path+"/ligand/equil-mdsim/npt_equil2"
    output:
        gro=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.cpt",
    threads: num_sim_threads
    shell:
        '''
            gmx grompp -f {params.run_dir}/npt_equil2.mdp -c {input.gro} -t {input.cpt} \
                    -p {input.top} -o {params.run_dir}/npt_equil2.tpr -maxwarn 2
            gmx mdrun -deffnm {params.run_dir}/npt_equil2 -ntmpi {threads}
        '''
