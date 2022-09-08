run_path = config["run_path"]

rule fep_gather_dGs:
    input:
        complex_var=run_path+"/complex/fep/ana/dg_results.csv",
        ligand_var=run_path+"/ligand/fep/ana/dg_results.csv"
    output:
        complex_dg=run_path+"/dg_complex_results.csv",
        ligand_dg=run_path+"/dg_ligand_results.csv"
    shell:
        '''
            cat {input.complex_var} > {output.complex_dg}
            cat {input.ligand_var} > {output.ligand_dg}
        '''
