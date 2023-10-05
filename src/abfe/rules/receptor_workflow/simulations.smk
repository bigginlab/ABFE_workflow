#Final Check Job
approach_path = config["out_approach_path"]

rule fep_ligand:
    input:
        in_ligand_gro= approach_path + "/{ligand_name}/input/ligand/ligand.gro",
        in_ligand_top= approach_path + "/{ligand_name}/input/ligand/ligand.top",
        output_dir= approach_path + "/{ligand_name}/{replica}/ligand"
    output:
        out_dg=approach_path+"/{ligand_name}/{replica}/ligand/fep/ana/dg_results.tsv"
    shell:
        """
            cd {input.output_dir}
            ./job_ligand.sh
        """

rule fep_complex:
    input:
        in_complex_gro = approach_path + "/{ligand_name}/input/complex/complex.gro",
        in_complex_top = approach_path + "/{ligand_name}/input/complex/complex.top",
        output_dir = approach_path + "/{ligand_name}/{replica}/complex"
    output:
        out_dg=approach_path + "/{ligand_name}/{replica}/complex/fep/ana/dg_results.tsv",

    shell:
        """
            cd {input.output_dir}
            ./job_complex.sh
        """

rule calculate_ABFE:
    input:
        out_dg_complex=approach_path + "/{ligand_name}/{replica}/complex/fep/ana/dg_results.tsv",
        out_dg_ligand=approach_path + "/{ligand_name}/{replica}/ligand/fep/ana/dg_results.tsv",
        output_dir = approach_path + "/{ligand_name}/{replica}"
    output:
        out_dg= approach_path + "/{ligand_name}/{replica}/dG_results.tsv"
    shell:
        """
            cd {input.output_dir}
            ./job.sh
        """