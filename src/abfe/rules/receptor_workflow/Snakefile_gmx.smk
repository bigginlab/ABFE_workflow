import abfe

approach_path = config["out_approach_path"]

# run Simulations:
include: 'simulations.smk'

# Gather Results
include: 'post_commands.smk'


rule final_result:
    input:
        result=approach_path+"/abfe_results.tsv"