import abfe

#Final Check Job
approach_path = config["approach_path"]

# System Generation
include: 'system_generation/build_system.smk'

# run Equilibration:
include: 'equilibration/Snakefile'

# run FEP
include: 'fep/Snakefile'

# Gather Results
include: 'post_commands.smk'


rule final_result:
    input:
        result=approach_path+"/abfe_final_result.csv",,