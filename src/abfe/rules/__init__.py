import os

root_path = os.path.dirname(__file__)

ligand_equilibration_workflow_path = root_path + "/ligand_equilibration_workflow/Snakefile.smk"
ligand_fep_workflow_path = root_path + "/ligand_fep_workflow/Snakefile.smk"

receptor_workflow_path = root_path + "/receptor_workflow/Snakefile.smk"
receptor_gmx_workflow_path = root_path + "/receptor_workflow/Snakefile_gmx.smk"
