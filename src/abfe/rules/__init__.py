import os

root_path = os.path.dirname(__file__)

ligand_equilibration_workflow_path = root_path + "/ligand_equilibration_workflow/Snakefile.smk.smk.smk"
ligand_fep_workflow_path = root_path + "/ligand_fep_workflow/Snakefile.smk.smk.smk"

ligand_equilibration_gpu_workflow_path = root_path + "/ligand_equilibration_workflow/Snakefile_GPU"
ligand_fep_gpu_workflow_path = root_path + "/ligand_fep_workflow/Snakefile_GPU"

receptor_workflow_path = root_path + "/receptor_workflow/Snakefile.smk.smk.smk"
receptor_gmx_workflow_path = root_path + "/receptor_workflow/Snakefile_gmx.smk"
