# ABFE_workflow
**WARNING**: **The repo is currently under development - Test systems do not work, code might be instabil :)**

**WARNING The repo is currently under heavy development!**

A snakemake based workflow for ABFE calculations using GMX. The workflow can be scaled on Slurm Queuing systems.
The here provided Cyclophilin D Test systems and experimental values originate from:
* [Alibay, I.; Magarkar, A.; Seeliger, D.; Biggin, P. C. Evaluating the use of absolute binding
free energy in the fragment optimisation process. Communications Chemistry 2022, 5,
105.](https://doi.org/10.1038/s42004-022-00721-4)
* [Grädler, U.; Schwarz, D.; Blaesse, M.; Leuthner, B.; Johnson, T. L.; Bernard, F.; Jiang, X.; Marx, A.; Gilardone, M.; Lemoine, H.; Roche, D.; Jorand-Lebrun, C. Discovery of novel Cyclophilin D inhibitors starting from three dimensional fragments with
millimolar potencies. Bioorganic Medicinal Chemistry Letters 2019, 29, 126717.](https://doi.org/10.1016/j.bmcl.2019.126717)

Here a visualization of the triggered process:
![](.img/full_snakemake_DAG.png)



## Usage: 
An example usage is provided with the `examples/example_execution.sh`, that uses the  `ABFE_Calculator.py` script.
If you remove the submit flag, you can a start a run, that only sets up the folder structure.

Additional script information is provided via:
```bash
  conda activate abfe

  ABFE_Calculator.py -h
```

Output:
```
> ABFE_Calculator.py -h

usage: ABFE_Calculator.py [-h] -p PROTEIN_PDB_PATH -l LIGAND_SDF_DIR -o OUTPUT_DIR_PATH [-c COFACTOR_SDF_PATH] [-nr NUMBER_OF_REPLICATES] [-njr NUMBER_OF_PARALLEL_RECEPTOR_JOBS] [-njl NUMBER_OF_PARALLEL_LIGAND_JOBS] [-ncl NUMBER_OF_CPUS_PER_LIGAND_JOB] [-submit]

optional arguments:
  -h, --help            show this help message and exit
  -p PROTEIN_PDB_PATH, --protein_pdb_path PROTEIN_PDB_PATH
                        Input protein pdb file path
  -l LIGAND_SDF_DIR, --ligand_sdf_dir LIGAND_SDF_DIR
                        Input ligand(s) sdf file path
  -o OUTPUT_DIR_PATH, --output_dir_path OUTPUT_DIR_PATH
                        Output approach folder
  -c COFACTOR_SDF_PATH, --cofactor_sdf_path COFACTOR_SDF_PATH
                        Input cofactor(s) sdf file path
  -nr NUMBER_OF_REPLICATES, --number_of_replicates NUMBER_OF_REPLICATES
                        Number of replicates
  -njr NUMBER_OF_PARALLEL_RECEPTOR_JOBS, --number_of_parallel_receptor_jobs NUMBER_OF_PARALLEL_RECEPTOR_JOBS
                        Number of jobs in parallel for receptor workflow
  -njl NUMBER_OF_PARALLEL_LIGAND_JOBS, --number_of_parallel_ligand_jobs NUMBER_OF_PARALLEL_LIGAND_JOBS
                        Number of jobs in parallel for ligand workflow
  -ncl NUMBER_OF_CPUS_PER_LIGAND_JOB, --number_of_cpus_per_ligand_job NUMBER_OF_CPUS_PER_LIGAND_JOB
                        Number of cpus per ligand job
  -submit               Will automatically submit the ABFE calculations

```

### Input
The input is suggested to be structured as follows for the commandline option:
  * Input
    * ligands
       * ligand1.sdf
       * ligand2.sdf
       * ligand3.sdf
       * ...
    * receptor.pdb

For the python call: 
 * ligand_sdfs:List[str] - paths to sdf files
 * protein_pdb_path: str - path to pdb file 

### Install:
The package can be installed like the following script:
```
  cd ABFE_workflow
  conda env create --file ./environment.yml
  conda activate abfe
  conda develop ${PWD}
```

### Running:
 if the input is set-up correctly and can be parsed, give it a run! (if you want to do the calculatios don't forget to `submit`)

Running an ABFE Campaign from Bash:
```bash
  conda activate abfe
  python ABFE_Calculator.py -p <path>/receptor.pdb \
                    -l <path>/myligands \
                    -o <path>/Out  \
                    -submit
```

Running an ABFE Campaign from Python
```python
#!/usr/bin/env python3

import glob
from abfe import calculate_abfe

ligand_sdfs = glob.glob("./myligands/*sdf")
receptor_pdb = "./receptor.pdb"
out_folder = "./Out"

calculate_abfe(protein_pdb_path=receptor_pdb, 
               ligand_sdf_path=ligand_sdfs, 
               out_root_folder_path=out_folder,
               submit=True
               )

```
