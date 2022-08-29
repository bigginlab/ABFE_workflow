# ABFE_workflow

A snakemake based workflow for ABFE calculations using GMX. The workflow can be scaled on Slurm Queuing systems.


## Install:
The package can be be used with the provided environment:

```
  cd ABFE_workflow
  conda env create -f ./environment.yml
  conda activate abfe
  conda develop ${PWD}
```

## Usage: 
An example usage is provided with the `calc_ABFE.py` package.

![](.img/dag-reduced.png)

