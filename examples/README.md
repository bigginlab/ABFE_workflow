Examples
=========

This directory contains example usages of our ABFE_Workflow from different starting points and from different layers (python or CLI).

The `data` folder contains 2 example systems CyclophilinD and HSP90, with different starting points: 
- from a simple protein `.pdb` and small molecule `.sdf` files
- from already prepared GROMACS systems with `.gro` and `.top` files.


Please find the following example scripts:

CLI:
- `example_executions.sh`
- `example_execution_gmx.sh`

Python:
- `example_execution.py`

  
We also provide the expected results in the `data` directories.
- `data/CyclophilinD_min/CyclophilinD_min_abfe_results.tsv`
- `data/HSP90_gmx/HSP90_gmx_abfe_results.tsv`
