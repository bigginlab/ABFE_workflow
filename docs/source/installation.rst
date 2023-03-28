Installation
============

Install conda dependencies
--------------------------

.. note::

  To deal with the complexity of the dependencies is advised to use ``mamba install`` instead of ``conda install``.
  The issues cam basically from ``biosimspace`` and ``snakemake``.

.. code-block:: bash

  mamba create -n abfe python=3.9 -y
  conda activate abfe
  mamba install -c conda-forge gromacs=2022.2 parmed pdbfixer openmm openff-toolkit -y
  mamba install -c openbiosim biosimspace -y
  mamba install -c bioconda snakemake -y

pip install from the repo
-------------------------

.. code-block:: bash

  pip install git+https://github.com/bigginlab/ABFE_workflow.git

If you want to modify the code and contribute, then:

.. code-block:: bash

    git clone https://github.com/bigginlab/ABFE_workflow.git
    cd ABFE_workflow 
    conda activate abfe
    python -m pip install -e .
