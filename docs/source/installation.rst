Installation
============

In the root of the project there is a ``environment.yml`` with the corresponded dependencies. You just need to copy this file to you local
machine, then


.. code-block:: bash

  cd ABFE_workflow
  conda env create --file environment.yml

In case that the conda takes forever during the step `Solving environment`, consider to use `mamba` instead to solve the environment.
For that you need to install mamba in the base environment like: `conda install mamba -n base -c conda-forge`. If the mamba Installation
on the base environment is also taking long time, then

.. code-block:: bash

  conda create -c conda-forge -n mamba  mamba -y
  conda activate mamba
  mamba env create --file environment.yml
  conda activate abfe
  pip install git+https://github.com/bigginlab/ABFE_workflow.git

If you want to modify the code and contribute, then:

.. code-block:: bash

    git clone https://github.com/bigginlab/ABFE_workflow.git
    cd ABFE_workflow 
    mamba env create --file environment.yml
    conda activate abfe
    pip install git+https://github.com/bigginlab/ABFE_workflow.git
    python -m pip install -e .
