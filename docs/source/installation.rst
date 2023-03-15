Installation
============

In the root of the project there is a ``environment.yml`` with the corresponded dependencies.


.. code-block:: bash

  cd ABFE_workflow
  conda env create --file environment.yml
  conda activate abfe
  conda develop ${PWD}

In case that the conda takes forever during the step `Solving environment`, consider to use `mamba` instead to solve the environment.
For that you need to install mamba in the base environment like: `conda install mamba -n base -c conda-forge`. If the mamba Installation
on the base environment is also taking long time, then

.. code-block:: bash

  conda create -c conda-forge -n mamba  mamba -y
  conda activate mamba
  mamba env create --file environment.yml
  conda activate abfe
  conda develop ${PWD}