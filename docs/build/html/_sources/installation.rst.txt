*********************************
Downloading and installing INGRID
*********************************

Requirements
============

To run INGRID on your machine, ``anaconda3`` and ``setuptools`` must be installed and up to date. ``anaconda3`` installers can be found `here <https://www.anaconda.com/products/individual>`_.

.. tip:: You can create a new conda environment with the command ``conda create --name myenv`` (replace ``myenv`` with the environment name).

Once the Anaconda package manager is installed, ``setuptools`` can be added to the conda environment by running:
::

    conda install -c anaconda setuptools

To update ``setuptools`` run:
::

    pip install setuptools --upgrade


Obtaining the code
==================

Clone the INGRID repo with the command:
::

    git clone https://github.com/LLNL/INGRID.git


Installing INGRID
=================

.. warning:: Users **not** on MacOS Mojave may skip this warning. **Read on otherwise**. MacOS Mojave has issues with certain backend libraries used in INGRID. These issues have been documented by Apple. As a workaround, a specific Conda evironment has been created and must be installed by Mojave users. Navigate into the cloned repo locate the file ``conda_env.yml``. Create the mentioned Conda environment by running ``conda env create -f conda_env.yml``. Activate the new Conda environment by running ``conda activate ingrid``. When active, the terminal prompt should begin with ``(ingrid)``. The ``ingrid`` Conda environment must be active for the next section.

The user will install INGRID with the ``setup.py`` file provided in the cloned repo. Installation begins by running: 
::

    python setup.py install --user

Contents
========
Within the cloned repo are a variety of directories containing source-code, drivers, example/template files for controlling INGRID (will be discussed later), and data that the provided example-files/demos use.

We will be utilizing the directory ``example_files`` in our tutorials, and we encourage you to utilize the items in directory ``template_files`` for your own INGRID usage.
