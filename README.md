# INGRID
Interactive Grid Generator for Tokamak Boundary Region

## Obtaining INGRID:
Clone the INGRID repo with the command:
```console
git clone username@https://github.com/LLNL/INGRID.git IngridDir
```
where ``IngridDir`` is the name of your clone destination.

## Installation prerequisites:
To run INGRID on your machine, you must have ``anaconda3`` and ``setuptools`` installed
and up to date. 

***MacOS Mojave users please read on. Otherwise, continue to the next section.***

MacOS Mojave has issues with certain backend libraries used in INGRID. This has been documented by Apple. To remedy this, a Conda evironment has been created and must be installed by the user. Inside the cloned repo should be the file ``conda_env.yml``. Creating the Conda environment can be done by running:
```console
conda env create -f conda_env.yml
```
Activate the created Conda environment by running:
```console
conda activate ingrid
```
When active, the terminal prompt should begin with ``(ingrid)``.

## Installing INGRID
The user can install INGRID via the ``setup.py`` file (***MacOS Mojave users: make sure the ``ingrid`` conda environment is activated***). Installation can be started by running: 
```console
python setup.py install --user
```

## Example cases
Included in the ``example_files`` directory are some test cases for Ingrid usuage.

### Running GUI from drivers
Navigate into the ``drivers`` directory.
```console
python StartGUI.py
```
From here a user can select an example case from ``example_files`` or create a new case by utilizing a template file.

### Running GUI from CL
Running the gui directly from a python session can be done as follows:
```python
from INGRID import ingrid
ingrid.QuickStart()
```