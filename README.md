# INGRID
Interactive Grid Generator for Tokamak Boundary Region

## Installing the INGRID Conda Environment:
Dependencies required to run this code are found in the "conda_env.yml" file. The first line of the yml file sets the new environment's name. This is set to _**ingrid**_ by default.
With the Conda package manager installed on your system, navigate to the cloned Ingrid root directory and run:
```console
conda env create -f conda_env.yml
```

If no edits were made to the provided .yml file, ctivate the newly created Conda environment by running:
```console
conda activate ingrid
```

Otherwise, run:
```console
conda activate MY_ENV_NAME
```

## Example cases
Included in the **Examples** directory are test cases for Ingrid usuage.

### Running via GUI
Navigate into the **Examples** directory. With the __ingrid__ Conda environment activated, run:
```console
python GUI_startup.py
```

### Running via command line
Navigate into the **Examples** directory. With the __ingrid__ Conda environment activated, run:
```console
python CL_startup.py
```
