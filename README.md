# INGRID
Interactive Grid Generator for Tokamak Boundary Region

Dependencies required to run this code can be found in the file "conda_env.yml." 
Installing the INGRID Conda Environment:
bash
>>> conda env create -f conda_env.yml

Note: The first line of the yml file sets the new environment's name. See Conda documentation for more details.

Activating the INGRID Conda Environment:
bash
>>> conda activate MY_ENV_NAME


Gfortranwas used in the developement of this code. To compile the fortran
scripts, call the makefile by typing make into the terminal while in
the src directory. This will generate the wrappers so the python codes
can call the fortran modules.



Running the GUI:
A demonstration file has been created to illustrate INGRID capabilities. Activate the above Conda Environment and run:
python
>>> python demo.py
