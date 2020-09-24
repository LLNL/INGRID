# INGRID
Interactive Grid Generator for Tokamak Boundary Region


## Installing INGRID:
After cloning the INGRID repo, the user can install INGRID via the ``setup.py`` file. This is done by running: 
```console
python setup.py install
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