# Import Ingrid class from ingrid.py
from INGRID.ingrid import Ingrid

"""
==================================================================
First ensure the .yml file used to generate the patch data
is on hand and populated as if operating in GUI mode.

Reminder: If 'dir_settings' are present in the .yml file,
the user must ensure that those path prefixes in combination with
file paths are valid relative where this script is executed.
==================================================================
"""
fname = './SNL/DIIID_SNL.yml'

"""
Create an Ingrid instance and load patch data.

No argument required for 'LoadPatches' since INGRID
will use the values provided from the parameter file.
"""
session = Ingrid(InputFile=fname)
session.LoadPatches()

loaded_settings = session.settings.copy()

"""
==================================================================
# Example 1: Global np = 2, nr = 3 per patch.
==================================================================
"""
updated_settings = {'np_default': 2, 'nr_default': 3}

# Overwrite all loaded settings.
session.settings['grid_settings']['grid_generation'] = updated_settings

"""
Repopulate any ommitted values in 'updated_settings' with default values

(e.g. 'poloidal_f_default', 'radial_f_default', 'DistortionCorrection')

Note we pass the entire 'settings' dictionary.
"""
session.PopulateSettings(session.settings)

# Create grid.
session.CreateSubgrid()
session.PlotSubgrid()

# Prepare gridue file and export (optional).
# session.PrepGridue()
# session.ExportGridue(fname='example_gridue_1')

"""
==================================================================
# Example 2: Update global np/nr, and refine target plates.
==================================================================
"""
updated_settings = {
                    'np_default': 5,
                    'nr_default': 5,
                    'np_A': 12,
                    'np_F': 12
                    }

# Overwrite current settings, and repopulate any missing entries.
session.settings['grid_settings']['grid_generation'] = updated_settings
session.PopulateSettings(session.settings)

# Create grid.
session.CreateSubgrid()
session.PlotSubgrid()

# Prepare gridue file and export (optional).
# session.PrepGridue()
# session.ExportGridue(fname='example_gridue_2')


"""
==================================================================
# Example 3: update settings with 'distortion_correction' globally
==================================================================
"""

# YAML entries interpreted as dictionaries. We nest accordingly.

distortion_correction_settings = {
  'all': {
        'active': True,
        'resolution': 1000,
        'theta_max': 140.0,
        'theta_min': 80.0
    }
}

# Same base use case as previous example, but with new key that
# accomadates our new 'distortion_correction' entry.
updated_settings = {
                    'np_default': 5,
                    'nr_default': 5,
                    'np_A': 12,
                    'np_F': 12,
                    'distortion_correction': None
                    }
updated_settings['distortion_correction'] = distortion_correction_settings

# Overwrite current settings, and repopulate any missing entries.
session.settings['grid_settings']['grid_generation'] = updated_settings
session.PopulateSettings(session.settings)

# Create grid.
session.CreateSubgrid()
session.PlotSubgrid()

# Prepare gridue file and export (optional).
# session.PrepGridue()
# session.ExportGridue(fname='example_gridue_3')

"""
==================================================================
# Example 4: Selectively apply 'distortion_correction' to SOL,
# and modify poloidal cells.
==================================================================
"""

# Patches of interest
desired_patches = ['A', 'C', 'E', 'F']

# Distortion correction values
resolution = 1000
theta_max = 140.0
theta_min = 100.0

# Poloidal cells
num_cells = 15

# Initial dics for us to populate.
# Important to initially turn OFF distortion correction! See below!
distortion_correction_settings = {'all': {'active': False}}
updated_settings = {'nr_default': 5, 'np_default': 5}

for tag in desired_patches:
    # Generate distortion correction entry.
    # We recall '2' corresponds to SOL in SNL cases.
    distortion_correction_settings[tag + '2'] = {
        'active': True,
        'resolution': resolution,
        'theta_max': theta_max,
        'theta_min': theta_min
    }

    # Generate poloidal cell entry
    updated_settings['np_' + tag] = num_cells

# Include distortion correction settings
updated_settings['distortion_correction'] = distortion_correction_settings

# Overwrite current settings, and repopulate any missing entries.
session.settings['grid_settings']['grid_generation'] = updated_settings
session.PopulateSettings(session.settings)

# Create grid.
session.CreateSubgrid()
session.PlotSubgrid()

# Prepare gridue file and export (optional).
# session.PrepGridue()
# session.ExportGridue(fname='example_gridue_4')
