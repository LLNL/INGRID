# Import Ingrid class from ingrid.py
from INGRID.ingrid import Ingrid

paths = {
    'eqdsk': '../data/SNL/DIII-D/',
    'target_plates': '../data/SNL/DIII-D/'
}

fnames = {
    'eqdsk': 'neqdsk',
    'plate_W1': 'd3d_itp.txt',
    'plate_E1': 'd3d_otp.txt',
}

# magx (r, z) can be obtained from neqdsk
grid_settings = {
    'rmagx': 1.5,
    'zmagx': 0.0,
    'rxpt': 1.0,
    'zxpt': -1.0,
    'psi_1': 1.06,  # SOL
    'psi_pf_1': 0.975,  # PF
    'psi_core': 0.95  # CORE
}


# Create an Ingrid instance
session = Ingrid()

"""
Preparing data for call to StartSetup()
"""

# Optional: Set paths to search for file names. Default = '.'
session.settings['dir_settings']['eqdsk'] = paths['eqdsk']
session.settings['eqdsk'] = fnames['eqdsk']

# SNL mode. Default value 1
grid_settings['num_xpt'] = 1

# Use target plates rather than limiter
(session.settings
    ['grid_settings']
    ['patch_generation']
    ['strike_pt_loc']) = 'target_plates'


'''
Provide target plate information.
    Lower Single Null:
      plate_W1 <--> inner target plate
      plate_E1 <--> outer target plate
'''
session.settings['dir_settings']['target_plates'] = paths['target_plates']
session.settings['target_plates']['plate_W1']['file'] = fnames['plate_W1']
session.settings['target_plates']['plate_E1']['file'] = fnames['plate_E1']

# zshift for target plate coordinates. useful for this case.
session.settings['target_plates']['plate_W1']['zshift'] = -1.6
session.settings['target_plates']['plate_E1']['zshift'] = -1.6

'''
Optional: rshift for target plates
    # session.settings['target_plates']['plate_W1']['rshift'] = 0.5
    # session.settings['target_plates']['plate_E1']['rshift'] = 0.0
'''

# Provide approximate (r, z) coordinate for magnetic axis
session.settings['grid_settings']['rmagx'] = grid_settings['rmagx']
session.settings['grid_settings']['zmagx'] = grid_settings['zmagx']

# Provide approximate (r, z) coordinate for primary x-point
session.settings['grid_settings']['rxpt'] = grid_settings['rxpt']
session.settings['grid_settings']['zxpt'] = grid_settings['zxpt']

# Provide psi-boundary values for SOL, PF, and CORE, respectively
session.settings['grid_settings']['psi_1'] = grid_settings['psi_1']
session.settings['grid_settings']['psi_pf_1'] = grid_settings['psi_pf_1']
session.settings['grid_settings']['psi_core'] = grid_settings['psi_core']

import pdb
pdb.set_trace()

session.PopulateSettings(session.settings)

# Loads EFIT, refines xpts and magx, calcs derivatives...
session.StartSetup()
session.SetGeometry({'limiter': session.settings['limiter']})
session.ShowSetup()

session.AnalyzeTopology()
session.ConstructPatches()
session.PlotPatches()

session.SaveSettingsFile('noninteractive_SNL_Example.yml')
session.SavePatches('noninteractive_SNL_patches.npy')
