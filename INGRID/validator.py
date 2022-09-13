import yaml
import argparse
import itertools
from schema import Or, And, Schema, Optional
from pathlib import Path

#
# Helper functions
#
validate_is_file = lambda p: Path(p).is_file()
validate_sympy_func = lambda x: True  # Defer this to INGRID core code
is_positive = lambda x: x > 0.0
is_real = Or(int, float)
is_positive_real = And(is_real, is_positive)
is_positive_int  = And(int, is_positive)

#
# Reference values for schema
#

#
# Supported Poloidal tags
#
PATCH_CHARS = list('ABCDEFGHI')

#
# All supported Patch Tags
#
PATCH_KEYS  = ['all'] + [
            ''.join(key) for key in list(
                itertools.product(
                    list(PATCH_CHARS), [str(i) for i in range(1, 4)]
                )
            )
        ]

#
# All supported radial 'nr' grid values
#
NR_KEYS = [
        '_'.join(key) 
        for key in itertools.product(
                ['nr'], [str(i) for i in range(1, 4)] + ['default']
            )
        ]
#
# All supported poloidal 'nr' grid values
#
NP_KEYS = [
        '_'.join(key) 
        for key in itertools.product(
                ['np'], PATCH_CHARS + ['default']
            )
        ]
#
# All supported radial distribution function keys
#
RADIAL_F_KEYS = [
        '_'.join(key) 
        for key in itertools.product(
                ['radial_f'], [str(i) for i in range(1, 4)] + ['default']
            )
        ]

#
# All supported poloidal distribution function keys
#
POLOIDAL_F_KEYS = [
        '_'.join(key) 
        for key in itertools.product(
                ['poloidal_f'], PATCH_CHARS + ['default']
            )
        ]

#
# All supported target plate keys
#
PLATE_KEYS = [
    '_'.join(key)
    for key in itertools.product(
        ['plate'], [''.join(subkey) for subkey in itertools.product('WE', '12')]
    )
]

#
# Default fields
#
default = {
    'DEBUG': {
        'verbose': {
            'grid_generation': False,
            'patch_generation': False,
            'target_plates': False
        },
        'visual': {
            'SF_analysis': False,
            'find_NSEW': False,
            'gridue': False,
            'patch_map': False,
            'subgrid': False
        }
    },
    'grid_settings': {
        'guard_cell_eps': 0.001,
        'nlevs': 30,
        'num_xpt': 1,
        'psi_1': 1.1,
        'psi_2': 1.1,
        'psi_core': 0.95,
        'psi_pf_1': 0.95,
        'psi_pf_2': 1.1,
        'rmagx': 0.0,
        'zmagx': 0.0,
        'rxpt': 0.0,
        'zxpt': 0.0,
        'rxpt2': 0.0,
        'zxpt2': 0.0,
        'view_mode': 'filled',
        'grid_generation': {
            'skewness_correction': {
                'all': {
                    'active': False,
                    'resolution': 1000,
                    'theta_max': 120,
                    'theta_min': 60
                }
            },
            'np_default': 3,
            'nr_default': 2,
            'poloidal_f_default': 'x, x',
            'radial_f_default': 'x, x'
        },
        'patch_generation': {
            'core_split_point_ratio': 0.5,
            'pf_split_point_ratio': 0.5,
            'strike_pt_loc': 'limiter',
            'rmagx_shift': 0.0,
            'zmagx_shift': 0.0,
            'magx_tilt_1': 0.0,
            'magx_tilt_2': 0.0,
            'use_xpt1_W': False,
            'use_xpt1_E': False,
            'use_xpt2_W': False,
            'use_xpt2_E': False,
            'xpt1_W_tilt': -0.785398,
            'xpt1_E_tilt': 0.785398,
            'xpt2_W_tilt': -0.785398,
            'xpt2_E_tilt': 0.785398,
        }
    },
    'integrator_settings': {
        'dt': 0.01,
        'eps': 5e-5,
        'first_step': 1e-5,
        'step_ratio': 0.02,
        'tol': 5e-3,
        'max_step': 0.064
    },
    'target_plates': {
        'plate_E1': {
            'file': '',
            'rshift': 0.0,
            'zshift': 0.0,
            'auto_order': True
        },
        'plate_E2': {
            'file': '',
            'rshift': 0.0,
            'zshift': 0.0,
            'auto_order': True
        },
        'plate_W1': {
            'file': '',
            'rshift': 0.0,
            'zshift': 0.0,
            'auto_order': True
        },
        'plate_W2': {
            'file': '',
            'rshift': 0.0,
            'zshift': 0.0,
            'auto_order': True
        },
    },
    'limiter': {
        'file': '',
        'use_efit_bounds': False,
        'efit_buffer_r': 1e-2,
        'efit_buffer_z': 1e-2,
        'rshift': 0.0,
        'zshift': 0.0
    },
    'patch_data': {
        'file': '',
        'use_file': False,
        'preferences': {
            'new_file': False,
            'new_fname': ''
        }
    }
}

#
# Shortcut references to the various levels of the YAML defaults
#
default_debug_settings      = default['DEBUG']
default_grid_settings       = default['grid_settings']
default_grid_generation     = default_grid_settings['grid_generation']
default_skewness_correction = default_grid_generation['skewness_correction']
default_patch_generation    = default_grid_settings['patch_generation']
default_integrator_settings = default['integrator_settings']
default_plate_settings      = default['target_plates']
default_limiter_settings    = default['limiter']
default_patch_data          = default['patch_data']

#
# Schema for 'DEBUG' level of YAML file validation
#
DEBUG_BODY = {
    Optional('verbose', default=default['DEBUG']['verbose']): {
        Optional('grid_generation', default=False):  Or(int, bool),
        Optional('patch_generation', default=False): Or(int, bool),
        Optional('target_plates', default=False):    Or(int, bool)
    },
    Optional('visual', default=default['DEBUG']['visual']): {
        Optional('SF_analysis', default=False): Or(int, bool),
        Optional('find_NSEW', default=False):   Or(int, bool),
        Optional('gridue', default=False):      Or(int, bool),
        Optional('patch_map', default=False):   Or(int, bool),
        Optional('subgrid', default=False):     Or(int, bool)
    }
}
#
# Schema for 'skewness_correction' level of YAML file validation
#
SKEWNESS_CORRECTION_BODY = {
    Optional('active', default=False): bool,
    Optional('resolution', default=1000): is_positive_int,
    Optional('theta_max', default=120.): And(is_positive_real, lambda x: x > 90., lambda x: x < 180.),
    Optional('theta_min', default=60.): And(is_positive_real, lambda x: x < 90.)
}
#
# Schema for 'grid_generation' level of YAML file validation
#
GRID_GENERATION_BODY = {
    Optional('skewness_correction', default=default_skewness_correction): {
        Optional(patch_key, default={**default_skewness_correction['all']}): SKEWNESS_CORRECTION_BODY
        for patch_key in PATCH_KEYS
    },
    **{
        Optional(nr_key, default=3): is_positive_int for nr_key in NR_KEYS
    },
    **{
        Optional(np_key, default=3): is_positive_int for np_key in NP_KEYS
    },
    **{
        Optional(poloidal_f, default='x, x'): validate_sympy_func for poloidal_f in POLOIDAL_F_KEYS
    },
    **{
        Optional(radial_f, default='x, x'): validate_sympy_func for radial_f in RADIAL_F_KEYS
    }
}
#
# Schema for 'patch_generation' level of YAML file validation
#
PATCH_GENERATION_BODY = {
    Optional('core_split_point_ratio', default=0.5): And(is_positive_real, lambda x: x < 1.0),
    Optional('pf_split_point_ratio', default=0.5): And(is_positive_real, lambda x: x < 1.0),
    Optional('strike_pt_loc', default='limiter'): lambda x: x in ['limiter', 'target_plates'],
    Optional('rmagx_shift', default=0.0): is_real,
    Optional('zmagx_shift', default=0.0): is_real,
    Optional('magx_tilt_1', default=0.0): is_real,
    Optional('magx_tilt_2', default=0.0): is_real,
    Optional('use_xpt1_W', default=False): Or(int, bool),
    Optional('use_xpt1_E', default=False): Or(int, bool),
    Optional('use_xpt2_W', default=False): Or(int, bool),
    Optional('use_xpt2_E', default=False): Or(int, bool),
    Optional('xpt1_W_tilt'): is_real,
    Optional('xpt1_E_tilt'): is_real,
    Optional('xpt2_W_tilt'): is_real,
    Optional('xpt2_E_tilt'): is_real,
}
#
# Schema for 'grid_settings' level of YAML file validation
#
GRID_SETTINGS_BODY = {
    Optional('guard_cell_eps', default=0.001): is_positive_real,
    Optional('nlevs', default=30): is_positive_int,
    'num_xpt': And(is_positive_int, lambda x: x <= 2),
    Optional('psi_1', default=1.1): And(is_positive_real, lambda x: x > 1.0),
    Optional('psi_2', default=1.1): And(is_positive_real, lambda x: x > 1.0),
    Optional('psi_core', default=0.95): And(is_positive_real, lambda x: x < 1.0),
    Optional('psi_pf_1', default=0.95): And(is_positive_real, lambda x: x < 1.0),
    Optional('psi_pf_2', default=1.1): And(is_positive_real, lambda x: x > 1.0),
    Optional('rmagx', default=0.0): is_real,
    Optional('zmagx', default=0.0): is_real,
    Optional('rxpt', default=0.0): is_real,
    Optional('zxpt', default=0.0): is_real,
    Optional('rxpt2', default=0.0): is_real,
    Optional('zxpt2', default=0.0): is_real,
    Optional('view_mode', default='filled'): And(str, lambda s: s in ['filled', 'lines']),
    Optional('grid_generation', default=default_grid_generation): 
        GRID_GENERATION_BODY,
    Optional('patch_generation', default=default_patch_generation): 
        PATCH_GENERATION_BODY
}
#
# Schema for 'integrator_settings' level of YAML file validation
#
INTEGRATOR_SETTINGS_BODY = {
    Optional('dt', default=0.01): is_positive_real,
    Optional('eps', default=5e-5): is_positive_real,
    Optional('first_step', default=1e-5): is_positive_real,
    Optional('step_ratio', default=0.02): is_positive_real,
    Optional('tol', default=5e-3): is_positive_real,
    Optional('max_step', default=0.064): is_positive_real
}
#
# Schema for target plate body
#
PLATE_BODY = {
    'file': validate_is_file,
    Optional('rshift', default=0.0): is_real,
    Optional('zshift', default=0.0): is_real,
    Optional('auto_order', default=True): bool
}
#
# Schema for 'target_plates' level of YAML file validation
#
TARGET_PLATE_BODY = {
    Optional(plate_key): PLATE_BODY for plate_key in PLATE_KEYS
}
#
# Schema for 'limiter' level of YAML file validation
#
LIMITER_BODY = {
    Optional('file'): validate_is_file,
    Optional('use_efit_bounds', default=False): bool,
    Optional('efit_buffer_r', default=1e-2): is_positive_real,
    Optional('efit_buffer_z', default=1e-2): is_positive_real,
    Optional('rshift', default=0.0): is_real,
    Optional('zshift', default=0.0): is_real
}
#
# Schema for 'preferences' in 'patch_data'
#
PATCH_DATA_PREFERENCES = {
    Optional('new_file'): bool,
    Optional('new_fname'): str
}
#
# Schema for 'patch_data' level of YAML file validation
#
PATCH_DATA_BODY = {
    Optional('file'): validate_is_file,
    Optional('use_file'): bool,
    Optional('preferences'): PATCH_DATA_PREFERENCES
}
#
# Final schema to use for validation
#
INGRID_SCHEMA = Schema(
    {
        'eqdsk': validate_is_file,
        Optional('DEBUG', default=default_debug_settings): DEBUG_BODY,
        Optional('grid_settings', default=default_grid_settings): GRID_SETTINGS_BODY,
        Optional('integrator_settings', default=default_integrator_settings): INTEGRATOR_SETTINGS_BODY,
        Optional('target_plates', default=default_plate_settings): TARGET_PLATE_BODY,
        Optional('limiter', default=default_limiter_settings): LIMITER_BODY,
        Optional('patch_data', default=default_patch_data): PATCH_DATA_BODY
    }
)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    args   = parser.parse_args()
    input_config = yaml.safe_load(open(args.config))
    validated    = INGRID_SCHEMA.validate(input_config)
    print(yaml.dump(validated))
