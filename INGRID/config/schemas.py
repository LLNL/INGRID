from concurrent.futures import process
from distutils.fancy_getopt import OptionDummy
import itertools
from reprlib import recursive_repr
from schema import Schema, And, Or, Use, Optional
from pathlib import Path
from INGRID.config.defaults import default_values_lookup

defaults = default_values_lookup

def optional_with_default(key, path='', delim='/'):
    """
    Helper for configuring an Optional object with
    a default matching a path
    """
    path = delim.join([path, key])
    split_list = path.split(delim)
    split_list = [c for c in split_list if c != '']
    data = {**defaults}
    for c in split_list:
        data = data[c]
    print(data)
    return Optional(c, default=data)
        
#
# General IO tests
#
def is_file(f):
    return Path(f).is_file()
def is_dir(d):
    return Path(d).is_dir()
def is_ftype(f, s):
    return Path(f).suffix == s
def is_npy(f):
    l = [is_ftype(f, suffix) for suffix in ['.np', '.npy', '.numpy']]
    return any(l)
def is_txt(f):
    l = [is_ftype(f, suffix) for suffix in ['.txt', '.text']]
    return any(l)

#
# Type tests, specialized IO tests
#
str_and_file  = And(str, is_file)
str_and_dir   = And(str, is_dir)
str_and_npy   = And(str, is_npy)
str_and_txt   = And(str, is_txt)
int_and_pos   = And(int, lambda x: x > 0)
int_and_non_neg  = And(int, lambda x: x >= 0)
float_and_pos = And(float, lambda x: x > 0.)
float_and_non_neg = And(float, lambda x: x >= 0.)
pos_real      = Or(int_and_pos, float_and_pos)
any_real      = Or(int, float)
non_neg_real  = Or(int_and_non_neg, float_and_non_neg)

#
# INGRID specific tests
#
def dist_func_check(f):
    nitems = len(f.split(',')) == 2
    pre, post = f.split(',')
    var_in_expression = pre in post
    return nitems and var_in_expression

dist_func_valid = And(str, dist_func_check)

def generate_skewness_correction_schema():

    path = 'grid_settings/grid_generation/skewness_correction/all'
    skewness_correction_body = {
        optional_with_default('theta_min', path): pos_real,
        optional_with_default('theta_max', path): pos_real,
        optional_with_default('resolution', path): int_and_pos,
        optional_with_default('active', path): bool
    }
    
    skewness_correction_schema = {
        Optional(''.join([c, i])): skewness_correction_body
        for c, i in itertools.product('ABCDEFGHI', '123')
    }

    skewness_correction_schema = {
        **{Optional('all'): skewness_correction_body},
        **skewness_correction_schema
    }
    return skewness_correction_schema

def generate_grid_generation_schema():

    path = 'grid_settings/grid_generation'
    np_nr_schema = {
        **{Optional('_'.join(['np', k])): int_and_pos for k in 'ABCDEFGHI'},
        **{Optional('_'.join(['nr', k])): int_and_pos for k in '123'},
        optional_with_default('np_default', path): int_and_pos,
        optional_with_default('nr_default', path): int_and_pos
    }

    poloidal_radial_f_schema = {
        **{
            optional_with_default(
                '_'.join(['poloidal_f', k]), path): dist_func_valid 
                for k in list('ABCDEFGHI')
        },
        **{
            optional_with_default(
                '_'.join(['radial_f', k]), path): dist_func_valid 
                for k in list('123')
        },
        optional_with_default('poloidal_f_default', path): dist_func_valid,
        optional_with_default('radial_f_default', path): dist_func_valid
    }
    grid_generation_schema = {
        **np_nr_schema,
        **poloidal_radial_f_schema,
        optional_with_default('skewness_correction', path): generate_skewness_correction_schema()
    }
    return grid_generation_schema

def generate_patch_generation_schema():

    path = 'grid_settings/patch_generation'
    patch_generation_schema = {
        'core_split_point_ratio': And(float_and_pos, lambda x: x <= 1.),
        'pf_split_point_ratio': And(float_and_pos, lambda x: x <= 1.),
        'strike_pt_loc': And(str, Use(str.lower), lambda v: v in ['limiter', 'target_plates']),
        'rmagx_shift': any_real,
        'zmagx_shift': any_real,
        'magx_tilt_1': any_real,
        'magx_tilt_2': any_real,
        'use_xpt1_W': bool,
        'use_xpt1_E': bool,
        'use_xpt2_W': bool,
        'use_xpt2_E': bool,
        'xpt1_W_tilt': any_real,
        'xpt1_E_tilt': any_real,
        'xpt2_W_tilt': any_real,
        'xpt2_E_tilt': any_real,
    }

    #
    # Wrap as optional and populate with defaults if missing
    #
    patch_generation_schema = {
        optional_with_default(k, path): v 
        for k, v in patch_generation_schema.items()
    }

    return patch_generation_schema

def generate_grid_settings_schema():
    path = 'grid_settings'
    grid_settings_schema = {
        'num_xpt': int_and_pos,
        'nlevs': int_and_pos,
        'view_mode': And(str, Use(str.lower), lambda v: v in ['filled', 'lines']),
        'psi_1': pos_real,
        'psi_core': pos_real,
        'psi_pf_1': pos_real,
        'psi_pf_2': pos_real,
        'psi_1': pos_real,
        'psi_2': pos_real,
        'rmagx': any_real,
        'zmagx': any_real,
        'rxpt': any_real,
        'zxpt': any_real,
        'rxpt2': any_real,
        'zxpt2': any_real,
        'guard_cell_eps': pos_real,
        'grid_generation': generate_grid_generation_schema(),
        'patch_generation': generate_patch_generation_schema()
    }
    #
    # Wrap as optional and populate with defaults if missing
    #
    grid_settings_schema = {
        optional_with_default(k, path): v 
        for k, v in grid_settings_schema.items()
    }
    return grid_settings_schema

def generate_integrator_settings_schema():
    path = 'integrator_settings'
    integrator_settings_schema = {
        'dt': pos_real,
        'eps': pos_real,
        'first_step': pos_real,
        'step_ratio': pos_real,
        'tol': pos_real,
        'max_step': pos_real
    }
    #
    # Wrap as optional and populate with defaults if missing
    #
    integrator_settings_schema = {
        optional_with_default(k, path): v 
        for k, v in integrator_settings_schema.items()
    }
    return integrator_settings_schema

def generate_target_plates_schema():
    path = 'target_plates/plate_E1'  # choose an arbitrary plate
    plate_body_schema = {
        'file': Or(str_and_file, str_and_dir),
        'rshift': any_real,
        'zshift': any_real,
        'auto_order': bool
    }
    #
    # Configure the defaults for a general plate entry
    #
    plate_body_schema = {
        optional_with_default(k, path): v
        for k, v in plate_body_schema.items()
    }

    #
    # Configure final 'target_plates' schema entry
    #
    path = 'target_plates'
    target_plates_schema = {
        optional_with_default(k, path): plate_body_schema
        for k in ['plate_E1', 'plate_E2', 'plate_W1', 'plate_W2']
    }

    return target_plates_schema
  

def generate_limiter_schema():

    path = 'limiter'
    limiter_schema = {
        'file': Or(str_and_file, str_and_dir),
        'use_efit_bounds': bool,
        'efit_buffer_r': pos_real,
        'efit_buffer_z': pos_real,
        'rshift': any_real,
        'zshift': any_real
    }
    #
    # Wrap as optional and populate with defaults
    #
    limiter_schema = {
        optional_with_default(k, path): v
        for k, v in limiter_schema.items()
    }
    return limiter_schema

def generate_patch_data_schema():
    #
    # Configure 'preferences' entry
    #
    path = 'patch_data/preferences'
    preferences_schema = {
        'new_file': bool,
        'new_fname': str
    }
    preferences_schema = {
        optional_with_default(k, path): v
        for k, v in preferences_schema.items()
    }

    #
    # Configure parent 'patch_data entry
    #
    path = 'patch_data'
    patch_data_schema = {
        'file': Or(str_and_file, str_and_dir),
        'use_file': bool,
        'preferences': preferences_schema
    }
    patch_data_schema = {
        optional_with_default(k, path): v
        for k, v in patch_data_schema.items()
    }
    return patch_data_schema

def generate_debug_settings_schema():
    #
    # Configure 'visual' schema
    #
    path = 'DEBUG/visual'
    visual_schema = {
            'find_NSEW': bool,
            'patch_map': bool,
            'subgrid': bool,
            'gridue': bool,
            'SF_analysis': bool
        }
    visual_schema = {
        optional_with_default(k, path): v
        for k, v in visual_schema.items()
    }

    #
    # Configure 'verbose' schema
    #
    path = 'DEBUG/verbose'
    verbose_schema = {
            'patch_generation': bool,
            'grid_generation': bool,
            'SF_analysis': bool
        }

    verbose_schema = {
        optional_with_default(k, path): v
        for k, v in verbose_schema.items()
    }

    #
    # Join together for parent 'DEBUG' schema
    #
    path = 'DEBUG'
    DEBUG_settings_schema = {
        optional_with_default('visual', path): visual_schema,
        optional_with_default('verbose', path): verbose_schema
    }
    return DEBUG_settings_schema

def generate_dir_settings_schema():
    path ='dir_settings'
    dir_settings_schema = {
        'eqdsk': str_and_dir,
        'limiter': str_and_dir,
        'patch_data': str_and_dir,
        'target_plates': str_and_dir
    }
    #
    # Wrap as optional and populate with defaults if missing
    #
    dir_settings_schema = {
        optional_with_default(k, path): v 
        for k, v in dir_settings_schema.items()
    }
    return dir_settings_schema

user_configuration_file_schema = {
    optional_with_default('eqdsk'): Or(str_and_file, str_and_dir), #Or(str_and_file, Use(resolve_path)),
    optional_with_default('dir_settings'): generate_dir_settings_schema(),
    optional_with_default('grid_settings'): generate_grid_settings_schema(),
    optional_with_default('integrator_settings'): generate_integrator_settings_schema(),
    optional_with_default('target_plates'): generate_target_plates_schema(),
    optional_with_default('limiter'):  generate_limiter_schema(),
    optional_with_default('patch_data'): generate_patch_data_schema(),
    optional_with_default('DEBUG'): generate_debug_settings_schema()
}
