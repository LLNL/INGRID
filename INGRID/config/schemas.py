import itertools
from schema import Schema, And, Or, Use, Optional
from pathlib import Path

def dist_func_check(f):
    nitems = len(f.split(',')) == 2
    pre, post = f.split(',')
    var_in_expression = pre in post
    return nitems and var_in_expression

#
# General IO check
#
is_file = lambda f: Path(f).is_file()
is_dir  = lambda d: Path(d).is_dir()
is_ftype = lambda f, s: Path(f).suffix == s
is_npy  = lambda f: any(is_ftype(f, '.np'), 
                        is_ftype(f, '.npy'), 
                        is_ftype(f, '.numpy'))
is_txt  = lambda f: any(is_ftype(f, '.txt'), is_ftype(f, '.text'))

#
# Type checks, specialized IO checks
#
str_and_file = And(str, is_file)
str_and_dir  = And(str, is_dir)
str_and_npy  = And(str, is_npy)
str_and_txt  = And(str, is_txt)
int_and_pos  = And(int, lambda x: x > 0)
float_and_pos = And(float, lambda x: x > 0.)
pos_real = Or(int_and_pos, float_and_pos)
any_real = Or(int, float)
#
# INGRID specifics
#
dist_func_valid = And(str, dist_func_check)

def generate_skewness_correction_schema():

    skewness_correction_body = {
        Optional('theta_min'): pos_real,
        Optional('theta_max'): pos_real,
        Optional('resolution'): int_and_pos,
        Optional('active'): bool
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
    np_nr_schema = {
        **{'_'.join(['np', k]): int_and_pos 
        for k in ['default'] + list('ABCDEFGHI')},
        **{'_'.join(['nr', k]): int_and_pos 
        for k in ['default'] + list('123')},
    }

    poloidal_radial_f_schema = {
        **{'_'.join(['poloidal_f', k]): dist_func_valid 
        for k in ['default'] + list('ABCDEFGHI')},
        **{'_'.join(['radial_f', k]): dist_func_valid 
        for k in ['default'] + list('123')},
    }

    grid_generation_schema = {
        **np_nr_schema,
        **poloidal_radial_f_schema,
        'skewness_correction': generate_skewness_correction_schema()
    }
    return grid_generation_schema

grid_settings = {
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
    'patch_generation': {
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
}

default_integrator_settings = {
    'dt': pos_real,
    'eps': pos_real,
    'first_step': pos_real,
    'step_ratio': pos_real,
    'tol': pos_real,
    'max_step': pos_real
}

default_target_plate_settings = {
    'plate_E1': {
        'file': str_and_file,
        'rshift': any_real,
        'zshift': any_real
    },

    'plate_E2': {
        'file': str_and_file,
        'rshift': any_real,
        'zshift': any_real
    },

    'plate_W1': {
        'file': str_and_file,
        'rshift': any_real,
        'zshift': any_real
    },

    'plate_W2': {
        'file': str_and_file,
        'rshift': any_real,
        'zshift': any_real
    },
}

default_limiter_settings = {
    'file': str_and_file,
    'use_efit_bounds': bool,
    'efit_buffer_r': pos_real,
    'efit_buffer_z': pos_real,
    'rshift': any_real,
    'zshift': any_real
}

default_patch_data_settings = {
    'file': str_and_file,
    'use_file': bool,
    'preferences': {
        'new_file': bool,
        'new_fname': str
    }
}

default_DEBUG_settings = {
    'visual': {
        'find_NSEW': bool,
        'patch_map': bool,
        'subgrid': bool,
        'gridue': bool,
        'SF_analysis': bool
    },

    'verbose': {
        'patch_generation': bool,
        'grid_generation': bool,
        'SF_analysis': bool
    }
}

default_dir_settings = {
    'eqdsk': str_and_dir,
    'limiter': str_and_dir,
    'patch_data': str_and_dir,
    'target_plates': str_and_dir
}

default_values_lookup = {
    'eqdsk': str_and_file,
    'dir_settings': default_dir_settings,
    'grid_settings': grid_settings,
    'integrator_settings': default_integrator_settings,
    'target_plates': default_target_plate_settings,
    'limiter':  default_limiter_settings,
    'patch_data': default_patch_data_settings,
    'DEBUG': default_DEBUG_settings
}

PlateData = {
    'plate_W1': {},
    'plate_E1': {},
    'plate_W2': {},
    'plate_E2': {}
}