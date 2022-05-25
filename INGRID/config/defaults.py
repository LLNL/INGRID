import itertools
def generate_skewness_correction_defaults():
    common_body = {
        'theta_min': 80.0,
        'theta_max': 120.0,
        'resolution': 1000,
        'active': False
    }

    p = itertools.product('ABCDEFGHI', '123')
    entries = ['all'] + [patch_tag + i for patch_tag, i in p]

    skewness_correction_defaults = {
        k: common_body for k in entries
    }
    return skewness_correction_defaults

def generate_np_nr_defaults():
    np_default = 2
    nr_default = 2

    np_entries = [
        '_'.join(['np', k]) for k in list('ABCDEFGHI') + ['default']
    ]
    nr_entries = [
        '_'.join(['nr', k]) for k in list('123') + ['default']
    ]

    np_nr_defaults = {
        **{k: np_default for k in np_entries},
        **{k: nr_default for k in nr_entries}
    }
    return np_nr_defaults

def generate_distribution_func_defaults():
    poloidal_f_default = 'x, x'
    radial_f_default = 'x, x'

    poloidal_f_entries = [
        '_'.join(['poloidal_f', k]) 
        for k in list('ABCDEFGHI') + ['default']
    ]
    radial_f_entries = [
        '_'.join(['radial_f', k]) 
        for k in list('123') + ['default']
    ]

    distribution_func_defaults = {
        **{k: poloidal_f_default for k in poloidal_f_entries},
        **{k: radial_f_default for k in radial_f_entries}
    }
    return distribution_func_defaults

default_grid_settings = {
    'num_xpt': 1,
    'nlevs': 30,
    'view_mode': 'filled',
    'psi_1': 0.01,
    'psi_core': 0.01,
    'psi_pf_1': 0.01,
    'psi_pf_2': 0.01,
    'psi_1': 0.01,
    'psi_2': 0.01,
    'rmagx': 0.0,
    'zmagx': 0.0,
    'rxpt': 0.0,
    'zxpt': 0.0,
    'rxpt2': 0.0,
    'zxpt2': 0.0,
    'guard_cell_eps': 1.0e-3,
    'grid_generation': {
        **generate_np_nr_defaults(),
        **generate_distribution_func_defaults(),
        'skewness_correction': generate_skewness_correction_defaults(),
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
        'xpt1_W_tilt': -0.785398,  # All values of pi / 4 radians.
        'xpt1_E_tilt': 0.785398,
        'xpt2_W_tilt': -0.785398,
        'xpt2_E_tilt': 0.785398,
    }
}

default_integrator_settings = {
    'dt': 0.01,
    'eps': 5e-5,
    'first_step': 1e-5,
    'step_ratio': 0.02,
    'tol': 5e-3,
    'max_step': 0.064
}

default_target_plate_settings = {
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
}

default_limiter_settings = {
    'file': '',
    'use_efit_bounds': False,
    'efit_buffer_r': 1e-2,
    'efit_buffer_z': 1e-2,
    'rshift': 0.0,
    'zshift': 0.0
}

default_patch_data_settings = {
    'file': '',
    'use_file': False,
    'preferences': {
        'new_file': False,
        'new_fname': ''
    }
}

default_DEBUG_settings = {
    'visual': {
        'find_NSEW': False,
        'patch_map': False,
        'subgrid': False,
        'gridue': False,
        'SF_analysis': False
    },

    'verbose': {
        'patch_generation': False,
        'grid_generation': False,
        'SF_analysis': False
    }
}

default_dir_settings = {
    'eqdsk': '.',
    'limiter': '.',
    'patch_data': '.',
    'target_plates': '.'
}

PlateData = {
    'plate_W1': {},
    'plate_E1': {},
    'plate_W2': {},
    'plate_E2': {}
}

default_values_lookup = {
    'eqdsk': '',
    'dir_settings': default_dir_settings,
    'grid_settings': default_grid_settings,
    'integrator_settings': default_integrator_settings,
    'target_plates': default_target_plate_settings,
    'limiter':  default_limiter_settings,
    'patch_data': default_patch_data_settings,
    'DEBUG': default_DEBUG_settings
}
