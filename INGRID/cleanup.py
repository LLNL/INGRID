import sys
sys.path.append('ingrid')
from ingrid import Ingrid
from pathlib import Path
import yaml as yml
import numpy as np

for p in Path('.').rglob('*.yml'):
    print(f'Checking file {p}...')
    d = Ingrid.ReadYamlFile(p)

    val = {'view_mode': 'filled'}

    if d.get('grid_settings') is None:
        continue

    d['grid_settings'].update(val)

    with open(p, 'w') as f:
        yml.dump(d, f)
        print(f'# Updating file: {p}')

    d = Ingrid.ReadYamlFile(p)

    val = {'all': {'theta_min': 80.0, 'theta_max': 120.0, 'resolution': 1000, 'active': False}}

    if d.get('grid_settings') is None:
        continue

    if d['grid_settings']['grid_generation'].get('CorrectDistortion') is not None:
        val = {k.lower(): v for k, v in d['grid_settings']['grid_generation'].items()}
        d['grid_settings'].pop('CorrectDistortion')

    if d['grid_settings']['grid_generation'].get('CorrectDistortion') is None:
        d['grid_settings']['grid_generation'].update({'distortion_correction': val})

        with open(p, 'w') as f:
            yml.dump(d, f)
            print(f'# Updating file: {p}')

    d = Ingrid.ReadYamlFile(p)

    val = 1.0e-3

    if d.get('grid_settings') is None:
        continue

    if d['grid_settings'].get('guard_cell_eps') is None:
        d['grid_settings'].update({'guard_cell_eps': val})

        with open(p, 'w') as f:
            yml.dump(d, f)
            print(f'# Updating file: {p}')

    d = Ingrid.ReadYamlFile(p)

    val = 'target_plates'

    if d.get('limiter') is not None \
            and d['limiter'].get('use_limiter') is not None:

        if d['limiter']['use_limiter'] is True:
            val = 'limiter'

        d['limiter'].pop('use_limiter')

    if d.get('grid_settings') is None:
        continue

    if d['grid_settings']['patch_generation'].get('strike_geometry') is None:
        d['grid_settings']['patch_generation'].update({'strike_geometry': val})

        with open(p, 'w') as f:
            yml.dump(d, f)
            print(f'# Updating file: {p}')

for p in Path('.').rglob('*.npy'):
    print(f'Checking file {p}...')
    data = np.load(p, allow_pickle=True)

    if data[3][-1].get('patchName') is None:
        continue

    for raw_patch in data[2:]:
        raw_patch[-1].update({'patch_name': raw_patch[-1]['patchName']})
        raw_patch[-1].pop('patchName')
        raw_patch[-1].update({'plate_patch': raw_patch[-1]['platePatch']})
        raw_patch[-1].pop('platePatch')
        raw_patch[-1].update({'plate_location': raw_patch[-1]['plateLocation']})
        raw_patch[-1].pop('plateLocation')

    np.save(p, data)
    print(f'# Updating file {p}')
