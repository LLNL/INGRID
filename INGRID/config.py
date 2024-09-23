from pydantic import BaseModel, Field, field_validator
from typing import Dict
import os
import re
from devtools import pprint

# Generic regex for capturing math equations.
generic_math_regex = r"([a-zA-Z])\s*,\s*([a-zA-Z0-9\.\+\-\*/\^\(\)\{\}\s]*\1[a-zA-Z0-9\.\+\-\*/\^\(\)\{\}\s]*)"

def validate_mathematical_expression(expression):
    pattern = re.compile(generic_math_regex)
    if not pattern.match(expression):
        raise ValueError("Invalid expression format")
    return expression

def validate_filepath(v):
    if not os.path.isfile(v):
        raise ValueError(f"The file {v} does not exist on the filesystem!")
    return v

class InputGeometry(BaseModel):
    file: str
    rshift: float = Field(default=0.0)
    zshift: float = Field(default=0.0)
    testing: float = Fiel
    
    @field_validator('file')
    def validate_filepaths(v):
        return validate_filepath(v)

class LimiterSettings(InputGeometry):
    use_efit_bounds: bool = Field(default=False)
    efit_buffer_r: float = Field(default=1e-2, gt=0.0)
    efit_buffer_z: float = Field(default=1e-2, gt=0.0)

class PlateSettings(InputGeometry):
    plate_E1: InputGeometry
    plate_W1: InputGeometry
    plate_E2: InputGeometry = Field(default_factory=dict)
    plate_W2: InputGeometry = Field(default_factory=dict)

class IntegratorSettings(BaseModel):
    dt: float = Field(default=0.01, gt=0.0)
    eps: float = Field(default=5e-5, gt=0.0)
    first_step: float = Field(default=1e-5, gt=0.0)
    step_ratio: float = Field(default=0.02, gt=0.0, lte=1.0)
    tol: float = Field(default=5e-3, gt=0.0)
    max_step: float = Field(default=0.064, gt=0.0, lte=1.0)

class DistortionCorrectionSettings(BaseModel):
    theta_min: float = Field(default=80.0, gt=0.0, lt=180.0)
    theta_max: float = Field(default=120.0, gt=0.0, lt=180.0)
    resolution: int = Field(default=1000, gte=1)
    activate: bool = Field(default=False)

class MeshSettings(BaseModel):
    np_default: int = Field(default=2, gte=1)
    nr_default: int = Field(default=2, gte=1)
    poloidal_f_default: str = Field(default='x, x')
    radial_f_default: str = Field(default='x, x')
    distortion_correction: Dict[str, DistortionCorrectionSettings] = Field(default_factory=lambda: {'all': DistortionCorrectionSettings()})

    @field_validator('poloidal_f_default', 'radial_f_default')
    def validate_expression(expression):
        return validate_mathematical_expression(expression)

class GenericGridSettings(BaseModel):
    psi_1: float = Field(default=1.1, gt=1.0)
    psi_2: float = Field(default=0.0)
    psi_core: float = Field(default=0.9, lt=1.0, gt=0.0)
    psi_pf_1: float = Field(default=0.9, lt=1.0, gt=0.0)
    psi_pf_2: float = Field(default=0.0)
    rmagx: float = Field(default=0.0)
    zmagx: float = Field(default=0.0)
    rxpt: float = Field(default=0.0)
    zxpt: float = Field(default=0.0)
    rxpt2: float = Field(default=0.0)
    zxpt2: float = Field(default=0.0)
    guard_cell_eps: float = Field(default=1e-3, gt=0.0)

class DefaultUserSettings(BaseModel):
    eqdsk: str
    generics: GenericGridSettings = Field(default_factory=GenericGridSettings)
    grid_generation: MeshSettings = Field(default_factory=MeshSettings)
    
    @field_validator('eqdsk')
    def validate_filepaths(fpath):
        return validate_filepath(fpath)
    

if __name__ == '__main__':
    defaults = DefaultUserSettings(eqdsk='~/clean.zsh')
    import json
    pprint(defaults.model_dump())