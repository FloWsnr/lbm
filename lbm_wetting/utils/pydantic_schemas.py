from pydantic import BaseModel
from typing import Literal, List


class InputOutput(BaseModel):
    file_name: str


class Crop(BaseModel):
    x1: int
    x2: int
    y1: int
    y2: int
    z1: int
    z2: int


class Geometry(BaseModel):
    crop: Crop


##############################
###### Domain Section ########
##############################
class DomainSize(BaseModel):
    nx: int
    ny: int
    nz: int


class PeriodicBoundary(BaseModel):
    x: bool
    y: bool
    z: bool


class Domain(BaseModel):
    periodic_boundary: PeriodicBoundary
    inlet_outlet_layers: int


##############################
###### Sim ########
##############################
class FluidInit(BaseModel):
    x1: int
    x2: int
    y1: int
    y2: int
    z1: int
    z2: int


class FluidData(BaseModel):
    Gc: float
    omega_f1: int
    omega_f2: int
    G_ads_f1_s1: float
    G_ads_f1_s2: float
    G_ads_f1_s3: float
    G_ads_f1_s4: float


class Simulation(BaseModel):
    restart_sim: bool
    rho_f1: int
    rho_f2: int
    force_f1: int
    force_f2: int
    pressure_bc: bool
    minimum_radius: int
    num_pressure_steps: int
    fluid_init: Literal["drainage", "custom"]
    fluid_1_init: FluidInit
    fluid_2_init: FluidInit
    fluid_data: FluidData
    convergence: float
    convergence_iter: int
    max_iterations: int
    save_sim: bool
    save_iter: int
    gif_iter: int
    vtk_iter: int
    rho_f2_vtk: bool
    print_geom: bool
    print_stl: bool


##############################
###### Main Config ###########
##############################
class Config(BaseModel):
    input_output: InputOutput
    geometry: Geometry
    domain: Domain
    materials: List[List[int]]
    simulation: Simulation
