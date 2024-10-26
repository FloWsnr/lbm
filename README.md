# Lattice Boltzmann - Wetting

Wetting of GDEs and other porous media using the Lattice Boltzmann method.
For the LB method, we use the [Palabos](https://gitlab.com/unigespc/palabos) library.

Additionally, we use [MPLBM-UT by Javier E. Santos](https://github.com/je-santos/MPLBM-UT) to run the simulations.
Palabos and MPLBM-UT are included as submodules in this repository.

## Installation

Clone the repository using

```
git clone git@git.rwth-aachen.de:avt.cvt/private/fluid_dynamics_models/lattice-boltzmann-method/lattice-boltzmann-wetting.git --recurse-submodules
```

Compile the Palabos and MPLBM-UT libraries using the install.sh script in the root directory of the repository.

```
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load CMake/3.24.3
```

followed by

```
chmod +x install.sh
./install.sh
```

To use the python code, you need to create a python virtual environment.
Download miniconda and install it.

Create a new environment using

```
conda create -n lbm
conda activate lbm
conda install ipykernel
conda install pydantic -c conda-forge
conda install scipy
conda install pyyaml
conda install pyvista
```
If you use another env name, make sure to change it in the `run.slrm` script.


## Usage

The simulation should be run using the `run.slrm` script. It loads the necessary modules. If run with the slurm command on the cluster, the script will automatically define the necessary variables (number of cores, etc.).
The actual simulation consists of two steps:
1. `2_phase_sim.py` creates the Palabos geometry and xml input file.
2. `twophase.cpp` runs the actual LBM simulation.

In the slrm script, specify the simulation name and directory, i.e. where the input file is located.
All other important parameters are given in the `twophase.yml` file. Change these parameters as needed.
Especially, specify the name of the input geometry file in `file_name`. Currently, numpy and vti files are supported.


## Notes

### Material encoding
Palabos and MPLBM-UT work with a different material encoding than most porous media simulations.
The materials 1,4,6 and 7 can be used for the solid phase, the pore space (wetting fluid, e.g. air) is represented by 0.
Material 2 is reserved for the internal solid material and material 3 for the non-wetting fluid.
Material 5 is used for the neutral-wet mesh, which is not used in this simulation.

Similarly, the contact angles of these materials are given by the adhesion parameters:

| Material | Adhesion parameter |
|----------|--------------------|
| 1        | G_ads_f1_s1        |
| 4        | G_ads_f1_s2        |
| 6        | G_ads_f1_s3        |
| 7        | G_ads_f1_s4        |

The adhesion parameters are given in lattice Boltzmann units, which are related to the surface tension $\sigma$ and the contact angle $\Theta$ by

$G_{ads} = \frac{2\cdot \sigma }{r}\cdot cos(\Theta )$

where $r$ is the radius of the pore space.
The **adhesion parameters** are taken from [Huang et al. 2007](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.76.066701).

Typical values (G_ads 1 = - G_ads 2) from Huang et al. 2007 are listed here:

| Case | G ads,2 | Contact angle computed from Eq | Contact angle measured |
|---|---|---|---|
| a | -0.4 | 156.4 | 158.3 |
| b | -0.3 | 133.4 | 135.1 |
| c | -0.2 | 117.3 | 117.0 |
| d | -0.1 | 103.2 | 103.2 |
| e | 0.1 | 76.8 | 75.3 |
| f | 0.2 | 62.7 | 59.5 |
| g | 0.3 | 46.6 | 40.6 |
| h | 0.4 | 23.6 | 18.9 |


### Rotation
Since MPLBM-UT simulates the drainage in the x-direction, the geometry is rotated before starting the palabos simulation.

### Run time
The runtime depends on the number of pressure steps, the number of cores and the resolution of the geometry.
For 48 cores and 12 pressure steps, a simulation with a domain of (200x200x280) takes approx. 5 hours per pressure step.

### Capillary pressure

MPLBM-UT calculates the **maximum capillary pressure** (i.e. contact angle and density)
from the selected adhesion parameters, whereby only G_ads_f1_s1 and the density of fluid 1 are used in the equation.
Additional parameters are the interparticle (cohesive) force, the minimum radius, the dissolved density and the interfacial tension.
Stable values are, $G_{c} = 0.9$, r = 3, $\rho _{d} = 0.06$ and $\sigma = 0.15$. The **unit conversion** can be done with the help of the
**Euler and Weber number**, detailed information can be found in the master thesis by Nicklas Bielfeldt 2023 titled "Multiphase simulation in porous electrodes".

$p_{c} = p_{nw} - p_{w} = \frac{\Delta \rho }{3} = \frac{2\cdot \sigma }{r}\cdot cos(\Theta )$

$cos(\Theta ) = \frac{4\cdot |G_{ads,s1}|}{G_{c}(\rho _{1} - \rho _{d} )}$