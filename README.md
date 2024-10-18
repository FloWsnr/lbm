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

Compile the Palabos library using the install.sh script in the root directory of the repository.

```
bash install.sh
```

To use the python code, you need to create a python virtual environment.
Download miniconda and install it.

Create a new environment using

```
conda create -n lbm
conda activate lbm
conda install ipykernel
conda install pydantic -c conda-forge
pip install -e mplbm-ut-mirror/src/python

```


## Notes

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



To install MPLBM-UT, first load the following modules using

```
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GCCcore/.11.3.0
module load Python/3.10.4
module load CMake/3.24.3
module load matplotlib/3.5.2
```

followed by

```
chmod +x install.sh
./install.sh
```

_Since changes to the RWTH cluster occur frequently, the commands **module unload, list, spider, save** and **restore** may be useful.
In case the installation is not sucsefull, try logging in to the FastX3 visual interface and check the **Run file as program** box.


## Simulation

In a first step, navigate e.g. to **/examples/unsteady_rel_perm/input** to add in the geometry as a NumPy array.
Then change the path and dimensions in the **input.yml** file, including (nx;ny;nz) if the geometry is to be cropped.
Since MPLBM-UT simulates the drainage in the x-direction **swap xz: True** should be set, in order to rotate the volume
(e.g. 500;500;280 results in 280;200;200). A detailed description of the general inputs can be found under [MPLBM-UT/examples](https://github.com/je-santos/MPLBM-UT/tree/master/examples).\
In order for the simulation to run stably, do not change any of the settings except **num procs**, **num pressure steps** and **G_ads_f1_s1 (...) s4**.
When using **48 cores and 12 pressure steps**, the computation time is approx. 5 h/pressure step. For more accurate results,
the number of pressure levels should be increased. The **adhesion parameters** are taken from [Huang et al. 2007](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.76.066701),
whereby negative values represent contact angles >90°. The table show the voxel labels and their corresponding value in the simulation.

| Voxel | Represents            |
|-------|-----------------------|
| 0     | Fluid 2 (wetting)     |
| 1     | G_ads_f1_s1 (solid)   |
| 2     | Inside solids         |
| 3     | Fluid 1 (non-wetting) |
| 4     | G_ads_f1_s2 (solid)   |
| 5     | Neutral-wet mesh      |
| 6     | G_ads_f1_s3 (solid)   |
| 7     | G_ads_f1_s4 (solid)   |

MPLBM-UT calculates the **maximum capillary pressure** (i.e. contact angle and density)
from the selected adhesion parameters, whereby only G_ads_f1_s1 and the density of fluid 1 are used in the equation.
Additional parameters are the interparticle (cohesive) force, the minimum radius, the dissolved density and the interfacial tension.
Stable values are, $G_{c} = 0.9$, r = 3, $\rho _{d} = 0.06$ and $\sigma = 0.15$. The **unit conversion** can be done with the help of the
**Euler and Weber number**, detailed information can be found in the master thesis by Nicklas Bielfeldt 2023 titled
"Multiphase simulation in porous electrodes".

$p_{c} = p_{nw} - p_{w} = \frac{\Delta \rho }{3} = \frac{2\cdot \sigma }{r}\cdot cos(\Theta )$

$cos(\Theta ) = \frac{4\cdot |G_{ads,s1}|}{G_{c}(\rho _{1} - \rho _{d} )}$

To simulate on the RWTH cluster, a batch script has to be submitted.
For general information, see [wiki](https://help.itc.rwth-aachen.de/service/rhr4fjjutttf/article/13ace46cfbb84e92a64c1361e0e4c104/).
A **run.slrm** may be structured as follows, and can be submitted with the command **sbatch run.slrm**.
Additional important commands are **squeue -u, r_wlm_usage -p** and **scancel <job ID>**.

```
#!/usr/bin/zsh

### Task name
#SBATCH --account=<foo>
#SBATCH --job-name=<foo>

### Output file
#SBATCH --output=<foo>.%J.out

### Maximum runtime
#SBATCH --time=24:00:00

### Ask for less than 4 GB memory per (CPU) task = MPI rank
#SBATCH --mem-per-cpu=2000

### Start a parallel job for a distributed-memory system on several nodes
#SBATCH --nodes=1

### For MPI, use one task per CPU
#SBATCH --cpus-per-task=1

### Number of processes per node (48 cores per node)
#SBATCH --ntasks-per-node=48

### Mail notification configuration
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<foo>@rwth-aachen.de

module purge
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GCCcore/.11.3.0
module load Python/3.10.4
module load CMake/3.24.3
module load matplotlib/3.5.2

python3 2_phase_sim.py
```

## Results

When the simulation is running, a **/tmp** folder is created in which the new structure, the .vti files, the capillary pressure
(lattice Boltzmann unit) and the saturation is stored. The Chan-Shen iterates for each pressure step until the density difference
reaches the convergence threshold.

## Hints

### Storage

The storage capacity on $WORK (250G) is not limited, but obtaining more space there is significantly more challenging and limited. More space available on $HPCWORK (1000G).
On the other hand, $HOME (150G) is reserved only for **truly** valuable data, which means any data that cannot be reproduced within a couple of weeks.
This is due to its higher cost compared to other storage options.

### Parallelization

Scalability and efficiency are not necessarily easy to achieve optimally, but in real-life scenarios,
they can be reduced to some **fundamental questions and rules** in the first approximation.

- [ ] Q1: Can the dataset be processed within the available memory? Or is it necessary to use multiple nodes to have enough RAM?
- [ ] Q2: Is the time to solution acceptable? (i.e., can one endure waiting for the results, as "faster is better" but usually more expensive).

The first rule of parallelization states that more cores lead to more overhead.
At a certain point, adding more resources may not result in any acceleration, and in fact, it could even slow down the computation.
The typical scalability curve is U-shaped. The second rule for codes that can run on multiple nodes is to use one node as long as the data fits within it.
Otherwise, try to minimize the number of nodes, as communication between nodes is slower than communication within a single node.
Keep in mind that larger jobs consume more resources; your quota is not unlimited, and a 4-node job running for 24 hours consumes 4600 core-hours or 10% of your quota in a single job.

**It is a good idea to do some performance testing**. Gradually increase the domain size and observe the scaling behavior.
If you double one dimension, does it take approximately 2x, 4x, 8x longer? This way, you can estimate how long the "large" dataset will run.
Also, test how scaling behaves concerning the number of cores. Choose a reasonably sized dataset and run it on, say, 12, 24, and 48 cores on a single node.
Observe how increasing the number of cores affects the computation. If the code scales well up to the full node, then compare runs with 1, 2, 3, and 4 nodes together.
Initially, try to stay within 1 node (to avoid network overhead) as long as the job runtimes remain acceptable
(a few days are acceptable; jobs can run up to 7 days, and for longer computations, consider dividing them into parts
using "restart files" and running them as [chained jobs](https://hpc-wiki.info/hpc/SLURM#Array_and_Chain_Jobs).