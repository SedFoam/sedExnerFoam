# SedExnerFoam

**SedExnerFoam** is an OpenFOAM solver in development aiming at the numerical study of various sediment transport and morphodynamics problems.


## Repository organization

The repository is organized as follows:
- **applications**: source code of the solvers.
- **src**: libraries required by the solver (bedload transport models, settling velocity models, morphodynamic utilities, etc.).
- **tests**: simplified test cases used to identify bugs and verify model behaviour.
- **tutorials**: example applications of `sedExnerFoam`. Tutorial cases are organized into **laminar**, **RAS** and idealized categories. The **PY** folder contains Python scripts and Jupyter notebooks for post-processing and comparison with the experimental datasets provided in the **DATA** folder.


## Installation

SedExnerFoam can be compiled with version *v2412* of [OpenFoam](https://www.openfoam.com/).

In the repository where *sedExnerFoam* code should be stored:

```bash
git clone --recurse-submodules 'https://github.com/Renaud-Matthias/sedExnerFoam/tree/main'
cd sedExnerFoam
./Allwmake
```


## Model description

The solver is derived from **pimpleFoam** and uses the **PIMPLE** algorithm to solve the incompressible Navier-Stokes equations.

Sediment transport is represented through two complementary modes:

- Suspended load, in which sediment particles are advected by the flow, settle under gravity, and are mixed by turbulence.
- Bedload transport, in which sediment particles move along the bed by rolling, sliding, and saltation.

An erosion/deposition formulation governs the exchange of sediment between the bed and the water column. The resulting sediment fluxes drive bed evolution through the **Exner** equation.

Morphological changes are handled using an Arbitrary Lagrangian-Eulerian (ALE) approach. In this framework, the sediment bed is represented as a moving boundary of the computational domain. The mesh is continuously deformed so that the bed boundary follows the evolving bed topography.

### Hydrodynamics

The hydrodynamics are described by the filtered Navier-Stokes equations for a homogeneous incompressible fluid:

$$ \frac{\partial u}{\partial t} + \nabla.(u \otimes u) = -\frac{1}{\rho_f}\nabla p + g + \nabla.(2\nu S + \tau_f) $$

$$ \nabla.u=0 $$

where $u$ is the fluid velocity, $p$ the pressure, $\rho_f$ is the fluid density, $\nu$ is the fluid kinematic viscosity, $g$ the gravitational acceleration, and $\tau_f$ is the filtering stress tensor. Depending on the selected turbulence modelling approach, $\tau_f$ may represent either the Reynolds stress tensor in Reynolds-Averaged Navier-Stokes (RANS) simulations or the subgrid-scale stress tensor in Large-Eddy Simulations (LES).

### suspended load transport:

The suspended load is described by a transport equation for the sediment volume fraction, $c_s$:

$$ \frac{\partial c_s}{\partial t} + \nabla.[(u+w_s)c_s] = \nabla.(\epsilon_s \nabla c_s) $$

where $w_s$ is the settling velocity vector, computed using empirical settling laws with optional hindered-settling corrections, and $\epsilon_s$ is the sediment turbulent diffusivity.
The exchange of sediment between the water column and the bed is represented through erosion and deposition fluxes applied as boundary conditions at the bed surface.

### bedload transport and morphodynamics

The bed elevation $z_b$, evolves in time under the influence of both bedload and suspended load transport:

$$ (1-\lambda_s)\frac{\partial z_b}{\partial t} + \nabla.q_b= D-E $$

where $\lambda_s$ is the granular material porosity, $q_b$ the bedload transport rate per unit width and $D$ and $E$ are the erosion and deposition per unit area, respectively.

The Exner equation is solved on the horizontal projection of the sediment bed using the **finite-area method**. The resulting bed elevation changes are then transferred to the three-dimensional computational mesh through an **Arbitrary Lagrangian-Eulerian (ALE)** formulation. In this approach, the sediment bed forms a moving boundary of the computational domain and the mesh is continuously deformed to follow the evolving bed topography.


## Developers

- Matthias Renaud
- Cyrille Bonamy
- Julien Chauchat
- Olivier Bertrand
- Arthur Lopez


## Acknowledgements

[OpenFoam](https://www.openfoam.com/) is the free, open source CFD software developed primarily by OpenCFD Ltd since 2004. The OpenFOAM trademark is owned by OpenCFD Ltd.
