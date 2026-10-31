# Solvers

Four solvers are currently available. Two are complete sediment transport solvers that couple hydrodynamics, suspended sediment transport, and bed morphodynamics through the Exner equation. The two others are idealized morphodynamic solvers that solve only the Exner equation without hydrodynamics.

In these idealized solvers, the bedload transport rate is prescribed as a function of the bed elevation, allowing the Exner equation to be solved independently from the flow equations.

The prefix **fa** indicates that the Exner equation is solved directly on the curved finite-area mesh using the native OpenFOAM finite-area operators (`fam::` and `fac::`). Solvers without the **fa** prefix solve the Exner equation on the horizontal projection of the finite-area mesh, which corresponds to the classical formulation based on horizontal coordinates.

## Complete Sediment Transport Solvers

- **sedExnerFoam**: fully coupled solver for hydrodynamics, suspended sediment transport, bedload transport, and bed evolution.

- **faSedExnerFoam**: finite-area version of `sedExnerFoam`, where the Exner equation is solved directly on the curved bed surface. This solver is currently under development and is not yet functional.

## Idealized Morphodynamic Solvers

- **exnerFoam**: simplified solver that solves the Exner equation on the horizontal projection of the bed using an analytical bedload transport formulation.

- **faExnerFoam**: finite-area counterpart of `exnerFoam`, solving the Exner equation directly on the curved bed surface using finite-area operators.