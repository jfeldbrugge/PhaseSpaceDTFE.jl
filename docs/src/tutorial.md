```@meta
CurrentModule = PhaseSpaceDTFE
```

# Tutorial

In this tutorial, we demonstrate the usage of the *PhaseSpaceDTFE* package to estimate the density and velocity fields from a *GADGET-4* simulation.

We start by importing the relevant libraries and loading the data:


The particle coordinates and velocities are `Float64` matrices of size `(N, 3)`. The particle mass `m` is a single `Float64` or a matrix of size `(N, 3)` for individual particle masses.

## Prequel: Delaunay Tesselation Field Estimator

Before going though through PS-DTFE method, we demonstrate the traditional DTFE method by calling the PS-DTFE code only on the final (*Eulerian*) particle positions. For details, see examples below.


## Phase-Space Delaunay Tessellation Field Estimator — basic implementation

We now demonstrate the use of the PS-DTFE method with the basic implementation that is suitable to simulations up to size 128^3 particles.

The first step is the construction of the estimator object from the initial (*Lagrangian*) and final (*Eulerian*) positions, `coords_q` and `coords_x`. This is only down once as a pre-computation step.


Note that `depth` specifies the simplex search tree depth in the estimator. Higher tree depths result faster field evaluations, but require longer construction times. It is recommended to start with `depth=5` and increase if required for high-resolution density fields.

The construction time should be of order 1-2 minutes for a 64^3 simulation at `depth=7`, or a 128^3 simulation at `depth=5`.

We now evaluate a density field with the `density()` function:


The corresponding number of streams is evaluated with the `numberOfStreams()` function:

Similarly, the velocity field is evaluated with the `velocity()`-function:

In multistream regions, the `velocity()`-function returns the velocities of the individual streams (or NaN if `single_stream=true` is set in the function). To obtain the stream-mass weighted summation of the velocities, call the `velocitySum()`-function (reducing to `velocity()` in single-stream regions):
