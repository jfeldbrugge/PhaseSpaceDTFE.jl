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



## The Phase-Space Delaunay Tessellation Field Estimator — subbox implementation

For the Phase-Space Delaunay Tessellation Field Estimator (PS-DTFE), use the same routine using both the initial and final positions and velocities of the $N$-body particles.

```julia
## construct estimators with velocities
ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box; N_target=32)

## construct estimator without velocities
# ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, m, depth, sim_box; N_target=32)

## it is recommended to save the estimator object (holding the subbox references) for further use
save("ps_dtfe_sb.jld2", "ps-dtfe-sb", ps_dtfe_sb)
ps_dtfe_sb = load("ps_dtfe_sb.jld2")["ps-dtfe-sb"]
nothing
```

We evaluate the density field 
```julia

coords_arr = [[L/2., y, z] for y in Range, z in Range]
density_field = density_subbox(coords_arr, ps_dtfe_sb)
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]") 
```
the number of streams
```julia
number_field = numberOfStreams_subbox(coords_arr, ps_dtfe_sb)
heatmap(Range, Range, log10.(number_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), xlabel="[Mpc]", ylabel="[Mpc]") 
```
and the mass weighted velocity field
```julia
velocitySum_field = velocitySum_subbox(coords_arr, ps_dtfe_sb)
```

Clear temporary files
```julia
rm("ps_dtfe", recursive=true)
rm("ps_dtfe_sb.jld2")
```