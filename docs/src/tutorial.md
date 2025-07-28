```@meta
CurrentModule = PhaseSpaceDTFE
```

# Tutorial

In this tutorial, we demonstrate the usage of the *PhaseSpaceDTFE* package to estimate the density and velocity fields from a *GADGET-4* simulation.

We start by importing the relevant libraries and loading the data:

```julia
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE

## set up simulation box
Ni = 64
L  = 100.
sim_box = SimBox(L, Ni)   ## need this for estimator creation

depth   = 5   # depth of estimator search tree

## load data 
function load_data(file)
    fid = h5open(file, "r")
    pos = convert.(Float64, read(fid["PartType1"]["Coordinates"]))
    vel = convert.(Float64, read(fid["PartType1"]["Velocities"]))
    ids = read(fid["PartType1"]["ParticleIDs"])
    time = read_attribute(fid["Header"], "Time")
    close(fid)

    ordering = sortperm(ids)
    return (copy(pos[:,ordering]'), copy(vel[:,ordering]'), time)
end

function load_mass(file)
    f = h5open(file, "r")
    read_attribute(f["Header"], "MassTable")[2]  # particle type 1
end

m = load_mass("../../test/data/snapshot_000.hdf5")
(coords_q, _, _) = load_data("../../test/data/snapshot_000.hdf5")
(coords_x, vels, _) = load_data("../../test/data/snapshot_002.hdf5")
```

The particle coordinates and velocities are `Float64` matrices of size `(N, 3)`. The particle mass `m` is a single `Float64` or a matrix of size `(N, 3)` for individual particle masses.

## Prequel: Delaunay Tesselation Field Estimator

Before going though through PS-DTFE method, we demonstrate the traditional DTFE method by calling the PS-DTFE code only on the final (*Eulerian*) particle positions. For details, see examples below.

```julia
## construct estimator
ps_dtfe = PS_DTFE_periodic(coords_x, coords_x, vels, m, depth, sim_box)

## evaluate density field
Range = 0:2.0:100.
density_field = [PhaseSpaceDTFE.density([L/2., y, z], ps_dtfe) for y in Range, z in Range]
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]")
```


## Phase-Space Delaunay Tessellation Field Estimator — basic implementation

We now demonstrate the use of the PS-DTFE method with the basic implementation that is suitable to simulations up to size 128^3 particles.

The first step is the construction of the estimator object from the initial (*Lagrangian*) and final (*Eulerian*) positions, `coords_q` and `coords_x`. This is only down once as a pre-computation step.

```julia
## construct estimator
ps_dtfe = PS_DTFE_periodic(coords_q, coords_x, vels, m, depth, sim_box)

## if want to ignore velocities
#ps_dtfe = PS_DTFE_periodic(coords_q, coords_x, m, depth, box)

## for further use without pre-computation, consider saving the estimator to file
#save("ps_dtfe.jld2, "ps-dtfe", ps_dtfe)
```

Note that `depth` specifies the simplex search tree depth in the estimator. Higher tree depths result faster field evaluations, but require longer construction times. It is recommended to start with `depth=5` and increase if required for high-resolution density fields.

The construction time should be of order 1-2 minutes for a 64^3 simulation at `depth=7`, or a 128^3 simulation at `depth=5`.

We now evaluate a density field with the `density()` function:

```julia
# evaluate density field
density_field = [PhaseSpaceDTFE.density([L/2., y, z], ps_dtfe) for y in Range, z in Range]
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]")
```

The corresponding number of streams is evaluated with the `numberOfStreams()` function:

```julia
nstreams_field = [numberOfStreams([L/2., y, z], ps_dtfe) for y in Range, z in Range]
heatmap(Range, Range, nstreams_field, aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), clim=(1, 7), xlabel="[Mpc]", ylabel="[Mpc]")
```

Similarly, the velocity field is evaluated with the `velocity()`-function:

```julia
vel_field = [velocity([L/2., y, z], ps_dtfe) for y in Range, z in Range]
```

In multistream regions, the `velocity()`-function returns the velocities of the individual streams (or NaN if `single_stream=true` is set in the function). To obtain the stream-mass weighted summation of the velocities, call the `velocitySum()`-function (reducing to `velocity()` in single-stream regions):


```julia
vel_field = [velocitySum([L/2., y, z], ps_dtfe) for y in Range, z in Range]
```

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