```@meta
CurrentModule = PhaseSpaceDTFE
```

# Tutorial

In this tutorial, we demonstrate the usage of the *PhaseSpaceDTFE* package to estimate the density and velocity fields from a *GADGET-4* simulation.

We start by importing the relevant libraries and loading the data:

```@example tutorial1
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

```@example tutorial1
## construct estimator
ps_dtfe = PS_DTFE_periodic(coords_x, coords_x, vels, m, depth, sim_box)

## evaluate density field
Range = 0:2.0:100.
density_field = [PhaseSpaceDTFE.density([L/2., y, z], ps_dtfe) for y in Range, z in Range]
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]")
```

## Phase-Space Delaunay Tessellation Field Estimator — basic implementation

We now demonstrate the use of the PS-DTFE method with the basic implementation that is suitable to simulations up to size 128^3 particles.

The first step is the construction of the estimator object from the initial (*Lagrangian*) and final (*Eulerian*) positions, `coords_q` and `coords_x`. This is done only once as a pre-computation step.

```@example tutorial1
## construct estimator
ps_dtfe = PS_DTFE_periodic(coords_q, coords_x, vels, m, depth, sim_box)

## if want to ignore velocities
#ps_dtfe = PS_DTFE_periodic(coords_q, coords_x, m, depth, box)

## for further use without pre-computation, consider saving the estimator to file
#save("ps_dtfe.jld2, "ps-dtfe", ps_dtfe)
nothing
```

The argument `depth` specifies the simplex search tree depth in the estimator. Higher tree depths result in faster field evaluations, but require longer construction times. We recommend to start with `depth=5` and increase this if required for high-resolution density fields.

The construction time should be of order 1-2 minutes for a 64^3 simulation at `depth=7`, or a 128^3 simulation at `depth=5` on a modern computer.

We now evaluate the density field with the `density()` function:

```@example tutorial1
# evaluate density field
density_field = [PhaseSpaceDTFE.density([L/2., y, z], ps_dtfe) for y in Range, z in Range]
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]")
```

The corresponding number of streams is evaluated with the `numberOfStreams()` function:

```@example tutorial1
nstreams_field = [numberOfStreams([L/2., y, z], ps_dtfe) for y in Range, z in Range]
heatmap(Range, Range, nstreams_field, aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), clim=(1, 7), xlabel="[Mpc]", ylabel="[Mpc]")
```

Similarly, the velocity field is evaluated with the `velocity()`-function:

```@example tutorial1
velocity_field = [velocity([L/2., y, z], ps_dtfe) for y in Range, z in Range]
```

In multistream regions, the `velocity()`-function returns the velocities of the individual streams (or `NaN` if `single_stream=true` is set in the function). To obtain the stream-mass weighted summation of the velocities, call the `velocitySum()`-function, which reduces to the `velocity()`-function in single-stream regions.

```julia
velocity_field = [velocitySum([L/2., y, z], ps_dtfe) for y in Range, z in Range]
```

## The Phase-Space Delaunay Tessellation Field Estimator — subbox implementation

We now demonstrate the use of the PS-DTFE method for simulations with more than 128^3 particles.

It is not feasible to directly apply the basic PS-DTFE implementation to high-resolution simulations, as the construction of the estimator's simplex search tree would require immense working memory (> 100 GB for 256^3 particles). To circumvent this, the subbox routine internally divides the simulation box into smaller subboxes, constructs an estimator for each of these and writes the estimator to file. The user constructs the `ps_dtfe_sb` object holding the subbox references as follows:

```@example tutorial1
## construct estimators with velocities
ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box; N_target=32)

## construct estimator without velocities
# ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, m, depth, sim_box; N_target=32)

## it is recommended to save the estimator object (holding the subbox references) for further use
save("ps_dtfe_sb.jld2", "ps-dtfe-sb", ps_dtfe_sb)
ps_dtfe_sb = load("ps_dtfe_sb.jld2")["ps-dtfe-sb"]
nothing
```

The keyword argument `N_target` specifies the particle number (`N_target`^3) of the subboxes. We recommend to use the default value `N_target=128`.

For a 256^3 simulation with 8 subboxes of size `N_target=128` at `depth=5-7`, the construction time should be of order 10-30 minutes. The estimator objects will require about 20-50 GB of storage space, which can be deleted after the field evaluations (see below).

For internal efficiency, the density field is evaluated by directly passing on the list of coordinates to the `density_subbox()`-function:

```@example tutorial1
coords_arr = [[L/2., y, z] for y in Range, z in Range]
density_field = density_subbox(coords_arr, ps_dtfe_sb)
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]") 
```

The number of streams follows analogously with `numberOfStreams_subbox()`-function:

```julia
nstreams_field = numberOfStreams_subbox(coords_arr, ps_dtfe_sb)
heatmap(Range, Range, nstreams_field, aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), clim=(1, 7), xlabel="[Mpc]", ylabel="[Mpc]") 
```

Finally, the velocities are evaluated with the `velocity()`- or `velocitySum()`-function:

```julia
velocitySum_field = velocitySum_subbox(coords_arr, ps_dtfe_sb)
```

We clear the temporary files here. The user might wish to consider storing the estimators for further use.

```@example tutorial1
rm("ps_dtfe", recursive=true)
rm("ps_dtfe_sb.jld2")
```