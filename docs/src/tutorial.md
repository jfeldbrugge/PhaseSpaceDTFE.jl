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
