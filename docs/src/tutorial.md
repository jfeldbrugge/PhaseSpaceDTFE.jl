```@meta
CurrentModule = PhaseSpaceDTFE
```

# Tutorial

We use the PhaseSpaceDTFE package to estimate the density and velocity fields of a GADGET-4 simulation. First, we load the data

```@example tutorial1
using JLD2, Plots, HDF5, ProgressMeter, PhaseSpaceDTFE

## set up simulation box
Ni = 64
L  = 100.

depth = 5
sim_box = SimBox(L, Ni)  # note that need this custom struct for subbox

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

## The Delaunay Tessellation Field Estimator
When setting the initial positions to the final positions of the $N$-body particles, we recover the Delaunay Tessellation Field Estimator (DTFE) method
```@example tutorial1
Range = 0.:0.2:100.
coords_arr  = [[L/2., y, z] for y in Range, z in Range]

dtfe_sb = ps_dtfe_subbox(coords_x, coords_x, vels, m, depth, sim_box; N_target=32)
density_field = density_subbox(coords_arr, dtfe_sb)
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]") 
```

## The Phase-Space Delaunay Tessellation Field Estimator
For the Phase-Space Delaunay Tessellation Field Estimator (PS-DTFE), use the same routine using both the initial and final positions and velocities of the $N$-body particles.

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

We evaluate the density field 
```@example tutorial1
density_field = density_subbox(coords_arr, ps_dtfe_sb)
heatmap(Range, Range, log10.(density_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), c=:grays, xlabel="[Mpc]", ylabel="[Mpc]") 
```
the number of streams
```@example tutorial1
number_field = numberOfStreams_subbox(coords_arr, ps_dtfe_sb)
heatmap(Range, Range, log10.(number_field), aspect_ratio=:equal, xlims=(0, L), ylims=(0, L), xlabel="[Mpc]", ylabel="[Mpc]") 
```
and the mass weighted velocity field
```@example tutorial1
Range = 0.:0.2:100.
coords_arr  = [[L/2., y, z] for y in Range, z in Range]
velocitySum_field = velocitySum_subbox(coords_arr, ps_dtfe_sb)
```

Clear temporary files
```@example tutorial1
rm("ps_dtfe", recursive=true)
rm("ps_dtfe_sb.jld2")
```