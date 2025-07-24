```@meta
CurrentModule = PSDTFE
```

# Tutorial

Load a GadGet-4 simulation.


```@example
using JLD2, Plots, HDF5, ProgressMeter
# , PSDTFE

# ## set up simulation box
# Ni = 64
# L  = 100.

# ## load data 
# function load_data(file)
#     fid = h5open(file, "r")
#     pos = convert.(Float64, read(fid["PartType1"]["Coordinates"]))
#     vel = convert.(Float64, read(fid["PartType1"]["Velocities"]))
#     ids = read(fid["PartType1"]["ParticleIDs"])
#     time = read_attribute(fid["Header"], "Time")
#     close(fid)

#     ordering = sortperm(ids)
#     return (copy(pos[:,ordering]'), copy(vel[:,ordering]'), time)
# end

# function load_mass(file)
#     f = h5open(file, "r")
#     read_attribute(f["Header"], "MassTable")[2]  # particle type 1
# end

# m = load_mass("../test/data/snapshot_000.hdf5")
# (coords_q, _, _) = load_data("../test/data/snapshot_000.hdf5")
# (coords_x, vels, _) = load_data("../test/data/snapshot_002.hdf5")
```