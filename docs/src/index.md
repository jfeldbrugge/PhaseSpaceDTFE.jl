```@meta
CurrentModule = PhaseSpaceDTFE
```

# Phase-Space DTFE

Documentation for [PhaseSpaceDTFE](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl).

![Density field](assets/figures/density.png)

## Installation

The PhaseSpaceDTFE package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl")
```

## Usage
Given the initial (`coords_q`) and final (`coords_x`) particle positions (and velocities (`vels`)) of an $N$-body simulation, we estimate the density, velocity and number of streams fields as follows:

```julia
using PhaseSpaceDTFE

m = 1.
depth = 5
sim_box = SimBox(L, Ni)

ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box; N_target=32)

Range = 0.:0.2:100.
coords_arr  = [[L/2., y, z] for y in Range, z in Range]
density_field = density_subbox(coords_arr, ps_dtfe_sb)
numberOfStreams_field = numberOfStreams_subbox(coords_arr, ps_dtfe_sb)
velocitySum_field = velocitySum_subbox(coords_arr, ps_dtfe_sb)
```

We strongly recommend going through the Tutorial page for more details.

## Contributors
This code was written by:
* Job Feldbrugge
* Benjamin Hertzsch

We thank:
* Bram Alferink