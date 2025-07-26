# PhaseSpaceDTFE
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl/dev/)
[![Build Status](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl/actions/workflows/CI.yml?query=branch%3Amain)

<figure>
<a href='docs/figures/density.png'><img src='docs/figures/density.png' width=100% /></a>
</figure>


The density and velocity fields of an N-body simulation are estimated with the Phase-Space Delaunay Tessellation Field Estimator implemented in Julia. This code accompanies the publication [Phase-Space Delaunay Tesselation Field Estimator](https://academic.oup.com/mnras/article/536/1/807/7915986). Please cite this publication when using the code.

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
Given the initial (coords_q) and final configuration (coords_x) (and velocities (vels)) of $N$-body particles, we estimate the density/velocity/number of stream fields using the code 
```julia
import PhaseSpaceDTFE

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

## Contributors
This code was written by:
* Job Feldbrugge
* Benjamin Hertzsch

We thank:
* Bram Alferink
