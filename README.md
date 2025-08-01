# PhaseSpaceDTFE
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl/dev/)
[![Build Status](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/github/jfeldbrugge/PhaseSpaceDTFE.jl/graph/badge.svg?token=QpOFnAiPby)](https://codecov.io/github/jfeldbrugge/PhaseSpaceDTFE.jl)


<figure>
<a href='docs/src/assets/figures/density.png'><img src='docs/src/assets/figures/density.png' width=100% /></a>
</figure>

The density and velocity fields of an N-body simulation are estimated with the Phase-Space Delaunay Tessellation Field Estimator implemented in Julia. This code accompanies the publication [Phase-Space Delaunay Tesselation Field Estimator](https://academic.oup.com/mnras/article/536/1/807/7915986). When using this code, please cite the publications listed in the citation section of the documentation.

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
Given the initial (`coords_q`) and final (`coords_x`) particle positions and velocities `vels` of an $N$-body simulation, we estimate the density, velocity and number of streams fields as follows:

```julia
import PhaseSpaceDTFE

m = 1.
depth = 5
sim_box = SimBox(L, Ni)

ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box)

Range = 0.:0.2:100.
coords_arr  = [[L/2., y, z] for y in Range, z in Range]
density_field = density_subbox(coords_arr, ps_dtfe_sb)
numberOfStreams_field = numberOfStreams_subbox(coords_arr, ps_dtfe_sb)
velocitySum_field = velocitySum_subbox(coords_arr, ps_dtfe_sb)
```

Please have a look at the Documentation for more details.

## Contributors
This code was written by:
* Job Feldbrugge ([job.feldbrugge@ed.ac.uk](mailto:job.feldbrugge@ed.ac.uk))
* Benjamin Hertzsch ([benjamin.hertzsch@ed.ac.uk](mailto:benjamin.hertzsch@ed.ac.uk))

We thank:
* Bram Alferink
