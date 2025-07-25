# PhaseSpaceDTFE
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl/dev/)
[![Build Status](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl/actions/workflows/CI.yml?query=branch%3Amain)

The density and velocity fields of an N-body simulation is estimated with the Phase-Space Delaunay Tessellation Field Estimator implemented in Julia. This code accompanies the publication [Phase-Space Delaunay Tesselation Field Estimator](https://academic.oup.com/mnras/article/536/1/807/7915986). Please cite this publication when using the code.

## Installation

The PhaseSpaceDTFE package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("https://github.com/jfeldbrugge/PhaseSpaceDTFE.jl")
```


## Contributors
This code was written by:
* Job Feldbrugge
* Benjamin Hertzsch

We thank:
* Bram Alferink
