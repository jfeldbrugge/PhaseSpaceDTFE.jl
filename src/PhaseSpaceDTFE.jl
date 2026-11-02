"""
    PhaseSpaceDTFE

The density, velocity and number-of-streams fields of a three-dimensional N-body simulation, estimated with the Phase-Space Delaunay Tessellation Field Estimator (Feldbrugge 2024, MNRAS 536, 807, https://academic.oup.com/mnras/article/536/1/807/7915986) on a tessellation of the initial (Lagrangian) particle positions: the periodic Delaunay tessellation, or the Kuhn split of the particle lattice.
"""
module PhaseSpaceDTFE

    using LinearAlgebra, Random, StaticArrays, TetGen

    include("lattice.jl")
    include("delaunay.jl")
    include("index.jl")
    include("estimator.jl")
    include("raster.jl")

    export SimBox, PS_DTFE, PS_DTFE_periodic, PS_DTFE_nonperiodic
    export density, numberOfStreams, velocity, velocitySum, inSimplices
    export density_grid, numberOfStreams_grid, ps_dtfe_grid
end
