"""
    The density and velocity fields of an N-body simulation are estimated with the Phase-Space Delaunay Tessellation Field Estimator implemented in Julia. This code accompanies the publication "Phase-Space Delaunay Tesselation Field Estimator", by Job Feldbrugge, 2024. (https://academic.oup.com/mnras/article/536/1/807/7915986).
"""
module PhaseSpaceDTFE

    using LinearAlgebra, StaticArrays, TetGen, Serialization, ProgressMeter

    include("bvh.jl")
    include("ps_dtfe.jl")
    include("ps_dtfe_subbox.jl")

    export SimBox, PS_DTFE, PS_DTFE_periodic, density, numberOfStreams, velocity, velocitySum, DTFE_periodic
    export PS_DTFE_subbox, ps_dtfe_subbox, density_subbox, numberOfStreams_subbox, velocity_subbox, velocitySum_subbox
end