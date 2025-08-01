using PhaseSpaceDTFE
using Test

@testset "PhaseSpaceDTFE.jl" begin
    # Write your tests here.
    @test π ≈ 3.14 atol=0.01 

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

    println(pwd())
    m = load_mass("data/snapshot_000.hdf5")
    (coords_q, _, _) = load_data("data/snapshot_000.hdf5")
    (coords_x, vels, _) = load_data("data/snapshot_002.hdf5")

    ## construct estimator
    ps_dtfe = PS_DTFE_periodic(coords_x, coords_x, vels, m, depth, sim_box)
    @test PhaseSpaceDTFE.density([L/2., L/2., L/2.], ps_dtfe) ≈ 6.626781014509928
    

    
end
