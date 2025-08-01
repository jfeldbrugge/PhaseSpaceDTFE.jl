using PhaseSpaceDTFE
using Test
using Suppressor

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

    m = load_mass("data/snapshot_000.hdf5")
    (coords_q, _, _) = load_data("data/snapshot_000.hdf5")
    (coords_x, vels, _) = load_data("data/snapshot_002.hdf5")


    ## Test the PS_DTFE_periodic estimator with velocity information
    ps_dtfe = PS_DTFE_periodic(coords_x, coords_x, vels, m, depth, sim_box)
    output = @capture_out PhaseSpaceDTFE.findBox([L/2. + 0.1, L/2. + 0.1, L/2. + 0.1], ps_dtfe.tree)
    @test output == "[50.0 53.125; 50.0 53.125; 50.0 53.125]\n"
    @test PhaseSpaceDTFE.density([L/2., L/2., L/2.], ps_dtfe) ≈ 6.626781014509928
    @test PhaseSpaceDTFE.numberOfStreams([L/2., L/2., L/2.], ps_dtfe) == 1
    #@test PhaseSpaceDTFE.inSimplices([L/4., L/4., L/4.], ps_dtfe) == [436370]
    println("inSimplices:", ps_dtfe.simplices[PhaseSpaceDTFE.inSimplices([L/4., L/4., L/4.], ps_dtfe)[1], :])
    @test ps_dtfe.simplices[PhaseSpaceDTFE.inSimplices([L/4., L/4., L/4.], ps_dtfe)[1], :] == [70672, 58387, 70608, 66511]
    # @test PhaseSpaceDTFE.inSimplices([L/2., L/2., L/2.], ps_dtfe) == [1686902]
    # @show PhaseSpaceDTFE.velocity([L/2. + 0.1, L/2. + 0.1, L/2. + 0.1], ps_dtfe)
    # @test PhaseSpaceDTFE.velocity([L/2. + 0.1, L/2. + 0.1, L/2. + 0.1], ps_dtfe) ≈ [-173.91874906668048 -315.74481895230974 418.73689416124455] 

    ## Test the PS_DTFE_periodic estimator without velocity information
    ps_dtfe = PS_DTFE_periodic(coords_x, coords_x, m, depth, sim_box)
    @test PhaseSpaceDTFE.density([L/2., L/2., L/2.], ps_dtfe) ≈ 6.626781014509928
    @test PhaseSpaceDTFE.numberOfStreams([L/2., L/2., L/2.], ps_dtfe) == 1


    ## Test the DTFE_periodic estimator with velocity information
    dtfe = DTFE_periodic(coords_x, vels, m, depth, sim_box)
    @test PhaseSpaceDTFE.density([L/2., L/2., L/2.], dtfe) ≈ 6.626781014509928
    @test PhaseSpaceDTFE.numberOfStreams([L/2., L/2., L/2.], dtfe) == 1

    ## Test the DTFE_periodic estimator without velocity information
    dtfe = DTFE_periodic(coords_x, m, depth, sim_box)
    @test PhaseSpaceDTFE.density([L/2., L/2., L/2.], dtfe) ≈ 6.626781014509928
    @test PhaseSpaceDTFE.numberOfStreams([L/2., L/2., L/2.], dtfe) == 1


    ## Test the PS_DTFE_subbox estimator with velocity information
    ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, vels, m, depth, sim_box; N_target=32)
    @test PhaseSpaceDTFE.density_subbox([[L/2., L/2., L/2.], [L/2., L/2., L/2.]], ps_dtfe_sb) ≈ [18.353994834770916, 18.353994834770916] 
    @test PhaseSpaceDTFE.numberOfStreams_subbox([[L/2., L/2., L/2.], [L/2., L/2., L/2.]], ps_dtfe_sb) == [3, 3]
    # @show PhaseSpaceDTFE.velocity_subbox([[L/2., L/2., L/2.], [L/2., L/2., L/2.]], ps_dtfe_sb) 
    # @test PhaseSpaceDTFE.velocitySum_subbox([[L/2., L/2., L/2.], [L/2., L/2., L/2.]], ps_dtfe_sb) ≈ [-142.50903658530848 -323.2234932699903 320.32784384240364; -142.50903658530848 -323.2234932699903 320.32784384240364]
    @test PhaseSpaceDTFE.get_subboxes(ps_dtfe_sb)[end] == [1, 1, 1]

    ## Test the PS_DTFE_subbox estimator without velocity information
    ps_dtfe_sb = ps_dtfe_subbox(coords_q, coords_x, m, depth, sim_box; N_target=32)
    @test PhaseSpaceDTFE.density_subbox([[L/2., L/2., L/2.], [L/2., L/2., L/2.]], ps_dtfe_sb) ≈ [18.353994834770916, 18.353994834770916] 
    @test PhaseSpaceDTFE.numberOfStreams_subbox([[L/2., L/2., L/2.], [L/2., L/2., L/2.]], ps_dtfe_sb) == [3, 3]

    ## Test the PS_DTFE_subbox estimator for a 3D density field
    Range  = 0:L/10:L
    coords = [[x, y, z] for x in Range, y in Range, z in Range]
    @test PhaseSpaceDTFE.density_subbox(coords, ps_dtfe_sb)[3, 7, 4] ≈ 3.1773005278948494

    rm("ps_dtfe", recursive=true)
end
