using PhaseSpaceDTFE
using Test, Random, HDF5, LinearAlgebra
const P = PhaseSpaceDTFE

"Integral of the estimator over the box, Σ_t V_t ⟨ρ⟩_t, over the stored originals of the tetrahedra."
integral(est) = sum(findall(est.original)) do t
    v = est.simplices[t, :]
    abs(det(P.tet_edges(est.positions, v..., est.L))) / 6 * sum(est.rho[v]) / 4
end

"Lattice of N³ particles in a box of size L, offset by `offset` lattice spacings, rows in random order."
function lattice(N, L; offset = 0.5, rng = MersenneTwister(1))
    h = L / N
    q = [((i - 1 + offset) * h, (j - 1 + offset) * h, (k - 1 + offset) * h) for i in 1:N, j in 1:N, k in 1:N]
    Q = permutedims(reduce(hcat, [collect(p) for p in vec(q)]))
    return Q[randperm(rng, N^3), :]
end

@testset "PhaseSpaceDTFE.jl" begin

    @testset "N-body simulation" begin
        ## set up simulation box
        Ni = 64
        L  = 100.
        sim_box = SimBox(L, Ni)

        ## load data
        function load_data(file)
            fid = h5open(file, "r")
            pos = convert.(Float64, read(fid["PartType1"]["Coordinates"]))
            vel = convert.(Float64, read(fid["PartType1"]["Velocities"]))
            ids = read(fid["PartType1"]["ParticleIDs"])
            close(fid)

            ordering = sortperm(ids)
            return (copy(pos[:, ordering]'), copy(vel[:, ordering]'))
        end

        function load_mass(file)
            f = h5open(file, "r")
            read_attribute(f["Header"], "MassTable")[2]  # particle type 1
        end

        data = joinpath(@__DIR__, "data")
        m = load_mass(joinpath(data, "snapshot_000.hdf5"))
        (coords_q, _) = load_data(joinpath(data, "snapshot_000.hdf5"))
        (coords_x, vels) = load_data(joinpath(data, "snapshot_002.hdf5"))
        M = m * Ni^3
        centre = [L/2., L/2., L/2.]

        ## DTFE: the estimator on the final positions alone
        dtfe = PS_DTFE(coords_x, coords_x, vels, m, sim_box)
        @test density(centre, dtfe) ≈ 6.626781014509928        # as PhaseSpaceDTFE v1
        @test numberOfStreams(centre, dtfe) == 1
        @test integral(dtfe) ≈ M

        ## PS-DTFE
        ps_dtfe = PS_DTFE(coords_q, coords_x, vels, m, sim_box)
        @test ps_dtfe.tessellation == :delaunay
        @test density([5., 7., 4.], ps_dtfe) ≈ 1.5741253121052377
        @test numberOfStreams([5., 7., 4.], ps_dtfe) == 1
        @test density(centre, ps_dtfe) ≈ 18.340397513345412
        @test numberOfStreams(centre, ps_dtfe) == 3
        @test size(velocity(centre, ps_dtfe)) == (3, 3)
        @test all(isnan, velocity(centre, ps_dtfe, true))
        @test vec(velocitySum(centre, ps_dtfe)) ≈ [14.209633852250965, -323.07651752343077, 190.00828416161292]
        @test length(inSimplices(centre, ps_dtfe)) == 3
        @test density(centre .+ [L, -2L, 3L], ps_dtfe) ≈ density(centre, ps_dtfe)   # points are taken modulo the box
        @test integral(ps_dtfe) ≈ M

        ## the call forms of PS_DTFE, PS_DTFE_periodic and PhaseSpaceDTFE v1 give the same estimator
        p = [L/4., L/4., L/4.]
        for est in (PS_DTFE_periodic(coords_q, coords_x, vels, m, sim_box),
                    PS_DTFE_periodic(coords_q, coords_x, vels, m, 5, sim_box),
                    PS_DTFE_periodic(coords_q, coords_x, m, 5, sim_box))
            @test density(p, est) ≈ density(p, ps_dtfe)
        end

        ## density and number of streams in a slice
        Range = 0:5.:100.
        density_field = [density([L/2., y, z], ps_dtfe) for y in Range, z in Range]
        number_field = [numberOfStreams([L/2., y, z], ps_dtfe) for y in Range, z in Range]
        @test all(>(0), density_field)
        @test all(isodd, number_field)

        ## density and number of streams on a grid that extends beyond the box equal the point queries
        Range = -10:5.:100.
        density_field = density_grid(Range, Range, Range, ps_dtfe)
        number_field = numberOfStreams_grid(Range, Range, Range, ps_dtfe)
        @test density_field ≈ [density([x, y, z], ps_dtfe) for x in Range, y in Range, z in Range]
        @test number_field == [numberOfStreams([x, y, z], ps_dtfe) for x in Range, y in Range, z in Range]
        @test density_field[length(Range) ÷ 2, :, :] ≈ [density([Range[length(Range) ÷ 2], y, z], ps_dtfe) for y in Range, z in Range]
        @test density_grid(Range, Range, Range, ps_dtfe; threads = false) ≈ density_field

        ## the same grid without constructing an estimator
        ρ, n = ps_dtfe_grid(coords_q, coords_x, m, sim_box, Range, Range, Range)
        @test ρ ≈ density_field && n == number_field

        ## the Kuhn split of the Lagrangian lattice
        kuhn = PS_DTFE_periodic(coords_q, coords_x, vels, m, sim_box; tessellation = :kuhn)
        @test kuhn.tessellation == :kuhn
        @test density(centre, kuhn) ≈ 18.890736143434097
        @test numberOfStreams(centre, kuhn) == 3
        @test integral(kuhn) ≈ M
        ρ, n = ps_dtfe_grid(coords_q, coords_x, m, sim_box, Range, Range, Range; tessellation = :kuhn)
        @test ρ ≈ density_grid(Range, Range, Range, kuhn) && n == numberOfStreams_grid(Range, Range, Range, kuhn)

        ## without periodic boundaries: equal to the periodic estimator away from the faces of the box
        nonperiodic = PS_DTFE(coords_q, coords_x, vels, m, sim_box; periodicity = :nonperiodic)
        for p in ([5., 7., 4.], centre, [L/4., L/4., L/4.])
            @test density(p, nonperiodic) ≈ density(p, ps_dtfe)
            @test numberOfStreams(p, nonperiodic) == numberOfStreams(p, ps_dtfe)
            @test velocitySum(p, nonperiodic) ≈ velocitySum(p, ps_dtfe)
        end
        @test !nonperiodic.periodic && ps_dtfe.periodic
        @test density(centre, PS_DTFE_nonperiodic(coords_q, coords_x, m, sim_box)) ≈ density(centre, nonperiodic)
        # grids equal the point queries, also where matter has crossed a face of the box
        @test density([95., 5., 5.], nonperiodic) == 0 < density([95., 5., 5.], ps_dtfe)
        @test density_grid(Range, Range, Range, nonperiodic) ≈ [density([x, y, z], nonperiodic) for x in Range, y in Range, z in Range]
        @test numberOfStreams_grid(Range, Range, Range, nonperiodic) == [numberOfStreams([x, y, z], nonperiodic) for x in Range, y in Range, z in Range]
        @test_throws ArgumentError PS_DTFE(coords_q, coords_x, m, sim_box; periodicity = :open)
        @test_throws ArgumentError PS_DTFE(coords_q, coords_x, m, sim_box; periodicity = :nonperiodic, tessellation = :kuhn)
    end

    N, L, m = 12, 1.0, 2.0
    box = SimBox(L, N)
    q = lattice(N, L)
    ρ̄ = m * N^3 / L^3
    pts = [[0.013, 0.977, 0.5013], [0.5013, 0.4987, 0.2511], [0.731, 0.264, 0.919], [0.998, 0.002, 0.4997], [1.37, -0.41, 2.08]]

    @testset "lattice and Kuhn triangulation" begin
        site = P.lattice_sites(q .+ 0.2 * (L / N) .* (rand(MersenneTwister(2), N^3, 3) .- 0.5), box)
        @test sort(vec(site)) == 1:N^3
        tets = P.lattice_tetrahedra(site)
        @test size(tets) == (6N^3, 4)
        @test all(count(==(i), tets) == 24 for i in 1:N^3)
        @test_throws ArgumentError P.lattice_sites(q[1:end - 1, :], box)
    end

    @testset "exact for translations and shears, also after several wraps" begin
        shift = [2.3, -3.6, 0.7] .* L
        maps = [("translation", q .+ shift'),
                ("shear 0.2 L", hcat(q[:, 1] .+ 0.2L .* sin.(2π .* q[:, 3] ./ L), q[:, 2], q[:, 3])),
                ("shear 0.6 L", hcat(q[:, 1], q[:, 2] .+ 0.6L .* sin.(2π .* q[:, 1] ./ L), q[:, 3]))]   # edges stay < L/2
        for (name, x) in maps, wrapped in (false, true)
            X = wrapped ? mod.(x, L) : x
            est = PS_DTFE_periodic(q, X, m, box; tessellation = :kuhn)
            @test all(density(p, est) ≈ ρ̄ for p in pts)
            @test all(numberOfStreams(p, est) == 1 for p in pts)
        end
    end

    @testset "velocities" begin
        v = hcat(q[:, 2], 2 .* q[:, 3], -q[:, 1])                 # linear velocity field (y, 2z, -x), single stream
        p = [0.3013, 0.4521, 0.6017]
        for tessellation in (:delaunay, :kuhn)
            est = PS_DTFE_periodic(q, q, v, m, box; tessellation)
            @test vec(velocity(p, est)) ≈ [p[2], 2p[3], -p[1]]
            @test vec(velocitySum(p, est)) ≈ [p[2], 2p[3], -p[1]]
            @test size(inSimplices(p, est)) == (1,)
            @test velocity(p, est, true) ≈ reshape([p[2], 2p[3], -p[1]], 1, 3)
        end
        @test_throws ArgumentError velocity(p, PS_DTFE_periodic(q, q, m, box))
    end

    @testset "multistream: grids equal point queries, mass is conserved" begin
        # shell crossing around x = L/2, where 1 + 2πA cos(2πq/L) < 0 for A > 1/2π
        x = hcat(q[:, 1] .+ 0.25L .* sin.(2π .* q[:, 1] ./ L), q[:, 2] .+ 0.05L .* sin.(2π .* q[:, 3] ./ L), q[:, 3])
        est = PS_DTFE_periodic(q, x, m, box; tessellation = :kuhn)
        @test maximum(numberOfStreams(p, est) for p in ([0.49, 0.3, 0.7], [0.51, 0.62, 0.2])) == 3
        xs, ys, zs = range(-0.11, 1.07, length = 23), range(0.003, 0.997, length = 17), range(-0.5, 0.6, length = 13)
        ρg, ng = P.rasterise(est, xs, ys, zs)
        @test ρg ≈ [density([a, b, c], est) for a in xs, b in ys, c in zs]
        @test ng == [numberOfStreams([a, b, c], est) for a in xs, b in ys, c in zs]
        @test density_grid(xs, ys, zs, est; threads = false) ≈ ρg
        zplane = range(0.37, 0.37, length = 1)                   # a single plane
        @test density_grid(xs, ys, zplane, est)[:, :, 1] ≈ [density([a, b, 0.37], est) for a in xs, b in ys]
        ρs, ns = ps_dtfe_grid(q, x, m, box, xs, ys, zs; tessellation = :kuhn)
        @test ρs ≈ ρg && ns == ng
        @test integral(est) ≈ m * N^3
        @test_throws ArgumentError density_grid(reverse(xs), ys, zs, est)
    end

    @testset "Delaunay tessellation, for glass initial conditions" begin
        rng = MersenneTwister(3)
        glass = q .+ 0.3 * (L / N) .* (rand(rng, N^3, 3) .- 0.5)            # irregular, non-lattice initial positions
        tets, VL = P.periodic_delaunay(glass, L)
        @test sum(VL) ≈ L^3                                                # the tetrahedra tile the periodic box
        @test sort(unique(vec(tets))) == 1:N^3
        @test sum(P.periodic_delaunay(q, L)[2]) ≈ L^3                      # degenerate case: an exact lattice
        for (name, x) in (("translation", glass .+ [2.3, -3.6, 0.7]' .* L), ("translation of the lattice", q .+ [0.4, 1.3, -2.2]' .* L))
            Q = startswith(name, "translation of") ? q : glass
            est = PS_DTFE_periodic(Q, mod.(x, L), m, box; vertex_density = :lagrangian)
            @test est.tessellation == :delaunay
            @test all(density(p, est) ≈ ρ̄ for p in pts)                    # exact for a uniform initial density
            @test all(numberOfStreams(p, est) == 1 for p in pts)
        end
        # with :mass the vertex densities scatter for irregular stars, but the total mass is conserved
        est = PS_DTFE_periodic(glass, glass, m, box)
        @test !all(density(p, est) ≈ ρ̄ for p in pts)
        @test integral(est) ≈ m * N^3
        # velocities, grids and the streaming grid with the Delaunay tessellation
        v = hcat(glass[:, 2], 2 .* glass[:, 3], -glass[:, 1])
        est = PS_DTFE_periodic(glass, glass, v, m, box; vertex_density = :lagrangian)
        p = [0.3013, 0.4521, 0.6017]
        @test vec(velocity(p, est)) ≈ [p[2], 2p[3], -p[1]]
        x = hcat(glass[:, 1] .+ 0.25L .* sin.(2π .* glass[:, 1] ./ L), glass[:, 2], glass[:, 3])
        est = PS_DTFE_periodic(glass, x, m, box; vertex_density = :lagrangian)
        xs, ys, zs = range(-0.11, 1.07, length = 19), range(0.003, 0.997, length = 13), range(-0.5, 0.6, length = 11)
        ρg, ng = P.rasterise(est, xs, ys, zs)
        @test ρg ≈ [density([a, b, c], est) for a in xs, b in ys, c in zs]
        @test ng == [numberOfStreams([a, b, c], est) for a in xs, b in ys, c in zs]
        ρs, ns = ps_dtfe_grid(glass, x, m, box, xs, ys, zs; vertex_density = :lagrangian)
        @test ρs ≈ ρg && ns == ng
        @test maximum(ng) >= 3
        @test_throws ArgumentError PS_DTFE_periodic(glass, x, m, box; tessellation = :voronoi)
    end

    @testset "without periodic boundaries" begin
        glass = q .+ 0.3 * (L / N) .* (rand(MersenneTwister(5), N^3, 3) .- 0.5)
        tets, VL = P.delaunay(glass)
        tets, VL, _ = P.drop_long_tets(glass, tets, VL, 3L / N)
        inner = [[0.5013, 0.4987, 0.2511], [0.431, 0.564, 0.619], [0.31, 0.72, 0.45]]       # away from the faces of the box
        # under a translation, :lagrangian gives the mean density over the volume covered by the tetrahedra
        est = PS_DTFE_nonperiodic(glass, glass .+ 0.01L, m, box; vertex_density = :lagrangian)
        @test all(density(p, est) ≈ m * N^3 / sum(VL) for p in inner)
        @test all(numberOfStreams(p, est) == 1 for p in inner)
        # displacements are taken modulo the box
        @test density(inner[1], PS_DTFE_nonperiodic(glass, mod.(glass .+ 0.01L, L), m, box)) ≈
              density(inner[1], PS_DTFE_nonperiodic(glass, glass .+ 0.01L, m, box))
        v = hcat(glass[:, 2], 2 .* glass[:, 3], -glass[:, 1])
        est = PS_DTFE(glass, glass, v, m, box; periodicity = :nonperiodic)
        @test vec(velocity(inner[2], est)) ≈ [inner[2][2], 2inner[2][3], -inner[2][1]]
        # grids equal the point queries for a flow that crosses the faces of the box
        x = hcat(glass[:, 1] .+ 0.25L .* sin.(2π .* glass[:, 1] ./ L) .+ 0.1L, glass[:, 2] .- 0.07L, glass[:, 3])
        est = PS_DTFE_nonperiodic(glass, x, m, box)
        xs, ys, zs = range(-0.11, 1.07, length = 19), range(0.003, 0.997, length = 13), range(-0.5, 0.6, length = 11)
        ρg, ng = P.rasterise(est, xs, ys, zs)
        @test ρg ≈ [density([a, b, c], est) for a in xs, b in ys, c in zs]
        @test ng == [numberOfStreams([a, b, c], est) for a in xs, b in ys, c in zs]
        @test maximum(ng) >= 3
    end

    @testset "particles of different mass" begin
        rng = MersenneTwister(4)
        glass = q .+ 0.3 * (L / N) .* (rand(rng, N^3, 3) .- 0.5)
        x = hcat(q[:, 1] .+ 0.25L .* sin.(2π .* q[:, 1] ./ L), q[:, 2], q[:, 3])
        xg = hcat(glass[:, 1] .+ 0.25L .* sin.(2π .* glass[:, 1] ./ L), glass[:, 2], glass[:, 3])
        masses = m .* (0.5 .+ rand(rng, N^3))
        for (Q, X, tess) in ((q, x, :kuhn), (glass, xg, :delaunay)), mode in (:mass, :lagrangian)
            est = PS_DTFE_periodic(Q, X, masses, box; tessellation = tess, vertex_density = mode)
            @test integral(est) ≈ sum(masses)                              # the total mass is conserved
        end
        # :mass is linear in the masses
        e1 = PS_DTFE_periodic(q, x, masses, box)
        e2 = PS_DTFE_periodic(q, x, 3 .* masses, box)
        @test all(density(p, e2) ≈ 3density(p, e1) for p in pts)
        # zoom-like initial conditions: masses proportional to the Lagrangian star volume represent a uniform density,
        # which :lagrangian recovers exactly under a translation
        tets, VL = P.periodic_delaunay(glass, L)
        W = zeros(N^3); for t in axes(tets, 1), c in tets[t, :]; W[c] += VL[t]; end
        zoom = sum(masses) .* W ./ sum(W)
        est = PS_DTFE_periodic(glass, glass .+ 1.3L, zoom, box; vertex_density = :lagrangian)
        @test all(density(p, est) ≈ sum(zoom) / L^3 for p in pts)
        xs, ys, zs = range(0, 1, length = 9), range(0, 1, length = 7), range(0, 1, length = 5)
        ρs, ns = ps_dtfe_grid(q, x, masses, box, xs, ys, zs)
        @test ρs ≈ density_grid(xs, ys, zs, e1)
        @test_throws ArgumentError PS_DTFE_periodic(q, x, masses[1:end - 1], box)
        @test_throws ArgumentError PS_DTFE_periodic(q, x, -1.0, box)
    end

    @testset "PhaseSpaceDTFE v1 call signatures" begin
        x = q .+ 0.01L
        @test density(pts[2], PS_DTFE_periodic(q, x, m, 5, box; tessellation = :kuhn)) ≈ ρ̄
        @test density(pts[2], PS_DTFE_periodic(q, x, zero(x), m, 5, box; tessellation = :kuhn)) ≈ ρ̄
        @test density(pts[2], PS_DTFE_periodic(q, x, m, box; tessellation = :kuhn, G = 5)) ≈ ρ̄
        @test density(pts[2], PS_DTFE_periodic(q, x, m, 5, box; vertex_density = :lagrangian)) ≈ ρ̄
    end
end
