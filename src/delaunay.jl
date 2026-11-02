"The 26 non-zero periodic shifts in units of the box size."
const PERIODIC_SHIFTS = [(a, b, c) for a in -1:1, b in -1:1, c in -1:1 if (a, b, c) != (0, 0, 0)]

"""
    periodic_delaunay(coords_q, L; pad = nothing, jitter = 1e-5)

Delaunay tessellation of the initial particle positions in a periodic box of size `L`, as used by `PhaseSpaceDTFE`, for lattice or glass initial conditions. The positions are tetrahedralised with TetGen together with their periodic images within a frame of width `pad * L` around the box (by default four mean particle spacings), and of every periodic tetrahedron the copy whose centroid lies in [c₀, c₀ + L)³ is kept. The reference cell is shifted by an irrational fraction c₀ of the mean particle spacing, so that no centroid of a (degenerate) lattice tessellation lies on its boundary. For degenerate configurations such as an exact lattice, where the Delaunay tessellation is not unique, the connectivity is computed from positions displaced by a fixed pseudo-random `jitter` (in mean particle spacings), the same for every periodic image of a particle, so that all images are split consistently; the volumes are computed from the exact positions. Returns the `(T, 4)` matrix of particle rows and the Lagrangian volumes of the tetrahedra. The kept tetrahedra must tile the box exactly; if their volumes do not add up to L³, the frame was too narrow and an error is thrown.
"""
function periodic_delaunay(coords_q::AbstractMatrix, L; pad = nothing, jitter = 1e-5)
    N = size(coords_q, 1)
    size(coords_q, 2) == 3 || throw(ArgumentError("coords_q must be an (N, 3) matrix"))
    q = mod.(Matrix{Float64}(coords_q), L)
    h = L / cbrt(N)
    w = pad === nothing ? 4h : pad * L
    0 < w < L / 2 || throw(ArgumentError("the frame width must lie between 0 and L/2, got $w"))

    # particles and their periodic images within the frame
    owner = collect(1:N)
    images = Vector{NTuple{3, Float64}}()
    @inbounds for p in 1:N, (a, b, c) in PERIODIC_SHIFTS
        img = (q[p, 1] + a * L, q[p, 2] + b * L, q[p, 3] + c * L)
        if all(-w <= y <= L + w for y in img)
            push!(images, img); push!(owner, p)
        end
    end
    pts = Matrix{Float64}(undef, 3, N + length(images))
    pts[:, 1:N] = q'
    @inbounds for (r, img) in enumerate(images)
        pts[1, N + r], pts[2, N + r], pts[3, N + r] = img
    end
    # the same small displacement for every image of a particle breaks the degeneracy of lattice-like configurations
    δ = jitter * h .* (rand(Random.Xoshiro(2024), 3, N) .- 0.5)

    io = TetGen.RawTetGenIO{Float64}()
    io.pointlist = pts .+ δ[:, owner]
    all_tets = tetrahedralize(io, "Q").tetrahedronlist             # (4, T), 1-based indices into pts

    c0 = (√2 - 1) / 10 * L / cbrt(N)
    keep = Int[]
    @inbounds for t in axes(all_tets, 2)
        if all(c0 <= sum(pts[d, all_tets[k, t]] for k in 1:4) / 4 < c0 + L for d in 1:3)
            push!(keep, t)
        end
    end
    tets = Matrix{Int32}(undef, length(keep), 4)
    VL = Vector{Float64}(undef, length(keep))
    @inbounds for (r, t) in enumerate(keep)
        v = ntuple(k -> all_tets[k, t], 4)
        for k in 1:4
            tets[r, k] = owner[v[k]]
        end
        E = SMatrix{3, 3, Float64}(ntuple(i -> pts[(i - 1) % 3 + 1, v[(i - 1) ÷ 3 + 2]] - pts[(i - 1) % 3 + 1, v[1]], 9))
        VL[r] = abs(det(E)) / 6
    end
    isapprox(sum(VL), L^3; rtol = 1e-8) ||
        error("the periodic Delaunay tetrahedra cover $(sum(VL) / L^3) of the box instead of 1: increase `pad`")
    return tets, VL
end



function delaunay(coords_q::AbstractMatrix; jitter = 1e-5)
    size(coords_q, 2) == 3 || throw(ArgumentError("coords_q must be an (N, 3) matrix"))
    N = size(coords_q, 1)
    N >= 4 || throw(ArgumentError("at least 4 points are needed, got $N"))
    q = Matrix{Float64}(coords_q)

    # mean particle spacing from the bounding box
    extent = [maximum(@view q[:, d]) - minimum(@view q[:, d]) for d in 1:3]
    h = cbrt(prod(extent) / N)
    h > 0 || throw(ArgumentError("the points are degenerate (zero bounding-box volume)"))

    pts = Matrix{Float64}(q')                                       # (3, N)
    # small fixed displacement breaks the degeneracy of lattice-like configurations
    δ = jitter * h .* (rand(Random.Xoshiro(2024), 3, N) .- 0.5)

    io = TetGen.RawTetGenIO{Float64}()
    io.pointlist = pts .+ δ
    all_tets = tetrahedralize(io, "Q").tetrahedronlist             # (4, T), 1-based indices into pts

    T = size(all_tets, 2)
    tets = Matrix{Int32}(undef, T, 4)
    VL = Vector{Float64}(undef, T)
    @inbounds for t in 1:T
        v = ntuple(k -> all_tets[k, t], 4)
        for k in 1:4
            tets[t, k] = v[k]
        end
        E = SMatrix{3, 3, Float64}(ntuple(i -> pts[(i - 1) % 3 + 1, v[(i - 1) ÷ 3 + 2]] - pts[(i - 1) % 3 + 1, v[1]], 9))
        VL[t] = abs(det(E)) / 6
    end

    # discard zero-volume slivers (exact positions can be coplanar even if the jittered ones are not)
    good = VL .> 1e-12 * h^3
    return tets[good, :], VL[good]
end