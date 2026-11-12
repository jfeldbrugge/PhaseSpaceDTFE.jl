"The 26 non-zero periodic shifts in units of the box size."
const PERIODIC_SHIFTS = [(a, b, c) for a in -1:1, b in -1:1, c in -1:1 if (a, b, c) != (0, 0, 0)]

"Positions 1:4 in increasing order of the values (a, b, c, d), by insertion sort without allocating."
@inline function sort4(a, b, c, d)
    v = (a, b, c, d)
    o = (1, 2, 3, 4)
    for i in 2:4
        for j in i:-1:2                                  # a separate loop: `break` leaves only this one
            v[o[j - 1]] <= v[o[j]] && break
            o = Base.setindex(Base.setindex(o, o[j], j - 1), o[j - 1], j)
        end
    end
    return o
end

"""
    canonical_tetrahedra(ids, pos)

Tetrahedra in a canonical form, independent of the order in which TetGen lists them: the corners of every tetrahedron
are sorted by particle id, the rows are sorted, and the Lagrangian volumes are computed from the sorted corners. `ids`
is a `(T, 4)` matrix of particle ids and `pos(t, k)` the position (an `SVector`) of corner k of tetrahedron t. Returns
the `(T, 4)` matrix of ids and the volumes.
"""
function canonical_tetrahedra(ids::AbstractMatrix, pos)
    T = size(ids, 1)
    rows = Vector{NTuple{4, Int32}}(undef, T)
    VL = Vector{Float64}(undef, T)
    Threads.@threads for t in 1:T
        @inbounds begin
            o = sort4(ids[t, 1], ids[t, 2], ids[t, 3], ids[t, 4])
            rows[t] = (Int32(ids[t, o[1]]), Int32(ids[t, o[2]]), Int32(ids[t, o[3]]), Int32(ids[t, o[4]]))
            p0 = pos(t, o[1])
            E = hcat(pos(t, o[2]) - p0, pos(t, o[3]) - p0, pos(t, o[4]) - p0)
            VL[t] = abs(det(E)) / 6
        end
    end
    perm = sortperm(rows)
    tets = Matrix{Int32}(undef, T, 4)
    @inbounds for (r, t) in enumerate(perm)
        tets[r, 1], tets[r, 2], tets[r, 3], tets[r, 4] = rows[t]
    end
    return tets, VL[perm]
end

"""
    periodic_delaunay(coords_q, L; pad = nothing, jitter = 1e-5, slabs = nothing)

Delaunay tessellation of the initial particle positions in a periodic box of size `L`, as used by `PhaseSpaceDTFE`, for lattice or glass initial conditions. The positions are tetrahedralised with TetGen together with their periodic images within a frame of width `pad * L` around the box (by default four mean particle spacings), and of every periodic tetrahedron the copy whose centroid lies in [c₀, c₀ + L)³ is kept. The reference cell is shifted by an irrational fraction c₀ of the mean particle spacing, so that no centroid of a (degenerate) lattice tessellation lies on its boundary. For degenerate configurations such as an exact lattice, where the Delaunay tessellation is not unique, the connectivity is computed from positions displaced by a fixed pseudo-random `jitter` (in mean particle spacings), the same for every periodic image of a particle, so that all images are split consistently; the volumes are computed from the exact positions.

The reference cell is split along z into `slabs` slabs (by default one per thread, at least two frame widths thick), which are tetrahedralised concurrently: each slab with the particles and periodic images within the frame width of it, keeping the tetrahedra whose centroid lies in the slab. The Delaunay tessellation is local, so the slabs together give the tetrahedra of a single tessellation of the box; `slabs = 1` uses one call of TetGen.

Returns the `(T, 4)` matrix of particle rows and the Lagrangian volumes of the tetrahedra, in the canonical order of `canonical_tetrahedra`. The kept tetrahedra must tile the box exactly; if their volumes do not add up to L³, the frame was too narrow and an error is thrown.
"""
function periodic_delaunay(coords_q::AbstractMatrix, L; pad = nothing, jitter = 1e-5, slabs = nothing)
    N = size(coords_q, 1)
    size(coords_q, 2) == 3 || throw(ArgumentError("coords_q must be an (N, 3) matrix"))
    q = mod.(Matrix{Float64}(coords_q), L)
    h = L / cbrt(N)
    w = pad === nothing ? 4h : pad * L
    0 < w < L / 2 || throw(ArgumentError("the frame width must lie between 0 and L/2, got $w"))
    ns = slabs === nothing ? clamp(min(Threads.nthreads(), floor(Int, L / (2w))), 1, N) : Int(slabs)
    ns >= 1 || throw(ArgumentError("slabs must be at least 1, got $slabs"))
    # the same small displacement for every image of a particle breaks the degeneracy of lattice-like configurations
    δ = jitter * h .* (rand(Random.Xoshiro(2024), 3, N) .- 0.5)
    c0 = (√2 - 1) / 10 * L / cbrt(N)
    zcut = [c0 + L * s / ns for s in 0:ns]                  # slab boundaries of the reference cell along z

    parts = Vector{Tuple{Matrix{Int32}, Vector{Float64}}}(undef, ns)
    @sync for s in 1:ns
        Threads.@spawn parts[s] = slab_delaunay(q, δ, L, w, c0, zcut[s], zcut[s + 1])
    end
    tets, VL = reduce(vcat, first.(parts)), reduce(vcat, last.(parts))
    perm = sortperm([(tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4]) for t in axes(tets, 1)])
    tets, VL = tets[perm, :], VL[perm]
    isapprox(sum(VL), L^3; rtol = 1e-8) ||
        error("the periodic Delaunay tetrahedra cover $(sum(VL) / L^3) of the box instead of 1: increase `pad`")
    return tets, VL
end

"""
    slab_delaunay(q, δ, L, w, c0, z0, z1)

The tetrahedra of the periodic Delaunay tessellation (see `periodic_delaunay`) whose centroid lies in the slab
[c0, c0 + L)² × [z0, z1): the particles `q` (in [0, L)³) and their periodic images within the frame width `w` of the
slab are tetrahedralised, displaced by `δ` of their particle. In canonical form, see `canonical_tetrahedra`.
"""
function slab_delaunay(q, δ, L, w, c0, z0, z1)
    N = size(q, 1)
    owner, pts = Int[], NTuple{3, Float64}[]
    shifts(y, lo, hi) = (s for s in (0, -1, 1) if lo <= y + s * L <= hi)
    @inbounds for p in 1:N
        for sz in shifts(q[p, 3], z0 - w, z1 + w), sy in shifts(q[p, 2], -w, L + w), sx in shifts(q[p, 1], -w, L + w)
            push!(owner, p)
            push!(pts, (q[p, 1] + sx * L, q[p, 2] + sy * L, q[p, 3] + sz * L))
        end
    end
    io = TetGen.RawTetGenIO{Float64}()
    io.pointlist = [pts[r][d] + δ[d, owner[r]] for d in 1:3, r in eachindex(pts)]
    all_tets = tetrahedralize(io, "Q").tetrahedronlist                # (4, T), 1-based indices into pts
    centre(t, d) = sum(pts[all_tets[k, t]][d] for k in 1:4) / 4
    keep = [t for t in axes(all_tets, 2) if c0 <= centre(t, 1) < c0 + L && c0 <= centre(t, 2) < c0 + L && z0 <= centre(t, 3) < z1]
    ids = [owner[all_tets[k, t]] for t in keep, k in 1:4]
    return canonical_tetrahedra(ids, (r, k) -> SVector(pts[all_tets[k, keep[r]]]))
end

"Pseudo-random number in [0, 1), fixed by the particle id and the axis d."
jitter_unit(id, d) = (hash((Int(id), d), 0x5053445446452024) >> 11) * 2.0^-53

"""
    delaunay(coords_q; jitter = 1e-5, ids = axes(coords_q, 1), spacing = nothing)

Delaunay tessellation of the points `coords_q` (no periodic images), in the canonical order of `canonical_tetrahedra`. For degenerate configurations such as a lattice, the connectivity is computed from positions displaced by `jitter` mean particle spacings in a pseudo-random direction fixed by the particle id `ids[r]` of every row, so that every copy of a particle (in another sub-box, or a periodic image in a frame) is displaced identically and degenerate cells are split consistently. `spacing` is the mean particle spacing (default: from the bounding box).
"""
function delaunay(coords_q::AbstractMatrix; jitter = 1e-5, ids = axes(coords_q, 1), spacing = nothing)
    size(coords_q, 2) == 3 || throw(ArgumentError("coords_q must be an (N, 3) matrix"))
    N = size(coords_q, 1)
    N >= 4 || throw(ArgumentError("at least 4 points are needed, got $N"))
    length(ids) == N || throw(ArgumentError("ids must have one entry per point"))
    q = Matrix{Float64}(coords_q)

    # mean particle spacing, by default from the bounding box
    extent = [maximum(@view q[:, d]) - minimum(@view q[:, d]) for d in 1:3]
    h = spacing === nothing ? cbrt(prod(extent) / N) : Float64(spacing)
    h > 0 || throw(ArgumentError("the points are degenerate (zero bounding-box volume)"))

    pts = Matrix{Float64}(q')                                       # (3, N)
    # small fixed displacement breaks the degeneracy of lattice-like configurations
    δ = [jitter * h * (jitter_unit(ids[r], d) - 0.5) for d in 1:3, r in 1:N]

    io = TetGen.RawTetGenIO{Float64}()
    io.pointlist = pts .+ δ
    all_tets = tetrahedralize(io, "Q").tetrahedronlist             # (4, T), 1-based indices into pts

    tets, VL = canonical_tetrahedra(permutedims(all_tets), (t, k) -> SVector(pts[1, all_tets[k, t]], pts[2, all_tets[k, t]], pts[3, all_tets[k, t]]))

    # discard zero-volume slivers (exact positions can be coplanar even if the jittered ones are not)
    good = VL .> 1e-12 * h^3
    return tets[good, :], VL[good]
end