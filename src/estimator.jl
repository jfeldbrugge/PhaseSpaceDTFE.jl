#using TimerOutputs
#const to = TimerOutput()

"""
    PS_DTFE

Phase-space DTFE estimator on a tessellation of the Lagrangian particle distribution: the periodic Delaunay
tessellation of the initial positions (default) or the Kuhn tetrahedralisation of the lattice. Every tetrahedron is
stored once, in the Morton order of its Eulerian centroid. Its properties are

    simplices    # (T, 4) particle rows (as in the input) of the corners of every tetrahedron
    positions    # (N, 3) final positions, rows as in the input
    geometry     # (12, T) first corner (at its image in [0, L)³) and inverse edge matrix (column-major, NaN for a
                 # degenerate tetrahedron); empty for `compact = true`, where they are recomputed from the positions

Fields are interpolated linearly with the barycentric coordinates λ = B d of a point in each tetrahedron that
contains it, where d is the separation of the point from the first corner, and summed over the streams. How a point
and a tetrahedron are related depends on `images`:

  * `:all` (periodic estimator, `PS_DTFE_periodic`): the point is taken modulo the box and its separation from the
    first corner at the nearest periodic image, so that every tetrahedron is found at all its periodic images.
  * `:modulo` (wrapped sub-box estimator): the point is taken modulo the box and the tetrahedra are taken where they
    are stored; the periodic images are stored as separate tetrahedra (a frame around the box).
  * `:none` (non-periodic estimator, `PS_DTFE_nonperiodic`): the point and the tetrahedra are taken as they are.

`estimator.periodicity` is `true` unless `images == :none`.
"""
struct PS_DTFE
    L::Float64
    Ni::Int
    tessellation::Symbol            # :kuhn or :delaunay
    rho::Vector{Float64}            # density at the particles
    xt::Matrix{Float64}             # final positions (3, N)
    velocities::Matrix{Float64}     # velocities (N, 3), rows as in the input; (0, 3) if not given
    corners::Matrix{Int32}          # (4, T) particle rows of the corners of every tetrahedron
    geometry::Matrix{Float64}       # (12, T), or (12, 0) for compact = true
    index::CellIndex
    images::Symbol                  # :all, :modulo or :none: how points relate to the tetrahedra
end

function Base.getproperty(est::PS_DTFE, s::Symbol)
    s === :simplices && return transpose(getfield(est, :corners))
    s === :positions && return transpose(getfield(est, :xt))
    s === :periodicity && return getfield(est, :images) != :none
    return getfield(est, s)
end
Base.propertynames(::PS_DTFE, private::Bool = false) = (fieldnames(PS_DTFE)..., :simplices, :positions, :periodicity)

"""
    assemble(L, Ni, tessellation, rho, x, velocities, tets, G, images; compact = false)

The estimator of the tetrahedra `tets` (a `(T, 4)` matrix of particle rows) on the final positions `x` with vertex
densities `rho`: the tetrahedra in Morton order, their geometry (unless `compact`) and the cell index, all with the
relation of points and tetrahedra given by `images` (see `PS_DTFE`).
"""
function assemble(L, Ni, tessellation, rho, x, velocities, tets, G, images; compact::Bool = false)
    images in (:all, :modulo, :none) || throw(ArgumentError("images must be :all, :modulo or :none, got :$images"))
    xt = permutedims(Matrix{Float64}(x))
    corners = morton_order(xt, tets, L)
    T = size(corners, 2)
    geometry = Matrix{Float64}(undef, 12, compact ? 0 : T)
    if !compact
        Threads.@threads for t in 1:T
            o, E = tet_frame(xt, corners, t, L, images)
            D, B = edge_inverse(E)
            D == 0 && (B = SMatrix{3, 3, Float64}(ntuple(_ -> NaN, 9)))
            @inbounds geometry[1:3, t] = o
            @inbounds geometry[4:12, t] = vec(B)
        end
    end
    index = build_index(T, L, G; periodic = images != :none) do t
        o, E = tet_frame(xt, corners, t, L, images)
        edge_inverse(E)[1] == 0 ? nothing : (o, o + E[:, 1], o + E[:, 2], o + E[:, 3])
    end
    return PS_DTFE(L, Ni, tessellation, rho, xt, Matrix{Float64}(velocities), corners, geometry, index, images)
end

"""
    tessellate(coords_q, sim_box, tessellation; pad = nothing, jitter = 1e-5)

The Lagrangian tetrahedra, as a `(T, 4)` matrix of particle rows, and their Lagrangian volumes: the Kuhn split of the lattice (`tessellation = :kuhn`), or the periodic Delaunay tessellation of the initial positions (`tessellation = :delaunay`, see `periodic_delaunay` for `pad` and `jitter`).
"""
function tessellate(coords_q, sim_box::SimBox, tessellation::Symbol; pad = nothing, jitter = 1e-5, slabs = nothing)
    if tessellation == :kuhn
        tets = lattice_tetrahedra(lattice_sites(coords_q, sim_box))
        return tets, fill((sim_box.L / sim_box.Ni)^3 / 6, size(tets, 1))
    elseif tessellation == :delaunay
        return periodic_delaunay(coords_q, sim_box.L; pad, jitter, slabs)
    elseif tessellation == :delaunay_nonperiodic
        return periodic_delaunay(coords_q, sim_box.L)
    end
    throw(ArgumentError("tessellation must be :kuhn or :delaunay, got :$tessellation"))
end

"""
    vertex_densities(x, tets, VL, m, mode)

Density at the particles of a non-periodic tessellation, as `vertex_densities_periodic` without periodic boundaries. For `mode = :lagrangian` the mean density is ρ̄ = M / Σₜ V_L(t), the total mass over the Lagrangian volume covered by the tetrahedra.
"""
function vertex_densities(x, tets, VL, m, mode::Symbol)
    mode in (:mass, :lagrangian) || throw(ArgumentError("vertex_density must be :mass or :lagrangian, got :$mode"))
    N = size(x, 1)
    check_masses(m, N)
    V, W = zeros(N), zeros(N)
    @inbounds for t in axes(tets, 1)
        v = (tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4])
        vol = abs(det(tet_edges(x, v...))) / 6
        for c in v
            V[c] += vol
            W[c] += VL[t]
        end
    end
    mode == :mass && return 4 .* m ./ V
    ρ̄ = (m isa Number ? m * N : sum(m)) / sum(VL)
    return ρ̄ .* W ./ V
end

"""
    vertex_densities_periodic(x, tets, VL, m, L, mode)

Density at the particles from the Eulerian volumes V of the tetrahedra around each particle (its star W):

  * `mode = :mass`: ρᵢ = 4 mᵢ / V(Wᵢ), as in `PhaseSpaceDTFE`. Valid for any particle masses: every particle's mass is shared by the tetrahedra of its star.
  * `mode = :lagrangian`: ρᵢ = ρ̄ V_L(Wᵢ) / V(Wᵢ), with the Lagrangian star volume V_L and the mean density ρ̄ = M/L³ (Feldbrugge 2024, eq. 26). It assumes a uniform initial density and is then exact, also when the stars differ in size (glass initial conditions, a Delaunay tessellation) or the masses differ with the local resolution of the initial conditions (zoom simulations).

Both conserve the total mass M, and they coincide for the Kuhn split of a lattice with equal masses. `m` is a single mass or a vector with the mass of every particle.
"""
function vertex_densities_periodic(x, tets, VL, m, L, mode::Symbol)
    mode in (:mass, :lagrangian) || throw(ArgumentError("vertex_density must be :mass or :lagrangian, got :$mode"))
    N = size(x, 1)
    check_masses(m, N)
    vol = Vector{Float64}(undef, size(tets, 1))
    Threads.@threads for t in axes(tets, 1)
        @inbounds vol[t] = abs(det(tet_edges(x, tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4], L))) / 6
    end
    V, W = zeros(N), zeros(N)
    @inbounds for t in axes(tets, 1), k in 1:4
        c = tets[t, k]
        V[c] += vol[t]
        W[c] += VL[t]
    end
    mode == :mass && return 4 .* m ./ V
    ρ̄ = (m isa Number ? m * N : sum(m)) / L^3
    return ρ̄ .* W ./ V
end

"Check that `m` is a single mass or one positive mass per particle."
function check_masses(m, N)
    m isa Number && return m > 0 || throw(ArgumentError("the particle mass must be positive"))
    length(m) == N || throw(ArgumentError("m must be a number or a vector of $N particle masses, got $(length(m)) masses"))
    all(>(0), m) || throw(ArgumentError("the particle masses must be positive"))
    return true
end


"""
    core_filter(q, tets, VL, (lo_c, hi_c))

The tetrahedra `tets`, and their Lagrangian volumes `VL`, whose centroid in the initial positions `q` lies in the core
[lo_c, hi_c) of a sub-box. Threaded, and without allocating per tetrahedron.
"""
function core_filter(q, tets, VL, (lo_c, hi_c))
    keep = Vector{Bool}(undef, size(tets, 1))
    Threads.@threads for t in axes(tets, 1)
        ok = true
        @inbounds for d in 1:3
            c = (q[tets[t, 1], d] + q[tets[t, 2], d] + q[tets[t, 3], d] + q[tets[t, 4], d]) / 4
            ok &= lo_c[d] <= c < hi_c[d]
        end
        @inbounds keep[t] = ok
    end
    return tets[keep, :], VL[keep]
end

"""
    drop_long_tets(q, tets, VL, hmax)

Remove the tetrahedra with two corners further apart than `hmax` along an axis in the initial positions `q`: the thin slivers on the convex hull of a non-periodic Delaunay tessellation. Returns the kept tetrahedra, their Lagrangian volumes and the number removed.
"""
function drop_long_tets(q, tets, VL, hmax)
    keep = trues(size(tets, 1))
    @inbounds for t in axes(tets, 1)
        for a in 1:3, b in a+1:4, d in 1:3
            abs(q[tets[t, a], d] - q[tets[t, b], d]) > hmax && (keep[t] = false)
        end
    end
    return tets[keep, :], VL[keep], count(!, keep)
end


"""
    PS_DTFE(coords_q, coords_x, velocities, m, sim_box; periodicity = :periodic, kwargs...)
    PS_DTFE(coords_q, coords_x, m, sim_box; ...)

Construct the phase-space DTFE estimator with (`periodicity = :periodic`, default) or without (`periodicity = :nonperiodic`) periodic boundary conditions; calls `PS_DTFE_periodic` or `PS_DTFE_nonperiodic` with the remaining keywords.
"""
function PS_DTFE(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                          periodicity::Symbol = :periodic, tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass, G::Integer = sim_box.Ni, pad = nothing,
                          jitter = 1e-5, drop_simplices=true, compact::Bool = false, slabs = nothing)
    if periodicity == :nonperiodic
        PS_DTFE_nonperiodic(coords_q, coords_x, velocities, m, sim_box;
                            tessellation, vertex_density, G, jitter, drop_simplices, compact)
    elseif periodicity == :periodic
        PS_DTFE_periodic(coords_q, coords_x, velocities, m, sim_box;
                        tessellation, vertex_density, G, pad, jitter, compact, slabs)
    else
        throw(ArgumentError("Invalid periodicity entered, must be 'periodic' or 'nonperiodic'.'"))
    end
end

PS_DTFE(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)


"""
    PS_DTFE_nonperiodic(coords_q, coords_x, velocities, m, sim_box; vertex_density = :mass, G, jitter, drop_simplices = true, compact = false)
    PS_DTFE_nonperiodic(coords_q, coords_x, m, sim_box; ...)

Construct the phase-space DTFE estimator without periodic boundaries, from the Delaunay tessellation of the initial
positions `coords_q` (only `tessellation = :delaunay`). The arguments are as for `PS_DTFE_periodic`; points are taken
as given (`images = :none`), not modulo the box.

  * The final positions are unwrapped to the periodic image nearest to the initial positions: every displacement
    `coords_x - coords_q` is taken in [-L/2, L/2]. This suits particles of a periodic simulation whose positions are
    wrapped into the box. A particle that truly moves more than L/2 along an axis is placed at the wrong image, L away.
  * With `drop_simplices = true` (default) the tetrahedra with two corners more than three lattice spacings L/Ni apart
    along an axis are removed: the thin slivers on the convex hull of the initial positions.
  * The tetrahedra cover the convex hull of the particles, which lies inside the box: points between the hull and the
    faces of the box lie in no tetrahedron and have density zero.
  * Near the hull the stars of the particles are incomplete. With `vertex_density = :mass` this biases the vertex
    densities high there (for a uniform glass the density within 1.5 spacings of a face exceeds the mean by about a
    third); `vertex_density = :lagrangian` avoids this, with the mean density ρ̄ = M / Σ V_L of the tetrahedra kept.
  * Both vertex densities conserve the total mass, and `:lagrangian` is exact for affine maps in the interior.
"""
function PS_DTFE_nonperiodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                             tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass, G::Integer = sim_box.Ni, jitter = 1e-5, drop_simplices=true,
                             core = nothing, verbose=false, ids = axes(coords_q, 1), compact::Bool = false)
    size(velocities) == size(coords_x) || isempty(velocities) ||
        throw(ArgumentError("velocities must have the same size as coords_x"))
    size(coords_q) == size(coords_x) || throw(ArgumentError("coords_q and coords_x must have the same size"))
    
    #coords_q = unwrap_x(coords_q, sim_box)
    #coords_x = unwrap_x(coords_x, sim_box)

    q = Matrix{Float64}(coords_q)
    x = Matrix{Float64}(coords_x)
    s = x .- q
    s .-= sim_box.L .* round.(s ./ sim_box.L)     # displacement in [-L/2, L/2]
    x = q .+ s 

    coords_q = q
    coords_x = x
    
    L = sim_box.L
    x = Matrix{Float64}(coords_x)

    if tessellation != :delaunay
        throw(ArgumentError("Non-periodic construction supports only Delaunay tessellation (under construction)."))
    end

    verbose ? println("doing tessellation") : nothing
    tets, VL = delaunay(coords_q; jitter, ids, spacing = L / sim_box.Ni)
    verbose ? println("got tessellation") : nothing

    if drop_simplices
        tets, VL, ndropped = drop_long_tets(Matrix{Float64}(coords_q), tets, VL, 3 * L / sim_box.Ni)   # 3 lattice spacings of the whole box, also for a sub-box
        verbose ? println("Removed $ndropped squeezed simplices from Delaunay tessellation.") : nothing
    end

    rho = vertex_densities(x, tets, VL, m, vertex_density)

    if core !== nothing
        tets, VL = core_filter(q, tets, VL, core)
    end

    #=
    h = L / cbrt(size(coords_q, 1))
    q = Matrix{Float64}(coords_q)
    longest(t) = maximum(abs(q[tets[t, a], d] - q[tets[t, b], d]) for a in 1:3 for b in a+1:4 for d in 1:3)
    nlong = count(t -> longest(t) > 3h, 1:size(tets, 1))
    @show nlong
    =#

    return assemble(L, sim_box.Ni, :delaunay, rho, x, velocities, tets, G, :none; compact)
end

PS_DTFE_nonperiodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE_nonperiodic(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)




"""
    PS_DTFE_periodic(coords_q, coords_x, velocities, m, sim_box; tessellation = :delaunay, vertex_density = :mass, G, pad, jitter, compact)
    PS_DTFE_periodic(coords_q, coords_x, m, sim_box; ...)

Construct the phase-space DTFE estimator from the initial positions `coords_q`, the final positions `coords_x` and (optionally) the velocities of an N-body simulation in a periodic box. These are `(N, 3)` matrices with one row per particle, in the same order. `m` is the particle mass: a number, or a vector with the mass of every particle for simulations whose particle masses differ.

  * `tessellation = :kuhn`: the Kuhn split of the Lagrangian lattice, for lattice initial conditions with N = Ni³; the lattice site of every particle is recovered from `coords_q`, so the rows may be in any order.
  * `tessellation = :delaunay` (default): the periodic Delaunay tessellation of `coords_q`, as in `PhaseSpaceDTFE`; for glass or other non-lattice initial conditions. `pad` is the width of the frame of periodic images in units of L (default: four mean particle spacings); `jitter` (default 1e-5 mean spacings) breaks the degeneracy of lattice-like initial positions consistently across periodic images. With `jitter = 0` the connectivity is exactly the Delaunay tessellation of `coords_q`, as in `PhaseSpaceDTFE`, which requires positions in general position (not an exact lattice).
  * `vertex_density = :mass` (default, as in `PhaseSpaceDTFE`) or `:lagrangian`; see `vertex_densities`. `:lagrangian` is recommended with `:delaunay`.
  * `G` sets the number of cells per side of the search index (default `sim_box.Ni`).
  * `slabs` sets the number of slabs of the periodic Delaunay tessellation that are tetrahedralised concurrently
    (default: one per thread); see `periodic_delaunay`.
  * `compact = true` does not store the first corner and inverse edge matrix of every tetrahedron, but recomputes
    them from the positions at every query, with identical results: the estimator needs about a third of the memory,
    and point queries are slower.

The forms with a `depth` argument after `m`, as in `PhaseSpaceDTFE`, are accepted; the depth is not used.
"""
function PS_DTFE_periodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                          tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass, G::Integer = sim_box.Ni, pad = nothing,
                          jitter = 1e-5, compact::Bool = false, slabs = nothing)
    size(velocities) == size(coords_x) || isempty(velocities) ||
        throw(ArgumentError("velocities must have the same size as coords_x"))
    size(coords_q) == size(coords_x) || throw(ArgumentError("coords_q and coords_x must have the same size"))
    L = sim_box.L
    x = Matrix{Float64}(coords_x)
    tets, VL = tessellate(coords_q, sim_box, tessellation; pad, jitter, slabs)
    rho = vertex_densities_periodic(x, tets, VL, m, L, vertex_density)

    #=
    q = Matrix{Float64}(coords_q)
    h = L / cbrt(size(q, 1))
    tets, VL = tessellate(coords_q, sim_box, tessellation; pad, jitter)

    function longest_edge(q, v, L)
        E = tet_edges(q, v[1], v[2], v[3], v[4], L)      # edges 1→2, 1→3, 1→4
        e = (E[:, 1], E[:, 2], E[:, 3], E[:, 2] - E[:, 1], E[:, 3] - E[:, 1], E[:, 3] - E[:, 2])
        return maximum(maximum(abs, ek) for ek in e)
    end

    lmax = [longest_edge(q, (tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4]), L) for t in axes(tets, 1)]
    @show count(>(3h), lmax) length(lmax) maximum(lmax) / h
    =#

    return assemble(L, sim_box.Ni, tessellation, rho, x, velocities, tets, G, :all; compact)
end

PS_DTFE_periodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE_periodic(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)





function PS_DTFE_wrapped(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                             tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass, G::Integer = sim_box.Ni, jitter = 1e-5, drop_simplices=true,
                             core = nothing, verbose=false, ids = axes(coords_q, 1), compact::Bool = false)
    size(velocities) == size(coords_x) || isempty(velocities) ||
        throw(ArgumentError("velocities must have the same size as coords_x"))
    size(coords_q) == size(coords_x) || throw(ArgumentError("coords_q and coords_x must have the same size"))
    
    #coords_q = unwrap_x(coords_q, sim_box)
    #coords_x = unwrap_x(coords_x, sim_box)

    q = Matrix{Float64}(coords_q)
    x = Matrix{Float64}(coords_x)
    s = x .- q
    s .-= sim_box.L .* round.(s ./ sim_box.L)     # displacement in [-L/2, L/2]
    x = q .+ s 

    coords_q = q
    coords_x = x
    
    L = sim_box.L
    x = Matrix{Float64}(coords_x)

    if tessellation != :delaunay
        throw(ArgumentError("Non-periodic construction supports only Delaunay tessellation (under construction)."))
    end

    tets, VL = delaunay(coords_q; jitter, ids, spacing = L / sim_box.Ni)

    if drop_simplices
        verbose ? println("dropping simplices") : nothing
        tets, VL, ndropped = drop_long_tets(Matrix{Float64}(coords_q), tets, VL, 3 * L / sim_box.Ni)   # 3 lattice spacings of the whole box, also for a sub-box
        verbose ? println("Removed $ndropped squeezed simplices from Delaunay tessellation.") : nothing
    end

    rho = vertex_densities(x, tets, VL, m, vertex_density)
    verbose ? println("got vertex densities") : nothing

    if core !== nothing
        tets, VL = core_filter(q, tets, VL, core)
    end

    #=
    h = L / cbrt(size(coords_q, 1))
    q = Matrix{Float64}(coords_q)
    longest(t) = maximum(abs(q[tets[t, a], d] - q[tets[t, b], d]) for a in 1:3 for b in a+1:4 for d in 1:3)
    nlong = count(t -> longest(t) > 3h, 1:size(tets, 1))
    @show nlong
    =#

    return assemble(L, sim_box.Ni, :delaunay, rho, x, velocities, tets, G, :modulo; compact)
end

PS_DTFE_wrapped(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE_wrapped(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)


qRange(n, L) = range(0, L, n + 1)[1:end - 1]

function unwrap_x(x::AbstractMatrix, simbox::SimBox)
    q_z = repeat(qRange(simbox.Ni,  simbox.L), 1, simbox.Ni,  simbox.Ni)
    q_y = permutedims(q_z, (3, 1, 2))
    q_x = permutedims(q_z, (2, 3, 1))
    q = [q_x[:] q_y[:] q_z[:]]
    return unwrap_s.(x - q,  simbox.L) + q
end


"The coordinates in which the estimator tests the point `p`: modulo the box, unless `images == :none`."
@inline query_point(est::PS_DTFE, p) = est.images == :none ? (Float64(p[1]), Float64(p[2]), Float64(p[3])) :
                                                             (mod(p[1], est.L), mod(p[2], est.L), mod(p[3], est.L))

"Separation of the coordinate `a` of a query point from the coordinate `o` of a first corner (see `PS_DTFE`)."
@inline separation(est::PS_DTFE, a, o) = est.images == :all ? unwrap_box(a - o, est.L) : a - o

"First corner and inverse edge matrix of tetrahedron t: stored, or recomputed for a compact estimator."
@inline function tet_inverse(est::PS_DTFE, t)
    g = est.geometry
    if size(g, 2) == 0
        o, E = tet_frame(est.xt, est.corners, t, est.L, est.images)
        return o, edge_inverse(E)[2]
    end
    @inbounds return SVector(g[1, t], g[2, t], g[3, t]), SMatrix{3, 3, Float64}(g[4, t], g[5, t], g[6, t], g[7, t], g[8, t],
                                                                                  g[9, t], g[10, t], g[11, t], g[12, t])
end

"Barycentric coordinates (λ₁, λ₂, λ₃) in tetrahedron t of the point (x, y, z) from `query_point`, and whether the point lies in it."
@inline function barycentric(est::PS_DTFE, t, x, y, z)
    o, B = tet_inverse(est, t)
    λ1, λ2, λ3 = barycentric_coordinates(B, separation(est, x, o[1]), separation(est, y, o[2]), separation(est, z, o[3]))
    return λ1, λ2, λ3, inside(λ1, λ2, λ3)
end

"""
Whether the barycentric coordinates describe a point of the closed tetrahedron (false for NaN). Evaluated without
branches: most candidate tetrahedra fail one of the tests at random, which a branch predictor cannot learn.
"""
@inline inside(λ1, λ2, λ3) = (λ1 >= 0) & (λ2 >= 0) & (λ3 >= 0) & (λ1 + λ2 + λ3 <= 1)

"Linear interpolation of the vertex values f on tetrahedron t."
@inline function interpolate(f, est::PS_DTFE, t, λ1, λ2, λ3)
    c = est.corners
    @inbounds f0 = f[c[1, t]]
    @inbounds return f0 + λ1 * (f[c[2, t]] - f0) + λ2 * (f[c[3, t]] - f0) + λ3 * (f[c[4, t]] - f0)
end

"""
    density(p, estimator::PS_DTFE)

Phase-space DTFE density at the point `p` (taken modulo the box), summed over all streams.
"""
function density(p::AbstractVector{<:Real}, est::PS_DTFE)
    x, y, z = query_point(est, p)
    ρ = 0.
    for t in cell_tetrahedra(est.index, x, y, z)
        λ1, λ2, λ3, in = barycentric(est, t, x, y, z)
        in && (ρ += interpolate(est.rho, est, t, λ1, λ2, λ3))
    end
    return ρ
end

"""
    numberOfStreams(p, estimator::PS_DTFE)

Number of streams at the point `p` (taken modulo the box).
"""
function numberOfStreams(p::AbstractVector{<:Real}, est::PS_DTFE)
    x, y, z = query_point(est, p)
    return count(t -> barycentric(est, t, x, y, z)[4], cell_tetrahedra(est.index, x, y, z))
end

"""
    inSimplices(p, estimator::PS_DTFE)

Indices of the tetrahedra that contain the point `p`; `estimator.simplices[i, :]` are the particle rows of their
corners.
"""
function inSimplices(p::AbstractVector{<:Real}, est::PS_DTFE)
    x, y, z = query_point(est, p)
    return [Int(t) for t in cell_tetrahedra(est.index, x, y, z) if barycentric(est, t, x, y, z)[4]]
end

"""
    velocity(p, estimator::PS_DTFE, single_stream = false)

Velocity of every stream at the point `p`, one row per stream (an `(n, 3)` matrix). With `single_stream = true`, `[NaN NaN NaN]` is returned in multistream regions.
"""
function velocity(p::AbstractVector{<:Real}, est::PS_DTFE, single_stream = false)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = query_point(est, p)
    vs = Vector{NTuple{3, Float64}}()
    V = est.velocities
    for t in cell_tetrahedra(est.index, x, y, z)
        λ1, λ2, λ3, in = barycentric(est, t, x, y, z)
        in || continue
        push!(vs, ntuple(d -> interpolate(@view(V[:, d]), est, t, λ1, λ2, λ3), 3))
    end
    single_stream && length(vs) != 1 && return [NaN NaN NaN]
    return isempty(vs) ? zeros(0, 3) : permutedims(reduce(hcat, [collect(v) for v in vs]))
end

"""
    velocitySum(p, estimator::PS_DTFE)

Stream-density weighted mean velocity at the point `p`, as a `(1, 3)` matrix.
"""
function velocitySum(p::AbstractVector{<:Real}, est::PS_DTFE)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = query_point(est, p)
    V = est.velocities
    ρsum, v = 0., zeros(1, 3)
    for t in cell_tetrahedra(est.index, x, y, z)
        λ1, λ2, λ3, in = barycentric(est, t, x, y, z)
        in || continue
        ρ = interpolate(est.rho, est, t, λ1, λ2, λ3)
        ρsum += ρ
        for d in 1:3
            v[1, d] += ρ * interpolate(@view(V[:, d]), est, t, λ1, λ2, λ3)
        end
    end
    return v ./ ρsum
end

# deformation of velocity ----------------------------------------------------------

"""
    velocity_gradient(V, v, B)

Velocity gradient J[i, j] = ∂vᵢ/∂xⱼ (constant) on the tetrahedron with corners `v`, from the `(N, 3)` velocities `V`
and the inverse edge matrix `B` (B[a, j] = ∂λₐ/∂xⱼ). The product is written out without fused multiply-adds (see
`edge_inverse`), so that point queries and grids obtain the same bits wherever it is inlined.
"""
@inline function velocity_gradient(V, v, B)
    @inbounds begin
        v0 = SVector(V[v[1], 1], V[v[1], 2], V[v[1], 3])
        dV = hcat(SVector(V[v[2], 1], V[v[2], 2], V[v[2], 3]) - v0,
                  SVector(V[v[3], 1], V[v[3], 2], V[v[3], 3]) - v0,
                  SVector(V[v[4], 1], V[v[4], 2], V[v[4], 3]) - v0)   # 3×3, columns = edges
        return SMatrix{3, 3, Float64}(ntuple(9) do k
            i, j = (k - 1) % 3 + 1, (k - 1) ÷ 3 + 1
            dV[i, 1] * B[1, j] + dV[i, 2] * B[2, j] + dV[i, 3] * B[3, j]
        end)
    end
end

"Velocity gradient J[i,j] = ∂vᵢ/∂xⱼ (constant) of stored tetrahedron t."
@inline function velocity_gradient(est::PS_DTFE, t)
    c = est.corners
    @inbounds v = (c[1, t], c[2, t], c[3, t], c[4, t])
    return velocity_gradient(est.velocities, v, tet_inverse(est, t)[2])
end

"""
    velocityDivergence(p, estimator::PS_DTFE)

∇·v of every stream at the point `p`, as a vector with one entry per stream.
"""
function divergence(p::AbstractVector{<:Real}, est::PS_DTFE)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = query_point(est, p)
    out = Float64[]
    for t in cell_tetrahedra(est.index, x, y, z)
        _, _, _, in = barycentric(est, t, x, y, z)
        in && push!(out, tr(velocity_gradient(est, t)))
    end
    return out
end

"Σρ and Σρ·(∇·v) at `p` for one estimator, not normalised, so they can be summed over blocks."
function divergenceSum_parts(p::AbstractVector{<:Real}, est::PS_DTFE)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = query_point(est, p)
    ρsum, num = 0., 0.
    for t in cell_tetrahedra(est.index, x, y, z)
        λ1, λ2, λ3, in = barycentric(est, t, x, y, z)
        in || continue
        ρ = interpolate(est.rho, est, t, λ1, λ2, λ3)
        ρsum += ρ
        num  += ρ * tr(velocity_gradient(est, t))
    end
    return ρsum, num
end

divergenceSum(p, est::PS_DTFE) = ((ρ, n) = divergenceSum_parts(p, est); n / ρ)

"Symmetric part S = (J + Jᵀ)/2 (rate-of-strain tensor)."
strain(J) = (J + J') / 2

"Traceless shear tensor σ = S − (tr S / 3) I."
shear_tensor(J) = (S = strain(J); S - (tr(S) / 3) * I)

"Scalar shear magnitude σ² = σ_ij σ_ij."
shear_scalar(J) = (σ = shear_tensor(J); sum(abs2, σ))

"Vorticity vector ω = ∇×v from J[i,j] = ∂vᵢ/∂xⱼ."
vorticity_vector(J) = SVector(J[3, 2] - J[2, 3],
                              J[1, 3] - J[3, 1],
                              J[2, 1] - J[1, 2])

"Apply f(J) to the velocity gradient of every stream at p; returns a Vector of the results."
function per_stream(f, p::AbstractVector{<:Real}, est::PS_DTFE)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = query_point(est, p)
    out = []
    for t in cell_tetrahedra(est.index, x, y, z)
        _, _, _, in = barycentric(est, t, x, y, z)
        in && push!(out, f(velocity_gradient(est, t)))
    end
    return out
end

"Σρ and Σρ·f(J) at p, not normalised (sum these over blocks, divide at the end)."
function weighted_parts(f, p::AbstractVector{<:Real}, est::PS_DTFE, zero_val)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = query_point(est, p)
    ρsum, num = 0.0, zero_val
    for t in cell_tetrahedra(est.index, x, y, z)
        λ1, λ2, λ3, in = barycentric(est, t, x, y, z)
        in || continue
        ρ = interpolate(est.rho, est, t, λ1, λ2, λ3)
        ρsum += ρ
        num   = num + ρ * f(velocity_gradient(est, t))
    end
    return ρsum, num
end


shear(p, est)     = per_stream(shear_tensor, p, est)
#velocityShearSq(p, est)   = per_stream(shear_scalar, p, est)
vorticity(p, est) = per_stream(vorticity_vector, p, est) 

shear_parts(p, est) = weighted_parts(shear_tensor, p, est, zero(SMatrix{3, 3, Float64}))
#shearSq_parts(p, est)   = weighted_parts(shear_scalar, p, est, 0.0)
vorticity_parts(p, est) = weighted_parts(vorticity_vector, p, est, zero(SVector{3, Float64}))

shearSum(p, est)    = ((ρ, n) = shear_parts(p, est); n / ρ)  
#shearSqSum(p, est)   = ((ρ, n) = shearSq_parts(p, est); n / ρ)
vorticitySum(p, est) = ((ρ, n) = vorticity_parts(p, est); n / ρ)