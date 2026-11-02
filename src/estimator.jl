"""
    PS_DTFE

Phase-space DTFE estimator on a tessellation of the Lagrangian particle distribution: the periodic Delaunay tessellation of the initial positions (default) or the Kuhn tetrahedralisation of the lattice. Every tetrahedron is stored at the periodic image of its first corner in [0, L)³, together with the periodic copies that overlap the box, as

    simplices::Matrix{Int32}     # (T, 4) particle rows (as in the input) of the corners of every stored tetrahedron
    origin::Matrix{Float64}      # (3, T) first corner at the stored periodic image
    Binv::Matrix{Float64}        # (9, T) inverse of the edge matrix, column-major; NaN for degenerate tetrahedra
    original::BitVector          # true for the image of the first corner in [0, L)³, false for a periodic copy

An estimator without periodic boundaries (`periodic == false`) stores every tetrahedron once, at the final positions of its corners, without periodic copies; points are still taken modulo the box, so a tetrahedron is only found at points inside [0, L)³.

Fields are interpolated linearly with the barycentric coordinates λ = B (x - origin) of the point in each tetrahedron that contains it, and summed over the streams. Construct it with `PS_DTFE_periodic`, `PS_DTFE_nonperiodic` or `PS_DTFE(...; periodicity)`.
"""
struct PS_DTFE
    L::Float64
    Ni::Int
    tessellation::Symbol            # :kuhn or :delaunay
    periodic::Bool                  # false for an estimator without periodic boundaries
    rho::Vector{Float64}            # density at the particles
    positions::Matrix{Float64}      # final positions (N, 3), rows as in the input
    velocities::Matrix{Float64}     # velocities (N, 3), rows as in the input; (0, 3) if not given
    simplices::Matrix{Int32}
    origin::Matrix{Float64}
    Binv::Matrix{Float64}
    original::BitVector
    index::CellIndex
end

"""
    tessellate(coords_q, sim_box, tessellation; pad = nothing, jitter = 1e-5)

The Lagrangian tetrahedra, as a `(T, 4)` matrix of particle rows, and their Lagrangian volumes: the Kuhn split of the lattice (`tessellation = :kuhn`), or the periodic Delaunay tessellation of the initial positions (`tessellation = :delaunay`, see `periodic_delaunay` for `pad` and `jitter`).
"""
function tessellate(coords_q, sim_box::SimBox, tessellation::Symbol; pad = nothing, jitter = 1e-5)
    if tessellation == :kuhn
        tets = lattice_tetrahedra(lattice_sites(coords_q, sim_box))
        return tets, fill((sim_box.L / sim_box.Ni)^3 / 6, size(tets, 1))
    elseif tessellation == :delaunay
        return periodic_delaunay(coords_q, sim_box.L; pad, jitter)
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

  * `mode = :mass`: ρᵢ = 4 mᵢ / V(Wᵢ), the classical DTFE estimate. Valid for any particle masses: every particle's mass is shared by the tetrahedra of its star.
  * `mode = :lagrangian`: ρᵢ = ρ̄ V_L(Wᵢ) / V(Wᵢ), with the Lagrangian star volume V_L and the mean density ρ̄ = M/L³ (Feldbrugge 2024, eq. 26). It assumes a uniform initial density and is then exact, also when the stars differ in size (glass initial conditions, a Delaunay tessellation) or the masses differ with the local resolution of the initial conditions (zoom simulations).

Both conserve the total mass M, and they coincide for the Kuhn split of a lattice with equal masses. `m` is a single mass or a vector with the mass of every particle.
"""
function vertex_densities_periodic(x, tets, VL, m, L, mode::Symbol)
    mode in (:mass, :lagrangian) || throw(ArgumentError("vertex_density must be :mass or :lagrangian, got :$mode"))
    N = size(x, 1)
    check_masses(m, N)
    V, W = zeros(N), zeros(N)
    @inbounds for t in axes(tets, 1)
        v = (tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4])
        vol = abs(det(tet_edges(x, v..., L))) / 6
        for c in v
            V[c] += vol
            W[c] += VL[t]
        end
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
                 periodicity::Symbol = :periodic, tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass,
                 G::Integer = sim_box.Ni, pad = nothing, jitter = 1e-5, drop_simplices = true)
    if periodicity == :nonperiodic
        PS_DTFE_nonperiodic(coords_q, coords_x, velocities, m, sim_box;
                            tessellation, vertex_density, G, jitter, drop_simplices)
    elseif periodicity == :periodic
        PS_DTFE_periodic(coords_q, coords_x, velocities, m, sim_box;
                         tessellation, vertex_density, G, pad, jitter)
    else
        throw(ArgumentError("periodicity must be :periodic or :nonperiodic, got :$periodicity"))
    end
end

PS_DTFE(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)

"""
    PS_DTFE_nonperiodic(coords_q, coords_x, velocities, m, sim_box; vertex_density = :mass, G, jitter, drop_simplices = true)
    PS_DTFE_nonperiodic(coords_q, coords_x, m, sim_box; ...)

Construct the phase-space DTFE estimator without periodic images of the tetrahedra, from the (non-periodic) Delaunay tessellation of the initial positions `coords_q`. The arguments are as for `PS_DTFE_periodic`; only `tessellation = :delaunay` is supported.

  * The final positions are unwrapped to the periodic image nearest to the initial position, i.e. the displacements `coords_x - coords_q` are taken in [-L/2, L/2].
  * With `drop_simplices = true` (default), the tetrahedra with two corners further apart than three mean particle spacings L/∛N along an axis are removed. These are the thin slivers that the Delaunay tessellation forms on the convex hull of the initial positions.
  * With `vertex_density = :lagrangian` the mean density is the total mass over the Lagrangian volume of the kept tetrahedra.
"""
function PS_DTFE_nonperiodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                             tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass, G::Integer = sim_box.Ni,
                             jitter = 1e-5, drop_simplices = true)
    size(velocities) == size(coords_x) || isempty(velocities) ||
        throw(ArgumentError("velocities must have the same size as coords_x"))
    size(coords_q) == size(coords_x) || throw(ArgumentError("coords_q and coords_x must have the same size"))
    tessellation == :delaunay ||
        throw(ArgumentError("the non-periodic estimator supports only tessellation = :delaunay, got :$tessellation"))
    L = sim_box.L
    q = Matrix{Float64}(coords_q)
    s = Matrix{Float64}(coords_x) .- q
    s .-= L .* round.(s ./ L)                       # displacement in [-L/2, L/2]
    x = q .+ s

    tets, VL = delaunay(q; jitter)
    if drop_simplices
        tets, VL, _ = drop_long_tets(q, tets, VL, 3 * L / cbrt(size(q, 1)))
    end
    rho = vertex_densities(x, tets, VL, m, vertex_density)

    # non-periodic: one stored image per tetrahedron, no shifts
    T = size(tets, 1)
    simplices, origin, Binv = Matrix{Int32}(undef, T, 4), Matrix{Float64}(undef, 3, T), Matrix{Float64}(undef, 9, T)
    lo, hi, original = Matrix{Float64}(undef, 3, T), Matrix{Float64}(undef, 3, T), trues(T)
    @inbounds for t in 1:T
        v = (tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4])
        simplices[t, 1], simplices[t, 2], simplices[t, 3], simplices[t, 4] = v
        o = SVector{3, Float64}(x[v[1], 1], x[v[1], 2], x[v[1], 3])
        E = tet_edges(x, v...)
        B = det(E) == 0 ? SMatrix{3, 3, Float64}(ntuple(_ -> NaN, 9)) : inv(E)
        origin[:, t] = o
        Binv[:, t] = vec(B)
        lo[:, t] = o .+ min.(0, E[:, 1], E[:, 2], E[:, 3])
        hi[:, t] = o .+ max.(0, E[:, 1], E[:, 2], E[:, 3])
    end

    return PS_DTFE(L, sim_box.Ni, :delaunay, false, rho, x, Matrix{Float64}(velocities), simplices, origin, Binv, original,
                   CellIndex(lo, hi, L, G))
end

PS_DTFE_nonperiodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE_nonperiodic(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)

"""
    PS_DTFE_periodic(coords_q, coords_x, velocities, m, sim_box; tessellation = :delaunay, vertex_density = :mass, G, pad, jitter)
    PS_DTFE_periodic(coords_q, coords_x, m, sim_box; ...)

Construct the phase-space DTFE estimator from the initial positions `coords_q`, the final positions `coords_x` and (optionally) the velocities of an N-body simulation in a periodic box. These are `(N, 3)` matrices with one row per particle, in the same order. `m` is the particle mass: a number, or a vector with the mass of every particle for simulations whose particle masses differ.

  * `tessellation = :delaunay` (default): the periodic Delaunay tessellation of `coords_q`, for lattice, glass or other initial conditions. `pad` is the width of the frame of periodic images in units of L (default: four mean particle spacings); `jitter` (default 1e-5 mean spacings) breaks the degeneracy of lattice-like initial positions consistently across periodic images. With `jitter = 0` the connectivity is exactly the Delaunay tessellation of `coords_q`, which requires positions in general position (not an exact lattice).
  * `tessellation = :kuhn`: the Kuhn split of the Lagrangian lattice, for lattice initial conditions with N = Ni³. Every particle has the same star of 24 tetrahedra, so the estimator is exact for affine maps. The lattice site of every particle is recovered from `coords_q`, so the rows may be in any order.
  * `vertex_density = :mass` (default) or `:lagrangian`; see `vertex_densities_periodic`. `:lagrangian` is recommended with `:delaunay`.
  * `G` sets the number of cells per side of the search index (default `sim_box.Ni`). A larger `G` gives faster queries at the cost of memory.

The forms with a `depth` argument after `m`, as in PhaseSpaceDTFE v1, are accepted; the depth is not used.
"""
function PS_DTFE_periodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                          tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass, G::Integer = sim_box.Ni, pad = nothing,
                          jitter = 1e-5)
    size(velocities) == size(coords_x) || isempty(velocities) ||
        throw(ArgumentError("velocities must have the same size as coords_x"))
    size(coords_q) == size(coords_x) || throw(ArgumentError("coords_q and coords_x must have the same size"))
    L = sim_box.L
    x = Matrix{Float64}(coords_x)
    tets, VL = tessellate(coords_q, sim_box, tessellation; pad, jitter)
    rho = vertex_densities_periodic(x, tets, VL, m, L, vertex_density)

    # stored images: the first-corner image of every tetrahedron and the periodic copies overlapping the box
    rows, shifts = Int32[], SVector{3, Float64}[]
    @inbounds for t in axes(tets, 1)
        o = first_corner(x, tets[t, 1], L)
        E = tet_edges(x, tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4], L)
        lo = o .+ min.(0, E[:, 1], E[:, 2], E[:, 3])
        hi = o .+ max.(0, E[:, 1], E[:, 2], E[:, 3])
        for sz in (-L, 0., L), sy in (-L, 0., L), sx in (-L, 0., L)
            s = SVector(sx, sy, sz)
            if all(hi .+ s .>= 0) && all(lo .+ s .<= L)
                push!(rows, t); push!(shifts, s)
            end
        end
    end

    T = length(rows)
    simplices, origin, Binv = Matrix{Int32}(undef, T, 4), Matrix{Float64}(undef, 3, T), Matrix{Float64}(undef, 9, T)
    lo, hi, original = Matrix{Float64}(undef, 3, T), Matrix{Float64}(undef, 3, T), falses(T)
    @inbounds for (r, (t, s)) in enumerate(zip(rows, shifts))
        v = (tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4])
        simplices[r, 1], simplices[r, 2], simplices[r, 3], simplices[r, 4] = v
        o = first_corner(x, v[1], L) + s
        E = tet_edges(x, v..., L)
        B = det(E) == 0 ? SMatrix{3, 3, Float64}(ntuple(_ -> NaN, 9)) : inv(E)
        origin[:, r] = o
        Binv[:, r] = vec(B)
        lo[:, r] = o .+ min.(0, E[:, 1], E[:, 2], E[:, 3])
        hi[:, r] = o .+ max.(0, E[:, 1], E[:, 2], E[:, 3])
        original[r] = iszero(s)
    end

    return PS_DTFE(L, sim_box.Ni, tessellation, true, rho, x, Matrix{Float64}(velocities), simplices, origin, Binv, original,
                   CellIndex(lo, hi, L, G))
end

PS_DTFE_periodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE_periodic(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)

# PhaseSpaceDTFE v1 call signatures with the depth of the search tree, which is not used
PS_DTFE_periodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, depth::Integer,
                 sim_box::SimBox; kwargs...) = PS_DTFE_periodic(coords_q, coords_x, velocities, m, sim_box; kwargs...)
PS_DTFE_periodic(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m::Union{Real, AbstractVector}, depth::Integer,
                 sim_box::SimBox; kwargs...) = PS_DTFE_periodic(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)


qRange(n, L) = range(0, L, n + 1)[1:end - 1]

function unwrap_x(x::AbstractMatrix, simbox::SimBox)
    q_z = repeat(qRange(simbox.Ni,  simbox.L), 1, simbox.Ni,  simbox.Ni)
    q_y = permutedims(q_z, (3, 1, 2))
    q_x = permutedims(q_z, (2, 3, 1))
    q = [q_x[:] q_y[:] q_z[:]]
    return unwrap_s.(x - q,  simbox.L) + q
end


"Barycentric coordinates (λ₁, λ₂, λ₃) of (x, y, z) in stored tetrahedron t, and whether the point lies in it."
@inline function barycentric(est::PS_DTFE, t, x, y, z)
    @inbounds begin
        dx, dy, dz = x - est.origin[1, t], y - est.origin[2, t], z - est.origin[3, t]
        B = @view est.Binv[:, t]
        λ1 = B[1] * dx + B[4] * dy + B[7] * dz
        λ2 = B[2] * dx + B[5] * dy + B[8] * dz
        λ3 = B[3] * dx + B[6] * dy + B[9] * dz
    end
    return λ1, λ2, λ3, inside(λ1, λ2, λ3)
end

"Whether the barycentric coordinates describe a point of the closed tetrahedron (false for NaN)."
@inline inside(λ1, λ2, λ3) = λ1 >= 0 && λ2 >= 0 && λ3 >= 0 && λ1 <= 1 && λ2 <= 1 && λ3 <= 1 && λ1 + λ2 + λ3 <= 1

"Linear interpolation of the vertex values f on stored tetrahedron t."
@inline function interpolate(f, est::PS_DTFE, t, λ1, λ2, λ3)
    @inbounds f0 = f[est.simplices[t, 1]]
    @inbounds return f0 + λ1 * (f[est.simplices[t, 2]] - f0) + λ2 * (f[est.simplices[t, 3]] - f0) + λ3 * (f[est.simplices[t, 4]] - f0)
end

"""
    density(p, estimator::PS_DTFE)

Phase-space DTFE density at the point `p` (taken modulo the box), summed over all streams.
"""
function density(p::AbstractVector{<:Real}, est::PS_DTFE)
    x, y, z = mod(p[1], est.L), mod(p[2], est.L), mod(p[3], est.L)
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
    x, y, z = mod(p[1], est.L), mod(p[2], est.L), mod(p[3], est.L)
    return count(t -> barycentric(est, t, x, y, z)[4], cell_tetrahedra(est.index, x, y, z))
end

"""
    inSimplices(p, estimator::PS_DTFE)

Indices of the stored tetrahedra that contain the point `p`; `estimator.simplices[i, :]` are the particle rows of their corners.
"""
function inSimplices(p::AbstractVector{<:Real}, est::PS_DTFE)
    x, y, z = mod(p[1], est.L), mod(p[2], est.L), mod(p[3], est.L)
    return [Int(t) for t in cell_tetrahedra(est.index, x, y, z) if barycentric(est, t, x, y, z)[4]]
end

"""
    velocity(p, estimator::PS_DTFE, single_stream = false)

Velocity of every stream at the point `p`, one row per stream (an `(n, 3)` matrix). With `single_stream = true`, `[NaN NaN NaN]` is returned in multistream regions.
"""
function velocity(p::AbstractVector{<:Real}, est::PS_DTFE, single_stream = false)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = mod(p[1], est.L), mod(p[2], est.L), mod(p[3], est.L)
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
    x, y, z = mod(p[1], est.L), mod(p[2], est.L), mod(p[3], est.L)
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
