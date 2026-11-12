const EDGES = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))

"Tolerance, in grid spacings, by which index ranges are widened so that points on a boundary are always tested."
const PAD = 1e-9

"First and last index of the grid range r (spacing Δ) inside [lo, hi], widened by PAD spacings."
@inline function index_range(lo, hi, r, Δ)
    i0 = ceil(Int, (lo - first(r)) / Δ - PAD) + 1
    i1 = floor(Int, (hi - first(r)) / Δ + PAD) + 1
    return max(1, i0), min(length(r), i1)
end

"Spacing of a grid range (1 for a single point, where it does not matter)."
gridstep(r) = length(r) > 1 ? step(r) : one(float(eltype(r)))

"Whether the grid range is usable: a single point or increasing."
valid_range(r) = length(r) <= 1 || step(r) > 0

"Range of integer periodic shifts k for which [lo, hi] + k L overlaps the grid range r."
shift_range(lo, hi, r, L) = ceil(Int, (first(r) - hi) / L):floor(Int, (last(r) - lo) / L)

"""
    section!(px, py, c, z)

Corners of the cross-section of the tetrahedron with corners `c` and the plane at height `z`, written to `px`, `py`; returns their number (0 if the plane misses the tetrahedron). The cross-section is a triangle or a quadrilateral.
"""
@inline function section!(px, py, c, z)
    n = 0
    @inbounds for (a, b) in EDGES
        za, zb = c[a][3] - z, c[b][3] - z
        if za * zb < 0
            s = za / (za - zb)
            n += 1
            px[n], py[n] = c[a][1] + s * (c[b][1] - c[a][1]), c[a][2] + s * (c[b][2] - c[a][2])
        end
    end
    @inbounds for a in 1:4
        if c[a][3] == z
            n += 1
            px[n], py[n] = c[a][1], c[a][2]
        end
    end
    return n
end

"Minimum and maximum of v[1:n]."
@inline function extent(v, n)
    lo, hi = Inf, -Inf
    @inbounds for a in 1:n
        lo, hi = min(lo, v[a]), max(hi, v[a])
    end
    return lo, hi
end

"Extent [xlo, xhi] of the convex polygon with corners (px, py)[1:n] along the line at height y (xlo > xhi if missed)."
@inline function row_extent(px, py, n, y)
    xlo, xhi = Inf, -Inf
    @inbounds for a in 1:n
        if py[a] == y
            xlo, xhi = min(xlo, px[a]), max(xhi, px[a])
        end
        for b in a + 1:n
            ya, yb = py[a] - y, py[b] - y
            if ya * yb < 0
                x = px[a] + ya / (ya - yb) * (px[b] - px[a])
                xlo, xhi = min(xlo, x), max(xhi, x)
            end
        end
    end
    return xlo, xhi
end

"""
    raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, B, f, L, coords, wrapped; images = :all)

Add one tetrahedron — first corner `o`, edge matrix `E` (columns), its inverse `B` and vertex values `f` — to the grid
arrays `ρ[ix, iy, iz]` and `ns[ix, iy, iz]` (either may be `nothing`) for the z-planes with indices in `zr`; see
`raster_visit!`.
"""
function raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, B, f, L, coords, wrapped; images::Symbol = :all)
    raster_visit!(xs, ys, zs, zr, o, E, B, L, coords, wrapped; images) do ix, iy, iz, λ1, λ2, λ3
        @inbounds begin
            ρ === nothing || (ρ[ix, iy, iz] += f[1] + λ1 * (f[2] - f[1]) + λ2 * (f[3] - f[1]) + λ3 * (f[4] - f[1]))
            ns === nothing || (ns[ix, iy, iz] += 1)
        end
    end
end

"""
    raster_visit!(visit, xs, ys, zs, zr, o, E, B, L, coords, wrapped; images = :all)

Call `visit(ix, iy, iz, λ1, λ2, λ3)` at every grid point (ix, iy, iz), with z-index in `zr`, that lies in the
tetrahedron with first corner `o`, edge matrix `E` and its inverse `B`, with the barycentric coordinates λ of the point.
Every plane cuts the tetrahedron, at every periodic image that overlaps the grid, in a convex polygon; for every grid
row the x-interval of the polygon is computed (scanline) and the grid points in it, with a tolerance of `PAD` grid
spacings, are tested with the barycentric coordinates. `coords` are the grid coordinates and `wrapped` the same modulo
L. The barycentric coordinates are computed exactly as in a point query of an estimator with the same `images` (see
`PS_DTFE`), so that the result equals point queries at the grid points:

  * `:all`: from the separation of the wrapped coordinate from the first corner at its nearest periodic image.
  * `:modulo`: from the plain separation of the wrapped coordinate; the image shifted by k L along an axis is only
    tested at the grid points g with `fld(g, L) == k`, the others cannot lie in it.
  * `:none`: from the plain separation of the grid coordinate; only the unshifted tetrahedron.
"""
function raster_visit!(visit, xs, ys, zs, zr, o, E, B, L, (xv, yv, zv), (xm, ym, zm); images::Symbol = :all)
    c = (o, o + E[:, 1], o + E[:, 2], o + E[:, 3])
    lo, hi = min.(c[1], c[2], c[3], c[4]), max.(c[1], c[2], c[3], c[4])
    Δx, Δy, Δz = gridstep(xs), gridstep(ys), gridstep(zs)
    px, py = MVector{10, Float64}(undef), MVector{10, Float64}(undef)
    shifts(a, b, r) = images == :none ? (0:0) : shift_range(a, b, r, L)
    modulo = images == :modulo
    # separation of grid coordinate a (wrapped: am) from the first-corner coordinate b, as in a point query
    sep(a, am, b) = images == :all ? unwrap_box(am - b, L) : images == :modulo ? am - b : a - b
    @inbounds for kz in shifts(lo[3], hi[3], zs)
        sz = kz * L
        iz0, iz1 = index_range(lo[3] + sz, hi[3] + sz, zs, Δz)
        iz0, iz1 = max(iz0, first(zr)), min(iz1, last(zr))
        for iz in iz0:iz1
            modulo && fld(zv[iz], L) != kz && continue
            zp = zv[iz] - sz                                   # plane in the coordinates of the tetrahedron
            dz = sep(zv[iz], zm[iz], o[3])
            n = section!(px, py, c, zp)
            n == 0 && continue
            ylo, yhi = extent(py, n)
            for ky in shifts(ylo, yhi, ys)
                sy = ky * L
                iy0, iy1 = index_range(ylo + sy, yhi + sy, ys, Δy)
                for iy in iy0:iy1
                    modulo && fld(yv[iy], L) != ky && continue
                    yp = yv[iy] - sy
                    dy = sep(yv[iy], ym[iy], o[2])
                    xlo, xhi = row_extent(px, py, n, yp)
                    if xlo > xhi                                    # row within PAD of a corner of the polygon
                        xlo, xhi = extent(px, n)
                    end
                    for kx in shifts(xlo, xhi, xs)
                        sx = kx * L
                        ix0, ix1 = index_range(xlo + sx, xhi + sx, xs, Δx)
                        for ix in ix0:ix1
                            modulo && fld(xv[ix], L) != kx && continue
                            λ1, λ2, λ3 = barycentric_coordinates(B, sep(xv[ix], xm[ix], o[1]), dy, dz)
                            inside(λ1, λ2, λ3) || continue
                            visit(ix, iy, iz, λ1, λ2, λ3)
                        end
                    end
                end
            end
        end
    end
    return nothing
end

"Whether the tetrahedron with z-extent [zlo, zhi] has an image overlapping the z-planes zs[zr]."
@inline function overlaps(zlo, zhi, zs, zr, L)
    z0, z1 = zs[first(zr)], zs[last(zr)]
    return !isempty(ceil(Int, (z0 - zhi) / L):floor(Int, (z1 - zlo) / L))
end

"""
    default_images(estimator)

The periodic images to rasterise so that grids equal the point queries of the estimator: its `images` (see `PS_DTFE`).
"""
default_images(est::PS_DTFE) = est.images

"""
    rasterise(estimator, xs, ys, zs; threads = true)

Density and number of streams on the regular grid xs × ys × zs, as arrays indexed `[ix, iy, iz]`, by rasterising the
tetrahedra of the estimator. The grid may extend beyond the box. With `threads = true` the z-planes are split over the
available threads. The result equals the point queries of the estimator at the grid points, bit for bit.
"""
function rasterise(est::PS_DTFE, xs::AbstractRange, ys::AbstractRange, zs::AbstractRange; kwargs...)
    ρ, ns = zeros(length(xs), length(ys), length(zs)), zeros(Int, length(xs), length(ys), length(zs))
    return rasterise!(ρ, ns, est, xs, ys, zs; kwargs...)
end

"""
    rasterise!(ρ, ns, estimator, xs, ys, zs; threads = true)

Add the density and the number of streams of the estimator on the grid xs × ys × zs to the arrays `ρ[ix, iy, iz]` and
`ns[ix, iy, iz]` (see `rasterise`), without allocating a grid. Either array may be `nothing` when it is not needed, and
the element types are free (e.g. `Float32` densities, `UInt16` stream counts).
"""
function rasterise!(ρ, ns, est::PS_DTFE, xs::AbstractRange, ys::AbstractRange, zs::AbstractRange; threads = true,
                    images::Symbol = default_images(est))
    dims = (length(xs), length(ys), length(zs))
    all(a -> a === nothing || size(a) == dims, (ρ, ns)) || throw(ArgumentError("the grid arrays must have size $dims"))
    foreach_raster_tet(est, xs, ys, zs; threads, images) do t, v, zr, o, E, B, coords, wrapped
        @inbounds f = (est.rho[v[1]], est.rho[v[2]], est.rho[v[3]], est.rho[v[4]])
        raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, B, f, est.L, coords, wrapped; images)
    end
    return ρ, ns
end

"""
    foreach_raster_tet(g, estimator, xs, ys, zs; threads = true, images = default_images(estimator))

Call `g(t, v, zr, o, E, B, coords, wrapped)` for every tetrahedron t of the estimator that reaches the grid
xs × ys × zs, with its corners `v`, the z-planes `zr` of the calling thread, its first corner `o`, edge matrix `E` and
inverse `B` (as in a point query, see `tet_frame`), and the grid coordinates `coords` and `wrapped` (modulo L); pass
these on to `raster_visit!`. With `threads = true` the z-planes are split over the threads, so that every grid point is
visited by one thread only.
"""
function foreach_raster_tet(g, est::PS_DTFE, xs::AbstractRange, ys::AbstractRange, zs::AbstractRange; threads = true,
                            images::Symbol = default_images(est))
    images in (:all, :modulo, :none) || throw(ArgumentError("images must be :all, :modulo or :none, got :$images"))
    any(isempty, (xs, ys, zs)) && return nothing
    all(valid_range, (xs, ys, zs)) || throw(ArgumentError("the grid ranges must be increasing"))
    L, xt, corners = est.L, est.xt, est.corners
    coords = (collect(xs), collect(ys), collect(zs))           # exact range elements, without range arithmetic
    wrapped = map(c -> mod.(c, L), coords)
    T = size(corners, 2)
    parts = threads ? chunks(length(zs)) : [1:length(zs)]
    # z-extent of every tetrahedron, computed once, whether it reaches the grid in x and y at all (so that a slice in
    # x or y skips the tetrahedra that miss it), and the tetrahedra that reach each slab of planes
    zext = Vector{NTuple{2, Float64}}(undef, T)
    reach = Vector{Bool}(undef, T)
    tol = 1e-9 * L                                               # never lose a grid point on a boundary to rounding
    Threads.@threads for t in 1:T
        o, E = tet_frame(xt, corners, t, L, images)
        lo, hi = o .+ min.(0, E[:, 1], E[:, 2], E[:, 3]), o .+ max.(0, E[:, 1], E[:, 2], E[:, 3])
        @inbounds zext[t] = (lo[3], hi[3])
        @inbounds reach[t] = overlaps(lo[1] - tol, hi[1] + tol, xs, 1:length(xs), L) &&
                             overlaps(lo[2] - tol, hi[2] + tol, ys, 1:length(ys), L)
    end
    buckets = [Int32[] for _ in parts]
    for (t, (zlo, zhi)) in enumerate(zext)
        @inbounds reach[t] || continue
        for (b, zr) in enumerate(parts)
            overlaps(zlo, zhi, zs, zr, L) && push!(buckets[b], t)
        end
    end
    Threads.@threads for b in eachindex(parts)
        zr = parts[b]
        @inbounds for t in buckets[b]
            o, E = tet_frame(xt, corners, t, L, images)
            D, B = edge_inverse(E)
            D == 0 && continue
            v = (corners[1, t], corners[2, t], corners[3, t], corners[4, t])
            g(t, v, zr, o, E, B, coords, wrapped)
        end
    end
    return nothing
end

"""
    density_grid(xs, ys, zs, estimator; threads = true)

Density on the regular grid xs × ys × zs, as an array indexed `[ix, iy, iz]`; equal to
`[density([x, y, z], estimator) for x in xs, y in ys, z in zs]`, but computed by rasterising the tetrahedra. """
density_grid(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, est::PS_DTFE; kwargs...) = rasterise(est, xs, ys, zs; kwargs...)[1]

"""
    numberOfStreams_grid(xs, ys, zs, estimator; threads = true)

Number of streams on the regular grid xs × ys × zs, as an array indexed `[ix, iy, iz]`; equal to `[numberOfStreams([x, y, z], estimator) for x in xs, y in ys, z in zs]`.
"""
numberOfStreams_grid(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, est::PS_DTFE; kwargs...) = rasterise(est, xs, ys, zs; kwargs...)[2]

"""
    StreamVelocityGrid

Velocities of all streams on the regular grid xs × ys × zs, stored compactly since the number of streams differs between grid points:

    nstreams::Array{Int32, 3}    # number of streams at every grid point [ix, iy, iz]
    offsets::Vector{Int}         # the streams of the grid point with linear index I are rows offsets[I]:offsets[I + 1] - 1
    velocity::Matrix{T}          # (S, 3): velocity of every stream
    density::Vector{T}           # (S,): density of every stream
    divergence::Vector{T}        # (S,): ∇·v of every stream                       (empty without `deformation`)
    shear::Matrix{T}             # (S, 6): shear tensor σ, components xx yy zz xy xz yz   (empty without `deformation`)
    vorticity::Matrix{T}         # (S, 3): vorticity ω = ∇×v                        (empty without `deformation`)

`G[ix, iy, iz]` is the `(n, 3)` view of the velocities of the n streams at the grid point, as `velocity([x, y, z], estimator)` (the streams in no particular order), and `stream_density(G, ix, iy, iz)` their densities. With `deformation = true`, `stream_divergence`, `stream_shear` and `stream_vorticity` give the velocity deformation of every stream, as `divergence`, `shear` and `vorticity`; `stream_average` gives density-weighted means over the streams. Construct it with `velocity_grid` or `velocity_grid_sb`.
"""
struct StreamVelocityGrid{T <: AbstractFloat, R <: AbstractRange}
    xs::R
    ys::R
    zs::R
    nstreams::Array{Int32, 3}
    offsets::Vector{Int}
    velocity::Matrix{T}
    density::Vector{T}
    divergence::Vector{T}
    shear::Matrix{T}
    vorticity::Matrix{T}
end

"Rows of the streams at grid point (ix, iy, iz)."
@inline stream_rows(G::StreamVelocityGrid, ix, iy, iz) =
    (I = LinearIndices(G.nstreams)[ix, iy, iz]; G.offsets[I]:G.offsets[I + 1] - 1)

Base.size(G::StreamVelocityGrid) = size(G.nstreams)
Base.getindex(G::StreamVelocityGrid, ix::Integer, iy::Integer, iz::Integer) = view(G.velocity, stream_rows(G, ix, iy, iz), :)
Base.getindex(G::StreamVelocityGrid, I::CartesianIndex{3}) = G[Tuple(I)...]
Base.CartesianIndices(G::StreamVelocityGrid) = CartesianIndices(G.nstreams)

"""
    stream_density(G::StreamVelocityGrid, ix, iy, iz)

Density of every stream at grid point (ix, iy, iz), in the order of the rows of `G[ix, iy, iz]`.
"""
stream_density(G::StreamVelocityGrid, ix::Integer, iy::Integer, iz::Integer) = view(G.density, stream_rows(G, ix, iy, iz))

"Whether the grid holds the velocity deformation of the streams."
has_deformation(G::StreamVelocityGrid) = size(G.vorticity, 1) == size(G.velocity, 1)

function check_deformation(G::StreamVelocityGrid)
    has_deformation(G) || throw(ArgumentError("the grid was constructed without `deformation = true`"))
end

"""
    stream_divergence(G::StreamVelocityGrid, ix, iy, iz)

Divergence ∇·v of every stream at grid point (ix, iy, iz), in the order of the rows of `G[ix, iy, iz]`; as
`divergence([x, y, z], estimator)`. Requires `deformation = true`.
"""
stream_divergence(G::StreamVelocityGrid, ix::Integer, iy::Integer, iz::Integer) =
    (check_deformation(G); view(G.divergence, stream_rows(G, ix, iy, iz)))

"""
    stream_shear(G::StreamVelocityGrid, ix, iy, iz)

Shear tensor σ = (J + Jᵀ)/2 - (∇·v / 3) I of every stream at grid point (ix, iy, iz), as a vector of symmetric 3 × 3
matrices in the order of the rows of `G[ix, iy, iz]`; as `shear([x, y, z], estimator)`. Requires `deformation = true`.
"""
function stream_shear(G::StreamVelocityGrid{T}, ix::Integer, iy::Integer, iz::Integer) where {T}
    check_deformation(G)
    σ = G.shear
    return [SMatrix{3, 3, T}(σ[r, 1], σ[r, 4], σ[r, 5], σ[r, 4], σ[r, 2], σ[r, 6], σ[r, 5], σ[r, 6], σ[r, 3])
            for r in stream_rows(G, ix, iy, iz)]
end

"""
    stream_vorticity(G::StreamVelocityGrid, ix, iy, iz)

Vorticity ω = ∇×v of every stream at grid point (ix, iy, iz), one row per stream (an `(n, 3)` view) in the order of the
rows of `G[ix, iy, iz]`; as `vorticity([x, y, z], estimator)`. Requires `deformation = true`.
"""
stream_vorticity(G::StreamVelocityGrid, ix::Integer, iy::Integer, iz::Integer) =
    (check_deformation(G); view(G.vorticity, stream_rows(G, ix, iy, iz), :))

"""
    stream_average(G::StreamVelocityGrid, field)

Density-weighted mean over the streams, Σᵢ ρᵢ fᵢ / Σᵢ ρᵢ, of `field` at every grid point; `NaN` where there is no
stream. `field` is `:velocity` (an `(nx, ny, nz, 3)` array, as `velocitySum`), `:divergence` (an `(nx, ny, nz)` array,
as `divergenceSum`), `:shear` (an `(nx, ny, nz, 6)` array of the components xx yy zz xy xz yz, as `shearSum`) or
`:vorticity` (an `(nx, ny, nz, 3)` array, as `vorticitySum`).
"""
function stream_average(G::StreamVelocityGrid{T}, field::Symbol) where {T}
    field in (:velocity, :divergence, :shear, :vorticity) ||
        throw(ArgumentError("field must be :velocity, :divergence, :shear or :vorticity, got :$field"))
    field == :velocity || check_deformation(G)
    data = field == :velocity ? G.velocity : field == :divergence ? reshape(G.divergence, :, 1) :
           field == :shear ? G.shear : G.vorticity
    k = size(data, 2)
    out = fill(T(NaN), size(G)..., k)
    Threads.@threads for I in CartesianIndices(G.nstreams)
        rows = stream_rows(G, Tuple(I)...)
        isempty(rows) && continue
        ρsum = sum(view(G.density, rows))
        for c in 1:k
            out[I, c] = sum(G.density[r] * data[r, c] for r in rows) / ρsum
        end
    end
    return field == :divergence ? dropdims(out; dims = 4) : out
end

"""
    StreamVelocityGrid(xs, ys, zs, nstreams, velocity_type; deformation = false)

Empty grid for the streams counted in `nstreams[ix, iy, iz]`: the rows are laid out, and `nstreams` is reset to zero, to count the rows filled by `fill_streams!`. With `deformation = true` it has room for the divergence, shear and vorticity of every stream.
"""
function StreamVelocityGrid(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, nstreams::AbstractArray{<:Integer, 3},
                            ::Type{T}; deformation::Bool = false) where {T}
    offsets = Vector{Int}(undef, length(nstreams) + 1)
    offsets[1] = 1
    for I in eachindex(nstreams)
        offsets[I + 1] = offsets[I] + nstreams[I]
    end
    S = offsets[end] - 1
    xs, ys, zs = promote(xs, ys, zs)
    Sd = deformation ? S : 0
    return StreamVelocityGrid(xs, ys, zs, zeros(Int32, size(nstreams)), offsets, Matrix{T}(undef, S, 3), Vector{T}(undef, S),
                              Vector{T}(undef, Sd), Matrix{T}(undef, Sd, 6), Matrix{T}(undef, Sd, 3))
end

"""
    fill_streams!(G, estimator; threads = true)

Write the density and the velocity of every stream of the estimator at the grid points of `G` into the next free rows of the grid point, counting them in `G.nstreams`; also its divergence, shear and vorticity, when `G` has room for them.
"""
function fill_streams!(G::StreamVelocityGrid, est::PS_DTFE; threads = true)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    ρ, V, L = est.rho, est.velocities, est.L
    lin = LinearIndices(G.nstreams)
    images = default_images(est)
    deform = has_deformation(G)
    foreach_raster_tet(est, G.xs, G.ys, G.zs; threads, images) do t, v, zr, o, E, B, coords, wrapped
        @inbounds fρ = (ρ[v[1]], ρ[v[2]], ρ[v[3]], ρ[v[4]])
        @inbounds fv = ntuple(d -> (V[v[1], d], V[v[2], d], V[v[3], d], V[v[4], d]), 3)
        # velocity deformation, constant on the tetrahedron, as in the point queries
        J = velocity_gradient(V, v, B)
        div, σ, ω = tr(J), shear_tensor(J), vorticity_vector(J)
        raster_visit!(G.xs, G.ys, G.zs, zr, o, E, B, L, coords, wrapped; images) do ix, iy, iz, λ1, λ2, λ3
            @inbounds begin
                I = lin[ix, iy, iz]
                row = G.offsets[I] + G.nstreams[I]
                row < G.offsets[I + 1] || error("more streams at grid point ($ix, $iy, $iz) than counted")
                G.nstreams[I] += 1
                interp(f) = f[1] + λ1 * (f[2] - f[1]) + λ2 * (f[3] - f[1]) + λ3 * (f[4] - f[1])
                G.density[row] = interp(fρ)
                G.velocity[row, 1], G.velocity[row, 2], G.velocity[row, 3] = interp(fv[1]), interp(fv[2]), interp(fv[3])
                if deform
                    G.divergence[row] = div
                    G.shear[row, 1], G.shear[row, 2], G.shear[row, 3] = σ[1, 1], σ[2, 2], σ[3, 3]
                    G.shear[row, 4], G.shear[row, 5], G.shear[row, 6] = σ[1, 2], σ[1, 3], σ[2, 3]
                    G.vorticity[row, 1], G.vorticity[row, 2], G.vorticity[row, 3] = ω[1], ω[2], ω[3]
                end
            end
        end
    end
    return G
end

"Check that every grid point received the streams counted when the grid was laid out."
function check_filled(G::StreamVelocityGrid)
    all(I -> G.nstreams[I] == G.offsets[I + 1] - G.offsets[I], eachindex(G.nstreams)) ||
        error("the streams found while filling the grid differ from those counted")
    return G
end

"""
    velocity_grid(xs, ys, zs, estimator; threads = true, velocity_type = Float64, deformation = false)

Velocities of all streams on the regular grid xs × ys × zs, as a `StreamVelocityGrid` `G`: `G[ix, iy, iz]` equals `velocity([x, y, z], estimator)` at the grid point, up to the order of the streams. With `deformation = true` the divergence, shear and vorticity of every stream are stored too (see `stream_divergence`, `stream_shear`, `stream_vorticity`, `stream_average`). Computed by rasterising the tetrahedra twice, to count the streams of every grid point and to fill in their values. The grid takes 12 bytes per grid point and `4 sizeof(velocity_type)` bytes per stream, `14 sizeof(velocity_type)` with `deformation`.
"""
function velocity_grid(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, est::PS_DTFE; threads = true,
                       velocity_type::Type{<:AbstractFloat} = Float64, deformation::Bool = false)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    ns = zeros(Int32, length(xs), length(ys), length(zs))
    rasterise!(nothing, ns, est, xs, ys, zs; threads)
    G = StreamVelocityGrid(xs, ys, zs, ns, velocity_type; deformation)
    return check_filled(fill_streams!(G, est; threads))
end

"""
    bucket_by_slab(extent, n, zs, parts, L)

The items 1:n whose z-extent `extent(i)` (a tuple `(zlo, zhi)`, or `nothing` to skip the item) has a periodic image
overlapping the z-planes `zs[zr]` of each part `zr` in `parts`: one list per part, in increasing order. Every item is
tested once, in parallel, so that the threads of a rasterisation need not each scan all items.
"""
function bucket_by_slab(extent, n, zs, parts, L)
    I = n <= typemax(UInt32) ? UInt32 : UInt64
    ranges = chunks(n, 8 * Threads.nthreads())
    lists = [[I[] for _ in parts] for _ in ranges]
    Threads.@threads for c in eachindex(ranges)
        @inbounds for i in ranges[c]
            e = extent(i)
            e === nothing && continue
            for (b, zr) in enumerate(parts)
                overlaps(e[1], e[2], zs, zr, L) && push!(lists[c][b], i)
            end
        end
    end
    buckets = Vector{Vector{I}}(undef, length(parts))
    for b in eachindex(parts)
        buckets[b] = reduce(vcat, (l[b] for l in lists); init = I[])
        foreach(l -> empty!(l[b]), lists)                    # release the partial lists
    end
    return buckets
end

"""
    ps_dtfe_grid(coords_q, coords_x, m, sim_box, xs, ys, zs; tessellation = :kuhn, vertex_density = :mass, pad, jitter, threads = true)

Density and number of streams on the regular grid xs × ys × zs directly from the particles, without constructing an estimator: the tetrahedra are unwrapped and rasterised one at a time. For the Kuhn split they are generated from the lattice on the fly, so the memory needed is that of the particles and the grid; for the Delaunay tessellation only its `(T, 4)` connectivity is kept. Equal to `density_grid` and `numberOfStreams_grid` of the estimator `PS_DTFE_periodic(coords_q, coords_x, m, sim_box; tessellation, vertex_density, pad)`. Use this for simulations whose estimator does not fit in memory.
"""
function ps_dtfe_grid(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox,
                      xs::AbstractRange, ys::AbstractRange, zs::AbstractRange;
                      tessellation::Symbol = :kuhn, vertex_density::Symbol = :mass, pad = nothing, jitter = 1e-5, threads = true,
                      slabs = nothing)
    all(valid_range, (xs, ys, zs)) || throw(ArgumentError("the grid ranges must be increasing"))
    L, N = sim_box.L, sim_box.Ni
    x = Matrix{Float64}(coords_x)
    ρ, ns = zeros(length(xs), length(ys), length(zs)), zeros(Int, length(xs), length(ys), length(zs))
    coords = (collect(xs), collect(ys), collect(zs))
    wrapped = map(c -> mod.(c, L), coords)
    parts = threads ? chunks(length(zs)) : [1:length(zs)]

    "Rasterise the tetrahedron with corners v, if it reaches the planes zr."
    @inline function add!(v, rho, zr)
        o = first_corner(x, v[1], L)
        E = tet_edges(x, v..., L)
        zlo, zhi = o[3] + min(0, E[3, 1], E[3, 2], E[3, 3]), o[3] + max(0, E[3, 1], E[3, 2], E[3, 3])
        overlaps(zlo, zhi, zs, zr, L) || return
        D, B = edge_inverse(E)
        D == 0 && return
        @inbounds raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, B, (rho[v[1]], rho[v[2]], rho[v[3]], rho[v[4]]), L, coords, wrapped)
    end

    if tessellation == :kuhn
        site = lattice_sites(coords_q, sim_box)
        check_masses(m, N^3)
        # every Kuhn star has Lagrangian volume 4h³, so :lagrangian is ρ̄ 4h³ / V(W), i.e. :mass with the mass ρ̄ h³
        vertex_density in (:mass, :lagrangian) || throw(ArgumentError("vertex_density must be :mass or :lagrangian"))
        ρ̄ = (m isa Number ? m * N^3 : sum(m)) / L^3
        rho = lattice_density(x, site, vertex_density == :mass ? m : ρ̄ * (L / N)^3, L)
        # the six Kuhn tetrahedra of a cube share its first corner, and lie within the z-extent of its eight corners
        cube(c) = ((c - 1) % N + 1, ((c - 1) ÷ N) % N + 1, (c - 1) ÷ N^2 + 1)
        function cube_extent(c)
            i, j, k = cube(c)
            @inbounds begin
                z0 = x[site[i, j, k], 3]
                o3 = z0 - fld(z0, L) * L
                lo = hi = 0.0
                for cz in 0:1, cy in 0:1, cx in 0:1
                    dz = unwrap_s(x[site[mod1(i + cx, N), mod1(j + cy, N), mod1(k + cz, N)], 3] - z0, L)
                    lo, hi = min(lo, dz), max(hi, dz)
                end
            end
            return (o3 + lo, o3 + hi)
        end
        buckets = bucket_by_slab(cube_extent, N^3, zs, parts, L)
        Threads.@threads for b in eachindex(parts)
            zr = parts[b]
            @inbounds for c in buckets[b], path in KUHN_PATHS
                i, j, k = cube(c)
                add!(kuhn_tetrahedron(site, i, j, k, path), rho, zr)
            end
        end
    else
        tets, VL = tessellate(coords_q, sim_box, tessellation; pad, jitter, slabs)
        rho = vertex_densities_periodic(x, tets, VL, m, L, vertex_density)
        function tet_extent(t)
            @inbounds v = (tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4])
            o, E = first_corner(x, v[1], L), tet_edges(x, v..., L)
            return (o[3] + min(0, E[3, 1], E[3, 2], E[3, 3]), o[3] + max(0, E[3, 1], E[3, 2], E[3, 3]))
        end
        buckets = bucket_by_slab(tet_extent, size(tets, 1), zs, parts, L)
        Threads.@threads for b in eachindex(parts)
            zr = parts[b]
            @inbounds for t in buckets[b]
                add!((tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4]), rho, zr)
            end
        end
    end
    return ρ, ns
end
