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
    raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, B, f, L; periodic = true)

Add one tetrahedron — first corner `o`, edge matrix `E` (columns), its inverse `B` and vertex values `f` — to the grid arrays `ρ[ix, iy, iz]` and `ns[ix, iy, iz]` at every periodic image that overlaps the grid, for the z-planes with indices in `zr`. Every plane cuts the tetrahedron in a convex polygon; for every grid row the x-interval of the polygon is computed (scanline) and the grid points in it, with a tolerance of `PAD` grid spacings, are tested with the barycentric coordinates, so that the result equals point queries at the grid points taken modulo the box. With `periodic = false` the image shifted by k L along an axis is only added at the grid points g with `fld(g, L) == k`, i.e. the tetrahedron is only found at the points whose image modulo the box lies in it, as for the point queries of a non-periodic estimator.
"""
function raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, B, f, L, (xv, yv, zv) = (xs, ys, zs); periodic = true)
    c = (o, o + E[:, 1], o + E[:, 2], o + E[:, 3])
    lo, hi = min.(c[1], c[2], c[3], c[4]), max.(c[1], c[2], c[3], c[4])
    Δx, Δy, Δz = gridstep(xs), gridstep(ys), gridstep(zs)
    px, py = MVector{10, Float64}(undef), MVector{10, Float64}(undef)
    @inbounds for kz in shift_range(lo[3], hi[3], zs, L)
        sz = kz * L
        iz0, iz1 = index_range(lo[3] + sz, hi[3] + sz, zs, Δz)
        iz0, iz1 = max(iz0, first(zr)), min(iz1, last(zr))
        for iz in iz0:iz1
            periodic || fld(zv[iz], L) == kz || continue
            zp = zv[iz] - sz                                  # plane in the coordinates of the tetrahedron
            n = section!(px, py, c, zp)
            n == 0 && continue
            ylo, yhi = extent(py, n)
            for ky in shift_range(ylo, yhi, ys, L)
                sy = ky * L
                iy0, iy1 = index_range(ylo + sy, yhi + sy, ys, Δy)
                for iy in iy0:iy1
                    periodic || fld(yv[iy], L) == ky || continue
                    yp = yv[iy] - sy
                    xlo, xhi = row_extent(px, py, n, yp)
                    if xlo > xhi                                    # row within PAD of a corner of the polygon
                        xlo, xhi = extent(px, n)
                    end
                    for kx in shift_range(xlo, xhi, xs, L)
                        sx = kx * L
                        ix0, ix1 = index_range(xlo + sx, xhi + sx, xs, Δx)
                        for ix in ix0:ix1
                            periodic || fld(xv[ix], L) == kx || continue
                            d =SVector(xv[ix] - sx - o[1], yp - o[2], zp - o[3])
                            λ = B * d
                            inside(λ[1], λ[2], λ[3]) || continue
                            ρ[ix, iy, iz] += f[1] + λ[1] * (f[2] - f[1]) + λ[2] * (f[3] - f[1]) + λ[3] * (f[4] - f[1])
                            ns[ix, iy, iz] += 1
                        end
                    end
                end
            end
        end
    end
    return nothing
end

"Split 1:n into one contiguous chunk per thread."
chunks(n, k = Threads.nthreads()) = [round(Int, (i - 1) * n / k) + 1:round(Int, i * n / k) for i in 1:k if round(Int, i * n / k) >= round(Int, (i - 1) * n / k) + 1]

"Whether the tetrahedron with z-extent [zlo, zhi] has an image overlapping the z-planes zs[zr]."
@inline function overlaps(zlo, zhi, zs, zr, L)
    z0, z1 = zs[first(zr)], zs[last(zr)]
    return !isempty(ceil(Int, (z0 - zhi) / L):floor(Int, (z1 - zlo) / L))
end

"""
    rasterise(estimator, xs, ys, zs; threads = true)

Density and number of streams on the regular grid xs × ys × zs, as arrays indexed `[ix, iy, iz]`, by rasterising the tetrahedra of the estimator. The grid may extend beyond the box. With `threads = true` the z-planes are split over the available threads.
"""
function rasterise(est::PS_DTFE, xs::AbstractRange, ys::AbstractRange, zs::AbstractRange; threads = true)
    ρ, ns = zeros(length(xs), length(ys), length(zs)), zeros(Int, length(xs), length(ys), length(zs))
    any(isempty, (xs, ys, zs)) && return ρ, ns
    all(valid_range, (xs, ys, zs)) || throw(ArgumentError("the grid ranges must be increasing"))
    L, x = est.L, est.positions
    coords = (collect(xs), collect(ys), collect(zs))           # exact range elements, without range arithmetic
    originals = findall(est.original)
    parts = threads ? chunks(length(zs)) : [1:length(zs)]
    # edges as stored: periodically unwrapped, or the plain differences of a non-periodic estimator
    edges(v) = est.periodic ? tet_edges(x, v..., L) : tet_edges(x, v...)
    # z-extent of every tetrahedron, computed once, and the tetrahedra that reach each slab of planes
    zext = Vector{NTuple{2, Float64}}(undef, length(originals))
    Threads.@threads for r in eachindex(originals)
        t = originals[r]
        @inbounds E = edges((est.simplices[t, 1], est.simplices[t, 2], est.simplices[t, 3], est.simplices[t, 4]))
        @inbounds zext[r] = (est.origin[3, t] + min(0, E[3, 1], E[3, 2], E[3, 3]), est.origin[3, t] + max(0, E[3, 1], E[3, 2], E[3, 3]))
    end
    buckets = [Int32[] for _ in parts]
    for (r, (zlo, zhi)) in enumerate(zext), (b, zr) in enumerate(parts)
        overlaps(zlo, zhi, zs, zr, L) && push!(buckets[b], r)
    end
    Threads.@threads for b in eachindex(parts)
        zr = parts[b]
        @inbounds for r in buckets[b]
            t = originals[r]
            isnan(est.Binv[1, t]) && continue
            v = (est.simplices[t, 1], est.simplices[t, 2], est.simplices[t, 3], est.simplices[t, 4])
            o = SVector(est.origin[1, t], est.origin[2, t], est.origin[3, t])
            E = edges(v)
            B = SMatrix{3, 3, Float64}(@view est.Binv[:, t])
            raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, B, (est.rho[v[1]], est.rho[v[2]], est.rho[v[3]], est.rho[v[4]]), L, coords;
                                periodic = est.periodic)
        end
    end
    return ρ, ns
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
    ps_dtfe_grid(coords_q, coords_x, m, sim_box, xs, ys, zs; tessellation = :delaunay, vertex_density = :mass, pad, jitter, threads = true)

Density and number of streams on the regular grid xs × ys × zs directly from the particles, without constructing an estimator: the tetrahedra are unwrapped and rasterised one at a time. For the Kuhn split they are generated from the lattice on the fly, so the memory needed is that of the particles and the grid; for the Delaunay tessellation only its `(T, 4)` connectivity is kept. Equal to `density_grid` and `numberOfStreams_grid` of the estimator `PS_DTFE_periodic(coords_q, coords_x, m, sim_box; tessellation, vertex_density, pad)`. Use this for simulations whose estimator does not fit in memory.
"""
function ps_dtfe_grid(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox,
                      xs::AbstractRange, ys::AbstractRange, zs::AbstractRange;
                      tessellation::Symbol = :delaunay, vertex_density::Symbol = :mass, pad = nothing, jitter = 1e-5, threads = true)
    all(valid_range, (xs, ys, zs)) || throw(ArgumentError("the grid ranges must be increasing"))
    L, N = sim_box.L, sim_box.Ni
    x = Matrix{Float64}(coords_x)
    ρ, ns = zeros(length(xs), length(ys), length(zs)), zeros(Int, length(xs), length(ys), length(zs))
    coords = (collect(xs), collect(ys), collect(zs))
    parts = threads ? chunks(length(zs)) : [1:length(zs)]

    "Rasterise the tetrahedron with corners v, if it reaches the planes zr."
    @inline function add!(v, rho, zr)
        o = first_corner(x, v[1], L)
        E = tet_edges(x, v..., L)
        zlo, zhi = o[3] + min(0, E[3, 1], E[3, 2], E[3, 3]), o[3] + max(0, E[3, 1], E[3, 2], E[3, 3])
        overlaps(zlo, zhi, zs, zr, L) || return
        det(E) == 0 && return
        @inbounds raster_tetrahedron!(ρ, ns, xs, ys, zs, zr, o, E, inv(E), (rho[v[1]], rho[v[2]], rho[v[3]], rho[v[4]]), L, coords)
    end

    if tessellation == :kuhn
        site = lattice_sites(coords_q, sim_box)
        check_masses(m, N^3)
        # every Kuhn star has Lagrangian volume 4h³, so :lagrangian is ρ̄ 4h³ / V(W), i.e. :mass with the mass ρ̄ h³
        vertex_density in (:mass, :lagrangian) || throw(ArgumentError("vertex_density must be :mass or :lagrangian"))
        ρ̄ = (m isa Number ? m * N^3 : sum(m)) / L^3
        rho = lattice_density(x, site, vertex_density == :mass ? m : ρ̄ * (L / N)^3, L)
        Threads.@threads for zr in parts
            @inbounds for k in 1:N, j in 1:N, i in 1:N, path in KUHN_PATHS
                add!(kuhn_tetrahedron(site, i, j, k, path), rho, zr)
            end
        end
    else
        tets, VL = tessellate(coords_q, sim_box, tessellation; pad, jitter)
        rho = vertex_densities_periodic(x, tets, VL, m, L, vertex_density)
        Threads.@threads for zr in parts
            @inbounds for t in axes(tets, 1)
                add!((tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4]), rho, zr)
            end
        end
    end
    return ρ, ns
end
