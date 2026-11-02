"""
    SimBox(L, Ni)

Holds the size `L` of the periodic simulation box and the number of particles `Ni` per side length.
# Examples
```julia-repl
julia> SimBox(100, 64)
SimBox(100.0, 64)
```
"""
struct SimBox
    L::Float64
    Ni::Int64
end

"""
    unwrap_s(s, L)

The periodic image of the displacement `s` in [-L/2, L/2); `s` itself if it is already in that interval.
"""
@inline unwrap_s(s, L) = s - L * floor(s / L + 0.5)

"""
    lattice_sites(coords_q, sim_box)

The Lagrangian lattice site of every particle, from its initial position: `site[i, j, k]` is the row of `coords_q` of the particle that started at lattice site (i, j, k). The lattice may be offset from the origin of the box; the initial displacements must be smaller than half a lattice spacing. The rows of `coords_q` may be in any order.
"""
function lattice_sites(coords_q, sim_box::SimBox)
    N, L = sim_box.Ni, sim_box.L
    h = L / N
    size(coords_q) == (N^3, 3) || throw(ArgumentError("coords_q must be an (Ni³, 3) matrix, got size $(size(coords_q))"))
    # offset of the lattice in units of the spacing, from the circular mean of the position within a cell
    offset = [angle(sum(cis(2π * coords_q[p, d] / h) for p in axes(coords_q, 1))) / 2π for d in 1:3]
    site = zeros(Int, N, N, N)
    @inbounds for p in axes(coords_q, 1)
        i, j, k = (mod(round(Int, coords_q[p, d] / h - offset[d]), N) + 1 for d in 1:3)
        site[i, j, k] == 0 || throw(ArgumentError("two particles start at lattice site ($i, $j, $k): the initial positions " *
                                                  "must lie within half a lattice spacing of a regular Ni³ lattice"))
        site[i, j, k] = p
    end
    return site
end

"The six orderings of the axes; each is a monotone path from (0, 0, 0) to (1, 1, 1), i.e. one Kuhn tetrahedron."
const KUHN_PATHS = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))

"""
    kuhn_tetrahedron(site, i, j, k, path)

Particle rows of the corners of the Kuhn tetrahedron of the lattice cube (i, j, k) that follows `path`: the first corner is the cube's origin and the last its opposite corner (periodic lattice).
"""
@inline function kuhn_tetrahedron(site, i, j, k, path)
    N = size(site, 1)
    c = (i, j, k)
    v0 = site[i, j, k]
    c = Base.setindex(c, c[path[1]] + 1, path[1]); v1 = site[mod1(c[1], N), mod1(c[2], N), mod1(c[3], N)]
    c = Base.setindex(c, c[path[2]] + 1, path[2]); v2 = site[mod1(c[1], N), mod1(c[2], N), mod1(c[3], N)]
    c = Base.setindex(c, c[path[3]] + 1, path[3]); v3 = site[mod1(c[1], N), mod1(c[2], N), mod1(c[3], N)]
    return (v0, v1, v2, v3)
end

"""
    lattice_tetrahedra(site)

Kuhn (Freudenthal) triangulation of the periodic Ni³ Lagrangian lattice: every cube is split into the six tetrahedra that share its main diagonal, so that every particle is a vertex of exactly 24 tetrahedra. Returns a (6 Ni³, 4) matrix of particle rows.
"""
function lattice_tetrahedra(site)
    N = size(site, 1)
    tets = Matrix{Int32}(undef, 6N^3, 4)
    t = 0
    @inbounds for k in 1:N, j in 1:N, i in 1:N, path in KUHN_PATHS
        t += 1
        tets[t, 1], tets[t, 2], tets[t, 3], tets[t, 4] = kuhn_tetrahedron(site, i, j, k, path)
    end
    return tets
end

"""
    tet_edges(x, v0, v1, v2, v3, L)

Edge matrix of a tetrahedron with corners at the rows v0…v3 of `x`: the columns are the edges from the first corner, each placed at the periodic image nearest to the first corner. Only assumes that the edges are shorter than L/2.
"""
@inline function tet_edges(x, v0, v1, v2, v3, L)
    @inbounds SMatrix{3, 3, Float64}(
        unwrap_s(x[v1, 1] - x[v0, 1], L), unwrap_s(x[v1, 2] - x[v0, 2], L), unwrap_s(x[v1, 3] - x[v0, 3], L),
        unwrap_s(x[v2, 1] - x[v0, 1], L), unwrap_s(x[v2, 2] - x[v0, 2], L), unwrap_s(x[v2, 3] - x[v0, 3], L),
        unwrap_s(x[v3, 1] - x[v0, 1], L), unwrap_s(x[v3, 2] - x[v0, 2], L), unwrap_s(x[v3, 3] - x[v0, 3], L))
end

"""
    tet_edges(x, v0, v1, v2, v3)

Edge matrix of a tetrahedron with corners at the rows v0…v3 of `x`: the columns are the edges from the first corner.
"""
@inline function tet_edges(x, v0, v1, v2, v3)
    @inbounds SMatrix{3, 3, Float64}(
        x[v1, 1] - x[v0, 1], x[v1, 2] - x[v0, 2], x[v1, 3] - x[v0, 3],
        x[v2, 1] - x[v0, 1], x[v2, 2] - x[v0, 2], x[v2, 3] - x[v0, 3],
        x[v3, 1] - x[v0, 1], x[v3, 2] - x[v0, 2], x[v3, 3] - x[v0, 3])
end

"First corner of the tetrahedron at the periodic image in [0, L)³."
@inline first_corner(x, v0, L) = @inbounds SVector(x[v0, 1] - fld(x[v0, 1], L) * L, x[v0, 2] - fld(x[v0, 2], L) * L,
                                                   x[v0, 3] - fld(x[v0, 3], L) * L)

"""
    lattice_density(x, site, m, L)

Density at the particles, ρ_i = 4 m_i / Σ_{t ∋ i} V_t, from the Eulerian volumes of the 24 tetrahedra around each particle. `m` is a single mass or a vector of particle masses.
"""
function lattice_density(x, site::Array{Int, 3}, m, L)
    N = size(site, 1)
    vol = zeros(size(x, 1))
    @inbounds for k in 1:N, j in 1:N, i in 1:N, path in KUHN_PATHS
        v = kuhn_tetrahedron(site, i, j, k, path)
        V = abs(det(tet_edges(x, v..., L))) / 6
        for c in v
            vol[c] += V
        end
    end
    return 4 .* m ./ vol
end
