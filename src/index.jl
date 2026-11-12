"""
    CellIndex

Uniform G³ cell index of a set of tetrahedra over the box [0, L)³, in compressed form: the tetrahedra that intersect
the cell k = i + (j - 1) G + (l - 1) G² are `tets[offsets[k]:offsets[k + 1] - 1]`, in increasing order. For a periodic
index the cells are taken at any periodic image; otherwise the cells at the faces extend to infinity, so that points
outside the box find the tetrahedra there.
"""
struct CellIndex
    G::Int
    L::Float64
    offsets::Vector{Int}
    tets::Vector{Int32}
end

"Margin, in cell sizes, by which the cells are widened in the overlap test, so that rounding never loses a tetrahedron."
const CELL_MARGIN = 1e-6

"Split 1:n into k contiguous, non-empty chunks (at most n)."
chunks(n, k = Threads.nthreads()) = [round(Int, (i - 1) * n / k) + 1:round(Int, i * n / k) for i in 1:k if round(Int, i * n / k) >= round(Int, (i - 1) * n / k) + 1]

"""
    separating_axes(c)

Candidate separating axes of the tetrahedron with corners `c` and an axis-aligned box, apart from the coordinate axes:
the four face normals and the cross products of the six edges with the coordinate axes.
"""
@inline function separating_axes(c)
    e = (c[2] - c[1], c[3] - c[1], c[4] - c[1], c[3] - c[2], c[4] - c[2], c[4] - c[3])
    faces = (cross(e[1], e[2]), cross(e[1], e[3]), cross(e[2], e[3]), cross(e[4], e[5]))
    return (faces..., map(a -> SVector(0.0, -a[3], a[2]), e)..., map(a -> SVector(a[3], 0.0, -a[1]), e)...,
            map(a -> SVector(-a[2], a[1], 0.0), e)...)
end

"""
    tet_cells!(keys, c, h, G, periodic = true)

Append to `keys` the cells (of size h, G per side) that the tetrahedron with corners `c` intersects, each once, and
return their number. Cells are tested exactly with the separating axis theorem, widened by `CELL_MARGIN`. With
`periodic = true` a cell outside the box is taken at its periodic image in the box; otherwise it is clamped to the
nearest cell at a face of the box.
"""
function tet_cells!(keys, c, h, G, periodic = true)
    lo, hi = min.(c[1], c[2], c[3], c[4]), max.(c[1], c[2], c[3], c[4])
    i0 = floor.(Int, lo ./ h .- CELL_MARGIN)
    i1 = floor.(Int, hi ./ h .+ CELL_MARGIN)
    ax = separating_axes(c)
    ivals = map(a -> extrema((dot(a, c[1]), dot(a, c[2]), dot(a, c[3]), dot(a, c[4]))), ax)
    half = h / 2 * (1 + 2CELL_MARGIN)
    radii = map(a -> half * (abs(a[1]) + abs(a[2]) + abs(a[3])), ax)
    n0 = length(keys)
    @inbounds for k in i0[3]:i1[3], j in i0[2]:i1[2], i in i0[1]:i1[1]
        m = SVector((i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h)
        separated = false
        for q in eachindex(ax)
            s = dot(ax[q], m)
            if s - radii[q] > ivals[q][2] || s + radii[q] < ivals[q][1]
                separated = true
                break
            end
        end
        separated && continue
        cell(a) = periodic ? mod(a, G) : clamp(a, 0, G - 1)
        push!(keys, UInt32(cell(i) + cell(j) * G + cell(k) * G^2 + 1))
    end
    # two periodic images of the tetrahedron may reach one cell, and clamping maps several cells to one
    if periodic ? any(i1 .- i0 .+ 1 .> G) : (any(i0 .< 0) || any(i1 .>= G))
        new = sort!(unique(keys[n0 + 1:end]))
        resize!(keys, n0)
        append!(keys, new)
    end
    return length(keys) - n0
end

"""
    build_index(corners, T, L, G; periodic = true)

Index of the tetrahedra 1:T over the box of size `L` with G³ cells, where `corners(t)` returns the four corners of
tetrahedron t (at any periodic image for a periodic index), or `nothing` for a tetrahedron that is not indexed. The
overlaps are computed in parallel.
"""
function build_index(corners, T::Integer, L, G::Integer; periodic::Bool = true)
    G >= 1 && G^3 < typemax(UInt32) || throw(ArgumentError("G must lie between 1 and 1625, got $G"))
    h = L / G
    parts = chunks(T, 8 * Threads.nthreads())
    keys = [UInt32[] for _ in parts]
    counts = [Vector{UInt32}(undef, length(r)) for r in parts]
    Threads.@threads for p in eachindex(parts)
        ks, cs = keys[p], counts[p]
        for (a, t) in enumerate(parts[p])
            c = corners(t)
            cs[a] = c === nothing ? 0 : tet_cells!(ks, c, h, G, periodic)
        end
    end
    offsets = zeros(Int, G^3 + 1)
    for ks in keys, k in ks
        @inbounds offsets[k + 1] += 1
    end
    offsets[1] = 1
    cumsum!(offsets, offsets)
    next = offsets[1:end - 1]
    tets = Vector{Int32}(undef, offsets[end] - 1)
    for p in eachindex(parts)
        e, ks = 0, keys[p]
        @inbounds for (a, t) in enumerate(parts[p]), _ in 1:counts[p][a]
            e += 1
            k = ks[e]
            tets[next[k]] = t
            next[k] += 1
        end
        keys[p] = UInt32[]                               # release the keys of this chunk
    end
    return CellIndex(G, Float64(L), offsets, tets)
end

"Tetrahedra that intersect the cell of the point (x, y, z) in [0, L]³."
@inline function cell_tetrahedra(index::CellIndex, x, y, z)
    h, G = index.L / index.G, index.G
    k = clamp(floor(Int, x / h) + 1, 1, G) + (clamp(floor(Int, y / h) + 1, 1, G) - 1) * G +
        (clamp(floor(Int, z / h) + 1, 1, G) - 1) * G^2
    return @inbounds @view index.tets[index.offsets[k]:index.offsets[k + 1] - 1]
end
