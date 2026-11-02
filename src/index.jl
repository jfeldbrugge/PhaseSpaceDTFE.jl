"""
    CellIndex

Uniform G³ cell index of a set of tetrahedra over the box [0, L)³, in compressed form: the tetrahedra whose bounding box overlaps the cell k = i + (j - 1) G + (l - 1) G² are `tets[offsets[k]:offsets[k + 1] - 1]`.
"""
struct CellIndex
    G::Int
    L::Float64
    offsets::Vector{Int}
    tets::Vector{Int32}
end

"""
    CellIndex(lo, hi, L, G)

Index of the tetrahedra with bounding boxes `lo[:, t]`…`hi[:, t]` (3 × T matrices of coordinates).
"""
function CellIndex(lo::AbstractMatrix, hi::AbstractMatrix, L, G)
    T, h = size(lo, 2), L / G
    cell(x) = clamp(floor(Int, x / h) + 1, 1, G)
    counts = zeros(Int, G^3)
    @inbounds for t in 1:T
        for l in cell(lo[3, t]):cell(hi[3, t]), j in cell(lo[2, t]):cell(hi[2, t]), i in cell(lo[1, t]):cell(hi[1, t])
            counts[i + (j - 1) * G + (l - 1) * G^2] += 1
        end
    end
    offsets = Vector{Int}(undef, G^3 + 1)
    offsets[1] = 1
    for k in 1:G^3
        offsets[k + 1] = offsets[k] + counts[k]
    end
    next = offsets[1:end - 1]
    tets = Vector{Int32}(undef, offsets[end] - 1)
    @inbounds for t in 1:T
        for l in cell(lo[3, t]):cell(hi[3, t]), j in cell(lo[2, t]):cell(hi[2, t]), i in cell(lo[1, t]):cell(hi[1, t])
            k = i + (j - 1) * G + (l - 1) * G^2
            tets[next[k]] = t
            next[k] += 1
        end
    end
    return CellIndex(G, Float64(L), offsets, tets)
end

"Tetrahedra whose bounding box overlaps the cell of the point (x, y, z) in [0, L)³."
@inline function cell_tetrahedra(index::CellIndex, x, y, z)
    h, G = index.L / index.G, index.G
    k = clamp(floor(Int, x / h) + 1, 1, G) + (clamp(floor(Int, y / h) + 1, 1, G) - 1) * G +
        (clamp(floor(Int, z / h) + 1, 1, G) - 1) * G^2
    return @inbounds @view index.tets[index.offsets[k]:index.offsets[k + 1] - 1]
end
