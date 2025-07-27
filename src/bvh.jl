"""
    BVH(data, box, depth, left, right)

A Bounding Volume Hierarchy (BVH) tree structure on a set of geometric objects.
"""
struct BVH
    data::Union{Vector{Int}, Nothing}
    box::Matrix{Float64}
    depth::Int
    left::Union{BVH, Nothing}
    right::Union{BVH, Nothing}
end

"""
    BVH(data, box::Matrix{Float64}, points::Matrix{Float64}, simplices, depth::Int)

Generates a Bounding Volume Hierarchy (BVH) tree structure.
"""
function BVH(data, box::Matrix{Float64}, points::Matrix{Float64}, simplices, depth::Int)
    if depth != 0
        dim = (depth % size(box, 1)) + 1
        @inbounds mins = vec(minimum(@view(points[@view(simplices[data, :]), dim]), dims = 2))
        @inbounds maxs = vec(maximum(@view(points[@view(simplices[data, :]), dim]), dims = 2))
            
        @inbounds div = sum(@view(box[dim,:])) / length(@view(box[dim,:]))
        L = mins .<= div
        R = maxs .>= div

        Lbox = copy(box)
        @inbounds Lbox[dim,:] = [Lbox[dim, 1], div]
        Rbox = copy(box)
        @inbounds Rbox[dim,:] = [div, Rbox[dim, 2]]

        @inbounds BVH(nothing, box, depth, BVH(data[L], Lbox, points, simplices, depth - 1), 
                                            BVH(data[R], Rbox, points, simplices, depth - 1))
    else 
        BVH(data, box, depth, nothing, nothing)
    end
end

"""
    findCandidateSimplices(p::Vector{Float64}, BVH_tree::BVH)

Given a Bounding Volume Hierarchy, findCandidateSimplices find a set simplices that may contain the point p. The function output indices of candidate simplices.
"""
function findCandidateSimplices(p::Vector{Float64}, BVH_tree::BVH)
    dim = (BVH_tree.depth % size(BVH_tree.box, 1)) + 1
    if BVH_tree.depth == 0
        return BVH_tree.data
    elseif p[dim] < BVH_tree.left.box[dim, 2]
        return findCandidateSimplices(p, BVH_tree.left)
    elseif p[dim] > BVH_tree.left.box[dim, 2]
        return findCandidateSimplices(p, BVH_tree.right)
    else
        return union(findCandidateSimplices(p, BVH_tree.left),
                    findCandidateSimplices(p, BVH_tree.right))
    end
end

"""
    findBox(p::Vector{Float64}, BVH_tree::BVH)

Find the box in the Bounding Volume Hierarchy that contains the point p.
"""
function findBox(p::Vector{Float64}, BVH_tree::BVH)
    dim = (BVH_tree.depth % size(BVH_tree.box, 1)) + 1
    if BVH_tree.depth == 0
        println(BVH_tree.box)
    elseif p[dim] < BVH_tree.left.box[dim, 2]
        return findBox(p, BVH_tree.left)
    elseif p[dim] > BVH_tree.left.box[dim, 2]
        return findBox(p, BVH_tree.right)
    else
        findBox(p, BVH_tree.left)
        findBox(p, BVH_tree.right)
    end
end

"""
    intersection(p::Vector{Float64}, simplex)

Check whether the point p is included in the simplex.
"""
function intersection(p::Vector{Float64}, simplex)
    @inbounds barry = inv(@view(simplex[2:end,:])' .- @view(simplex[1,:])) * (p .- @view(simplex[1,:]))
    return all(barry .>= 0) & all(barry .<= 1) & (sum(barry) <= 1.)
end

"""
    findIntersections(p::Vector{Float64}, BVH_tree::BVH, points, simplices)

Find the indices of the simplices that intersect the point p.
"""
function findIntersections(p::Vector{Float64}, BVH_tree::BVH, points, simplices)
    candidates = findCandidateSimplices(p, BVH_tree)
    @inbounds filter(i -> intersection(p, @view(points[@view(simplices[i,:]),:])), candidates)
end

"""
    volume(sim, points)

Evaluate the volume of the simplex.
"""
volume(sim, points) = @inbounds abs(det(@view(points[@view(sim[2:end]),:])' .- @view(points[sim[1],:]))) /  factorial(size(points, 2))