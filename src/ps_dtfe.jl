"""
    SimBox(L, Ni)

Sets the size and resolution of the simulation box.

# Examples
```julia-repl
L = 100
Ni = 64
julia> SimBox(L, Ni)
```
"""
struct SimBox
    L::Float64
    Ni::Int64
end

function unwrap_x_(q, x, L)
    unwrap_s.(x - q, L) + q
end

unwrap_s(s, L) = mod((s + L / 2), L) - L / 2

function translate(coords_q, coords_x, L)
    coords_q_ = copy(coords_q)
    coords_x_ = copy(coords_x)

    let index_x = coords_x .< 0.
        coords_q_[index_x] .+= L
        coords_x_[index_x] .+= L
    end
        
    let index_x = coords_x .> L
        coords_q_[index_x] .-= L
        coords_x_[index_x] .-= L
    end
    return coords_q_, coords_x_
end

function frame(coords_q, coords_x, L, pad=0.05)
    for d in 1:3
        indexL = coords_x[:,d] .> (1. - pad) * L
        indexR = coords_x[:,d] .< pad * L

        newPointsL_q = coords_q[indexL, :]
        newPointsL_x = coords_x[indexL, :]
        newPointsL_q[:,d] .-= L 
        newPointsL_x[:,d] .-= L 
        
        newPointsR_q = coords_q[indexR, :]
        newPointsR_x = coords_x[indexR, :]
        newPointsR_q[:,d] .+= L 
        newPointsR_x[:,d] .+= L 

        coords_q = vcat(coords_q, newPointsL_q)        
        coords_q = vcat(coords_q, newPointsR_q)
        coords_x = vcat(coords_x, newPointsL_x)        
        coords_x = vcat(coords_x, newPointsR_x)
    end

    return coords_q, coords_x
end

function frame_velocities(coords_x, velocities, L, pad=0.05)
    for d in 1:3
        indexL = coords_x[:,d] .> (1. - pad) * L
        indexR = coords_x[:,d] .< pad * L

        newVelocitiesL = velocities[indexL, :]
        newPointsL_x   = coords_x[indexL, :]
        newPointsL_x[:,d] .-= L 
        
        newVelocitiesR = velocities[indexR, :]
        newPointsR_x   = coords_x[indexR, :]
        newPointsR_x[:,d] .+= L 

        coords_x = vcat(coords_x, newPointsL_x)        
        coords_x = vcat(coords_x, newPointsR_x)
        velocities = vcat(velocities,  newVelocitiesL)        
        velocities = vcat(velocities,  newVelocitiesR)        
    end

    return velocities
end



struct PS_DTFE
    rho::Vector{Float64}
    Drho::Matrix{Float64}
    Dv::Array{Float64}
    tree::BVH
    simplices::Matrix{Int}
    positions::Matrix{Float64}
    velocities::Matrix{Float64}
    positions_initial::Matrix{Float64}

    function PS_DTFE(positions_initial, positions, velocities, m, depth, box)
            
        inputTetGen = TetGen.RawTetGenIO{Float64}()
        inputTetGen.pointlist = copy(positions_initial')
        simplices = Int64.(tetrahedralize(inputTetGen,"Q").tetrahedronlist')
        inputTetGen = nothing
            
        dim = size(box, 1)
        
        BVH_tree = BVH(1:size(simplices, 1), box, positions, simplices, depth * dim)
        
        rho = zeros(size(positions,1))
        @inbounds for i in axes(simplices, 1)
            vol = volume(@view(simplices[i,:]), positions)
            @inbounds for index in @view(simplices[i,:])
                rho[index] += vol
            end
        end
        rho = (1. + dim) * m ./ rho
        
        Drho = zeros(size(simplices, 1), dim)
        Dv   = zeros(size(simplices, 1), dim, dim)
        
        @inbounds for i in axes(simplices, 1)
            sim = @view(simplices[i,:])
            p = @view(positions[sim,:])
            r = @view(rho[sim,:])
            v = @view(velocities[sim,:])
            A_inv = inv(@view(p[2:end,:])' .- @view(p[1,:]))'
            Drho[i,:] = A_inv * (@view(r[2:end]) .- r[1])
            Dv[i,:,:] = A_inv * ((@view(v[2:end,:])' .- @view(v[1,:]))')
        end
        
        return new(rho, Drho, Dv, BVH_tree, simplices, positions, velocities, positions_initial)
    end
end

function PS_DTFE_periodic(coords_q, coords_x, velocities, m, depth, sim_box; pad=0.05)
    coords_x = unwrap_x_(coords_q, coords_x, sim_box.L);
    coords_q_, coords_x_ = translate(coords_q, coords_x, sim_box.L)
    velocities_          = frame_velocities(coords_x_, velocities, sim_box.L, pad)
    coords_q_, coords_x_ = frame(coords_q_, coords_x_, sim_box.L, pad)

    box = [0 sim_box.L; 0 sim_box.L; 0 sim_box.L]
    PS_DTFE(coords_q_, coords_x_, velocities_, m, depth, box)
end

function PS_DTFE_periodic(coords_q, coords_x, m, depth, sim_box; pad=0.05)
    coords_x = unwrap_x_(coords_q, coords_x, sim_box.L);
    coords_q_, coords_x_ = translate(coords_q, coords_x, sim_box.L)
    coords_q_, coords_x_ = frame(coords_q_, coords_x_, sim_box.L, pad)
    Ni_ = size(coords_q_)[1]

    box = [0 sim_box.L; 0 sim_box.L; 0 sim_box.L]
    PS_DTFE(coords_q_, coords_x_, zeros(Ni_, 3), m, depth, box)
end

function density(p::Vector{Float64}, estimator::PS_DTFE)
    simplexIndices = findIntersections(p, estimator.tree, estimator.positions, estimator.simplices)

    dens = 0.
    @inbounds for simplexIndex in simplexIndices
        pointIndex = estimator.simplices[simplexIndex,1]
        dens += estimator.rho[pointIndex] + @view(estimator.Drho[simplexIndex,:])' * (p .- @view(estimator.positions[pointIndex,:]))
    end
    return dens
end

function numberOfStreams(p::Vector{Float64}, estimator::PS_DTFE)
    simplexIndices = findIntersections(p, estimator.tree, estimator.positions, estimator.simplices)
    return length(simplexIndices)
end

function velocity(p::Vector{Float64}, estimator::PS_DTFE, single_stream=false)
    simplexIndices = findIntersections(p, estimator.tree, estimator.positions, estimator.simplices)

    vs = zeros(length(simplexIndices), length(p))

    @inbounds for (i, simplexIndex) in pairs(simplexIndices)
        pointIndex = estimator.simplices[simplexIndex, 1]
        vs[i,:] = @view(estimator.velocities[pointIndex,:]) + @view(estimator.Dv[simplexIndex,:,:]) * (p .- @view(estimator.positions[pointIndex,:]))
    end

    if single_stream && size(vs) != (1,3)  # in multistream region
        return [NaN NaN NaN]
    else # in single-stream region
        return vs
    end
end

function velocitySum(p::Vector{Float64}, estimator::PS_DTFE)
    simplexIndices = findIntersections(p, estimator.tree, estimator.positions, estimator.simplices)

    vs = zeros(length(simplexIndices), length(p))

    @inbounds for (i, simplexIndex) in pairs(simplexIndices)
        pointIndex = estimator.simplices[simplexIndex, 1]
        vs[i,:] = @view(estimator.velocities[pointIndex,:]) + @view(estimator.Dv[simplexIndex,:,:]) * (p .- @view(estimator.positions[pointIndex,:]))
    end

    if size(vs) == (1,3)  # single-stream region
        return vs
    else   # multistream region
        densities = zeros(length(simplexIndices))
        for (i, simplexIndex) in pairs(simplexIndices)
            pointIndex   = estimator.simplices[simplexIndex,1]
            densities[i] = estimator.rho[pointIndex] + @view(estimator.Drho[simplexIndex,:])' * (p .- @view(estimator.positions[pointIndex,:]))
        end

        #stream-density weighted sum of velocities
        return sum(densities .* vs, dims=1) / sum(densities)
    end
end

function inSimplices(p::Vector{Float64}, estimator::PS_DTFE)
    simplexIndices = findIntersections(p, estimator.tree, estimator.positions, estimator.simplices)
end