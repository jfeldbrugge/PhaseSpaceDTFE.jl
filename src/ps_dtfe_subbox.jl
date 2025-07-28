"""
    PS_DTFE_subbox

Phase-Space DTFE subbox object containing

    N_sub::Int64
    N_target::Int64
    m::Float64
    depth::Int64
    dir::String
    L::Float64
    Ni::Int64
"""
struct PS_DTFE_subbox
    N_sub::Int64
    N_target::Int64
    m::Float64
    depth::Int64
    dir::String
    L::Float64
    Ni::Int64
end

"""
    ps_dtfe_subbox(coords_q, coords_x, m, depth, sim_box::SimBox; N_target=128, pad=0.05, dir="./ps_dtfe")

Generates an Phase-Space DTFE subbox given the initial positions and the final positions of the N-body particles. The Boundary Volume Hirarchy goes depth levels deep.
"""
function ps_dtfe_subbox(coords_q, coords_x, m, depth, sim_box::SimBox; N_target=128, pad=0.05, dir="./ps_dtfe")

    mkpath(dir)  # succeeds even if already exists
    rm(dir, force=true, recursive=true)
    mkdir(dir)

    # subboxes per side length
    N_sub = sim_box.Ni ÷ N_target
    println("Calculating $(N_sub^3) subbox estimators with $N_target^3 particles each (total simulation $(sim_box.Ni)^3 particles).")

    # pre-process coordinates: unwrapping, translating and padding
    coords_x = unwrap_x_(coords_q, coords_x, sim_box.L);
    coords_q_, coords_x_ = translate(coords_q, coords_x, sim_box.L)
    coords_q_, coords_x_ = frame(coords_q_, coords_x_, sim_box.L, pad);

    # iterate over subboxes and create subbox ps_dtfe
    idx_list = [[i, j, k] for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1]
    @showprogress for (i, j, k) in idx_list
        println("subbox i, j, k = " * string([i, j, k]))
        get_subbox_estimator(coords_q_, coords_x_, [i, j, k], N_sub, m, depth, sim_box; pad=pad, dir=dir)
    end

    return PS_DTFE_subbox(N_sub, N_target, m, depth, dir, sim_box.L, sim_box.Ni)
end

"""
    ps_dtfe_subbox(coords_q, coords_x, velocities, m, depth, sim_box::SimBox; N_target=128, pad=0.05, dir="./ps_dtfe")

Generates an Phase-Space DTFE subbox given the initial positions and the final positions and final velocities of the N-body particles. The Boundary Volume Hirarchy goes depth levels deep.
"""
function ps_dtfe_subbox(coords_q, coords_x, velocities, m, depth, sim_box::SimBox; N_target=128, pad=0.05, dir="./ps_dtfe")

    mkpath(dir)  # succeeds even if already exists
    rm(dir, force=true, recursive=true)
    mkdir(dir)

    # subboxes per side length
    N_sub = sim_box.Ni ÷ N_target
    println("Calculating $(N_sub^3) subbox estimators with $N_target^3 particles each (total simulation $(sim_box.Ni)^3 particles).")

    # pre-process coordinates: unwrapping, translating and padding
    coords_x = unwrap_x_(coords_q, coords_x, sim_box.L);
    coords_q_, coords_x_ = translate(coords_q, coords_x, sim_box.L)
    velocities_          = frame_velocities(coords_x_, velocities, sim_box.L, pad)
    coords_q_, coords_x_ = frame(coords_q_, coords_x_, sim_box.L, pad);

    # iterate over subboxes and create subbox ps_dtfe
    idx_list = [[i, j, k] for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1]
    @showprogress for (i, j, k) in idx_list
        println("subbox i, j, k = " * string([i, j, k]))
        get_subbox_estimator(coords_q_, coords_x_, velocities_, [i, j, k], N_sub, m, depth, sim_box; pad=pad, dir=dir)
    end

    return PS_DTFE_subbox(N_sub, N_target, m, depth, dir, sim_box.L, sim_box.Ni)
end

"""
    get_subbox_estimator(coords_q, coords_x, idx, N_sub, m, depth, sim_box::SimBox; pad=0.05, dir="./ps_dtfe")

Obtain the subbox estimator associated to a given position.
"""
function get_subbox_estimator(coords_q, coords_x, idx, N_sub, m, depth, sim_box::SimBox; pad=0.05, dir="./ps_dtfe")
    box = 1. / N_sub * sim_box.L .* hcat(idx, idx .+ 1)
        
    idx_x = (box[1,1] - pad*sim_box.L .< coords_x[:,1] .< box[1,2] + pad*sim_box.L) .& 
            (box[2,1] - pad*sim_box.L .< coords_x[:,2] .< box[2,2] + pad*sim_box.L) .& 
            (box[3,1] - pad*sim_box.L .< coords_x[:,3] .< box[3,2] + pad*sim_box.L) 
    Ni_ = size(coords_q[idx_x, :])[1]

    ps_dtfe = PS_DTFE(coords_q[idx_x, :], coords_x[idx_x, :], zeros(Ni_, 3), m, depth, box)
    i, j, k = idx
    serialize(dir * "/box_" * string(i) * "_" * string(j) * "_" * string(k), ps_dtfe)
end

"""
    get_subbox_estimator(coords_q, coords_x, velocities, idx, N_sub, m, depth, sim_box::SimBox; pad=0.05, dir="./ps_dtfe")

Obtain the subbox estimator associated to a given position.
"""
function get_subbox_estimator(coords_q, coords_x, velocities, idx, N_sub, m, depth, sim_box::SimBox; pad=0.05, dir="./ps_dtfe")
    box = 1. / N_sub * sim_box.L .* hcat(idx, idx .+ 1)
        
    idx_x = (box[1,1] - pad*sim_box.L .< coords_x[:,1] .< box[1,2] + pad*sim_box.L) .& 
            (box[2,1] - pad*sim_box.L .< coords_x[:,2] .< box[2,2] + pad*sim_box.L) .& 
            (box[3,1] - pad*sim_box.L .< coords_x[:,3] .< box[3,2] + pad*sim_box.L) 
    Ni_ = size(coords_q[idx_x, :])[1]

    ps_dtfe = PS_DTFE(coords_q[idx_x, :], coords_x[idx_x, :], velocities[idx_x, :], m, depth, box)
    i, j, k = idx
    serialize(dir * "/box_" * string(i) * "_" * string(j) * "_" * string(k), ps_dtfe)
end

"""
    get_coords_in_subbox(coords, idx, N_sub, L)

Obtain the positions in a given subbox.
"""
function get_coords_in_subbox(coords, idx, N_sub, L)
    box = 1. / N_sub * L .* hcat(idx, idx .+ 1)

    coord_indices  = CartesianIndices(coords)
    in_subbox      = CartesianIndex[]

    for coord_idx in coord_indices
        p = coords[coord_idx]

        ## if box has index such that it contains a "right-hand boundary" (x = L or y = L or z = L),
        ## need to consider for subbox range

        if all(idx .<  N_sub-1) ## not at "right-hand boundary"
            if (box[1,1] <= p[1] < box[1,2]) & (box[2,1] <= p[2] < box[2,2]) & (box[3,1] <= p[3] < box[3,2])
                push!(in_subbox, coord_idx)
            end

        else  ## at least one "right-hand boundary" -> distinguish cases

            if (idx[1] == N_sub-1) && !(idx[2] == N_sub-1) && !(idx[3] == N_sub-1)
                if (box[1,1] <= p[1] <= box[1,2]) & (box[2,1] <= p[2] < box[2,2]) & (box[3,1] <= p[3] < box[3,2])
                    push!(in_subbox, coord_idx)
                end

            elseif !(idx[1] == N_sub-1) && (idx[2] == N_sub-1) && !(idx[3] == N_sub-1)
                if (box[1,1] <= p[1] < box[1,2]) & (box[2,1] <= p[2] <= box[2,2]) & (box[3,1] <= p[3] < box[3,2])
                    push!(in_subbox, coord_idx)
                end

            elseif !(idx[1] == N_sub-1) && !(idx[2] == N_sub-1) && (idx[3] == N_sub-1)
                if (box[1,1] <= p[1] < box[1,2]) & (box[2,1] <= p[2] < box[2,2]) & (box[3,1] <= p[3] <= box[3,2])
                    push!(in_subbox, coord_idx)
                end

            elseif (idx[1] == N_sub-1) && (idx[2] == N_sub-1) && !(idx[3] == N_sub-1)
                if (box[1,1] <= p[1] <= box[1,2]) & (box[2,1] <= p[2] <= box[2,2]) & (box[3,1] <= p[3] < box[3,2])
                    push!(in_subbox, coord_idx)
                end

            elseif (idx[1] == N_sub-1) && !(idx[2] == N_sub-1) && (idx[3] == N_sub-1)
                if (box[1,1] <= p[1] <= box[1,2]) & (box[2,1] <= p[2] < box[2,2]) & (box[3,1] <= p[3] <= box[3,2])
                    push!(in_subbox, coord_idx)
                end

            elseif !(idx[1] == N_sub-1) && (idx[2] == N_sub-1) && (idx[3] == N_sub-1)
                if (box[1,1] <= p[1] < box[1,2]) & (box[2,1] <= p[2] <= box[2,2]) & (box[3,1] <= p[3] <= box[3,2])
                    push!(in_subbox, coord_idx)
                end

            elseif (idx[1] == N_sub-1) && (idx[2] == N_sub-1) && (idx[3] == N_sub-1)
                if (box[1,1] <= p[1] <= box[1,2]) & (box[2,1] <= p[2] <= box[2,2]) & (box[3,1] <= p[3] <= box[3,2])
                    push!(in_subbox, coord_idx)
                end
                
            else ## shouldn't appear if have considered all cases
                println("Warning: unexpected exception in subbox density calulation.")
            end
        end
    end
    return in_subbox
end


## field calculation from subbox estimators -----------------------------------
"""
    density_subbox(coords_arr, ps_dtfe_sb)

Reconstruct the Phase-Space DTFE density in the point p.
"""
function density_subbox(coords_arr, ps_dtfe_sb)
    N_sub = ps_dtfe_sb.N_sub
    dir   = ps_dtfe_sb.dir

    density_arr  = zeros(Float64, size(coords_arr)...)

    box_idx_list = [[i, j, k] for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1]

    for (i, j, k) in box_idx_list
        println("subbox i, j, k = " * string([i, j, k]))
        in_subbox = get_coords_in_subbox(coords_arr, [i, j, k], N_sub, ps_dtfe_sb.L)

        println("   -> computing " * string(length(in_subbox)) * " elements")
        if length(in_subbox) > 0  ## only load ps-dtfe if necessary
            subbox_ps_dtfe = deserialize(dir * "/box_" * string(i) * "_" * string(j) * "_" * string(k))
            @showprogress for idx in in_subbox
                density_arr[idx] = density(coords_arr[idx], subbox_ps_dtfe)
            end
        end
    end
    return density_arr
end

"""
    numberOfStreams_subbox(coords_arr, ps_dtfe_sb)

Evaluates the number of streams the point p is in.
"""
function numberOfStreams_subbox(coords_arr, ps_dtfe_sb)
    N_sub = ps_dtfe_sb.N_sub
    dir   = ps_dtfe_sb.dir

    nstreams_arr  = zeros(Float64, size(coords_arr)...)

    box_idx_list = [[i, j, k] for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1]

    for (i, j, k) in box_idx_list
        println("subbox i, j, k = " * string([i, j, k]))
        in_subbox = get_coords_in_subbox(coords_arr, [i, j, k], N_sub, ps_dtfe_sb.L)

        println("   -> computing " * string(length(in_subbox)) * " elements")
        if length(in_subbox) > 0  ## only load ps-dtfe if necessary
            subbox_ps_dtfe = deserialize(dir * "/box_" * string(i) * "_" * string(j) * "_" * string(k))
            @showprogress for idx in in_subbox
                nstreams_arr[idx] = numberOfStreams(coords_arr[idx], subbox_ps_dtfe)
            end
        end
    end
    return nstreams_arr
end

"""
    velocity_subbox(coords_arr, ps_dtfe_sb)

Reconstructs the Phase-Space DTFE velocities in an array of points.
"""
function velocity_subbox(coords_arr, ps_dtfe_sb)
    N_sub = ps_dtfe_sb.N_sub
    dir   = ps_dtfe_sb.dir

    velocity_arr  = zeros(Float64, size(coords_arr)..., 3)

    box_idx_list = [[i, j, k] for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1]

    for (i, j, k) in box_idx_list
        println("subbox i, j, k = " * string([i, j, k]))
        in_subbox = get_coords_in_subbox(coords_arr, [i, j, k], N_sub, ps_dtfe_sb.L)

        println("   -> computing " * string(length(in_subbox)) * " elements")
        if length(in_subbox) > 0  ## only load ps-dtfe if necessary
            subbox_ps_dtfe = deserialize(dir * "/box_" * string(i) * "_" * string(j) * "_" * string(k))
            @showprogress for idx in in_subbox
                velocity_arr[idx, :] = velocity(coords_arr[idx], subbox_ps_dtfe, true)[:]
            end
        end
    end
    return velocity_arr
end

"""
    velocitySum_subbox(coords_arr, ps_dtfe_sb)

Reconstructs the mass weighted sum of the Phase-Space DTFE velocities in an array of points.
"""
function velocitySum_subbox(coords_arr, ps_dtfe_sb)
    N_sub = ps_dtfe_sb.N_sub
    dir   = ps_dtfe_sb.dir

    velocity_arr  = zeros(Float64, size(coords_arr)..., 3)

    box_idx_list = [[i, j, k] for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1]

    for (i, j, k) in box_idx_list
        println("subbox i, j, k = " * string([i, j, k]))
        in_subbox = get_coords_in_subbox(coords_arr, [i, j, k], N_sub, ps_dtfe_sb.L)

        println("   -> computing " * string(length(in_subbox)) * " elements")
        if length(in_subbox) > 0  ## only load ps-dtfe if necessary
            subbox_ps_dtfe = deserialize(dir * "/box_" * string(i) * "_" * string(j) * "_" * string(k))
            @showprogress for idx in in_subbox
                velocity_arr[idx, :] = velocitySum(coords_arr[idx], subbox_ps_dtfe)
            end
        end
    end
    return velocity_arr
end

## modular subbox calculations: single subbox ---------------------------------
"""
    get_subboxes(ps_dtfe_sub::PS_DTFE_subbox)
"""
function get_subboxes(ps_dtfe_sub::PS_DTFE_subbox)
    [[i, j, k] for i in 0:ps_dtfe_sub.N_sub-1, j in 0:ps_dtfe_sub.N_sub-1, k in 0:ps_dtfe_sub.N_sub-1][:]
end

"""
    get_coords_chunk(coords, m, m_idx, random=false)
"""
function get_coords_chunk(coords, m, m_idx, random=false)
    indices = CartesianIndices(size(coords))
        
    flat_coords  = reshape(coords, :)
    flat_indices = reshape(collect(indices), :)
    println("flat_coords = ", flat_coords[1], " ", flat_coords[end])
    println("flat_indices = ", flat_indices[1], " ", flat_indices[end])

    if random  # permute coordinates for approx. equal computational resources across jobs
        rng = MersenneTwister(1)   # seed fixed for same permutation on each job
        perm = randperm(rng, length(flat_coords))

        flat_coords  = flat_coords[perm]
        flat_indices = flat_indices[perm]
    end

    println("flat_coords = ", flat_coords[1], " ", flat_coords[end])
    println("flat_indices = ", flat_indices[1], " ", flat_indices[end])

    coords_chunks = chunks(flat_coords; n=m)
    index_chunks  = chunks(flat_indices; n=m)

    collect(coords_chunks[m_idx]), collect(index_chunks[m_idx])
end