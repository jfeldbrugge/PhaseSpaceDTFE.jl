using ProgressMeter, Serialization

###############################################################################
# particles of a sub-box
###############################################################################

"""
    subbox_range(a, N_sub, L, pad, h)

Range of the initial positions along one axis of the particles of sub-box `a` (0-based): the core [a/N_sub, (a+1)/N_sub) L widened by `pad` lattice spacings `h`, open towards the faces of the box for the first and last sub-box.
"""
subbox_range(a, N_sub, L, pad, h) = (a == 0 ? -Inf : a / N_sub * L - pad * h, a == N_sub - 1 ? Inf : (a + 1) / N_sub * L + pad * h)

"""
    subbox_particles(coords_q, coords_x, sim_box, (i, j, k); N_sub, pad = 5, periodicity = true)

The particles of sub-box (i, j, k), selected from the whole simulation one particle at a time, so that only the particles of the sub-box are copied. Returns their initial and final positions, as `(n, 3)` matrices, and the row of `coords_q` of every selected particle.

With `periodicity = true` the particles are first wrapped as for a periodic box: the final position is unwrapped to the image nearest to the initial position (displacement in [-L/2, L/2)), particles with a final position outside the box are translated by ±L (initial and final position alike), and particles within `pad` lattice spacings of a face are copied to the opposite side of the box, as a frame of periodic images. A particle may thus be selected several times, with different shifts. With `periodicity = false` the positions are taken as given. Particles are selected when their (shifted) initial position lies in the `subbox_range` of the sub-box along every axis.
"""
function subbox_particles(coords_q::AbstractMatrix, coords_x::AbstractMatrix, sim_box::SimBox, (i, j, k);
                          N_sub, pad = 5, periodicity = true)
    L, h, w = sim_box.L, sim_box.L / sim_box.Ni, pad / sim_box.Ni
    ranges = (subbox_range(i, N_sub, L, pad, h), subbox_range(j, N_sub, L, pad, h), subbox_range(k, N_sub, L, pad, h))
    parts = chunks(size(coords_q, 1))
    found = [(Int[], NTuple{3, Float64}[], NTuple{3, Float64}[]) for _ in parts]
    Threads.@threads for c in eachindex(parts)
        local rows, qs, xs = found[c]                   # local: the loop body is a closure over the function's variables
        cand = MMatrix{3, 4, Float64}(undef)            # per axis: up to two candidates (q, x), as columns (q₁, x₁, q₂, x₂)
        ncand = MVector{3, Int}(undef)
        @inbounds for r in parts[c]
            ok = true
            for d in 1:3
                local q, x = Float64(coords_q[r, d]), Float64(coords_x[r, d])
                if periodicity
                    x = unwrap_s(x - q, L) + q                          # final position nearest to the initial one
                    if x < 0                                            # translate into the box
                        q, x = q + L, x + L
                    elseif x > L
                        q, x = q - L, x - L
                    end
                end
                lo, hi = ranges[d]
                n = 0
                if lo < q < hi
                    n += 1; cand[d, 2n - 1], cand[d, 2n] = q, x
                end
                if periodicity                                          # frame: copy to the opposite side of the box
                    if x > (1. - w) * L
                        qc, xc = q - L, x - L
                    elseif x < w * L
                        qc, xc = q + L, x + L
                    else
                        qc = NaN
                    end
                    if lo < qc < hi
                        n += 1; cand[d, 2n - 1], cand[d, 2n] = qc, xc
                    end
                end
                ncand[d] = n
                n == 0 && (ok = false; break)
            end
            ok || continue
            for a in 1:ncand[1], b in 1:ncand[2], e in 1:ncand[3]
                push!(rows, r)
                push!(qs, (cand[1, 2a - 1], cand[2, 2b - 1], cand[3, 2e - 1]))
                push!(xs, (cand[1, 2a], cand[2, 2b], cand[3, 2e]))
            end
        end
    end
    sel_rows = reduce(vcat, (f[1] for f in found))
    sel_q = [t[d] for t in reduce(vcat, (f[2] for f in found)), d in 1:3]
    sel_x = [t[d] for t in reduce(vcat, (f[3] for f in found)), d in 1:3]
    return sel_q, sel_x, sel_rows
end



###############################################################################
# subbox estimators
###############################################################################

"""
    get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, (i, j, k); N_sub, pad = 5, periodicity = true)

The estimator of sub-box (i, j, k): the Delaunay tessellation of the particles of the sub-box (see `subbox_particles`), of which the tetrahedra with their Lagrangian centroid in the core of the sub-box are kept.
"""
function get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, (i, j, k); N_sub, pad = 5, periodicity = true)
    q, x, rows = subbox_particles(coords_q, coords_x, sim_box, (i, j, k); N_sub, pad, periodicity)
    v = isempty(velocities) ? velocities : velocities[rows, :]
    m_sb = m isa Number ? m : m[rows]

    # subbox boundaries
    L = sim_box.L
    lo_c = [i == 0 ? -Inf : i/N_sub*L, j == 0 ? -Inf : j/N_sub*L, k == 0 ? -Inf : k/N_sub*L]
    hi_c = [i == N_sub-1 ? Inf : (i+1)/N_sub*L, j == N_sub-1 ? Inf : (j+1)/N_sub*L, k == N_sub-1 ? Inf : (k+1)/N_sub*L]

    construct = periodicity ? PS_DTFE_wrapped : PS_DTFE_nonperiodic
    construct(q, x, v, m_sb, sim_box; core = (lo_c, hi_c), ids = rows)
end

"Whether the sub-box routines show progress bars: unless the environment variable `JULIA_PROGRESS` is `\"false\"` (as when building the documentation)."
show_progress() = get(ENV, "JULIA_PROGRESS", "true") != "false"

"The sub-boxes (i, j, k), in the order of the loops over i, j and k (k fastest)."
subbox_list(N_sub) = vec([(i, j, k) for k in 0:N_sub-1, j in 0:N_sub-1, i in 0:N_sub-1])

"""
    foreach_batch(f, make, boxes, parallel, prog)

For the sub-boxes `boxes`, in batches of `parallel`: make the estimators of a batch concurrently with `make(ijk)`, then
call `f(est, ijk)` on each of them in the order of `boxes`, and release them before the next batch. At most `parallel`
estimators are held in memory at a time; the results do not depend on `parallel`, since `f` is called in order.
"""
function foreach_batch(f, make, boxes, parallel::Integer, prog)
    parallel >= 1 || throw(ArgumentError("parallel must be at least 1, got $parallel"))
    n = 0
    for batch in Iterators.partition(boxes, parallel)
        update!(prog, n; showvalues = [(:subboxes, "$(n + 1)–$(n + length(batch)) / $(length(boxes))"), (:stage, "constructing estimators")])
        ests = Vector{Any}(undef, length(batch))
        if length(batch) == 1
            ests[1] = make(batch[1])
        else
            @sync for (b, ijk) in enumerate(batch)
                Threads.@spawn ests[b] = make(ijk)
            end
        end
        for (b, ijk) in enumerate(batch)
            n += 1
            update!(prog, n - 1; showvalues = [(:subbox, "$n / $(length(boxes)) $ijk"), (:stage, "evaluating")])
            f(ests[b], ijk)
            ests[b] = nothing
        end
        GC.gc()
    end
    finish!(prog)
end

struct PS_DTFE_subbox
    sim_box::SimBox
    N_sub::Int64
    dir::String
end

"""
    PS_DTFE_subbox(coords_q, coords_x, velocities, m, sim_box; N_sub = 2, pad = 5, periodicity = true, dir = "./ps_dtfe", parallel = 1)
    PS_DTFE_subbox(coords_q, coords_x, m, sim_box; ...)

Construct the N_sub³ sub-box estimators and save them in `dir` (replacing the estimators of an earlier run, and no other
files), for repeated use with `density_grid_sb`, `numberOfStreams_grid_sb` and `velocity_grid_sb`. With `parallel = k`,
k estimators are constructed concurrently, which needs k times the memory of one.
"""
function PS_DTFE_subbox(coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                 N_sub=2, pad=5, periodicity=true, dir="./ps_dtfe", parallel::Integer = 1)

    # remove only the estimators of an earlier run, never other files in `dir`
    mkpath(dir)
    for f in readdir(dir)
        startswith(f, "est_") && rm(joinpath(dir, f))
    end

    prog = Progress(N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())
    make(ijk) = get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, ijk; N_sub, pad, periodicity)
    foreach_batch(make, subbox_list(N_sub), parallel, prog) do est, (i, j, k)
        serialize(joinpath(dir, "est_$(i)_$(j)_$(k)"), est)
    end

    return PS_DTFE_subbox(sim_box, N_sub, dir)
end

PS_DTFE_subbox(coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox; kwargs...) =
    PS_DTFE_subbox(coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)



#= Commented out: sub-box routines that loop over the points, evaluating every point in every sub-box.
   The grid routines below (density_grid_sb, numberOfStreams_grid_sb, velocity_grid_sb) loop over the
   tetrahedra instead and give the same values at the grid points.

###############################################################################
# on-the-fly methods
###############################################################################

function density_sb(points::AbstractArray{<:AbstractVector}, coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                 N_sub=2, pad=5, periodicity=true)


    field = zeros(Float64, size(points))
    prog  = Progress(N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())

    n = 0
    for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1
        n += 1
        status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(N_sub^3) ($i, $j, $k)"), (:stage, s)])

        status("constructing estimator")
        est     = get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, (i, j, k); N_sub, pad, periodicity)
        
        status("evaluating density")
        field .+= density.(points, Ref(est))

        est = nothing
        GC.gc()
        #@info "block $n" live_GB = Base.gc_live_bytes()/1e9 rss_GB = Sys.maxrss()/1e9

    end
    finish!(prog)
    return field
end

function density_sb(points::AbstractArray{<:AbstractVector}, coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox;
                     kwargs...)
    density_sb(points, coords_q, coords_x, zeros(0,3), m, sim_box; kwargs...)
end

function numberOfStreams_sb(points::AbstractArray{<:AbstractVector}, coords_q::AbstractMatrix, coords_x::AbstractMatrix, velocities::AbstractMatrix, m, sim_box::SimBox;
                 N_sub=2, pad=5, periodicity=true)


    field = zeros(Float64, size(points))
    prog  = Progress(N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())

    n = 0
    for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1
        n += 1
        status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(N_sub^3) ($i, $j, $k)"), (:stage, s)])

        status("constructing estimator")
        est     = get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, (i, j, k); N_sub, pad, periodicity)
        
        status("evaluating number of streams")
        field .+= numberOfStreams.(points, Ref(est))

        est = nothing
        GC.gc()
        #@info "block $n" live_GB = Base.gc_live_bytes()/1e9 rss_GB = Sys.maxrss()/1e9

    end
    finish!(prog)
    return field
end

function numberOfStreams_sb(points::AbstractArray{<:AbstractVector}, coords_q::AbstractMatrix, coords_x::AbstractMatrix, m, sim_box::SimBox;
                     kwargs...)
    numberOfStreams_sb(points, coords_q, coords_x, zeros(0,3), m, sim_box; kwargs...)
end

function velocity_sb(points::AbstractArray{<:AbstractVector}, coords_q::AbstractMatrix, coords_x::AbstractMatrix,
                  velocities::AbstractMatrix, m, sim_box::SimBox; N_sub = 2, pad = 5, periodicity=true)
    

    field = [zeros(0, 3) for _ in points]
    prog  = Progress(N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())

    n = 0
    for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1
        n += 1
        status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(N_sub^3) ($i, $j, $k)"), (:stage, s)])

        status("constructing estimator")
        est = get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, (i, j, k); N_sub, pad, periodicity)
        
        status("evaluating velocity")
        for idx in eachindex(points)
            vs = velocity(points[idx], est)          # (streams in this block, 3), maybe (0, 3)
            isempty(vs) || (field[idx] = vcat(field[idx], vs))
        end

        est = nothing
        GC.gc()
        #@info "block $n" live_GB = Base.gc_live_bytes()/1e9 rss_GB = Sys.maxrss()/1e9
    end
    finish!(prog)
    return field
end

function velocitySum_parts(p::AbstractVector{<:Real}, est::PS_DTFE)
    isempty(est.velocities) && throw(ArgumentError("the estimator was constructed without velocities"))
    x, y, z = query_point(est, p)
    V = est.velocities
    ρsum = 0.
    num  = zeros(3)
    for t in cell_tetrahedra(est.index, x, y, z)
        λ1, λ2, λ3, in = barycentric(est, t, x, y, z)
        in || continue
        ρ = interpolate(est.rho, est, t, λ1, λ2, λ3)
        ρsum += ρ
        for d in 1:3
            num[d] += ρ * interpolate(@view(V[:, d]), est, t, λ1, λ2, λ3)
        end
    end
    return ρsum, num
end

function velocitySum_sb(points::AbstractArray{<:AbstractVector}, coords_q::AbstractMatrix, coords_x::AbstractMatrix,
                        velocities::AbstractMatrix, m, sim_box::SimBox; N_sub = 2, pad = 5, periodicity=true)
    
    
    ρsum = zeros(Float64, size(points))               # denominator: Σρ over all streams and blocks
    vel  = zeros(Float64, size(points)..., 3)         # numerator:   Σρ·v

    prog  = Progress(N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())
    n = 0
    for i in 0:N_sub-1, j in 0:N_sub-1, k in 0:N_sub-1
        n += 1
        status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(N_sub^3) ($i, $j, $k)"), (:stage, s)])

        status("constructing estimator")
        est = get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, (i, j, k); N_sub, pad, periodicity)
        
        status("evaluating velocity")
        for idx in CartesianIndices(points)
            ρ, vp = velocitySum_parts(points[idx], est)
            ρ == 0 && continue                        # no stream of this block at the point
            ρsum[idx] += ρ
            vel[idx, :] .+= vp
        end

        est = nothing                                 # release the estimator before building the next one
        GC.gc()
    end
    finish!(prog)

    for n in CartesianIndices(points)
        vel[n, :] ./= ρsum[n]                         # 0/0 = NaN where no block has a stream
    end
    return vel                                        # (size(points)..., 3)
end
=#



# deformation of velocity ----------------------------------------------------------

"""
    foreach_subbox(f, coords_q, coords_x, velocities, m, sim_box; N_sub = 2, pad = 5, periodicity = true, parallel = 1)
    foreach_subbox(f, ps_dtfe_sb::PS_DTFE_subbox; parallel = 1)

Call `f(est)` on the estimator of every sub-box, in a fixed order, constructing it (or loading it from disk) and
releasing it again, so that at most `parallel` estimators are held in memory at a time. With `parallel = k`, k
estimators are constructed (or loaded) concurrently.
"""
function foreach_subbox(f, coords_q, coords_x, velocities, m, sim_box::SimBox;
                        N_sub = 2, pad = 5, periodicity = true, parallel::Integer = 1)
    prog = Progress(N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())
    make(ijk) = get_subbox_estimator(coords_q, coords_x, velocities, m, sim_box, ijk; N_sub, pad, periodicity)
    foreach_batch((est, _) -> f(est), make, subbox_list(N_sub), parallel, prog)
end

function foreach_subbox(f, sb::PS_DTFE_subbox; parallel::Integer = 1)
    prog = Progress(sb.N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())
    make((i, j, k)) = deserialize(joinpath(sb.dir, "est_$(i)_$(j)_$(k)"))
    foreach_batch((est, _) -> f(est), make, subbox_list(sb.N_sub), parallel, prog)
end

#= Commented out: sub-box routines that loop over the points, evaluating every point in every sub-box.
   The grid routines below (density_grid_sb, numberOfStreams_grid_sb, velocity_grid_sb) loop over the
   tetrahedra instead and give the same values at the grid points.

# ---- generic per-stream field: concatenates the streams of all blocks ----
# eval(p, est) must return a Vector (one entry per stream), e.g. velocityDivergence
function perstream_sb(eval, points, args...; kwargs...)
    field = [Any[] for _ in points]
    foreach_subbox(args...; kwargs...) do est
        for idx in eachindex(points)
            append!(field[idx], eval(points[idx], est))
        end
    end
    return field
end

# ---- generic density-weighted field: Σρ·f / Σρ over streams and blocks ----
# parts(p, est) must return (ρsum, num), e.g. shear_parts
function weighted_sb(parts, zero_val, points, args...; kwargs...)
    ρsum = zeros(Float64, size(points))
    num  = fill(zero_val, size(points))              # same shape as points, entries of type typeof(zero_val)
    foreach_subbox(args...; kwargs...) do est
        for idx in eachindex(points)
            ρ, nm = parts(points[idx], est)
            ρ == 0 && continue
            ρsum[idx] += ρ
            num[idx]   = num[idx] + nm
        end
    end
    return num ./ ρsum                               # NaN where no block has a stream
end

divergence_sb(points, cq, cx, v, m, sb::SimBox; kw...) = perstream_sb(divergence, points, cq, cx, v, m, sb; kw...)
shear_sb(points, cq, cx, v, m, sb::SimBox; kw...)      = perstream_sb(shear,      points, cq, cx, v, m, sb; kw...)
#velocityShearSq_sb(points, cq, cx, v, m, sb::SimBox; kw...)    = perstream_sb(shearSq,    points, cq, cx, v, m, sb; kw...)
vorticity_sb(points, cq, cx, v, m, sb::SimBox; kw...)  = perstream_sb(vorticity,  points, cq, cx, v, m, sb; kw...)

# density-weighted
divergenceSum_sb(points, cq, cx, v, m, sb::SimBox; kw...) = weighted_sb(divergenceSum_parts, 0.0,                         points, cq, cx, v, m, sb; kw...)
shearSum_sb(points, cq, cx, v, m, sb::SimBox; kw...)      = weighted_sb(shear_parts,         zero(SMatrix{3,3,Float64}),  points, cq, cx, v, m, sb; kw...)
#shearSqSum_sb(points, cq, cx, v, m, sb::SimBox; kw...)    = weighted_sb(shearSq_parts,       0.0,                         points, cq, cx, v, m, sb; kw...)
vorticitySum_sb(points, cq, cx, v, m, sb::SimBox; kw...)  = weighted_sb(vorticity_parts,     zero(SVector{3,Float64}),    points, cq, cx, v, m, sb; kw...)











###############################################################################
# from-estimator methods
###############################################################################

function density_sb(points::AbstractArray{<:AbstractVector}, ps_dtfe_sb::PS_DTFE_subbox)

    field = zeros(Float64, size(points))
    prog  = Progress(ps_dtfe_sb.N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())

    n = 0
    for i in 0:ps_dtfe_sb.N_sub-1, j in 0:ps_dtfe_sb.N_sub-1, k in 0:ps_dtfe_sb.N_sub-1
        n += 1
        update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3)  ($i, $j, $k)")])
        #status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3) ($i, $j, $k)"), (:stage, s)])

        #status("get estimator")
        est     = deserialize(ps_dtfe_sb.dir * "/est_$(i)_$(j)_$(k)")
        
        #status("evaluating density")
        field .+= density.(points, Ref(est))

        est = nothing
        GC.gc()
        Sys.islinux() && ccall(:malloc_trim, Cint, (Cint,), 0)
    end
    finish!(prog)
    return field
end

function numberOfStreams_sb(points::AbstractArray{<:AbstractVector}, ps_dtfe_sb::PS_DTFE_subbox)

    field = zeros(Float64, size(points))
    prog  = Progress(ps_dtfe_sb.N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())

    n = 0
    for i in 0:ps_dtfe_sb.N_sub-1, j in 0:ps_dtfe_sb.N_sub-1, k in 0:ps_dtfe_sb.N_sub-1
        n += 1
        update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3)  ($i, $j, $k)")])
        #status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3) ($i, $j, $k)"), (:stage, s)])

        #status("get estimator")
        est     = deserialize(ps_dtfe_sb.dir * "/est_$(i)_$(j)_$(k)")
        
        #status("evaluating density")
        field .+= numberOfStreams.(points, Ref(est))

        est = nothing
        GC.gc()
        Sys.islinux() && ccall(:malloc_trim, Cint, (Cint,), 0)
    end
    finish!(prog)
    return field
end

function velocity_sb(points::AbstractArray{<:AbstractVector}, ps_dtfe_sb::PS_DTFE_subbox)
    field = [zeros(0, 3) for _ in points] 
    prog  = Progress(ps_dtfe_sb.N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())

    n = 0
    for i in 0:ps_dtfe_sb.N_sub-1, j in 0:ps_dtfe_sb.N_sub-1, k in 0:ps_dtfe_sb.N_sub-1
        n += 1
        update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3)  ($i, $j, $k)")])
        #status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3) ($i, $j, $k)"), (:stage, s)])

        #status("get estimator")
        est     = deserialize(ps_dtfe_sb.dir * "/est_$(i)_$(j)_$(k)")
        
        #status("evaluating density")
        for p_idx in eachindex(points)
            v = velocity(points[p_idx], est)
            isempty(v) || (field[p_idx] = vcat(field[p_idx], v))
        end

        est = nothing
        GC.gc()
    end
    finish!(prog)
    return field
end

function velocitySum_sb(points::AbstractArray{<:AbstractVector}, ps_dtfe_sb::PS_DTFE_subbox)
    ρsum = zeros(Float64, size(points))               # denominator: Σρ over all streams and blocks
    vel  = zeros(Float64, size(points)..., 3)         # numerator:   Σρ·v

    prog  = Progress(ps_dtfe_sb.N_sub^3; desc = "subboxes: ", showspeed = true, enabled = show_progress())

    n = 0
    for i in 0:ps_dtfe_sb.N_sub-1, j in 0:ps_dtfe_sb.N_sub-1, k in 0:ps_dtfe_sb.N_sub-1
        n += 1
        update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3)  ($i, $j, $k)")])
        #status(s) = update!(prog, n - 1; showvalues = [(:subbox, "$n / $(ps_dtfe_sb.N_sub^3) ($i, $j, $k)"), (:stage, s)])

        #status("get estimator")
        est     = deserialize(ps_dtfe_sb.dir * "/est_$(i)_$(j)_$(k)")
        
        #status("evaluating density")
        for n in CartesianIndices(points)
            ρ, v = velocitySum_parts(points[n], est)
            ρ == 0 && continue                        # no stream of this block at the point
            ρsum[n] += ρ
            vel[n, :] .+= v
        end

        est = nothing                                 # release the estimator before loading the next one
        GC.gc()
    end
    finish!(prog)

    for n in CartesianIndices(points)
        vel[n, :] ./= ρsum[n]                         # 0/0 = NaN where no block has a stream
    end
    return vel                                        # (size(points)..., 3)
end


#velocity_sb(points, sb::PS_DTFE_subbox)          = perstream_sb(velocityStreams,   points, sb)
divergence_sb(points, sb::PS_DTFE_subbox) = perstream_sb(divergence, points, sb)
shear_sb(points, sb::PS_DTFE_subbox)     = perstream_sb(shear,     points, sb)
#shearSq_sb(points, sb::PS_DTFE_subbox)   = perstream_sb(shearSq,   points, sb)
vorticity_sb(points, sb::PS_DTFE_subbox) = perstream_sb(vorticity, points, sb)

# ---------- density-weighted: Σρ·f / Σρ over streams and blocks ----------
#velocitySum_sb(points, sb::PS_DTFE_subbox)   = weighted_sb(velocitySum_parts_s, zero(SVector{3, Float64}),   points, sb)
divergenceSum_sb(points, sb::PS_DTFE_subbox) = weighted_sb(divergenceSum_parts, 0.0,                          points, sb)
shearSum_sb(points, sb::PS_DTFE_subbox)      = weighted_sb(shear_parts,         zero(SMatrix{3, 3, Float64}), points, sb)
#shearSqSum_sb(points, sb::PS_DTFE_subbox)    = weighted_sb(shearSq_parts,       0.0,                          points, sb)
vorticitySum_sb(points, sb::PS_DTFE_subbox)  = weighted_sb(vorticity_parts,     zero(SVector{3, Float64}),    points, sb)
=#

###############################################################################
# fields on regular grids
###############################################################################

"""
    rasterise_sb(xs, ys, zs, coords_q, coords_x, m, sim_box; N_sub = 2, pad = 5, periodicity = true, threads = true,
                 density_type = Float64, streams_type = Int, parallel = 1)
    rasterise_sb(xs, ys, zs, ps_dtfe_sb::PS_DTFE_subbox; threads = true, density_type = Float64, streams_type = Int, parallel = 1)

Density and number of streams on the regular grid xs × ys × zs, as arrays indexed `[ix, iy, iz]`, from the sub-box estimators: the tetrahedra of every sub-box are rasterised into the same grid (see `rasterise!`), holding one estimator in memory at a time. Equal at every grid point to the sum over the sub-boxes of the point queries `density` and `numberOfStreams` of the sub-box estimators.

The grid takes `sizeof(density_type) + sizeof(streams_type)` bytes per point; `density_type = Float32` and `streams_type = UInt16` (at most 65535 streams) reduce it from 16 to 6 bytes. Pass `nothing` for an array that is not needed. With `parallel = k`, k sub-box estimators are constructed (or loaded) concurrently, which needs k times the memory of one; the result does not depend on k. The same keywords apply to `density_grid_sb`, `numberOfStreams_grid_sb` and `velocity_grid_sb`.
"""
function rasterise_sb(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, args...; threads = true,
                      density_type = Float64, streams_type = Int, kwargs...)
    dims = (length(xs), length(ys), length(zs))
    ρ  = density_type === nothing ? nothing : zeros(density_type, dims)
    ns = streams_type === nothing ? nothing : zeros(streams_type, dims)
    foreach_subbox(args...; kwargs...) do est
        rasterise!(ρ, ns, est, xs, ys, zs; threads)
    end
    return ρ, ns
end

rasterise_sb(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, coords_q::AbstractMatrix, coords_x::AbstractMatrix,
             m, sim_box::SimBox; kwargs...) =
    rasterise_sb(xs, ys, zs, coords_q, coords_x, zeros(0, 3), m, sim_box; kwargs...)

"""
    density_grid_sb(xs, ys, zs, coords_q, coords_x, m, sim_box; N_sub = 2, pad = 5, periodicity = true, threads = true, density_type = Float64)
    density_grid_sb(xs, ys, zs, ps_dtfe_sb::PS_DTFE_subbox; threads = true, density_type = Float64)

Density on the regular grid xs × ys × zs, as an array indexed `[ix, iy, iz]`; equal at every grid point to the sum of `density` over the sub-box estimators, but computed by rasterising the tetrahedra of every sub-box. Only the density grid is allocated.
"""
density_grid_sb(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, args...; kwargs...) =
    rasterise_sb(xs, ys, zs, args...; streams_type = nothing, kwargs...)[1]

"""
    numberOfStreams_grid_sb(xs, ys, zs, coords_q, coords_x, m, sim_box; N_sub = 2, pad = 5, periodicity = true, threads = true, streams_type = Int)
    numberOfStreams_grid_sb(xs, ys, zs, ps_dtfe_sb::PS_DTFE_subbox; threads = true, streams_type = Int)

Number of streams on the regular grid xs × ys × zs, as an array indexed `[ix, iy, iz]`; equal at every grid point to the sum of `numberOfStreams` over the sub-box estimators. Only the stream grid is allocated.
"""
numberOfStreams_grid_sb(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, args...; kwargs...) =
    rasterise_sb(xs, ys, zs, args...; density_type = nothing, kwargs...)[2]

"""
    velocity_grid_sb(xs, ys, zs, coords_q, coords_x, velocities, m, sim_box; N_sub = 2, pad = 5, periodicity = true,
                     threads = true, velocity_type = Float64, deformation = false)
    velocity_grid_sb(xs, ys, zs, ps_dtfe_sb::PS_DTFE_subbox; threads = true, velocity_type = Float64, deformation = false)

Velocities of all streams on the regular grid xs × ys × zs from the sub-box estimators, as a `StreamVelocityGrid` `G`: `G[ix, iy, iz]` holds the rows of `velocity` of every sub-box estimator at the grid point, up to the order of the streams. With `deformation = true` the divergence, shear and vorticity of every stream are stored too, as in `velocity_grid`. The sub-boxes are rasterised twice, holding one estimator at a time: first to count the streams of every grid point, then to fill in their velocities. On the fly, every sub-box estimator is therefore constructed twice; with the estimators saved by `PS_DTFE_subbox`, they are only read twice.
"""
function velocity_grid_sb(xs::AbstractRange, ys::AbstractRange, zs::AbstractRange, args...; threads = true,
                          velocity_type::Type{<:AbstractFloat} = Float64, deformation::Bool = false, kwargs...)
    ns = numberOfStreams_grid_sb(xs, ys, zs, args...; threads, streams_type = Int32, kwargs...)
    G = StreamVelocityGrid(xs, ys, zs, ns, velocity_type; deformation)
    foreach_subbox(args...; kwargs...) do est
        fill_streams!(G, est; threads)
    end
    return check_filled(G)
end
