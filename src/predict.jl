
function assign(library, targetdata, targetind, nbcprob, popDict)

    # Split into worker threads
    NCHUNKS = nthreads()
    chunks = LocalAncestry.vecsplit(targetind, NCHUNKS)
    npopulations = length(keys(popDict))
    blocks = sort(UnitRange.(keys(library)))
    inddict = Dict{String,Int}(targetind .=> 1:length(targetind))
    default_probabilities = repeat([1], npopulations) ./ npopulations
    out = DataFrame(individual=Vector{String}(),
        haplotype=Vector{Int64}(),
        block=Vector{UnitRange}(),
        ancestry=Vector{String}())
    # Locks
    writelock = ReentrantLock()

    # Do work
    @threads for c in 1:NCHUNKS
        # Internal initialization
        probabilities = zeros(Float64, length(keys(popDict)), length(blocks))
        ancestry = zeros(Int8, length(blocks))
        internal_ancestries = zeros(Int8, length(blocks), PLOIDITY * length(chunks[c]))
        haplotyperow::Int = 0
        internal_haplotypecol::Int = 0

        for (i, ind) in enumerate(chunks[c])
            for h in 1:PLOIDITY
                ancestry .= 0
                probabilities .= 0.0
                haplotyperow = PLOIDITY * inddict[ind] - abs(h - 2)
                internal_haplotypecol = PLOIDITY * i - abs(h - 2)

                for (ib, b) in enumerate(blocks)
                    probabilities[:, ib] = get(library[b], targetdata[haplotyperow, b], default_probabilities)
                end

                for p in 1:npopulations
                    ancestry[probabilities[p, :].>=nbcprob] .= p
                end

                # Assign missing
                if any(ancestry .> 0)
                    assign_missing!(probabilities, ancestry)
                end

                ## Allocate to internal memory
                internal_ancestries[:, internal_haplotypecol] = ancestry

            end
        end
        # Allocate to common memory
        @lock writelock append!(out, anc_reduce(internal_ancestries, chunks[c], popDict, blocks))
    end
    return out
end

function anc_reduce(a::Matrix{Int8}, i::Vector{String}, popDict, blocks)
    p = string.(keys(popDict))
    t = LocalAncestry.vec_reduce.(eachcol(a))
    o = DataFrame(individual=LocalAncestry.breapeat(repeat(i, inner=2), length.(t)),
        haplotype=LocalAncestry.breapeat(repeat(collect(1:PLOIDITY), outer=length(i)), length.(t)),
        block=vcat(t...),
        ancestry="")
    for j in axes(o, 1)
        jcol = 2 * findfirst(o.individual[j] .== i) - abs(o.haplotype[j] - 2)
        jrow = first(o.block[j])
        o.ancestry[j] = a[jrow, jcol] > 0 ? p[a[jrow, jcol]] : ""
    end
    o.block = UnitRange.(first.(blocks[first.(o.block)]), last.(blocks[last.(o.block)]))
    return o
end

function vec_reduce(v::AbstractVector)

    o = Vector{UnitRange}()
    s::Union{Nothing,Int} = 1
    e::Union{Nothing,Int} = 1
    n = length(v)

    while ~isnothing(e)
        e = findnext(x -> x != v[s], v, s)
        isnothing(e) ? push!(o, s:n) : push!(o, s:(e-1))
        s = e
    end

    return o
end

function hmm!(ancestry, probabilities, s, e)

    workrange = (s-1):e
    workrange2 = workrange[2:(end-1)]
    looprange = (workrange2) .- first(workrange) .+ 1
    ib::Int = 0
    # Setup
    forward = copy(probabilities[[ancestry[first(workrange)], ancestry[last(workrange)]], workrange])
    backward = copy(forward)

    for (j, i) in enumerate(looprange)
        # Forward
        @views forward[:, i] = forward[:, i] .* (forward[:, i-1] .* (1 - HMM_STATECHANGE_PROB) + forward[2:-1:1, i-1] .* HMM_STATECHANGE_PROB)
        # Backward
        ib = reverse(looprange)[j]
        @views backward[:, ib] = backward[:, ib] .* (backward[:, ib+1] .* (1 - HMM_STATECHANGE_PROB) + backward[2:-1:1, ib+1] .* HMM_STATECHANGE_PROB)
    end

    # Viterby
    forward = forward .* backward

    ancestry[workrange2] = ancestry[[s - 1, e]][last.(findmax.(eachcol(forward)[2:(end-1)]))]
    return nothing
end

function assign_missing!(probabilities, ancestry)

    s::Union{Nothing,Int} = 1
    e::Union{Nothing,Int} = 1
    n::Int = length(ancestry)

    while (~isnothing(e) & ~isnothing(s)) && (s <= n)

        e = findnext(x -> x > 0, ancestry, s)
        if ~isnothing(e)
            probabilities[:, e] .= 0.0
            probabilities[ancestry[e], e] = 1.0
            if (s == 1) & (e > 0)
                ancestry[s:e] .= ancestry[e]
            elseif ancestry[s] == ancestry[e]
                ancestry[s:e] .= ancestry[e]
            else
                hmm!(ancestry, probabilities, s, e)
            end
            s = findnext(x -> x == 0, ancestry, e)
        else
            ancestry[(s):n] .= ancestry[s-1]
        end
    end
end
