function assign(library, targetdata, targetsamples, popDict; totcm::Float64=100.0)

    # Split into worker threads
    targetind = unique(targetsamples.individual)
    NCHUNKS = nthreads()
    chunks = vecsplit(targetind, NCHUNKS)
    npopulations = length(keys(popDict))
    blocks = sort(UnitRange.(keys(library)))
    inddict = Dict{String,Int}(targetind .=> 1:length(targetind))
    default_probabilities = fill(1.0 / npopulations, npopulations)
    out = DataFrame(individual=Vector{String}(),
        haplotype=Vector{Int64}(),
        block=Vector{UnitRange}(),
        ancestry=Vector{String}())

    # Locks
    writelock = ReentrantLock()

    # Do work
    @threads for c in 1:NCHUNKS
        # Internal initialization
        percm = totcm/(size(targetdata, 2)-1)
        pcross = (1/2) * (1 - exp(-2 * percm/100))
        pdif = (npopulations-1)/npopulations
        r = pdif * pcross
        pstay = 1 - npopulations*r
        pswitch = r
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
                    probabilities[:, ib] = 0.9 .* get(library[b], targetdata[haplotyperow, b], default_probabilities) .+ (0.1 * 1/npopulations)
                end

                hmm_all!(ancestry, probabilities, pstay, pswitch)

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
    t = vec_reduce.(eachcol(a))
    o = DataFrame(individual=brepeat(repeat(i, inner=2), length.(t)),
        haplotype=brepeat(repeat(collect(1:PLOIDITY), outer=length(i)), length.(t)),
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

    while !isnothing(e)
        e = findnext(x -> x != v[s], v, s)
        isnothing(e) ? push!(o, s:n) : push!(o, s:(e-1))
        s = e
    end

    return o
end

function hmm_all!(ancestry, probabilities, pstay, pswitch)

    # Setup
    ib::Int = 0
    n = length(ancestry)
    forward = deepcopy(probabilities)
    backward = deepcopy(probabilities)

    for i in axes(forward, 2)[2:end]
        # Forward
        @views forward[:, i] = forward[:, i] .* (forward[:, i-1] .* pstay .+ (pswitch*sum(forward[:, i-1])))
        forward[:, i] = forward[:, i]/sum(forward[:, i])

        # Backward
        ib = n - i + 1
        @views backward[:, ib] = backward[:, ib] .* (backward[:, ib+1] .* pstay .+ (pswitch*sum(backward[:, ib+1])))
        backward[:, ib] = backward[:, ib]/sum(backward[:, ib])
    end

    # Posterior decoding (forward-backward)
    forward = forward .* backward

    ancestry[:] = last.(findmax.(eachcol(forward)))
    return nothing
end
