
function get_unique_haplotypes!(data::Matrix{Int8}, m::Int, v::Vector{Int}, offset::Int)
    for c in 1:m
        hapDict = Dict{Tuple{Int8,Int},Int}()
        if c == 1
            v .= data[:, c+offset] .+ 1
        end
        j::Int = 1
        for i in eachindex(v)
            if haskey(hapDict, (data[i, c+offset], v[i]))
                v[i] = hapDict[(data[i, c+offset], v[i])]
            else
                get!(hapDict, (data[i, c+offset], v[i]), j)
                v[i] = j
                j = j + 1
            end
        end
    end
end

function get_haplotype_frequencies(data, wv, p, bstart, bend)

    nunique = maximum(wv)
    tmpdict = Dict{Array{Int8},Int}(eachrow(data[[findfirst(wv .== i) for i in 1:nunique], bstart:bend]) .=> 1:nunique)
    countmat = zeros(Float64, nunique, length(keys(p))) .+ NEARZERO_FLOAT
    # Calculate IA
    for (popi, pop) in enumerate(keys(p))
        for h in 1:nunique
            countmat[h, popi] = sum(wv[p[pop]] .== h)
        end
        countmat[:, popi] = countmat[:, popi] ./ length(p[pop])
    end

    return countmat, tmpdict

end

function get_haplotype_library(refdata::Matrix{Int8}, popDict::Dict{String, UnitRange{Int}}, threshold::Float64, refloci::DataFrame)

    # Initialize
    ## Chunks
    NCHUNKS = nthreads()

    nloci::Int = size(refdata, 2)
    chunks::Vector{UnitRange} = LocalAncestry.vecsplit(1:nloci, NCHUNKS)
    blocks = Vector{UnitRange}()

    ## Locks
    writelock = ReentrantLock()
    println("   Block ranges")
    # Work
    @threads for i in 1:NCHUNKS

        internal_blocks = LocalAncestry.get_haplotype_blocks(refdata, chunks[i], popDict, threshold, refloci)
        @lock writelock push!(blocks, internal_blocks...)

    end

    # Create 
    println("   Haplotype frequencies")
    haploLib = Dict{UnitRange,Dict{Vector{Int8},Vector{Float64}}}()
    blockchunks = LocalAncestry.vecsplit(blocks, NCHUNKS)

    @threads for i in 1:NCHUNKS
        blocks =  blockchunks[i]
        internal_haploLib = get_haplolib_from_blocks(blocks, refdata, popDict)
        @lock writelock merge!(haploLib, internal_haploLib)
    end

    # Get 
    return haploLib
end

function compute_IA(wv::Vector{Int}, popDict, countmat, p_bar_v, n, npopulations; outtype="all")
    countmat[1:n, :] .= NEARZERO_FLOAT
    for (popi, pop) in enumerate(keys(popDict))
        for i in values(popDict[pop])
            countmat[wv[i], popi] += 1
        end
        countmat[1:n, popi] = countmat[1:n, popi] ./ length(popDict[pop])
    end
    p_bar_v[1:n] .= LocalAncestry.mean.(eachrow(countmat[1:n, :]))

    if outtype == "all"
        return IAall(countmat, p_bar_v, npopulations, n)
    elseif outtype == "min"
        IA = log(size(countmat, 2)) + 1
        for i in 1:npopulations, j in 1:npopulations
            if i < j
                IA = min(IAsome(countmat, [i,j], n), IA)
            end
        end
        return IA
    else
        error("IA specification not valid")
    end
end

function IAall(countmat, p_bar_v, npopulations, n)
    return sum(countmat[1:n, :] .* log.(countmat[1:n, :])) / npopulations - sum(p_bar_v[1:n] .* log.(p_bar_v[1:n]))
end

function IAsome(countmat, these, n)
    npopulations = length(these)
    p_bar_v = LocalAncestry.mean.(eachrow(countmat[1:n, these]))
    return sum(countmat[1:n, these] .* log.(countmat[1:n, these])) / npopulations - sum(p_bar_v .* log.(p_bar_v))
end

function get_haplotype_blocks(refdata, c, p, threshold, refloci)
    # Initialize
    ## Amounts
    nhaplotypes = size(refdata, 1)
    npopulations = length(keys(p))
    nloci = nrow(refloci)

    ## Chunk-related
    n = length(c)
    firsti = first(c)

    ## Containers
    v = zeros(Int, nhaplotypes)
    countmat = zeros(Float64, nhaplotypes, npopulations)
    p_bar_v = zeros(Float64, nhaplotypes)
    o = Vector{UnitRange}(undef, nloci)
    IA1::Float64 = 0.0
    IA2::Float64 = 0.0
    hapDict = Dict{Tuple{Int8,Int},Int}()
    sizehint!(hapDict, nhaplotypes)
    j::Int = 1
    thisi = firsti
    oj = 1

    maxn = Int(ceil(MAX_BLOCK_FRACTION*nloci))

    # Do work: 
    for l in 1:n
        j = 1
        # If first locus, insert haplotypes
        if thisi == (firsti + l - 1)
            v .= refdata[:, thisi] .+ 1
            IA1 = compute_IA(v, p, countmat, p_bar_v, 2, npopulations) 
        else
            for i in eachindex(v)
                if haskey(hapDict, (refdata[i, l+firsti-1], v[i]))
                    v[i] = hapDict[(refdata[i, l+firsti-1], v[i])]
                else
                    get!(hapDict, (refdata[i, l+firsti-1], v[i]), j)
                    v[i] = j
                    j = j + 1
                end
            end

            IA2 = compute_IA(v, p, countmat, p_bar_v, maximum(values(hapDict)), npopulations) 
            empty!(hapDict)

            if (IA2 <= threshold) & (length(thisi:(firsti+l-1)) < maxn)
                IA1 = IA2
            else
                o[oj] = thisi:(firsti+l-1)
                thisi = firsti + l
                oj = oj + 1
            end
        end
    end

    # If no block was assigned
    if oj == 1
        o[1] = firsti:(firsti+n-1)
        oj = 2
    end

    # If last locus wasnt assigned to a block
    if (last(o[oj-1]) < (n + firsti - 1))
        o[oj] = (last(o[oj-1])+1):(n+firsti-1)
        oj = oj + 1
    end

    return o[1:(oj-1)]
end


function get_haplolib_from_blocks(blocks, refdata, popDict)

    npopulations = length(keys(popDict))
    outlib = Dict{UnitRange,Dict{Vector{Int8},Vector{Float64}}}()
    sizehint!(outlib, length(blocks))

    # 
    for block in blocks
    tmplib = Dict{Vector{Int8},Vector{Float64}}()
        for (ip, p) in enumerate(keys(popDict))
            for r in popDict[p]

                t = @view refdata[r, block]
                if haskey(tmplib, t)
                    tmplib[t][ip] += 1
                else
                    get!(tmplib, t, zeros(Float64, npopulations))
                    tmplib[t][ip] += 1
                end
            end
        end

        # Standardize
        for k in keys(tmplib)
            tmplib[k] = tmplib[k] ./ sum(tmplib[k])
        end
                 outlib[block] = tmplib
    end

    return outlib
end

