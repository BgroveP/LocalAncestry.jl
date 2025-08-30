
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
    countmat = zeros(Float64, nunique, length(keys(p))) .+ 0.000000001
    # Calculate IA
    for (popi, pop) in enumerate(keys(p))
        for h in 1:nunique
            countmat[h, popi] = sum(wv[p[pop]] .== h)
        end
        countmat[:, popi] = countmat[:, popi] ./ length(p[pop])
    end

    return countmat, tmpdict

end

function get_haplolib_from_blocks(block, refdata, popDict)

    refview = @view refdata[:, block]
    V = unique(refview, dims=1)
    nunique = size(V, 1)
    tmpdict = Dict{Array{Int8},Int}(eachrow(V) .=> 1:nunique)
    countmat = zeros(Float64, nunique, length(keys(popDict))) .+ 0.000000001

    # Calculate IA
    for (popi, pop) in enumerate(keys(popDict))
        for h in 1:nunique
            countmat[h, popi] = sum(eachrow(refview[popDict[pop], :]) .== [V[h, :]])
        end
        countmat[:, popi] = countmat[:, popi] ./ length(popDict[pop])
    end

    return Dict{Vector{Int8},Vector{Float64}}(keys(tmpdict) .=> eachrow(countmat) ./ sum.(eachrow(countmat)))
end

function get_haplotype_library(refdata::Matrix{Int8}, popDict::Dict{String,Vector{Int}}, threshold::Float64)

    # Initialize
    ## Chunks
    nloci::Int = size(refdata, 2)
    chunks::Vector{UnitRange} = vecsplit(1:nloci, NCHUNKS)
    nhaplotypes::Int = size(refdata, 1)
    nhaplotypesperblock::Vector{Int} = length.(values(popDict))
    npopulations::Int = length(keys(popDict))
    blocks = Vector{UnitRange}()

    ## Locks
    writelock = ReentrantLock()
    readlock = ReentrantLock()

    # Work
    @threads for i in 1:NCHUNKS

        # Internal initialize
        wv = zeros(Int, nhaplotypes)
        countmat = zeros(Float64, nhaplotypes, npopulations)
        p_bar_v = zeros(Float64, nhaplotypes)

        # Do work until there are no chunks left
        chunk::UnitRange = chunks[i]

        internal_blocks = get_haplotype_blocks(refdata, length(chunk), first(chunk), wv, popDict, countmat, nhaplotypesperblock, p_bar_v, npopulations, threshold)

        @lock writelock push!(blocks, internal_blocks...)

    end

    # Create 
    haploLib = Dict{UnitRange,Dict{Vector{Int8},Vector{Float64}}}()
    blockchunks = vecsplit(blocks, NCHUNKS)

    @threads for i in 1:NCHUNKS
        blockchunk::Vector{UnitRange} = [0:0]
        blockchunk = blockchunks[i]

        internal_haploLib = Dict{UnitRange,Dict{Vector{Int8},Vector{Float64}}}(block => get_haplolib_from_blocks(block, refdata, popDict) for block in blockchunk)
        @lock writelock merge!(haploLib, internal_haploLib)


    end

    # Get 
    return haploLib
end

function compute_IA(wv::Vector{Int}, p, countmat, nhaplotypes)

    # Get maximum number
    n = maximum(wv)
    for (popi, pop) in enumerate(keys(p))
        for i in 1:n
            countmat[i, popi] = sum(wv[p[pop]] .== i) .+ 0.00001
        end
        countmat[1:n, popi] = countmat[1:n, popi] ./ nhaplotypes[pop]
    end

    p_bar_v = mean.(eachrow(countmat[1:n, :]))

    IA = sum(sum(countmat[1:n, :] .* log.(countmat[1:n, :])) ./ length(keys(p))) - sum(p_bar_v .* log.(p_bar_v))
    return IA
end

function compute_IA2(wv::Vector{Int}, popDict, countmat, nhaplotypesperblock, p_bar_v, npopulations, n)
    countmat[1:n, :] .= NEARZERO_FLOAT
    for (popi, pop) in enumerate(keys(popDict))
        for i in values(popDict[pop])
            countmat[wv[i], popi] += 1
        end
        countmat[1:n, popi] = countmat[1:n, popi] ./ nhaplotypesperblock[popi]
    end
    p_bar_v[1:n] .= mean.(eachrow(countmat[1:n, :]))
    return sum(countmat[1:n, :] .* log.(countmat[1:n, :])) / npopulations - sum(p_bar_v[1:n] .* log.(p_bar_v[1:n]))
end


function get_haplotype_dictionary(countmat, tmpdict)

    countmat = countmat ./ sum.(eachrow(countmat))

    return Dict{Vector{Int8},Vector{Float64}}(keys(tmpdict) .=> eachrow(countmat))
end

function get_haplotype_blocks(refdata::Matrix{Int8}, n::Int, firsti::Int, v::Vector{Int}, p, countmat, nhaplotypesperblock, p_bar_v, npopulations, threshold)
    # Initialize
    o = Vector{UnitRange}(undef, n)
    IA1::Float64 = 0.0
    IA2::Float64 = 0.0
    hapDict = Dict{Tuple{Int8,Int},Int}()
    sizehint!(hapDict, n)
    j::Int = 1
    thisi = firsti
    oj = 1

    # Do work: 
    for l in 1:n
        j = 1
        # If first locus, insert haplotypes
        if thisi == (firsti + l - 1)

            v .= refdata[:, thisi] .+ 1
            IA1 = compute_IA2(v, p, countmat, nhaplotypesperblock, p_bar_v, npopulations, 2)

            if IA1 < NEARZERO_FLOAT
                thisi = firsti + l
            end
        else
            for i in eachindex(v)
                if haskey(hapDict, (refdata[i, l+firsti], v[i]))
                    v[i] = hapDict[(refdata[i, l+firsti], v[i])]
                else
                    get!(hapDict, (refdata[i, l+firsti], v[i]), j)
                    v[i] = j
                    j = j + 1
                end
            end
            IA2 = compute_IA2(v, p, countmat, nhaplotypesperblock, p_bar_v, npopulations, maximum(values(hapDict)))
            empty!(hapDict)
            if (log(IA2 / IA1) > threshold)
                IA1 = IA2
            else
                o[oj] = thisi:(firsti+l-1)
                thisi = firsti + l
                oj = oj + 1
            end
        end
    end

    if last(o[oj-1]) < (n + firsti - 1)
        o[oj] = (last(o[oj-1])+1):(n+firsti-1)
        oj = oj + 1
    end

    return o[1:(oj-1)]
end