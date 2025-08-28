
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


function get_haplotype_library(data::Matrix{Int8}, popDict::Dict{String,Vector{Int}}, threshold)

    # Initialize
    maxchunk = 0.02
    nloci = size(data, 2)
    lociperchunk = Int(ceil(nloci * maxchunk))
    nhaplotypes = Dict{String,Int}(key => length(value) for (key, value) in popDict)
    writelock = ReentrantLock()
    claimlock = ReentrantLock()
    chunks = [[(1+(i-1)*lociperchunk):min(nloci, i * lociperchunk) for i in 1:50]; repeat([0:0], inner = Threads.nthreads())]

    # Output container
    haploLib = Dict{UnitRange{Int},Dict{Vector{Int8},Vector{Float64}}}()
    haploLib.slots = zeros(Int8, Int(ceil(nloci / 10)))

    Threads.@threads for _ in 1:Threads.nthreads()
        chunk = 0:0
        @lock claimlock chunk = popfirst!(chunks)
        if chunk == 0:0
            break
        end

        tmpdict = Dict{UnitRange{Int},Dict{Vector{Int8},Vector{Float64}}}()
        bstart = first(chunk)
        bend = bstart

        IA = zeros(Float64, 2)
        wv = zeros(Int, sum(values(nhaplotypes)))
        countmat = zeros(Float64, sum(values(nhaplotypes)), length(keys(popDict)))

        while bend <= last(chunk)
            get_unique_haplotypes!(data, bend - bstart + 1, wv, bstart - 1)

            if bstart == bend
                IA[1] = 0.00000000001
            end

            IA[2] = bend == last(chunk) ? 0.0 : compute_IA(wv, popDict, countmat, nhaplotypes)
            
            if (IA[1] == IA[2]) | (log(IA[2] / IA[1]) < threshold) | (bend == (last(chunk)))
                get_unique_haplotypes!(data, max(1, bend - bstart), wv, bstart - 1)
                counts, dict = get_haplotype_frequencies(data, wv, popDict, bstart, max(bstart, bend - 1))
                get!(tmpdict, bstart:max(bstart, bend - 1), get_haplotype_dictionary(counts, dict))
                bstart = max(bend, bstart + 1)
                bend = bstart
            else
                    IA[1] = IA[2]
                bend = bend + 1
            end

        end

        @lock writelock merge!(haploLib, tmpdict)
    end
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

    p_bar_v = mean.(eachrow(countmat[1:n,:]))

    IA = sum(sum(countmat[1:n,:] .* log.(countmat[1:n,:])) ./ length(keys(p))) - sum(p_bar_v .* log.(p_bar_v))
    return IA
end

function get_haplotype_dictionary(countmat, tmpdict)

    countmat = countmat ./ sum.(eachrow(countmat))

    return Dict{Vector{Int8}, Vector{Float64}}(keys(tmpdict) .=> eachrow(countmat))
end

