
function computeIA(h, b, p, n)

    countmat, _ = haplotypeFrequencies(h,b,p,n)
    #this computes freq as mean of freqs
    p_bar_v = mean.(eachrow(countmat))

    IA = sum(sum(countmat .* log.(countmat)) ./ length(n)) - sum(p_bar_v .* log.(p_bar_v))
    return IA
end

function haplotypeFrequencies(h::Matrix{Int8}, b::UnitRange{Int}, p::Dict{String, Vector{Int}}, n::Vector{Int})
    # Get unique haplotypes
    V = unique(h[:, b]; dims=1)
    tmpdict = OrderedDict{Array{Int8},Int}(eachrow(V) .=> 1:size(V, 1))
    countmat = zeros(Float64, size(V, 1), length(n)) .+ 0.00000001

    # Calculate IA
    for (popi, pop) in enumerate(keys(p))
        for i in p[pop]
            countmat[tmpdict[h[i, b]], popi] += 1
        end
        countmat[:, popi] = countmat[:, popi] ./ n[popi]
    end
    
    return countmat, tmpdict
end 

function haplotypeDict(countmat, tmpdict)

    countmat = countmat ./ sum.(eachrow(countmat))

    return Dict{Vector{Int8}, Vector{Float64}}(keys(tmpdict) .=> eachrow(countmat))
end

