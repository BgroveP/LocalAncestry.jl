
#gets unique haplotypes from the data of breeds for a given block
#n or N refers to haplotype number
function computeIA(h, p)

    # Get unique haplotypes
    V = unique(h; dims=1)

    # Number of haplotypes per population
    n = zeros(Int64, length(p))
    for (i, k) in enumerate(keys(p))
        n[i] = length(p[k])
    end

    # Calculate IA
    IA = 0
    for v in eachrow(V)
        count_v = []
        p_v = [] #within pop freq
        for (i, k) in enumerate(keys(p))
            count_b_v = 0
            for r in p[k]
                if all(h[r, :] .== v)
                    count_b_v += 1
                else
                    nothing # Why nothing?
                end
            end
            push!(count_v, count_b_v)
            push!(p_v, count_b_v / n[i])
        end
        #this computes freq as mean of freqs
        p_bar_v = mean(p_v)

        p_v_n0 = p_v[p_v.!=0]

        IA += (-1.0 * p_bar_v * log(p_bar_v) + sum(p_v_n0 .* log.(p_v_n0)) / length(n))
    end
    return IA, V
end

#gets unique haplotypes from the data of breeds for a given block
#n or N refers to haplotype number
function computeIA2(h, p, n)

    # Get unique haplotypes
    V = unique(h; dims=1)

    # Calculate IA
    IA = 0.0
    p_v = zeros(Float64, length(p)) #within pop freq
    count_b_v::Int = 0
    p_bar_v = 0.0
    for v in axes(V, 1)
        for (i, k) in enumerate(keys(p))
            count_b_v = 0
            for r in p[k]
                if all(h[r, :] .== V[v, :])
                    count_b_v += 1
                end
            end
            p_v[i] = count_b_v / n[i]
        end
        #this computes freq as mean of freqs
        p_bar_v = mean(p_v)

        p_v_n0 = p_v[p_v.!=0]

        IA += (-1.0 * p_bar_v * log(p_bar_v) + sum(p_v_n0 .* log.(p_v_n0)) / length(n))
    end
    return IA, V
end


#gets unique haplotypes from the data of breeds for a given block
#n or N refers to haplotype number
function computeIA3(h, b, p, n)

    # Get unique haplotypes
    V = unique(h[:, b]; dims=1)
    tmpdict = Dict{Array{Int8},Int}(eachrow(V) .=> 1:size(V, 1))
    countmat = zeros(Float64, size(V, 1), length(n)) .+ 0.00000001

    # Calculate IA
    for (popi, pop) in enumerate(keys(p))
        for i in p[pop]
            countmat[tmpdict[h[i, b]], popi] += 1
        end
        countmat[:, popi] = countmat[:, popi] ./ n[popi]
    end
    #this computes freq as mean of freqs
    p_bar_v = mean.(eachrow(countmat))

    IA = sum(sum(countmat .* log.(countmat)) ./ length(n)) - sum(p_bar_v .* log.(p_bar_v))
    return IA, V, countmat
end



function computeIA4(h, b, p, n)

    # Get unique haplotypes
    V = unique(h[:, b]; dims=1)
    tmpdict = Dict{Array{Int8},Int}(eachrow(V) .=> 1:size(V, 1))
    countmat = zeros(Float64, size(V, 1), length(n)) .+ 0.00000001

    # Calculate IA
    for (popi, pop) in enumerate(keys(p))
        for i in p[pop]
            countmat[tmpdict[h[i, b]], popi] += 1
        end
        countmat[:, popi] = countmat[:, popi] ./ n[popi]
    end
    
    #this computes freq as mean of freqs
    p_bar_v = mean.(eachrow(countmat))

    IA = sum(sum(countmat .* log.(countmat)) ./ length(n)) - sum(p_bar_v .* log.(p_bar_v))
    return IA
end

function computeIA5(h, b, p, n)

    countmat, _ = haplotypeFrequencies(h,b,p,n)
    #this computes freq as mean of freqs
    p_bar_v = mean.(eachrow(countmat))

    IA = sum(sum(countmat .* log.(countmat)) ./ length(n)) - sum(p_bar_v .* log.(p_bar_v))
    return IA
end

function haplotypeFrequencies(h::Matrix{Int8}, b::UnitRange{Int}, p::Dict{String, Vector{Int}}, n::Vector{Int})
    # Get unique haplotypes
    V = unique(h[:, b]; dims=1)
    tmpdict = Dict{Array{Int8},Int}(eachrow(V) .=> 1:size(V, 1))
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

    return Dict{Vector{Int8}, Vector{Float64}}(keys(tmpdict) .=> eachrow(countmat))
end

