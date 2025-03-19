

#gets unique haplotypes from the data of breeds for a given block
#n or N refers to haplotype number
function computeIA(h, p)

    # Get unique haplotypes
    V = unique(h,dims=1)

    # Number of haplotypes per population
    n = zeros(Int64, length(p))
    for (i, k) in enumerate(keys(p))
        n[i] = length(p[k])
    end

    # Calculate IA
    IA = 0
    for v in eachrow(V)
        count_v = []
        p_v  = [] #within pop freq
        for (i, k) in enumerate(keys(p))
            count_b_v = 0
            for r in p[k]
                if all(h[r,:] .== v)
                        count_b_v += 1
                else nothing # Why nothing?
                end
            end
            push!(count_v, count_b_v)
            push!(p_v, count_b_v / n[i])
        end
        #this computes freq as mean of freqs
        p_bar_v = mean(p_v)

        p_v_n0 = p_v[p_v.!=0]

        IA  += (-1.0*p_bar_v*log(p_bar_v) + sum(p_v_n0 .* log.(p_v_n0))/length(n))
    end
    return IA,V
end
