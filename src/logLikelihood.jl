#returns LL for each haplo within a region
function getLL(region,haplotype,data, popDict)
    #store LL
    LogL = Dict{Vector{Int8}, Dict{String, Float64}}()
    #iterate for each region
#        println(region," ",size(Haplo,2)) #key is region
        nUniqueHaplo = size(haplotype,1)
        #total haplotypes in this breed + total unique haploptypes in the population
        #total haplotypes is equal to total n for this breed
        for (i,p) in enumerate(keys(popDict))
            sum_H_bh = length(popDict[p]) + size(nUniqueHaplo,1) # Why this??? size(nUniqueHaplo,1) = 1??
            for h in eachrow(haplotype)
                #how many times we observe a specific haplotype in each breed + 1
                H_bh = countThisHaploNumber(data[popDict[p],region], h) + 1
                logH_bh = log(H_bh/sum_H_bh)
                haskey(LogL, h) ? LogL[h][p] = logH_bh : LogL[h]=Dict(p=> logH_bh)
            end
        end
        return LogL
end
