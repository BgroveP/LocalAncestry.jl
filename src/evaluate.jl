
function evaluate(predla,truela)

    nhaplotypes = nrow(unique(predla[:,["individual", "haplotype"]]))
    out = DataFrame(locus = 1:maximum(truela.last), ncorrect = 0)
    for i in axes(truela, 1)
        
        tmp = predla[(predla.individual .== truela.individual[i]) .& 
        (predla.haplotype .== truela.haplotype[i]) .& 
        (length.(intersect.(UnitRange.(truela.first[i],truela.last[i]), predla.block)) .> 0) .& 
        (predla.ancestry .== truela.ancestry[i]), :] 
        for j in axes(tmp,1)
            out.ncorrect[tmp.block[j]] .+= 1
        end
    end
    out.pcorrect = out.ncorrect/nhaplotypes
    return out
end
