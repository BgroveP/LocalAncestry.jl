
function check_vcfinfo(v)

    for f in names(v)[(length(VCF_HEADER) + 1):end]

        if any(.~v[:,f])
            error("One or more loci is not present in all files")
        end
    end
    return nothing
end

function check_individuals(i)
    if length(unique(vcat(i...))) != sum(length.(unique.(i)))
        error("One or more individuals are present in more than one file")
    end
    return nothing
end