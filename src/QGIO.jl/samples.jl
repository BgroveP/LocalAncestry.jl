
function samples(path)
    # Open file and initialize string
    file = open_vcf(path)
    s::String = "1234"
    while ~eof(file)
        s = readline(file)
        if s[1:7] == "#CHROM\t"
            close(file)
            return string.(split(replace(s, '\n' => ""), '\t')[10:end])
        end
    end
    close(file)
    error("Reached end of file without finding samples")
end

function _mergesamples(samples, ancestries, omits; removeomits = true)

    # Setup vcf individuals
    sdf = DataFrame(individual=repeat(samples, inner=PLOIDITY),
        haplotype=repeat(1:PLOIDITY, outer=length(samples)),
        invcf=true,
        index=1:(2*length(samples)))

    # Setup omitted individuals
    odf = deepcopy(omits)
    odf[:, "omit"] .= true

    # Join information
    leftjoin!(sdf, odf, on=["individual", "haplotype"])
    leftjoin!(sdf, ancestries, on="individual")

    # Convert missing values
    sdf.omit = coalesce.(sdf.omit, false)
    sdf.population = coalesce.(sdf.population, "unknown")

    # Remove omitted individuals
    if removeomits
        deleteat!(sdf, findall(sdf.omit))
    end
    # Return
    return sdf
end
