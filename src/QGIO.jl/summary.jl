function allelefreq!(loci,
    path;
    ancestries::DataFrame=DataFrame(individual=String[], population=String[]),
    omits::DataFrame=DataFrame(individual=String[], haplotype=Int[]))

    # 
    println("Calculating allele counts (allelecount) and frequencies (allelefreq)")

    sdf = _mergesamples(path, ancestries, omits)
    
    #
    pops = unique(sdf.population)
    ipops = [sdf.index[sdf.population.==p] for p in pops]
    locus_dictionary = Dict(collect(zip(loci.chromosome, loci.identifier)) .=> 1:nrow(loci))
    for p in pops
        loci[:, "allelecount_"*p] .= 0
    end

    # Read
    file = QGIO.open_vcf(path)
    buffer = QGIO.create_buffer()
    haplotypes = zeros(Int8, PLOIDITY * length(samples))
    locusentry = ("0", "0")
    while QGIO.readline!(file, buffer)
        if buffer[1] != UInt8('#')
            locusentry = (_buffer_chromosome(buffer), _buffer_identifier(buffer))
            _buffer_haplotypes!(haplotypes, buffer)

            # For each population

            for (ip, p) in enumerate(pops)
                @views loci[locus_dictionary[locusentry], "allelecount_"*p] = sum(haplotypes[ipops[ip]])
            end
        end
    end

    # 
    loci[:,"allelecount_across"] .= 0
    for (ip,p) in enumerate(pops)
        loci[:,"allelecount_across"] += loci[:,"allelecount_$(p)"]
        loci[:,"allelefreq_$(p)"] = loci[:,"allelecount_$(p)"] ./ length(ipops[ip])
    end
    loci[:,"allelefreq_across"] =  loci[:,"allelecount_across"] ./ sum(length.(ipops))

    return nothing
end

function inform_for_assign!(loci; pops = String[], mode = "standard", scale = "asis")

    s = "allelefreq_"
    o = "infoforassign_"
    if pops == String[]
        pops = replace.([i for i in names(loci) if occursin(s, i) & (i != "$(s)across")], "allelefreq_" => "")
    end
    pops = sort(pops)
    println("populations: ", join(pops, ", "))
    pbar = sum([loci[:,s * p] for p in pops]) ./ length(pops)

    for p in pops
        loci[:,o * p] = (loci[:,s * p] .+ NEARZERO_FLOAT) .* log.((loci[:,s * p] .+ NEARZERO_FLOAT)) / length(pops) - (pbar .+ NEARZERO_FLOAT) .* log.((pbar .+ NEARZERO_FLOAT))
    end
    loci[:,o * "across"] .= 0.0
    for p in pops
        loci[:,o * "across"] += loci[:,o * p]
    end
end