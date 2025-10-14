function allelefreq!(loci,
    path;
    ancestries::DataFrame=DataFrame(individual=String[], population=String[]),
    omits::DataFrame=DataFrame(individual=String[], haplotype=Int[]))

    # Initialize
    samples = QGIO.samples(path)
    locus_dictionary = Dict(collect(zip(loci.chromosome, loci.position, loci.identifier)) .=> 1:nrow(loci))

    # Merge information
    sdf = QGIO._mergesamples(samples, ancestries, omits)

    # Get populations
    pops = sort(unique(sdf.population))
    ipops = [sdf.index[sdf.population.==p] for p in pops]

    if length(pops) > 1
        outcols1 = "allelecount_" .* pops
    else
        outcols1 = ["allelecount"]
    end

    for oc in outcols1
        loci[:, oc] .= 0
    end

    # Read
    file = QGIO.open_vcf(path)
    buffer = QGIO.create_buffer()
    haplotypes = zeros(Int8, PLOIDITY * length(samples))
    locusentry = 0
    while QGIO.readline!(file, buffer)
        if buffer[1] != UInt8('#')
            locusentry = locus_dictionary[(QGIO._buffer_chromosome(buffer), QGIO._buffer_position(buffer), QGIO._buffer_identifier(buffer))]
            QGIO._buffer_haplotypes!(haplotypes, buffer)

            # For each population

            for ip in eachindex(pops)
                @views loci[locusentry, outcols1[ip]] = sum(haplotypes[ipops[ip]])
            end
        end
    end

    # Calculate across counts and frequencies
    outcols2 = replace.(outcols1, "allelecount" => "allelefreq")
    for (ioc, oc) in enumerate(outcols1)
        loci[:, outcols2[ioc]] = loci[:, oc] ./ length(ipops[ioc])
    end
    if length(outcols1) > 1
        loci[:, "allelecount"] = sum.(eachrow(loci[:, outcols1]))
        loci[:, "allelefreq"] = loci[:, "allelecount"] ./ sum(length.(ipops))
    end

    return nothing
end

function inform_for_assign!(loci; mode="standard", scale="asis")

    s = "allelefreq_"
    o = "infoforassign_"
    pops = replace.([i for i in names(loci) if occursin(s, i) & (i != "$(s)across")], "allelefreq_" => "")
    pops = sort(pops)
    pbar = loci.allelefreq

    if mode == "standard"
        for p in pops
            loci[:, o*p] = (loci[:, s*p] .+ NEARZERO_FLOAT) .* log.((loci[:, s*p] .+ NEARZERO_FLOAT)) / length(pops) - (pbar .+ NEARZERO_FLOAT) .* log.((pbar .+ NEARZERO_FLOAT))
        end
        loci[:, "infoforassign"] .= 0.0
        for p in pops
            loci[:, "infoforassign"] += loci[:, o*p]
        end
    else
        loci[:, "infoforassign"] .= log(2)
        pmat = loci[:, s.*pops] |> Matrix
        pmat = max.(pmat, NEARZERO_FLOAT)
        for i in eachindex(pops), j in eachindex(pops)
            if i < j
                pbar = sum.(eachrow(pmat[:,[i,j]])) ./ 2 
                loci[:, "infoforassign"] = min.(loci[:, "infoforassign"], sum.(eachrow(pmat[:, [i,j]] .* log.(pmat[:, [i,j]]))) ./ 2 .- sum.(eachrow(pbar .* log.(pbar))))
            end
        end
    end
    return nothing
end
