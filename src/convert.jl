
function read_rfmix(file::String, loci)

    # Constants
    outcols = [:individual, :chromosome, :haplotype, :basepairs, :ancestry]

    # Reads
    ancestry = CSV.read(file, DataFrame, header=2, types=String)
    ancestry = stack(ancestry, 7:ncol(ancestry))
    populations = replace.(split(replace(readline(file), r".+: " => ""), "\t"), r"=.+" .=> [""])
    enddict = Dict{String,Int}(string.(loci.position[2:end]) .=> loci.position[1:(end-1)])

    ancestry[:, "chromosome"] = ancestry[:, "#chm"]
    ancestry[:, "individual"] = replace.(ancestry.variable, r"\..+" => "")
    ancestry[:, "haplotype"] = parse.(Int, replace.(ancestry.variable, r".+\." => "")) .+ 1
    ancestry[:, "basepairs"] = UnitRange.(parse.(Int, ancestry.spos), get.(Ref(enddict), ancestry.epos, missing))
    ancestry[:, "ancestry"] = get.(Ref(populations), parse.(Int, ancestry.value) .+ 1, missing)

    # Remove consecutive samples
    ancestry.startrow = axes(ancestry, 1)
    ancestry.endrow = axes(ancestry, 1)
    ancestry.discard .= false

    for i in axes(ancestry, 1)[2:end]
        if ancestry[i-1, outcols[[1, 2, 3, 5]]] == ancestry[i, outcols[[1, 2, 3, 5]]]

            ancestry.discard[i] = true
            ancestry.startrow[i] = ancestry.startrow[i-1]
            ancestry.endrow[ancestry.startrow[i]] = ancestry.endrow[i]
        end
    end
    
    
    ancestry.basepairs = UnitRange.(first.(ancestry.basepairs[ancestry.startrow]), last.(ancestry.basepairs[ancestry.endrow]))
    deleteat!(ancestry, ancestry.discard)

    # 
    return ancestry[:, outcols]
end


function write_flare(la, hapfile::String, outfile::String)

    # Read input
    populations = sort(unique(la.ancestry))
    header =
        [replace.(LocalAncestry.QGIO.header(hapfile), '\n' => "");
            ["##FORMAT=<ID=AN1,Number=1,Type=Integer,Description=\"Ancestry of first haplotype\">",
                "##FORMAT=<ID=AN2,Number=1,Type=Integer,Description=\"Ancestry of second haplotype\">",
                "##ANCESTRY=<$(join(populations .* "=" .* string.(collect(0:(length(populations)-1))), ","))>"
            ]
        ]

    loci = LocalAncestry.QGIO.loci(hapfile)
    haplotypes, samples = LocalAncestry.QGIO.haplotypes(hapfile)
    loci.format = loci.format .* ":AN1:AN2"
    popdict = Dict(populations .=> string.(0:(length(populations)-1)))
    indsort = [findfirst(unique(la.individual) .== i) for i in unique(samples.individual)]
   
    # Write
    ofile = GZip.open(outfile, "w")

    ## Header
    println.(ofile, header)

    ## Samples
    println(ofile, join(["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"; unique(samples.individual)], '\t'))

    # Filter loci
    posstart = minimum(first.(la.basepairs))
    posend = maximum(last.(la.basepairs))

    ## Body
    for r in axes(loci,1) 
        if (loci.position[r] >= posstart) & (loci.position[r] <= posend)
            a1 = get.(Ref(popdict), get_ancestry(la, 1, loci.position[r]), missing)
            a2 = get.(Ref(popdict), get_ancestry(la, 2, loci.position[r]), missing)
            if length(a1) == length(a2) == Int(nrow(samples)/2)
                println(ofile, join(loci[r,:], '\t') * '\t' * join(string.(haplotypes[1:2:end,r]) .* "|" .* string.(haplotypes[2:2:end,r]) .* ":" .* a1[indsort]  .* ":" .* a2[indsort], '\t'))
            end
        end
    end
    close(ofile)

    return nothing
end

