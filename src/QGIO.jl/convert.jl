function string2UInt8(s::AbstractString)
    o = Vector{UInt8}(undef, length(s))
    for (i, c) in enumerate(s)
        o[i] = UInt8(c)
    end
    return o
end

function convert_chromosome(c, loci)
    chrom =  ["", "chr"] .* string.(repeat([c], 2))

    i = [any(i .== loci.chromosome) for i in chrom]
    return chrom[i][1]
end
