function translate_chromosome(c::Union{String,Int}, vcfinfo)
    ic = Vector{UInt8}()
    # Translate chromosome input to internal format
    if c == ""
        ic = string2UInt8(vcfinfo.chromosome[1])
    elseif any(c .== ("chr" * string(c)))
        ic = string2UInt8("chr" * string(c))
    elseif any(c .== (string(c)))
        ic = string2UInt8(string(c))
    else
        error("focal chromosome $(c) not found in the vcf file")
    end

    return ic
end

function string2UInt8(s::AbstractString)
    o = Vector{UInt8}(undef, length(s))
    for (i, c) in enumerate(s)
        o[i] = UInt8(c)
    end
    return o
end

