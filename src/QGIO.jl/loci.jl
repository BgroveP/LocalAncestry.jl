
function loci(path)
    
    # Initialize
    buffer = create_buffer()
    out = DataFrame(
        chromosome = Vector{String}(), 
        position = Vector{String}(), 
        identifier = Vector{String}(), 
        reference = Vector{String}(), 
        alternative = Vector{String}(), 
        quality = Vector{String}(), 
        info = Vector{String}(), 
        filter = Vector{String}(), 
        format = Vector{String}() 
    )

    # Get loci
    i1 = 1
    i2 = 1
    file = open_vcf(path)
    while readline!(file, buffer)
        if buffer[1] != UInt('#')
            resize!(out, nrow(out) + 1)
            for c in names(out)
                i2 = findnext((x -> x == UInt8('\t')), buffer, i1)
                out[nrow(out),c] = String(buffer[i1:(i2-1)])
                i1 = i2 + 1
            end
            i1 = 1
        end
    end
    close(file)
    out.position = parse.(Int, out.position)
    return out
end
