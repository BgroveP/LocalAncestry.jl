
function skipline!(s::GZip.GZipStream, buf::Vector{UInt8})
    while !eof(s)
        # Read a chunk
        ptr = gzgets(s, buf)
        ptr == C_NULL && error("Error reading from GZipStream")
        # Find the first newline in the chunk
        idx = findfirst(==(UInt8('\n')), buf)
        if idx !== nothing
            # Found a newline; skip to the next line
            return true
        end
    end
    # Reached EOF without finding a newline
    return false
end

function skipline!(s::IOStream)
    _ = readline(s)
end
