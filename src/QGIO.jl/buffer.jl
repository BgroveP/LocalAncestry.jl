
function create_buffer()
    return zeros(UInt8, READLINE_BUFFER_SIZE)
end

function print_buffer(b)
    i = findnext((x -> x == 0x00), b, 1) - 1
    if i > 0
        println(join(Char.(b[1:i])))
    end
    return nothing
end

function _reset_buffer!(b)
    b .= 0x00
    return nothing
end

function _buffer_entry(b, n)
    i = 1
    if n > 1
        for _ in 1:(n-1)
            i = findnext(x -> x == UInt8('\t'), b, i) + 1
        end
    end
    return String(b[i:(findnext(x -> x == UInt8('\t'), b, i)-1)])
end

# Simpler interfaces
function _buffer_chromosome(b)
    return _buffer_entry(b, 1)
end

function _buffer_identifier(b)
    return _buffer_entry(b, 3)
end

function _buffer_field(b)
    return _buffer_entry(b, 9)
end

function _buffer_haplotypes!(c, b)
    # Initialize
    bufferposition = 1
    # Get necessary information
    fields = _buffer_field(b)
    nfields = sum(fields .== ':') + 1
    gtfield = _get_entry("GT", fields, ':')
    # Scroll to data (one column before data )
    for _ in 1:8
        bufferposition = findnext(x -> x == UInt8('\t'), b, bufferposition) + 1
    end
    # Get haplotypes1
    if gtfield == 1
        for i in eachindex(c)[1:2:end]
            bufferposition = findnext(x -> x == UInt8('\t'), b, bufferposition) + 1
            @views c[i] = b[bufferposition] - 48
            @views c[i+1] = b[bufferposition+2] - 48
        end
    elseif gtfield > 1
        error("Not implemented yet that GT field is not first")
    else
        error("Did not find GT field in buffer")
    end
end

function _get_entry(p, v, d)

    o = 0
    n = length(p)
    m = length(v)
    for i in eachindex(v)
        if v[i] == d
            o += 1
        elseif v[i:min(m, i + n - 1)] == p
            if ((i == (m - n + 1)) || (v[i+1] == d)) & ((i == 1) || (v[i-1] == d))
                o += 1
                break
            end
        end
        if i == m
            o = 0
        end
    end
    return o
end
