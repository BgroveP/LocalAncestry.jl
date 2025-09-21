function mean(x)
    return sum(x) / length(x)
end



function vecsplit(x::Union{UnitRange{Int}, AbstractVector}, n::Int)
    y = Vector{typeof(x)}(undef, n)
    z::Int = length(x)
    r::Int = z
    for i in eachindex(y)
        thistake = ((z-r + 1):z)[1:(Int(ceil(r / (n+1-i))))]
        y[i] = x[thistake]
        r = r - length(thistake)
    end
    return y
end

# Taken from StatsBase (renamed. Thank you!)
function breapeat(vals::AbstractVector{T}, lens::AbstractVector{<:Integer}) where T
    m = length(vals)
    mlens = length(lens)
    mlens == m || throw(DimensionMismatch(
                        "number of vals ($m) does not match the number of lens ($mlens)"))
    n = sum(lens)
    n >= 0 || throw(ArgumentError("lengths must be non-negative"))

    r = Vector{T}(undef, n)
    p = 0
    @inbounds for i = 1 : m
        j = lens[i]
        j >= 0 || throw(ArgumentError("lengths must be non-negative"))
        v = vals[i]
        while j > 0
            r[p+=1] = v
            j -=1
        end
    end
    return r
end

function pretty!(assignments, chromosome, loci)
    assignments.chromosome .= chromosome
    assignments.basepairs = UnitRange.(loci.position[first.(assignments.block)], loci.position[last.(assignments.block)])
end