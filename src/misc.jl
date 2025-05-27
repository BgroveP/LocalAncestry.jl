
function vvToV(x)
    # Initialize an empty Vector{Int64} to store the result
    flattened_vector = Int64[]

    # Iterate over each vector of ranges
    for range_vector in x
        # Iterate over each range in the current vector
        for range in range_vector
            # Append each integer in the range to the flattened vector
            append!(flattened_vector, collect(range))
        end
    end

    return flattened_vector
end


function haplotypeOrigins(i::Vector{String}, o::DataFrames.DataFrame)
    x = DataFrames.DataFrame(individual=repeat(i, inner=2))
    y = string.(DataFrames.leftjoin(x, o, on="individual")[!, "population"])
    
    # Asserts
    assertVector(y, "The vector with origins of reference individuals")
    
    # Return
    return y
end


rangeChange(rObject; firstInc=0, firstDec=0, lastInc=0, lastDec=0) = (first(rObject)+firstInc-firstDec):(last(rObject)+lastInc-lastDec)

function countThisHaploNumber(X, t)
    count = 0
    for x in eachrow(X)
        if all(x .== t)
            count += 1
        end
    end #Loop over training set
    return count
end

function getDictionaryCrosssection(x, y)
    return getindex.(collect(values.(Ref(x))), y)
end

function searchForward(est0_ind, prob0_ind, pos)
    countForward = 1
    while ismissing(est0_ind[pos+countForward])
        countForward += 1
        if (pos + countForward) > length(est0_ind)
            countForward = 0
            break
        end
    end
    forward = getDictionaryCrosssection(prob0_ind, pos + countForward)
    return forward, countForward
end

function searchBackwards(est0_ind, prob0_ind, pos)
    countBackwards = 1
    while ismissing(est0_ind[pos-countBackwards])
        countBackwards += 1
        if (pos - countBackwards) < 1
            countBackwards = 0
            break
        end
    end
    #    println("counted $countBackwards steps backwards for pos $pos")
    backwards = getDictionaryCrosssection(prob0_ind, pos - countBackwards)
    return backwards, countBackwards
end


"""
    rangeFromString(x)

...
# Arguments
- `x::String`: The range that should be converted from string.
...

Converts a integer range from string to actual range.
Can only convert to ranges with increments of one.

# Examples
```julia-repl
julia> BoA.rangeFromString("1:2")
1:2

julia> BoA.rangeFromString("1:10")
1:10

julia> typeof(BoA.rangeFromString("1:10"))
UnitRange{Int64}
```
"""
function rangeFromString(x)
    startAndStop = parse.(Int, split(x, ":"))
    range = startAndStop[1]:startAndStop[2]
    return range
end

function getPopulationDictionary(x)
    y = Dict{String,Vector{Int64}}()
    for i in unique(x)
        y[i] = findall(i .== x)
    end
    return y
end

function getPopulations(x::Vector{String})::Vector{String}
    return unique(x)
end

function alleleFrequencies(x, y)
    pops = getPopulations(y)
    p = zeros(Float32, size(x, 2), length(pops))
    for (j, pop) in enumerate(pops)
        rows = findall(pop .== y)
        p[:, j] = Statistics.mean(x[rows, :], dims=1)
    end

    return p
end

function instantiateOutput()

    postClass = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    return postProb, postClass
end


function calculateBlockFrequencies(haplotypeLibrary, referenceData, popDict)
    LL = OrderedDict()
    for (region, Haplo) in haplotypeLibrary
        LL[region] = getLL(region, Haplo, referenceData, popDict)
    end

    return LL
end

function mean(x)
    return sum(x) / length(x)
end
