
# Used
function haplotype_ancestries(i::Vector{String}, o::DataFrame)
    x = DataFrames.DataFrame(; individual=repeat(i; inner=2))
    return string.(DataFrames.leftjoin(x, o; on="individual")[!, "population"])
end


function get_pop_dict(x)
    y = Dict{String,Vector{Int64}}()
    for i in unique(x)
        y[i] = findall(i .== x)
    end
    return y
end

function mean(x)
    return sum(x) / length(x)
end


function string2UInt8(s::AbstractString)
    o = Vector{UInt8}(undef, length(s))
    for (i, c) in enumerate(s)
        o[i] = UInt8(c)
    end
    return o
end


function vecsplit(x::Vector, n::Int)
    y = Vector{typeof(x)}(undef, n)
    for i in 1:length(x)
        if i <= n
            y[i] = [x[i]]
        else
            z = (i % n) == 0 ? n : i % n 
            push!(y[z],x[i])
        end
    end
    return y
end

# Not used

function rangeChange(rObject; firstInc=0, firstDec=0, lastInc=0, lastDec=0)
    (first(rObject)+firstInc-firstDec):(last(rObject)+lastInc-lastDec)
end

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

function rangeFromString(x)
    startAndStop = parse.(Int, split(x, ":"))
    range = startAndStop[1]:startAndStop[2]
    return range
end

function alleleFrequencies(x, y)
    pops = getPopulations(y)
    p = zeros(Float32, size(x, 2), length(pops))
    for (j, pop) in enumerate(pops)
        rows = findall(pop .== y)
        p[:, j] = Statistics.mean(x[rows, :]; dims=1)
    end

    return p
end

function instantiateOutput()
    postClass = OrderedDict(
        zip(
            [i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
            [
                Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for
                l in 1:length(targetIndividuals) for h in 1:ploidity
            ],
        ),
    )

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
