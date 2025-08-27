
function getHaploBlocks(initSize::Int, stepSize::Int, threshold::Float64, data::Matrix{Int8}, popDict::Dict{String, Vector{Int}})

    nloci = size(data, 2)
    haploLib = OrderedDict{UnitRange{Int},Dict{Vector{Int8}, Vector{Float64}}}()

    n = zeros(Int, length(popDict))
    for (i, k) in enumerate(keys(popDict))
        n[i] = length(popDict[k])
    end
    oldBlock::UnitRange{Int} = 1:(initSize)
    newBlock::UnitRange{Int} = 1:(initSize+stepSize)
    notdone = true
    save = false
    oldIA::Float64 = LocalAncestry.computeIA(data, oldBlock, popDict, n)
    newIA::Float64 = 0.0

    while notdone
        # General case
        if last(newBlock) <= (nloci - initSize - stepSize)
            newIA = LocalAncestry.computeIA(data, newBlock, popDict, n)

            if (newIA - oldIA) > threshold
                oldBlock = newBlock
                oldIA = newIA
            else
                save = true
            end
        else
            notdone = false
            save = true
            oldBlock = first(oldBlock):nloci
        end

        if save
            counts, dict = LocalAncestry.haplotypeFrequencies(data, oldBlock, popDict, n)
            haploLib[oldBlock] = LocalAncestry.haplotypeDict(counts,dict)
            
            save = false
            oldBlock = (last(oldBlock) + 1):(min(last(oldBlock) + initSize, nloci))
            oldIA = LocalAncestry.computeIA(data, oldBlock, popDict, n)
        end
        newBlock = first(oldBlock):(last(oldBlock)+stepSize)
    end
    return haploLib
end

function getHaploBlocks2(initSize::Int, stepSize::Int, threshold::Float64, data::Matrix{Int8}, popDict::Dict{String, Vector{Int}})

    nloci = size(data, 2)
    locusmap = Matrix{Vector{Int}}(undef, nloci, 2)
    map!(x -> x .== 0, locusmap[:,1], eachcol(data) )
    haploLib = OrderedDict{UnitRange{Int},Dict{Vector{Int8}, Vector{Float64}}}()

    n = zeros(Int, length(popDict))
    for (i, k) in enumerate(keys(popDict))
        n[i] = length(popDict[k])
    end
    oldBlock::UnitRange{Int} = 1:(initSize)
    newBlock::UnitRange{Int} = 1:(initSize+stepSize)
    notdone = true
    save = false
    oldIA::Float64 = LocalAncestry.computeIA(data, oldBlock, popDict, n)
    newIA::Float64 = 0.0

    while notdone
        # General case
        if last(newBlock) <= (nloci - initSize - stepSize)
            newIA = LocalAncestry.computeIA(data, newBlock, popDict, n)

            if (newIA - oldIA) > threshold
                oldBlock = newBlock
                oldIA = newIA
            else
                save = true
            end
        else
            notdone = false
            save = true
            oldBlock = first(oldBlock):nloci
        end

        if save
            counts, dict = LocalAncestry.haplotypeFrequencies(data, oldBlock, popDict, n)
            haploLib[oldBlock] = LocalAncestry.haplotypeDict(counts,dict)
            
            save = false
            oldBlock = (last(oldBlock) + 1):(min(last(oldBlock) + initSize, nloci))
            oldIA = LocalAncestry.computeIA(data, oldBlock, popDict, n)
        end
        newBlock = first(oldBlock):(last(oldBlock)+stepSize)
    end
    return haploLib
end
