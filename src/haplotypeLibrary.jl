function getHaploBlocks(initSize, stepSize, threshold, data, popDict, startFrom=1)
    nCol = size(data, 2)
    haploLib = OrderedDict()
    while in(startFrom, 1:nCol)
        if (startFrom == nCol) || (startFrom + stepSize + initSize >= nCol)
            thisBlock = startFrom:nCol
            ia, haplotypes = computeIA(data[:, thisBlock], popDict)
            haploLib[thisBlock] = haplotypes
            break
        else
            key, value = haploSearch(
                initSize, stepSize, threshold, data, startFrom, popDict
            )
            haploLib[key] = value
            startFrom = last(key) + 1
        end
    end
    return haploLib, length(haploLib)
end


function getHaploBlocks2(initSize::Int, stepSize::Int, threshold::Float64, data::Matrix{Int8}, popDict::Dict{String, Vector{Int}})

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
    oldIA::Float64 = LocalAncestry.computeIA2(data, oldBlock, popDict, n)
    newIA::Float64 = 0.0

    while notdone
        # General case
        if last(newBlock) <= (nloci - initSize - stepSize)
            newIA = LocalAncestry.computeIA2(data, newBlock, popDict, n)

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
            oldBlock = (last(oldBlock) + 1):(min((last(oldBlock) + 1) + initSize, nloci))
            oldIA = LocalAncestry.computeIA2(data, oldBlock, popDict, n)
        end
        newBlock = first(oldBlock):(last(oldBlock)+stepSize)
    end
    return haploLib
end
