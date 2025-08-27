
function predict(haplotypeLibrary, targetData, targetIndividuals, nPopulations::Int, ploidity)

    nBlocks = length(haplotypeLibrary)
    nIndividuals = length(targetIndividuals)
    nHaplotypes = ploidity * nIndividuals
    nProbabilities = nPopulations * nHaplotypes
    postProb = zeros(Float16, nProbabilities, nBlocks)
    probDict = Dict{String,UnitRange{Int}}(targetIndividuals .=> UnitRange.(1:(nPopulations*ploidity):nProbabilities, (nPopulations*ploidity):(nPopulations*ploidity):nProbabilities))
    rowDict = Dict{Int,UnitRange{Int}}(1:nHaplotypes .=> UnitRange.(1:(nPopulations):nProbabilities, (nPopulations):(nPopulations):nProbabilities))

    # Predict
    for (j, b) in enumerate(keys(haplotypeLibrary))
        for i in 1:nHaplotypes
            thisRange = (nPopulations*(i-1)+1):(nPopulations*(i-1)+nPopulations)
            if haskey(haplotypeLibrary[b], targetData[i, b])
                postProb[thisRange, j] = haplotypeLibrary[b][targetData[i, b]]
            else
                postProb[thisRange, j] .= 1.0 / nPopulations
            end
        end
    end
    return postProb, probDict, rowDict
end

function predict2(haplotypeLibrary, targetData, targetIndividuals, nPopulations::Int, ploidity)

    nBlocks = length(haplotypeLibrary)
    nIndividuals = length(targetIndividuals)
    nHaplotypes = ploidity * nIndividuals
    nProbabilities = nPopulations * nHaplotypes
    postProb = zeros(Float16, nProbabilities, nBlocks)
    probDict = Dict{String,UnitRange{Int}}(targetIndividuals .=> UnitRange.(1:(nPopulations*ploidity):nProbabilities, (nPopulations*ploidity):(nPopulations*ploidity):nProbabilities))
    rowDict = Dict{Int,UnitRange{Int}}(1:nHaplotypes .=> UnitRange.(1:(nPopulations):nProbabilities, (nPopulations):(nPopulations):nProbabilities))

    # Predict
    write_lock = ReentrantLock()
    hapdict = Dict(1:length(haplotypeLibrary) .=> keys(haplotypeLibrary))
    Threads.@threads for j in 1:length(haplotypeLibrary)
        b = hapdict[j]
        for i in 1:nHaplotypes
            thisRange = (nPopulations*(i-1)+1):(nPopulations*(i-1)+nPopulations)
            if haskey(haplotypeLibrary[b], targetData[i, b])
                @lock write_lock postProb[thisRange, j] = haplotypeLibrary[b][targetData[i, b]]
            else
                @lock write_lock postProb[thisRange, j] .= 1.0 / nPopulations
            end
        end
    end
    return postProb, probDict, rowDict
end



function assignCertain(postProb, nPopulations, ploidity, targetIndividuals, minNBCProb)
    nHaplotypes::Int = size(postProb, 1) / nPopulations
    nBlocks = size(postProb, 2)
    postClass = zeros(Int8, nHaplotypes, nBlocks)
    classDict = Dict{String,UnitRange{Int}}(targetIndividuals .=> UnitRange.(1:ploidity:nHaplotypes, ploidity:ploidity:nHaplotypes))
    for h in 1:nHaplotypes
        thisRange = (nPopulations*(h-1)+1):(nPopulations*(h-1)+nPopulations)
        postClass[h, :] = [first(i) > minNBCProb ? last(i) : 0 for i in findmax.(eachcol(postProb[thisRange, :]))]
    end

    return postClass, classDict
end
