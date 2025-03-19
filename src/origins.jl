

function origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)
    # Function

    ## Derived information
    populations = unique(referenceOrigins[!, "population"])

    ## Data
    referenceData, referenceIndividuals = readVCF(reference_path, chromosome)
    targetData, targetIndividuals = readVCF(target_path, chromosome)
    referenceOriginsVector = haplotypeOrigins(referenceIndividuals, referenceOrigins)

    # reference populations 
    popDict = Dict{String,Vector{Int64}}()
    for i in unique(referenceOriginsVector)
        popDict[i] = findall(i .== referenceOriginsVector)
    end

    ## Priors
    priorProb = makePriors(populations, targetIndividuals, originPriors)

    ## Haplotype library
    haplotypeLibrary = getHaploBlocks(minHaploSize, incHaploSize, haploCrit, referenceData, popDict, 1)
    nHaplotypeBlocks = length(haplotypeLibrary)

    #Store, donor haplo is free of ID, already removed above
    postProb = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [OrderedDict(populations .=> [Vector{Union{Missing,Float64}}(missing, length(haplotypeLibrary)) for i in 1:length(populations)]) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    postClass = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for l in 1:length(targetIndividuals) for h in 1:ploidity]))


    # Get log-likelihoods
    LL = OrderedDict()
    for (region, Haplo) in haplotypeLibrary
        LL[region] = getLL(region, Haplo, referenceData, popDict)
    end

    # predict
    for (i, id) in enumerate(targetIndividuals)
        for h in 1:ploidity
            idname = id * "_hap" * string(h)

            postProb[idname], postClass[idname] = predInd(priorProb[id], targetData[2*i+h-2, :], LL, populations, nHaplotypeBlocks, minProb)
            # Assign missing
            postClass[idname] = refineBoA(postClass[idname], postProb[idname])
        end
    end

    return postProb, postClass
end


