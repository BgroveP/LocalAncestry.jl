

function origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)
    # Function

    ## Derived information
    populations = unique(referenceOrigins[!, "population"])

    ## Data
    referenceData, referenceIndividuals = ARV.readVCF(reference_path, chromosome)
    targetData, targetIndividuals = ARV.readVCF(target_path, chromosome)
    referenceOriginsVector = ARV.haplotypeOrigins(referenceIndividuals, referenceOrigins)

    # reference populations 
    popDict = Dict{String,Vector{Int64}}()
    for i in unique(referenceOriginsVector)
        popDict[i] = findall(i .== referenceOriginsVector)
    end

    ## Priors
    if originPriors == "flat"
        priorProb = priorsFlat(populations, targetIndividuals)
    elseif originPriors == "CGR"
        priorProb = priorsCGR(referenceData, targetData, targetIndividuals, referenceOriginsVector, "auto")
    elseif originPriors == "CGRfull"
        priorProb = priorsCGR(referenceData, targetData, targetIndividuals, referenceOriginsVector, "full")
    elseif originPriors == "CGRsqrt"
        priorProb = priorsCGR(referenceData, targetData, targetIndividuals, referenceOriginsVector, "autosqrt")
    else
        throw(DomainError(originPriors, "Expected 'flat' or 'CGR'"))
    end

    ## Haplotype library
    haplotypeLibrary = ARV.getHaploBlocks(minHaploSize, incHaploSize, haploCrit, referenceData, popDict, 1)
    nHaplotypeBlocks = length(haplotypeLibrary)

    #Store, donor haplo is free of ID, already removed above
    postProb = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [OrderedDict(populations .=> [Vector{Union{Missing,Float64}}(missing, length(haplotypeLibrary)) for i in 1:length(populations)]) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    postClass = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for l in 1:length(targetIndividuals) for h in 1:ploidity]))


    # Get log-likelihoods
    LL = OrderedDict()
    for (region, Haplo) in haplotypeLibrary
        LL[region] = ARV.getLL(region, Haplo, referenceData, popDict)
    end

    # predict
    for (i, id) in enumerate(targetIndividuals)
        for h in 1:ploidity
            idname = id * "_hap" * string(h)

            postProb[idname], postClass[idname] = ARV.predInd(priorProb[idname], targetData[2*i+h-2, :], LL, populations, nHaplotypeBlocks, minProb)
            # Assign missing
            postClass[idname] = ARV.refineBoA(postClass[idname], postProb[idname])
        end
    end

    return postProb, postClass, haplotypeLibrary
end


