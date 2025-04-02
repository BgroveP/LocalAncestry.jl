

function origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb; assignType="probonly", predictType="Naive Bayes", probStayState=0.9)


    ## Read haplotype data
    referenceData, referenceIndividuals = ARV.readVCF(reference_path, chromosome)
    targetData, targetIndividuals = ARV.readVCF(target_path, chromosome)
    referenceOriginsVector = haplotypeOrigins(referenceIndividuals, referenceOrigins)

    # Get population information
    populations = getPopulations(referenceOriginsVector)
    popDict = getPopulationDictionary(referenceOriginsVector)

    # Get haplotype library
    haplotypeLibrary, nHaplotypeBlocks = getHaploBlocks(minHaploSize, incHaploSize, haploCrit, referenceData, popDict, 1)

    # Get priors
    priorProb, priorLevel = getPriors(referenceOriginsVector, referenceData, targetIndividuals, targetData, originPriors)

    # Get log-likelihoods
    LL = calculateBlockFrequencies(haplotypeLibrary, referenceData, popDict)

    # predict
    postProb = getProbabilities(predictType, targetIndividuals, ploidity, LL, populations, nHaplotypeBlocks, minProb, priorProb, priorLevel, probStayState, targetData)
    postClass = getAssignments(assignType, postProb, populations, LL, minProb, targetIndividuals, ploidity, haplotypeLibrary)

    return postProb, postClass, haplotypeLibrary
end

function getPriors(referenceOriginsVector, referenceData, targetIndividuals, targetData, originPriors)

    if originPriors == "flat"
        x, n = priorsFlat(referenceOriginsVector, targetIndividuals)
    elseif originPriors[1:3] == "CGR"
        x, n = priorsCGR(referenceData, targetData, targetIndividuals, referenceOriginsVector, originPriors)
    else
        throw(DomainError(originPriors, "Expected 'flat' or 'CGR'"))
    end

    return x, n
end

function getProbabilities(predictType, targetIndividuals, ploidity, LL, populations, nHaplotypeBlocks, minProb, priorProb, priorLevel, probStayState, targetData)

    # Instantiate output
    postProb = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [OrderedDict(populations .=> [Vector{Union{Missing,Float64}}(missing, length(haplotypeLibrary)) for i in 1:length(populations)]) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    if predictType == "Naive Bayes"
        predictNaiveBayes!(postProb, targetIndividuals, ploidity, LL, populations, nHaplotypeBlocks, priorProb, priorLevel)
    elseif predictType == "Hidden Markov"
        predictHiddenMarkov!(postProb, LL, populations, targetIndividuals, probStayState, targetData, ploidity=ploidity)
    else
        throw(DomainError(predictType, "Expected 'Naive Bayes' or 'Hidden Markov'"))
    end
    return postProb
end

function getAssignments(assignType, postProb, populations, LL, minProb, targetIndividuals, ploidity, haplotypeLibrary)

    # Instantiate output
    postClass = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    # Assign certain blocks
    assignCertain!(postClass, postProb, populations, LL, minProb)

    # Assign the rest
    if assignType == "none"
    elseif assignType == "probonly"
        assignFirst!(postClass, postProb)
    elseif assignType == "Hidden Markov"
        assignHMM!(postClass, postProb)
    else
        throw(DomainError(assignType, "Expected 'probonly' or 'Hidden Markov'"))
    end
    return postClass
end
