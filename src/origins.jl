"""
    getLocalAncestry(chromosome, referenceVCF, targetVCF, referenceAncestries; priorsMethod = "flat", minBlockSize = 5, incrBlockSize = 1, blockCrit = 0.2, minNBCProb = 0.95)

# Purpose
This function infers local ancestries. It is meant as a one-function interface to the entire inference process. 

# Arguments
- `chromosome::Int64`: The focal chromosome.
- `referenceVCF::String`: The relative path to .vcf file with phased genotypes of reference individuals.
- `targetVCF::String`: The relative path to .vcf file with phased genotypes of target individuals.
- `referenceAncestries::DataFrame`: Two-column (["individual", "ancestry"]) DataFrame with ancestries of reference individuals.
- `priorsMethod::String`: The method for calculating priors for the Naive Bayes Classification step (flat, CGR). We recommend flat priors for now.
- `minBlockSize::Int64`: The minimal size of haplotype blocks.
- `incrBlockSize::Int64`: The incremental size of haplotype blocks.
- `blockCrit::Float64`: The stopping criterion for building haplotype blocks. Smaller values provide larger haplotype blocks.
- `minNBCProb::Float64`: The lower threshold for posterior probabilities. Posterior probabilities above this threshold is assigned with the Naive Bayes Classification step, while those below the threshold will be assigned with the Hidden Markov step. 

# Returns
- `postProb::OrderedDict{String, Vector{OrderedDict{String, Float64}}}`: The area of the rectangle.
- `postClass::OrderedDict{String, Vector{String}}`: The area of the rectangle.
- `haplotypeLibrary::OrderedDict{}`: The area of the rectangle.

"""
function getLocalAncestry(chromosome, referenceVCF, targetVCF, referenceAncestries; priorsMethod = "flat", minBlockSize = 5, incrBlockSize = 1, blockCrit = 0.2, minNBCProb = 0.95)

    # Constants
    ploidity = 2
    predictType = "probonly"
    assignType = "Hidden Markov"

    ## Read haplotype data
    referenceData, referenceIndividuals = readVCF(referenceVCF, chromosome)
    targetData, targetIndividuals = readVCF(targetVCF, chromosome)
    referenceAncestriesVector = haplotypeOrigins(referenceIndividuals, referenceAncestries)

    # Get population information
    populations = unique(referenceAncestriesVector)
    popDict = getPopulationDictionary(referenceAncestriesVector)

    # Get haplotype library
    haplotypeLibrary, nHaplotypeBlocks = getHaploBlocks(minBlockSize, incrBlockSize, blockCrit, referenceData, popDict, 1)

    # Get priors
    priorProb, priorLevel = getPriors(referenceAncestriesVector, referenceData, targetIndividuals, targetData, priorsMethod)

    # Get log-likelihoods
    LL = calculateBlockFrequencies(haplotypeLibrary, referenceData, popDict)

    # predict
    postProb = getProbabilities(predictType, targetIndividuals, ploidity, LL, populations, nHaplotypeBlocks, priorProb, priorLevel, probStayState, targetData)
    postClass = getAssignments(assignType, postProb, populations, LL, minNBCProb, targetIndividuals, ploidity, haplotypeLibrary)

    return postProb, postClass, haplotypeLibrary
end

function getPriors(referenceAncestriesVector, referenceData, targetIndividuals, targetData, priorsMethod)

    if priorsMethod == "flat"
        x, n = priorsFlat(referenceAncestriesVector, targetIndividuals)
    elseif priorsMethod[1:3] == "CGR"
        x, n = priorsCGR(referenceData, targetData, targetIndividuals, referenceAncestriesVector, priorsMethod)
    else
        throw(DomainError(priorsMethod, "Expected 'flat' or 'CGR'"))
    end

    return x, n
end

function getProbabilities(predictType, targetIndividuals, ploidity, LL, populations, nHaplotypeBlocks, priorProb, priorLevel, probStayState, targetData)

    # Instantiate output
    postProb = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [OrderedDict(populations .=> [Vector{Union{Missing,Float64}}(missing, nHaplotypeBlocks) for i in 1:length(populations)]) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    if predictType == "Naive Bayes"
        predictNaiveBayes!(postProb, targetData, targetIndividuals, ploidity, LL, populations, nHaplotypeBlocks, priorProb, priorLevel)
    elseif predictType == "Hidden Markov"
        predictHiddenMarkov!(postProb, LL, populations, targetIndividuals, probStayState, targetData, ploidity=ploidity)
    else
        throw(DomainError(predictType, "Expected 'Naive Bayes' or 'Hidden Markov'"))
    end
    return postProb
end

function getAssignments(assignType, postProb, populations, LL, minNBCProb, targetIndividuals, ploidity, haplotypeLibrary)

    # Instantiate output
    postClass = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    # Assign certain blocks
    assignCertain!(postClass, postProb, populations, LL, minNBCProb)

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
