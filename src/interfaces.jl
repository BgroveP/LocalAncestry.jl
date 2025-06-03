"""
    get_local_ancestries(chromosome::Union{Int64,String}, 
                          referenceVCF::String, 
                          targetVCF::String, 
                          referenceAncestries::DataFrame; 
                          priorsMethod::String = "flat", 
                          minBlockSize::Int64 = 5, 
                          incrBlockSize::Int64 = 1, 
                          blockCrit::Float64 = 0.2, 
                          minNBCProb::Float64 = 0.95)
# Purpose
This function infers local ancestries. It is meant as a one-function interface to the entire inference process. 

# Arguments
- `chromosome::Union{Int64, String}`: The focal chromosome. Autosomal chromosomes can be denoted by their number as e.g.: 1, "1", or "chr1".
- `referenceVCF::String`: The relative path to .vcf file with phased genotypes of reference individuals.
- `targetVCF::String`: The relative path to .vcf file with phased genotypes of target individuals.
- `referenceAncestries::DataFrame`: Two-column (["individual", "ancestry"]) DataFrame with ancestries of reference individuals.
- `priorsMethod::String`: The method for calculating priors for the Naive Bayes Classification step (flat, CGR). We recommend flat priors for now.
- `minBlockSize::Int64`: The minimal size of haplotype blocks.
- `incrBlockSize::Int64`: The incremental size of haplotype blocks.
- `blockCrit::Float64`: The stopping criterion for building haplotype blocks. Smaller values provide larger haplotype blocks.
- `minNBCProb::Float64`: The lower threshold for posterior probabilities. Posterior probabilities above this threshold is assigned with the Naive Bayes Classification step, while those below the threshold will be assigned with the Hidden Markov step. 

# Returns
- `postProb::OrderedDict{String, Vector{OrderedDict{String, Float64}}}`: The posterior probabilities from the Naive Bayes step.
- `postClass::OrderedDict{String, Vector{String}}`: The assigned populations after the Hidden Markov model step.
- `haplotypeLibrary::OrderedDict{}`: The library of haplotype blocks.

"""
function get_local_ancestries(
    chromosome::Union{Int64,String},
    referenceVCF::String,
    targetVCF::String,
    referenceAncestries::DataFrame;
    priorsMethod::String="flat",
    minBlockSize::Int64=5,
    incrBlockSize::Int64=1,
    blockCrit::Float64=0.2,
    minNBCProb::Float64=0.95,
)

    # Constants
    ploidity = 2
    predictType = "Naive Bayes"
    assignType = "Hidden Markov"
    probStayState = 0.99

    ## Read haplotype data
    referenceData, referenceIndividuals = LocalAncestry.readVCF(referenceVCF, chromosome)
    targetData, targetIndividuals = LocalAncestry.readVCF(targetVCF, chromosome)
    referenceAncestriesVector = LocalAncestry.haplotypeOrigins(
        referenceIndividuals, referenceAncestries
    )

    # Get population information
    populations = unique(referenceAncestriesVector)
    popDict = LocalAncestry.getPopulationDictionary(referenceAncestriesVector)

    # Get haplotype library
    haplotypeLibrary, nHaplotypeBlocks = LocalAncestry.getHaploBlocks(
        minBlockSize, incrBlockSize, blockCrit, referenceData, popDict, 1
    )

    # Get priors
    priorProb, priorLevel = LocalAncestry.getPriors(
        referenceAncestriesVector,
        referenceData,
        targetIndividuals,
        targetData,
        priorsMethod,
    )

    # Get log-likelihoods
    LL = LocalAncestry.calculateBlockFrequencies(haplotypeLibrary, referenceData, popDict)

    # predict
    postProb = LocalAncestry.getProbabilities(
        predictType,
        targetIndividuals,
        ploidity,
        LL,
        populations,
        nHaplotypeBlocks,
        priorProb,
        priorLevel,
        probStayState,
        targetData,
    )
    postClass = LocalAncestry.getAssignments(
        assignType,
        postProb,
        populations,
        LL,
        minNBCProb,
        targetIndividuals,
        ploidity,
        haplotypeLibrary,
    )

    return postProb, postClass, haplotypeLibrary
end

function getPriors(
    referenceAncestriesVector, referenceData, targetIndividuals, targetData, priorsMethod
)
    if priorsMethod == "flat"
        x, n = priorsFlat(referenceAncestriesVector, targetIndividuals)
    elseif priorsMethod[1:3] == "CGR"
        x, n = priorsCGR(
            referenceData,
            targetData,
            targetIndividuals,
            referenceAncestriesVector,
            priorsMethod,
        )
    else
        throw(DomainError(priorsMethod, "Expected 'flat' or 'CGR'"))
    end

    return x, n
end

function getProbabilities(
    predictType,
    targetIndividuals,
    ploidity,
    LL,
    populations,
    nHaplotypeBlocks,
    priorProb,
    priorLevel,
    probStayState,
    targetData,
)

    # Instantiate output
    postProb = OrderedDict(
        zip(
            [i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
            [
                OrderedDict(
                    populations .=> [
                        Vector{Union{Missing,Float64}}(missing, nHaplotypeBlocks) for
                        i in 1:length(populations)
                    ],
                ) for l in 1:length(targetIndividuals) for h in 1:ploidity
            ],
        ),
    )

    if predictType == "Naive Bayes"
        predictNaiveBayes!(
            postProb,
            targetData,
            targetIndividuals,
            ploidity,
            LL,
            populations,
            nHaplotypeBlocks,
            priorProb,
            priorLevel,
        )
    elseif predictType == "Hidden Markov"
        predictHiddenMarkov!(
            postProb,
            LL,
            populations,
            targetIndividuals,
            probStayState,
            targetData;
            ploidity=ploidity,
        )
    else
        throw(DomainError(predictType, "Expected 'Naive Bayes' or 'Hidden Markov'"))
    end
    return postProb
end

function getAssignments(
    assignType,
    postProb,
    populations,
    LL,
    minNBCProb,
    targetIndividuals,
    ploidity,
    haplotypeLibrary,
)

    # Instantiate output
    postClass = OrderedDict(
        zip(
            [i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
            [
                Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for
                l in 1:length(targetIndividuals) for h in 1:ploidity
            ],
        ),
    )

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
