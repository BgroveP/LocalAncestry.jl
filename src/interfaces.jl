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
    minBlockSize::Int64=5,
    incrBlockSize::Int64=1,
    blockCrit::Float64=0.01,
    minNBCProb::Float64=0.95,
)

    # Constants
    ploidity = 2

    ## Read reference haplotype data
    println("Reading reference data")
    referenceData, referenceIndividuals = LocalAncestry.readVCF(referenceVCF, chromosome)
    referenceAncestriesVector = LocalAncestry.haplotypeOrigins(
        referenceIndividuals, referenceAncestries
    )

    # Get population information
    println("Getting population information")
    populations = unique(referenceAncestriesVector)
    popDict = LocalAncestry.getPopulationDictionary(referenceAncestriesVector)

    # Get haplotype library
    println("Constructing haplotype blocks")
    haplotypeLibrary = LocalAncestry.getHaploBlocks(
        minBlockSize, incrBlockSize, blockCrit, referenceData, popDict)

    # Load target haplotype data
    println("Reading target data")
    targetData, targetIndividuals = LocalAncestry.readVCF(targetVCF, chromosome)

    # Estimating Local ancestries
    println("Estimating local ancestries")
    postProb, postClass, classDict, probDict = lai(haplotypeLibrary, targetData, targetIndividuals, minNBCProb, length(popDict)::Int, ploidity)

    # Return
    return postProb, postClass, haplotypeLibrary, keys(popDict), classDict, probDict
end

function lai(haplotypeLibrary, targetData, targetIndividuals, minNBCProb, nPopulations::Int, ploidity)

    # Predict
    println(" - Predict initial ancestries")
    postProb, probDict, probRowdict = predict(haplotypeLibrary, targetData, targetIndividuals, nPopulations, ploidity)

    # Assign certain
    println(" - Assign certain ancestries")
    postClass, classDict = assignCertain(postProb, nPopulations, ploidity, targetIndividuals, minNBCProb)

    # Assign2
    println(" - Predict final ancestries")
    assign_missing!(postProb, postClass, nPopulations, probRowdict, ploidity)

    return postProb, postClass, classDict, probDict
end


function get_local_ancestries2(
    chromosome::Union{Int64,String},
    referenceVCF::String,
    targetVCF::String,
    referenceAncestries::DataFrame;
    minBlockSize::Int64=5,
    incrBlockSize::Int64=1,
    blockCrit::Float64=0.01,
    minNBCProb::Float64=0.95,
)

    # Constants
    ploidity = 2

    ## Read reference haplotype data
    println("Reading reference data")
    referenceData, _, referenceIndividuals = LocalAncestry.readVCF2(referenceVCF, chromosome)
    referenceAncestriesVector = LocalAncestry.haplotypeOrigins(
        string.(referenceIndividuals), referenceAncestries
    )

    # Get population information
    println("Getting population information")
    populations = unique(referenceAncestriesVector)
    popDict = LocalAncestry.getPopulationDictionary(referenceAncestriesVector)


    # Get haplotype library
    println("Constructing haplotype blocks")
    haplotypeLibrary = LocalAncestry.getHaploBlocks(
        minBlockSize, incrBlockSize, blockCrit, referenceData, popDict)

    # Load target haplotype data
    println("Reading target data")
    targetData, _, targetIndividuals = LocalAncestry.readVCF2(targetVCF, chromosome)

    # Estimating Local ancestries
    println("Estimating local ancestries")
    postProb, postClass, classDict, probDict = LocalAncestry.lai(haplotypeLibrary, targetData, targetIndividuals, minNBCProb, length(popDict)::Int, ploidity)

    # Return
    return postProb, postClass, haplotypeLibrary, keys(popDict), classDict, probDict
end
