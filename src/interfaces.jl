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
    chromosome::Union{Int,AbstractString},
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestries::DataFrame;
    threshold::Float64=0.01,
    nbcprob::Float64=0.95,
)
    # Read input
    println("Reading .vcf files")
    t1 = @spawn LocalAncestry.readVCF(referencepath, chromosome)
    t2 = @spawn LocalAncestry.readVCF(targetpath, chromosome)
    refdata, refloci, refind = fetch(t1)
    refancestries = LocalAncestry.haplotype_ancestries(refind, ancestries)
    
    # Checks and prints
    
    # Get population information
    println("Getting population information")
    popDict = LocalAncestry.get_pop_dict(refancestries)
    
    # Get haplotype library
    println("Constructing haplotype blocks")
    library = LocalAncestry.get_haplotype_library(refdata, popDict, threshold)
    
    # Debug prints
    # Estimating Local ancestries
    println("Estimating local ancestries")
    targetdata, targetloci, targetind = fetch(t2)
    return assign(library, targetdata, targetind, nbcprob, popDict)
end