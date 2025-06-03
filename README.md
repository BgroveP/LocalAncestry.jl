# LocalAncestry.jl: Local Ancestry Inference using a combination of Naive Bayes Classification and Hidden Markov Model assignment

[![Build Status](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://github.com/BgroveP/ARV.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/BgroveP/ARV.jl/actions/workflows/documentation.yml?query=branch%3Amain)

The aim of this package is to provide the users with an accurate, fast, and accessible means of inferring local ancestries. 
The package is set apart from existing local ancestry software by allowing for varying window sizes across the chromosome, and by using a combination of Naive Bayes Classification and Hidden Markov Models to assign local ancestries.

**Local Ancestry Inference:** Local Ancestry Inference is the prediction of ancestry for each combination of single nucleotide polymorphism and individual. The contrast is Global Ancestry Inference, which is the prediction of ancestral proportions at the individual level.

## Installation
The package can be installed through github:
```julia
using Pkg
Pkg.add(url = "https://github.com/BgroveP/LocalAncestry.jl", rev="0.1.0")
```

## Required inputs
1. A Variant Call Format file with phased genotypes for reference individuals. 
2. A Variant Call Format file with phased genotypes for target individuals. 
3. A DataFrame object with ancestries of reference individuals (one ancestry per individual). 

## Public API
The user of this package only need to run one function to estimate local ancestries:
```julia
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

```
This function loads phased genotypes for both reference individuals and target individuals, 
builds the library of haplotype blocks, 
predicts local ancestry using Naive Bayes Classification, 
and lastly it reinforces the assignment using Hidden Markov models.

The walkthrough of the function is based on the toy data in the *test/data/* folder of the repository.