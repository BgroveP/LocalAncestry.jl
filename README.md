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
getLocalAncestries(chromosome::Union{Int64,String}, 
                          referenceVCF::String, 
                          targetVCF::String, 
                          referenceAncestries::DataFrame; 
                          priorsMethod::String = "flat", 
                          minBlockSize::Int64 = 5, 
                          incrBlockSize::Int64 = 1, 
                          blockCrit::Float64 = 0.2, 
                          minNBCProb::Float64 = 0.95)
```
This function loads phased genotypes for both reference individuals and target individuals, 
builds the library of haplotype blocks, 
predicts local ancestry using Naive Bayes Classification, 
and lastly it reinforces the assignment using Hidden Markov models.

The walkthrough of the function is based on the toy data in the *test/data/* folder of the repository.