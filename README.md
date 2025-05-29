# LocalAncestry.jl: Local Ancestry Inference using a combination of Naive Bayes Classification and Hidden Markov Model assignment


[![Build Status](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml?query=branch%3Amain)


The aim of this package is to provide the users with an accurate, fast, and accessible means of inferring local ancestries. The current version of the inference approach is set apart from existing local ancestry software by allowing for varying window sizes across the chromosome, and by using both Naive Bayes Classification and Hidden Markov Models for assign local ancestries.

**Local Ancestry Inference:** Local Ancestry Inference is the prediction of ancestry for each combination of single nucleotide polymorphism and individual. The contrast is Global Ancestry Inference, which is the prediction of ancestral proportions at the individual level.

## Current state of the package
We aim to have the first stable version of the package by October 1st, 2025.
The plan is to only change the under-the-hood functions, such that the calls for the interface function (LocalAncestry.origins()) stays the same.

**Planned changes:**
- Include approaches for creating priors
- Parallelize computationally intensive steps.
- Expand the number of methods for calculating priors for Naive Bayes Classification.
- Remove steps that don't consistently improve the accuracy of the local ancestry assignment. 

## Installation
The package can be installed through github:
```julia
using Pkg
Pkg.add(url = "https://github.com/BgroveP/LocalAncestry.jl", rev="0.1.0")
```

## Input
- The path to a Variant Call Format (.vcf) file with phased genotypes of reference individuals
- The path to a .vcf file with phased genotypes of target individuals
- A DataFrames.DataFrame object with ancestry of reference individuals (one ancestry per individual)
- An integer object that declares the focal chromosome.
- A string object that declares the approach for calculating priors for the Naive Bayes Classification step.
- An integer object that declares the minimum size of haplotype blocks.
- A float object that declares the stopping criterion for building haplotype blocks.
- An integer object that declares the ploidity.
- A float object that declares the minimal threshold for the Naive Bayes Classification step.

## Local ancestry inference with the LocalAncestry.origins function
The LocalAncestry.origins function provides a high-level interface, and enables the user to perform the entire analysis with one function.  

## Local ancestry inference with manual calls of subfunctions

## Comparison with other softwares
Work in progress.
