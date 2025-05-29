# LocalAncestry.jl: Local Ancestry Inference using a combination of Naive Bayes Classification and Hidden Markov Model assignment

The aim of this package is to provide the users with an accurate, fast, and accessible means of inferring local ancestries. The current version of the inference approach is set apart from existing local ancestry software by allowing for varying window sizes across the chromosome, and by using both Naive Bayes Classification and Hidden Markov Models for assign local ancestries.

**Local Ancestry Inference:** Local Ancestry Inference is the prediction of ancestry for each combination of single nucleotide polymorphism and individual. The contrast is Global Ancestry Inference, which is the prediction of ancestral proportions at the individual level.

## Current state of the package
We aim to have the first stable version of the package by October 1st, 2025.
The plan is to only change the under-the-hood functions, such that the calls for the interface function (LocalAncestry.origins()) stays the same.
The current version is fully functional.

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
## Top-level tutorial

```julia
```