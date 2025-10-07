# LocalAncestry.jl: Fast and accurate local ancestry inference using Bayesian Classification and Hidden Markov Model

[![Build Status](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://github.com/BgroveP/ARV.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/BgroveP/ARV.jl/actions/workflows/documentation.yml?query=branch%3Amain)

The purpose of this package is to provide an accurate, fast, and accessible means of inferring local ancestries. 
The package is set apart from existing local ancestry software by allowing for varying window sizes across the chromosome, and by using a combination of Bayesian Classification and Hidden Markov Models to assign local ancestries.

**Local Ancestry Inference:** Local Ancestry Inference is the prediction of ancestry for each combination of single nucleotide polymorphism and individual. The contrast is Global Ancestry Inference, which is the prediction of ancestral proportions at the individual level.

## Installation
The package can be installed through github:
```julia
using Pkg
Pkg.add("LocalAncestry")
```

## Required inputs
1. A Variant Call Format file with phased genotypes for reference individuals. 
2. A Variant Call Format file with phased genotypes for target individuals. 
3. A delimited file with ancestries of reference individuals (one ancestry per individual). 

## Public API
The user of this package only need to run one function to estimate local ancestries:
```julia
function get_local_ancestries(
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestrypath::String;
    omitpath::String="",
    chromosome::Union{Int,AbstractString}="",
    threshold::Float64=0.13,
    nbcprob::Float64=0.95,
    maf::Float64=0.0001,
    printlevel::String="standard"
)

```
### Purpose
This function infers local ancestries. It is meant as a one-function interface to the entire inference process. 

### Arguments

- `referencepath::AbstractString`: The relative path to .vcf file with phased genotypes of reference individuals.
- `targetpath::AbstractString`: The relative path to .vcf file with phased genotypes of target individuals.
- `ancestrypath::String`: The relative path to a delimited file with ancestries of reference individuals with two columns: individual and population.
- `omitpath::String`: The relative path to a delimited file with omitted reference individuals with two columns: individual and haplotype.
- `chromosome::Union{Int,AbstractString}`: The focal chromosome as either integer or string: 1, "1", or "chr1".
- `threshold::Float64`: The lower threshold for the informativeness for assignment statistic when building haplotype blocks.
- `nbcprob::Float64`: The lower limit for posterior probabilities for the Bayesian Classification step.
- `maf::Float64`: The lower limit for the minor allele frequency among reference individuals (omission of loci).
- `printlevel::String`: The output level. 

### Returns
- `x::DataFrame`: A DataFrame object with columns ....

## Planned changes
- Expanded support for input and output to make incorporation into existing pipelines easier.
- Make the code adhere to the blue style for the Julia language.
