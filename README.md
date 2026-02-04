# LocalAncestry.jl: A Novel Tool for Local Ancestry Inference Using Genome Regions With Variable Lengths

[![Build Status](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BgroveP/ARV.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://github.com/BgroveP/ARV.jl/actions/workflows/documentation.yml/badge.svg?branch=main)](https://github.com/BgroveP/ARV.jl/actions/workflows/documentation.yml?query=branch%3Amain)

This Julia package provides an accurate, fast, and accessible means of inferring local ancestries. 
The package is set apart from existing local ancestry software by being a region-based local ancestry inference tool that allows for regions of variying sizes across the chromosome.

## Installation
The package is publicly available through the Julia Ecosystem:
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
get_local_ancestries(
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestrypath::String;
    ,
    ,
    nbcprob::Float64=0.99,
    maf::Float64=0.0001,
    printlevel::String="standard"
)
```

### Optional parameters
The optional parameters for the get_local_ancestries function are:
1. *omitpath::String=""* which is an absolute or relative path to a delimited file with information that denotes the haplotypes that should be excluded from the reference panel. The delimited file contains two columns: *individual* and *haplotype*.
2. *chromosome::Union{Int,AbstractString}=""* which is the focal chromosome for the local ancestry inference. The default behaviour is to analyze the chromosome that is encountered first in the Variant Call Format file for reference individuals.
3. *threshold::Float64=0.66* is the stopping criterion for building chromosomal regions using the informativeness for assignment statistic.
4. *maf::Float64=0.0001* is the across-population minor allele frequencies required among reference individuals for a locus to be included in the analysis.

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
