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
    ...)
```

### Optional parameters
The optional parameters for the get_local_ancestries function are:
1. *omitpath::String=""* which is an absolute or relative path to a delimited file with information that denotes the haplotypes that should be excluded from the reference panel. The delimited file contains two columns: *individual* and *haplotype*.
2. *chromosome::Union{Int,AbstractString}=""* which is the focal chromosome for the local ancestry inference. The default behaviour is to analyze the chromosome that is encountered first in the Variant Call Format file for reference individuals.
3. *threshold::Float64=0.66* is the stopping criterion for building chromosomal regions using the informativeness for assignment statistic.
4. *maf::Float64=0.0001* is the across-population minor allele frequencies required among reference individuals for a locus to be included in the analysis.

### Returns
- `x::DataFrame`: A DataFrame object with columns ....


