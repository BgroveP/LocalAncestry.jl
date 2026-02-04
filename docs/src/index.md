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

```@meta
CurrentModule = LocalAncestry
```

```@docs
get_local_ancestries(
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestrypath::String;
    omitpath::String="",
    chromosome::Union{Int,AbstractString}="",
    threshold::Float64=0.66,
    maf::Float64=0.0001,
    printlevel::String="standard"
    )
```

### Returns
- `x::DataFrame`: A DataFrame object with columns ....

