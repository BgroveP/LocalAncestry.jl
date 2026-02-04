# LocalAncestry.jl: Local Ancestry Inference using a combination of Naive Bayes Classification and Hidden Markov Model assignment

The aim of this package is to provide the users with an accurate, fast, and accessible means of inferring local ancestries. 
The package is set apart from existing local ancestry software by allowing for varying window sizes across the chromosome, and by using a combination of Naive Bayes Classification and Hidden Markov Models to assign local ancestries.

**Local Ancestry Inference:** Local Ancestry Inference is the prediction of ancestry for each combination of single nucleotide polymorphism and individual. The contrast is Global Ancestry Inference, which is the prediction of ancestral proportions at the individual level.

## Installation
The package can be installed through github:
```julia
using Pkg
Pkg.add(url = "https://github.com/BgroveP/LocalAncestry.jl", rev="0.1.0")
```

## Public API
```@meta
CurrentModule = LocalAncestry
```

```@docs
    get_local_ancestries(chromosome::Union{Int64,String}, 
                          referenceVCF::String, 
                          targetVCF::String, 
                          referenceAncestries::DataFrame; 
                          priorsMethod::String = "flat", 
                          minBlockSize::Int64 = 5, 
                          incrBlockSize::Int64 = 1, 
                          blockCrit::Float64 = 0.2, 
                          minNBCProb::Float64 = 0.95)
                          
```
