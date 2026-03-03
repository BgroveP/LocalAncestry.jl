# LocalAncestry.jl: Local Ancestry Inference Using Genome Regions With Variable Lengths

This Julia package provides an accurate, fast, and accessible means of inferring local ancestries. 
The package is mainly set apart from existing local ancestry software by being a region-based local ancestry inference tool that allows for regions of variying sizes across the chromosome.

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

```@autodocs
Modules = [LocalAncestry]
```

### Returns
- `x::DataFrame`: A DataFrame object with columns *individual* with String elements, *chromosome* with String elements, *haplotype* with Int elements, *basepairs* with UnitRange elements, and *ancestry* with String elements.

## Citation
Please cite the package if you use it for scientific publications. We are working on a scientific article that later can be used to refer to the package. 

## License
