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
localancestry(
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

## Citation
Please cite the package if you use it for scientific publications. We are working on a scientific article that later can be used to refer to the package. 

## Example using human data simulated from The 1000 Genomes project and with the HapMap3 genomic map

There is only one function in LocalAncestry.jl that is exposed for public use. This function is the *localancestry* function, and it is a wrapper function that reads the reference information, creates a library of haplotype blocks, and infers local ancestries for target individuals using the library and a Hidden Markov model. This function can be used on the data that is provided with the LocalAncestry.jl package:

```julia
using LocalAncestry

packageroot = dirname(dirname(pathof(LocalAncestry)))
x = localancestry("$(packageroot)/data/reference.vcf.gz", 
                  "$(packageroot)/data/target.vcf.gz", 
                  "$(packageroot)/data/globalancestries.csv")

```
which returns:
```
LocalAncestry.jl v1.0.0
Started at 12:18:11 on Mar 3, 2026

Files
   Reference        .../data/reference.vcf.gz
   Target           .../data/target.vcf.gz
   Ancestry         .../data/globalancestries.csv

Parameters
   Chromosome
   MAF threshold    0.0001
   IA threshold     0.66
   N threads        1

Reading the locus information
   Chromosomes        chr22
   Loci               20129
   Loci w. GT field   20129

Reading the reference ancestries
   Populations     AFR, AFRanc, EAS, EASanc, EUR, EURanc, admixed
   N individuals   10000, 384, 10000, 411, 10000, 323, 10000

Number of haplotypes
 population  N      unique  omit  use   
────────────────────────────────────────
 EUR          4000    4000     0   4000
 EAS          4000    4000     0   4000
 AFR          4000    4000     0   4000
 Total       12000   12000     0  12000

Subsetting loci
   Loci on focal chromosome:   20129/20129
   Loci with MAF above limit:  10621/20129
   Loci for later analysis:    10621/20129

Reading the reference haplotypes

Mapping reference ancestries

Getting haplotype library
   Block ranges
   Block frequencies
   Haplotype blocks: 180

Reading the target haplotypes

Assigning local ancestries
```

To understand what LocalAncestry.jl does, we will explain the output in segments. First, some standard information and the input is printed to the screen for later reference:

```
LocalAncestry.jl v0.9.0
Started at 12:18:11 on Mar 3, 2026

Files
   Reference        .../LocalAncestry.jl/data/reference.vcf.gz
   Target           .../data/target.vcf.gz
   Ancestry         .../data/globalancestries.csv

Parameters
   Chromosome
   MAF threshold    0.0001
   IA threshold     0.66
   N threads        1
```

This output contains the version of the package, when the analysis was started, the paths to input files, and parameters used for the analysis. In the above, the parameters are chromosome "" which states that the focal chromosome is the first encountered chromosome of the reference VCF file, the average minor allele frequency across ancestral population for inclusion into the analysis is 0.0001, the informativeness for assignment of genome regions should be larger than 0.66, and the number of threads used is one. A higher value for the IA threshold results in larger but fewer genome regions across the chromosome.

The next step is to get the locus information from the reference VCF file. That is, no genotypes are read yet. 
```
Reading the locus information
   Chromosomes        chr22
   Loci               20129
   Loci w. GT field   20129
```
This states that the reference VCF file contains information on chromosome 22, there are 20129 loci, and all of these loci has a genotype field. This is important because only the genotype field is used for the local ancestry inference.

The next step is to read the global ancestries from the global ancestry file:
```
Reading the reference ancestries
   Populations     AFR, AFRanc, EAS, EASanc, EUR, EURanc, admixed
   N individuals   10000, 384, 10000, 411, 10000, 323, 10000
```
Where the output states that we found 7 different ancestries in our global ancestry file, most of these ancestries have 10K individuals, while others have fewer. By combining this information with the samples in the vcf file, the *localancestry* function identifies 4K haplotypes (2K individuals) for the three populations EUR, EAS, and AFR in the reference VCF file.
```
Number of haplotypes
 population  N      unique  omit  use   
────────────────────────────────────────
 EUR          4000    4000     0   4000
 EAS          4000    4000     0   4000
 AFR          4000    4000     0   4000
 Total       12000   12000     0  12000
```
These are the reference individuals that will be used for further analysis. Optionally, individual haplotypes from reference individuals could be omitted by supplying a path to the omit parameter.

The next step is to determine which loci in the VCF file that will be used for further analysis. Currently, the only omission factors for loci are whether they are on the focal chromosome and whether their minor allele frequencies is above the lower limit. 
```
Subsetting loci
   Loci on focal chromosome:   20129/20129
   Loci with MAF above limit:  10621/20129
   Loci for later analysis:    10621/20129
```
In the output above, it is shown that the focal chromosome had genotype information for 20129 loci of which 10621 had minor allele frequencies above the lower limit. Now we know the number of loci and the number of haplotypes, and we can read the reference haplotypes:
```
Reading the reference haplotypes

Mapping reference ancestries
```
and construct the library of haplotype blocks:
```
Getting haplotype library
   Block ranges
   Block frequencies
   Genome regions: 180
```
which states that the library contains haplotype blocks in 180 genome regions. Now the target haplotypes are read and analyzed:
```
Reading the target haplotypes

Assigning local ancestries
```

The output is a DataFrame object with 5 columns:
```
18161×5 DataFrame
   Row │ individual      chromosome  haplotype  basepairs          ancestry 
       │ String          String      Int64      UnitRang…          String
───────┼────────────────────────────────────────────────────────────────────
     1 │ admixed_126793  chr22               1  16396408:17217190  EAS
     2 │ admixed_126793  chr22               1  17221888:17382205  EUR
     3 │ admixed_126793  chr22               1  17382326:18144133  EAS
   ⋮   │       ⋮             ⋮           ⋮              ⋮             ⋮
 18160 │ admixed_138779  chr22               2  31133786:44265564  EUR
 18161 │ admixed_138779  chr22               2  44270505:50780578  EAS
```
where elements in the *individual* column are the names of individuals as stated in the VCF file, the elements of the *chromosome* column is the focal chromosome, the elements of the *haplotype* column are the haplotype number (LocalAncestry currently only works for diploid organisms), elements of the *basepairs* column are the basepair range, and the elements in the *ancestry* are the predicted ancestries. This output format was chosen because it compressed the output to some extend (block-based output rather than locus-based), while preserving human readability. 

