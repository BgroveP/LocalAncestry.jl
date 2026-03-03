
# Example using human data simulated from The 1000 Genomes project and with the HapMap3 genomic map

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