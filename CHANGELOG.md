## LocalAncestry.jl v0.3.0
### Breaking Changes
- The estimating function changed from *get_local_ancestries* to *localancestry*.
- The mandatory and ordered input changed from *chromosome, referenceVCF, targetVCF, and referenceAncestries* to *referencepath, targetpath, and ancestrypath*.
- The function arguments changed to be more general, rather than related to VCF files.
- The argument for ancestries for reference individuals was changed from a DataFrame object to the path for a delimited file with similar information.
- The arguments related to block sizes have been removed: minBlockSize, incrBlockSize, blockCrit.
- The priorsMethod function has been removed because preliminary analyses showed no improvement for any of the investigated approaches.
- The minNBCProb argument has been removed, because there no longer is a cut-off for posterior probabilities for when they are used in the second step of the estimation.
- The output is no longer a tuple of posterior probabilities, labels, and the haplotype library.

### Internal changes
- Changed dependency on external libraries such as OrderedCollections and VariantCallFormat.
- Changed the internal algorithm for genetic distances such that they are based on Haldanes function and the number of loci.
- The maximum block size is fixed at 10 % of the loci.

### New Features
- It is now possible to omit haplotypes from individuals using the omitpath argument.
- The tool can omit loci based on minor allele frequencies using the maf argument to the estimating function.
- The tool now prints summary information to the console when the function is called.
- The calculation speed has improved drastically.
- The new stopping criterion *threshold* is still based on the *informativeness for assignment* 
