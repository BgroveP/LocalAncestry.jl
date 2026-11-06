
# Release notes

## Version 0.3.2: November 6th, 2026

* The maximum length of a genome region is no longer fixed at 10 % of the loci; it can be set with the new `maxblockfraction` keyword argument to the `localancestry` function.
* The genetic length of the focal chromosome is no longer fixed at 100 cM; it can be set with the new `totcm` keyword argument to the `localancestry` function.
* The `read_rfmix` and `write_flare` functions are now exported and documented.
* The `printlevel="debug"` option works again because the `Sockets` standard library is now a declared dependency, and the printed package version is now read from the package metadata.
* Documentation fixes and an expanded test suite.

## Version 0.3.0: September 22nd, 2026

* Maintenance release with no changes to the public API.
* The repository is now also available at the QGG-Computational-Genomics GitHub organization.

## Version 0.2.0: March 3rd, 2026

* The *get\_local\_ancestry* function was renamed to *localancestry*, and the function arguments changed from being VCF-specific to general file paths.
* The reference ancestries are now read from a delimited file rather than a DataFrame object.
* The arguments related to block sizes were removed in favour of the informativeness-for-assignment-based stopping criterion.
* It is now possible to omit haplotypes with the `omitpath` argument and to omit loci based on minor allele frequencies with the `maf` argument.
* Package dependencies for *OrderedCollections*, *VariantCallFormat*, and *Tables* were removed in favour of the internal QGIO module.
* The internal algorithm for genetic distances changed to be based on Haldane's function and the number of loci.
* The calculation speed improved drastically.
* The tests now include accuracies when analyzing AdmixPop data simulated from The 1000 Genomes data.
* The evaluate function was updated.
* The printing under *Getting haplotype library* was changed.
