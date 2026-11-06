# Changelog

All notable changes to this project are documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.2] - 2026-11-06

### Added
- The maximum length of a genome region can now be set with the new `maxblockfraction` keyword argument to `localancestry` (default 0.1, i.e. 10% of the loci on the focal chromosome). It was previously fixed by the internal constant `MAX_BLOCK_FRACTION`, which has been removed.
- The genetic length of the focal chromosome can now be set with the new `totcm` keyword argument to `localancestry` (default 100.0 centiMorgans). It is used for the state-transition probabilities of the hidden Markov model and was previously fixed at 100 cM.
- `read_rfmix` and `write_flare` are now exported and documented. They convert between the output of `localancestry` and the RFMix and FLARE formats.

### Fixed
- `localancestry` no longer errors when run with `printlevel="debug"`: the `Sockets` standard library is now a declared dependency.
- The printed package version is now read from the package metadata instead of a hard-coded string.
- Documentation: the badges now point to the LocalAncestry.jl repository, stale references to the removed `get_local_ancestries` function were replaced with `localancestry`, the output columns are documented, and the `totcm` and `printlevel` parameters are described.

### Changed
- Replaced the package-internal `mean` function with the one from the `Statistics` standard library.
- Renamed the internal `breapeat` function to `brepeat` and removed unused internal functions.
- The test suite now covers the internal helper functions, the hidden Markov model, the haplotype library, and the optional arguments of `localancestry`.

## [0.3.0] - 2026-09-22

Maintenance release with no changes to the public API.

- The repository is now also available at the QGG-Computational-Genomics GitHub organization.
- The continuous-integration and documentation workflows were updated.

## [0.2.0] - 2026-03-03

### Breaking changes
- The estimating function changed from `get_local_ancestries` to `localancestry`.
- The mandatory and ordered input changed from `chromosome`, `referenceVCF`, `targetVCF`, and `referenceAncestries` to `referencepath`, `targetpath`, and `ancestrypath`.
- The function arguments changed to be more general, rather than related to VCF files.
- The argument for ancestries for reference individuals changed from a DataFrame object to the path for a delimited file with similar information.
- The arguments related to block sizes were removed: `minBlockSize`, `incrBlockSize`, `blockCrit`.
- The `priorsMethod` function was removed because preliminary analyses showed no improvement for any of the investigated approaches.
- The `minNBCProb` argument was removed because there no longer is a cut-off for posterior probabilities for when they are used in the second step of the estimation.
- The output is no longer a tuple of posterior probabilities, labels, and the haplotype library; it is a DataFrame with the estimated local ancestries.

### Added
- It is now possible to omit haplotypes from reference individuals using the `omitpath` argument.
- The tool can omit loci based on minor allele frequencies using the `maf` argument to the estimating function.
- The tool now prints summary information to the console when the function is called.

### Changed
- Replaced the dependencies on external libraries such as OrderedCollections and VariantCallFormat with the internal QGIO module.
- Changed the internal algorithm for genetic distances such that they are based on Haldane's function and the number of loci.
- The maximum block size is fixed at 10 % of the loci.
- The stopping criterion for building genome regions is the informativeness for assignment statistic: a region is extended until the statistic drops below the `threshold` argument (default 0.66).
- The calculation speed has improved drastically.

[Unreleased]: https://github.com/BgroveP/LocalAncestry.jl/compare/v0.3.2...HEAD
[0.3.2]: https://github.com/BgroveP/LocalAncestry.jl/compare/v0.3.0...v0.3.2
[0.3.0]: https://github.com/BgroveP/LocalAncestry.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/BgroveP/LocalAncestry.jl/compare/v0.1.0...v0.2.0
