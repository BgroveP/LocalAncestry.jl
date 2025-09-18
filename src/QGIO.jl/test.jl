
using CSV
using DataFrames
using LocalAncestry
include("src/QGIO.jl/QGIO.jl")

#
chromosome = ""
referencepath = "test/data/AdmixPop_human1_chr22_ancestral_all.vcf.gz"
targetpath = "test/data/AdmixPop_human1_chr22_admixed_all.vcf.gz"
ancestrypath = "test/data/AdmixPop_human1_ancestral.csv"
omitpath = "test/data/omitancestral.csv"
threshold = 0.66
nbcprob = 0.95
maf = 0.001
printlevel = "debug"
PLOIDITY = 2
NEARZERO_FLOAT = 0.0000000001
#

test = get_local_ancestries(referencepath, targetpath, ancestrypath, omitpath = omitpath)


header = QGIO.header(path)
samples = QGIO.samples(path)
loci = QGIO.loci(path)

omits = CSV.read("test/data/omitancestral.csv", DataFrame)
ancestries = CSV.read("test/data/AdmixPop_human1_ancestral.csv", DataFrame)

QGIO.allelefreq!(loci, path, ancestries = ancestries, omits = omits)
QGIO.inform_for_assign!(loci)


refdata = QGIO.haplotypes(path)

@benchmark QGIO.haplotypes(path)
@benchmark LocalAncestry.readVCF(path, 22)
