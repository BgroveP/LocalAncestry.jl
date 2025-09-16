
using CSV
using DataFrames
include("src/QGIO.jl/QGIO.jl")

path = "test/data/AdmixPop_human1_chr22_ancestral_all.vcf.gz"
header = QGIO.header(path)
samples = QGIO.samples(path)
loci = QGIO.loci(path)

omits = CSV.read("test/data/omitancestral.csv", DataFrame)
ancestries = CSV.read("test/data/AdmixPop_human1_ancestral.csv", DataFrame)

QGIO.allelefreq!(loci, path, ancestries = ancestries, omits = omits)
QGIO.inform_for_assign!(loci)

