
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


test = get_local_ancestries(referencepath, targetpath, ancestrypath, printlevel="debug")
truela = CSV.read("test/data/AdmixPop_human1_admixed.csv", DataFrame)
LocalAncestry.evaluate(test, truela)

truela.basepairs = UnitRange.(refloci.position[truela.first], refloci.position[truela.last])
refloci = QGIO.loci("test/data/AdmixPop_human1_chr22_ancestral_all.vcf.gz")
gmap = refloci[:, ["chromosome", "position"]]
gmap[:,"position(bp)"] = gmap[:, "position"]
function evaluate(predla, truela, gmap)

    nhaplotypes = nrow(unique(predla[:, ["individual", "haplotype"]]))
    outdict = Dict(gmap[!, "position"] .=> 1:nrow(gmap))
    out = deepcopy(gmap)
    out.ncorrect .= 0

    for i in axes(truela, 1)
        tmp = predla[(predla.individual.==truela.individual[i]).&(predla.haplotype.==truela.haplotype[i]).&(length.(intersect.([truela.basepairs[i]], predla.basepairs)).>0).&(predla.ancestry.==truela.ancestry[i]), :]

        if nrow(tmp) > 0
            tstart = gmap[searchsortedfirst(gmap[!, "position(bp)"], first(truela.basepairs[i])), "position(bp)"]
            tend = gmap[searchsortedfirst(gmap[!, "position(bp)"], last(truela.basepairs[i]))-1, "position(bp)"]
            for j in axes(tmp, 1)
                start = outdict[max(first(tmp.basepairs[j]), tstart)]
                stop = outdict[min(last(tmp.basepairs[j]), tend)]
                out.ncorrect[start:stop] .+= 1
            end
        end
    end
    out.pcorrect = out.ncorrect ./ nhaplotypes
    return out
end

LocalAncestry.mean(out[out.pcorrect .> NEARZERO_FLOAT, :pcorrect])