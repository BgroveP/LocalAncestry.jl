
using CSV
using DataFrames
using LocalAncestry

refpath = "test/data/AdmixPop_human1_chr22_ancestral_all.vcf"
refpathgz = "test/data/AdmixPop_human1_chr22_ancestral_all.vcf.gz"
targetpathgz = "test/data/AdmixPop_human1_chr22_admixed_all.vcf.gz"
targetpath = "test/data/AdmixPop_human1_chr22_admixed_all.vcf"
gapath = "test/data/AdmixPop_human1_ancestral.csv"
true_origins = "test/data/AdmixPop_human1_admixed.csv"
globalancestries = CSV.read(gapath, DataFrame, types=String)

# 
chromosome = 22

output = get_local_ancestries2(chromosome, refpath, targetpath, globalancestries)

probabilities, ela, library, populations, laDict, probDict = output
tla = CSV.read(true_origins, DataFrame, types=[String, Int, Int, String])
tla.accuracy .= 0.0
populations = string.(populations)
libraryregions = UnitRange.(keys(library))
nregions = length(libraryregions)
writelock = ReentrantLock()

tla.h .= 1
tla.done .= false

for i in 2:nrow(tla)
    if i == nrow(tla)
        tla.h[i] = tla.h[i-1]
    else
        if tla.first[i+1] < tla.first[i]
            tla.h[i] = (tla.h[i-1] == 1) ? 2 : 1
        else
            tla.h[i] = tla.h[i-1]
        end
    end
end

Threads.@threads for tlapos in 1:nrow(tla)

    elapos = findfirst((tla.first[tlapos] .<= first.(libraryregions)) .& (first.(libraryregions) .<= tla.last[tlapos]))
    bstart = tla.first[tlapos]
    bend = tla.last[tlapos]
    counter = 0
    h = 1
    while (length(intersect(bstart:bend, libraryregions[elapos])) > 0) 
        counter = counter + (tla.ancestry[tlapos] == populations[ela[laDict[tla.individual[tlapos]], elapos][tla.h[tlapos]]]) * length(intersect(bstart:bend, libraryregions[elapos]))
        elapos = elapos + 1
        if elapos == nregions
            break
        end
    end

    if bstart >= tla.first[tlapos+1]
        elapos = 1
    end

    @lock writelock tla.accuracy[tlapos] = counter / (bend - bstart + 1)
    @lock writelock tla.done[tlapos] = true
    @lock writelock print(tla[tlapos,:])
end


using CSV, DataFrames, Plots
data = CSV.read("C:/Users/au488376/Downloads/AdmixPop_human1.csv", DataFrame)
scatter(data.locus, data.LocalAncestry, label = "LocalAncestry")
scatter!(data.locus, data.Flare, label = "Flare")


dataf = data[data.Flare .> 0.2, :]

scatter(dataf.Flare, dataf.LocalAncestry)

using StatsBase

mean(dataf.Flare)
mean(dataf.LocalAncestry)