
using VariantCallFormat, BenchmarkTools, CSV, DataFrames, OrderedCollections

# 
reference_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/vcf/rotationalcattle_any_purebredparents_replicate2.vcf"
target_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/vcf/rotationalcattle_any_F1_replicate2.vcf"
origins_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/population_origins_of_individuals/rotationalcattle_purebredparents_replicate2.csv"

# Input
chromosome = 1
referenceOrigins = CSV.read(origins_path, DataFrames.DataFrame, header=["individual", "population"], types=String)
originPriors = []
minHaploSize = 1
incHaploSize = 1
haploCrit = 0.3
ploidity = 2
minProb = 0.9

function origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)
    # Function

    ## Derived information
    populations = unique(referenceOrigins[!, "population"])

    ## Data
    referenceData, referenceIndividuals = readVCF(reference_path, chromosome)
    targetData, targetIndividuals = readVCF(target_path, chromosome)
    referenceOriginsVector = haplotypeOrigins(referenceIndividuals, referenceOrigins)

    # reference populations 
    popDict = Dict{String,Vector{Int64}}()
    for i in unique(referenceOriginsVector)
        popDict[i] = findall(i .== referenceOriginsVector)
    end

    ## Priors
    priorProb = makePriors(populations, targetIndividuals, originPriors)

    ## Haplotype library
    haplotypeLibrary = getHaploBlocks(minHaploSize, incHaploSize, haploCrit, referenceData, popDict, 1)
    nHaplotypeBlocks = length(haplotypeLibrary)

    #Store, donor haplo is free of ID, already removed above
    postProb = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [OrderedDict(populations .=> [Vector{Union{Missing,Float64}}(missing, length(haplotypeLibrary)) for i in 1:length(populations)]) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    postClass = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for l in 1:length(targetIndividuals) for h in 1:ploidity]))


    # Get log-likelihoods
    LL = OrderedDict()
    for (region, Haplo) in haplotypeLibrary
        LL[region] = getLL(region, Haplo, referenceData, popDict)
    end

    # predict
    for (i, id) in enumerate(targetIndividuals)
        for h in 1:ploidity
            idname = id * "_hap" * string(h)

            postProb[idname], postClass[idname] = predInd(priorProb[id], targetData[2*i+h-2, :], LL, breeds, nHaplotypeBlocks, minProb)
            # Assign missing
            postClass[idname] = refineBoA(postClass[idname], postProb[idname])
        end
    end

    return postProb, postClass
end

@benchmark origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb);


