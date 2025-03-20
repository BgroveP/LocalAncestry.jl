
using AlleleOrigins, DataFrames, CSV, BenchmarkTools
 
reference_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/vcf/rotationalcattle_purebredparents_replicate2.vcf"
target_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/vcf/rotationalcattle_F1_replicate2.vcf"
origins_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/population_origins_of_individuals/rotationalcattle_purebredparents_replicate2.csv"
true_origins = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/haplotype_boa/rotational_crossbred_replicate2.csv"

# Input
chromosome = 1
referenceOrigins = CSV.read(origins_path, DataFrames.DataFrame, header=["individual", "population"], types=String)
originPriors = []
minHaploSize = 1
incHaploSize = 1
haploCrit = 0.3
ploidity = 2
minProb = 0.9


BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
@benchmark AlleleOrigins.origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)
probs, classes = AlleleOrigins.origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)

trueO = CSV.read(true_origins, DataFrame, types = Int32) |> Matrix

tI = trueO[:, 1]
trueO = trueO[:, 2:end]

# Evaluate
individuals_with_predictions = string.(keys(classes))
individuals_with_true = repeat(string.(tI), inner = 2) .* "_hap" .* repeat(string.(1:2), outer = size(tI, 1))

individuals_with_both = intersect(individuals_with_true, individuals_with_predictions)

# Print if some animals not present
if length(individuals_with_both) != length(individuals_with_predictions)
    println(string(length(individuals_with_both) - length(individuals_with_predictions)) * " individuals with predictions not in reference!")
end

out = zeros(Float16, size(individuals_with_both, 1), size(trueO, 2))
[findfirst(i .== individuals_with_true) for i in individuals_with_both]
for (i, ind) in enumerate(individuals_with_both)
    classes[ind]
end

popDict = Dict{String, Int}("holstein" => 1, "jersey" => 2, "reddairy" => 3)

[popDict[j] for j in classes[ind]]