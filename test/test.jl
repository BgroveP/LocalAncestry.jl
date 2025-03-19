
using AlleleOrigins, DataFrames, CSV, BenchmarkTools
 
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


BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
@benchmark AlleleOrigins.origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)
AlleleOrigins.origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)
