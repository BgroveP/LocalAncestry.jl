
using AlleleOrigins, DataFrames, CSV, BenchmarkTools



for r in 1:10
    println("replicate: ", r)
    for c in 1:5
        reference_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/vcf/rotationalcattle_purebredparents_replicate$r.vcf"
        target_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/vcf/rotationalcattle_F1_replicate$r.vcf"
        origins_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/population_origins_of_individuals/rotationalcattle_purebredparents_replicate$r.csv"
        true_origins = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/haplotype_boa/rotationalcattle_crossbred_replicate$r.csv"
        map_path = "/usr/home/qgg/bgrovep/projects/AdmixPop/data/genomic_map/rotationalcattle.csv"

        # Input
        chromosome = c
        referenceOrigins = CSV.read(origins_path, DataFrames.DataFrame, header=["individual", "population"], types=String)
        originPriors = []
        minHaploSize = 1
        incHaploSize = 1
        haploCrit = 0.3
        ploidity = 2
        minProb = 0.8


        probs, classes, library = AlleleOrigins.origins(chromosome, reference_path, target_path, referenceOrigins, originPriors, minHaploSize, incHaploSize, haploCrit, ploidity, minProb)
        overall, per_locus = AlleleOrigins.evaluate(classes, map_path, true_origins, library, chromosome)
        
        println("- chromosome $c: ", round(overall, sigdigits=2))

    end
end



