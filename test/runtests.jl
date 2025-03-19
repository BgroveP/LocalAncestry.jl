using AlleleOrigins
using Test
using CSV
using Tables

@testset "Reading" begin

    # Read test
    ## from VCF files
    for c in 1:2
    haplotypes, individuals = AlleleOrigins.readVCF("../data/readvcfin.vcf", c)
    control = CSV.read("../data/readvcfout" * string(c) * ".csv", Tables.matrix, types = Int8, header = false)
    @test all(haplotypes .== control)
    end


end

@testset "Priors" begin
    alphabet = 'a':'z'
    for p in 1:10
        tmp = AlleleOrigins.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1/p)) 
    end
end

