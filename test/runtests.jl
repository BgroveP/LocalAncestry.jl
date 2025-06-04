using CSV
using LocalAncestry
using Tables
using Test
using DataFrames

# Full
printstyled("Testing the interface\n"; color=:blue)
ra = DataFrame(individual="individual" .* string.([1, 2, 3]), population=["a", "b", "c"])
@testset verbose = true "Interface" begin
    @testset verbose = true "big blocks" begin
        # Read
        for c in 1:2
            postProb, postClass, haplotypeLibrary = get_local_ancestries(c, "data/reference.vcf", "data/target.vcf", ra)


            ancestryControl = CSV.read("data/targetancestries.csv", DataFrame)
            for (j, i) in enumerate(ancestryControl.individual), h in 1:2
                name = i * "_hap" * string(h)

                ancestryVector = postClass[name]
                for k in ancestryVector
                    @test k == ancestryControl.population[j]
                end
            end
        end
    end
    # Read
    @testset verbose = true "small blocks" begin
        for c in 1:2
            postProb, postClass, haplotypeLibrary = get_local_ancestries(c, "data/reference.vcf", "data/target.vcf", ra, minBlockSize=1, blockCrit = 0.2)


            ancestryControl = CSV.read("data/targetancestries.csv", DataFrame)
            for (j, i) in enumerate(ancestryControl.individual), h in 1:2
                name = i * "_hap" * string(h)

                ancestryVector = postClass[name]
                for k in ancestryVector
                    @test k == ancestryControl.population[j]
                end
            end
        end
    end
end

# Read function
printstyled("\n\nTesting the read function\n"; color=:blue)
@testset verbose = true "Read .vcf files" begin
    @testset "individuals" begin
        # Read
        _, individualsVCF = LocalAncestry.readVCF("data/reference.vcf", 1)
        individualsControl = CSV.read("data/referenceindividuals.csv", DataFrame).individuals
        @test all(individualsVCF .== individualsControl)

    end

    @testset "haplotypes" begin
        for c in 1:2
            # Read
            haplotypesVCF, _ = LocalAncestry.readVCF("data/reference.vcf", c)

            # Test haplotypes
            haplotypesControl = CSV.read("data/referencehaplotypes$(c).csv", Tables.matrix, header=false)
            @test all(haplotypesVCF .== haplotypesControl)
        end
    end
end

# Make test for haplotype library search with different initial sizes and incremental sizes