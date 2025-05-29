using CSV
using LocalAncestry
using Tables
using Test
using DataFrames

# Test of read functions
include("testobjects.jl")
@testset verbose = true "Package" begin
    @testset verbose = true "IO" begin
        @testset "vcf" begin
            for c in 1:2
                    # Read
                    haplotypesVCF, individualsVCF = LocalAncestry.readVCF("data/reference.vcf", c)
                    
                    # Test individuals
                    individualsControl = CSV.read("data/referenceindividuals.csv", DataFrame).individuals
                    @test all(individualsVCF .== individualsControl)
                    
                    # Test haplotypes
                    haplotypesControl = CSV.read("data/referencehaplotypes$(c).csv", Tables.matrix, header = false)
                    @test all(haplotypesVCF .== haplotypesControl)

                end
            end
        end
    end

    # Test of function that constructs dictionary of priors
    #populations = 'a':'z'
    #individuals = collect(1:100)
    #@testset verbose = true "Priors" begin
    #    @testset "Same" begin
    #        for (i, p) in enumerate(populations)
    #            pop_subset = collect('a':p)
    #            tmp = makePriors(pop_subset, individuals, [])
    #            combined_vector = [v for v in values(tmp[i]) for k in pop_subset for i in individuals]
    #            @test all(combined_vector .== log(1 / i))
    #        end
    #    end
    #end


    # Test of function that constructs haplotypes library
    @testset verbose = true "Library" begin
        @testset "Overall" begin
            for (i, c) in enumerate(criteria)
                x, _ = NBLA.getHaploBlocks(1, 1, c, haplotypes, popDict, 1)
                @test expectedOutputLibraryOverall[i] == x
            end
        end

        @testset verbose = true "Search" begin
            @testset "Key" begin
                for (i, c) in enumerate(criteria)
                    for k in keys(expectedOutputLibraryOverall[i])
                        if (k[1] + 1 < 6)
                            tmpkey, tmpdata = NBLA.haploSearch(1, 1, c, haplotypes, k[1], popDict)
                            @test tmpkey == k
                        end
                    end
                end
            end
            @testset verbose = true "Haplotypes" begin
                for (i, c) in enumerate(criteria)
                    for k in keys(expectedOutputLibraryOverall[i])
                        if (k[1] + 1 < 6)
                            tmpkey, tmpdata = NBLA.haploSearch(1, 1, c, haplotypes, k[1], popDict)
                            @test tmpdata == expectedOutputLibraryOverall[i][k]
                        end
                    end
                end
            end
        end
        @testset verbose = true "Informativeness" begin end
    end
end