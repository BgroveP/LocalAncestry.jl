using ARV
using CSV
using Tables
using Test

# Test of read functions such as readVCF 
include("test/testobjects.jl")

@testset verbose = true "IO" begin
    @testset "vcf" begin
        for c in 1:2
            mktempdir() do temp_dir
                temp_file_path = joinpath(temp_dir, "test.vcf")
                write(temp_file_path, referencevcf)
                result, _ = ARV.readVCF(temp_file_path, c)
                @test all(result .== referencehaplotypes[c])
            end
        end
    end
end

# Test of function that constructs dictionary of priors
populations = 'a':'z'
individuals = collect(1:100)
@testset verbose = true "Priors" begin
    @testset "Same" begin
        for (i, p) in enumerate(populations)
            pop_subset = collect('a':p)
            tmp = ARV.makePriors(pop_subset, individuals, [])
            combined_vector = [v for v in values(tmp[i]) for k in pop_subset for i in individuals]
            @test all(combined_vector .== log(1 / i))
        end
    end
end


# Test of function that constructs haplotypes library
@testset verbose = true "Library" begin
    @testset "Overall" begin
        for (i, c) in enumerate(criteria)
            @test expectedOutputLibraryOverall[i] == ARV.getHaploBlocks(1, 1, c, haplotypes, popDict, 1)
        end
    end

    @testset verbose = true "Search" begin
        @testset "Key" begin
            for (i, c) in enumerate(criteria)
                for k in keys(expectedOutputLibraryOverall[i])
                    if (k[1] + 1 < 6)
                    tmpkey, tmpdata = ARV.haploSearch(1,1, c, haplotypes, k[1], popDict);
                   @test tmpkey == k
                end
            end
        end 
    end
    @testset verbose = true "Haplotypes" begin
        for (i, c) in enumerate(criteria)
            for k in keys(expectedOutputLibraryOverall[i])
                if (k[1] + 1 < 6)
                tmpkey, tmpdata = ARV.haploSearch(1,1, c, haplotypes, k[1], popDict);
                @test tmpdata == expectedOutputLibraryOverall[i][k]
            end
        end
    end
end
end
@testset verbose = true "Informativeness" begin
end
end


# Test function that compur
@testset "Likelihood" begin
    popDict = Dict{String,Vector{Int64}}("a" => [1, 2, 3, 4], "b")
    for p in 1:10
        tmp = ARV.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1 / p))
    end
end




@testset "Prediction" begin
    popDict = Dict{String,Vector{Int64}}("a" => [1, 2, 3, 4], "b")
    for p in 1:10
        tmp = ARV.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1 / p))
    end
end



@testset "Refinement" begin
    popDict = Dict{String,Vector{Int64}}("a" => [1, 2, 3, 4], "b")
    for p in 1:10
        tmp = ARV.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1 / p))
    end
end

