using AlleleOrigins
using Test
using CSV
using Tables
using AlleleOrigins

# File contents

vcfcontents = """##fileformat=VCFv4.2
                    ##FILTER=<ID=PASS,Description="All filters passed">
                    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                    ##contig=<ID=1,length=40>
                    ##bcftools_mergeVersion=1.11+htslib-1.11
                    ##bcftools_mergeCommand=merge -o data/vcf/smallcattle_reference_all_replicate1.vcf.gz -O z data/vcf/smallcattle_reference_holstein_replicate1.vcf.gz data/vcf/smallcattle_reference_jersey_replicate1.vcf.gz data/vcf/smallcattle_reference_reddairy_replicate1.vcf.gz; Date=Fri Feb 14 12:40:10 2025
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tholstein1\tholstein2\tholstein3
                    1\t1\tlocus1\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    1\t10\tlocus2\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    1\t20\tlocus3\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    1\t30\tlocus4\tA\tG\t.\t.\t.\tGT\t1|0\t0/1\t1|1
                    1\t40\tlocus5\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t0|1
                    1\t50\tlocus6\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    1\t60\tlocus7\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    1\t70\tlocus8\tA\tG\t.\t.\t.\tGT\t0|1\t1/0\t1|0
                    1\t80\tlocus9\tA\tG\t.\t.\t.\tGT\t1|0\t1/1\t0|1
                    1\t100\tlocus10\tA\tG\t.\t.\t.\tGT\t1|1\t0/1\t1|1
                    2\t1\tlocus11\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    2\t10\tlocus12\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    2\t20\tlocus13\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    2\t30\tlocus14\tA\tG\t.\t.\t.\tGT\t1|0\t0/1\t1|1
                    2\t40\tlocus15\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t0|1
                    2\t50\tlocus16\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    2\t60\tlocus17\tA\tG\t.\t.\t.\tGT\t1|1\t1/1\t1|1
                    2\t70\tlocus18\tA\tG\t.\t.\t.\tGT\t0|1\t1/0\t1|0
                    2\t80\tlocus19\tA\tG\t.\t.\t.\tGT\t1|0\t1/1\t0|1
                    2\t100\tlocus20\tA\tG\t.\t.\t.\tGT\t1|1\t0/1\t1|1

                    """
vcfexpected = [[1 1 1 1 1 1 1 0 1 1;
        1 1 1 0 1 1 1 1 0 1;
        1 1 1 0 1 1 1 1 1 0;
        1 1 1 1 1 1 1 0 1 1;
        1 1 1 1 0 1 1 1 0 1;
        1 1 1 1 1 1 1 0 1 1],
    [1 1 1 1 1 1 1 0 1 1;
        1 1 1 0 1 1 1 1 0 1;
        1 1 1 0 1 1 1 1 1 0;
        1 1 1 1 1 1 1 0 1 1;
        1 1 1 1 0 1 1 1 0 1;
        1 1 1 1 1 1 1 0 1 1]
]

@testset verbose = true "IO" begin
    @testset "vcf" begin
        for c in 1:2
        mktempdir() do temp_dir
            temp_file_path = joinpath(temp_dir, "test.vcf")
            write(temp_file_path, vcfcontents)
                result, _ = AlleleOrigins.readVCF(temp_file_path, c)
                @test all(result .== vcfexpected[c])
            end
        end
    end
end

@testset verbose = true "Priors" begin
    @testset "Same" begin
        alphabet = 'a':'z'
        for p in 1:10
            tmp = AlleleOrigins.makePriors(alphabet[1:p], string.(collect(1:p)), [])
            combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
            @test all(combined_vector .== log(1 / p))
        end
    end
end

@testset "Library" begin
    popDict = Dict{String,Vector{Int64}}("a" => [1, 2, 3, 4], "b")
    for p in 1:10
        tmp = AlleleOrigins.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1 / p))
    end
end

@testset "Likelihood" begin
    popDict = Dict{String,Vector{Int64}}("a" => [1, 2, 3, 4], "b")
    for p in 1:10
        tmp = AlleleOrigins.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1 / p))
    end
end




@testset "Prediction" begin
    popDict = Dict{String,Vector{Int64}}("a" => [1, 2, 3, 4], "b")
    for p in 1:10
        tmp = AlleleOrigins.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1 / p))
    end
end



@testset "Refinement" begin
    popDict = Dict{String,Vector{Int64}}("a" => [1, 2, 3, 4], "b")
    for p in 1:10
        tmp = AlleleOrigins.makePriors(alphabet[1:p], string.(collect(1:p)), [])
        combined_vector = [v for v in values(tmp["1"]) for k in alphabet[1:p]]
        @test all(combined_vector .== log(1 / p))
    end
end

