using CSV
using DataFrames
using LocalAncestry
using Test
using Base.Threads

packageroot = dirname(dirname(pathof(LocalAncestry)))
referencepath = "$(packageroot)/data/reference.vcf.gz"
targetpath = "$(packageroot)/data/target.vcf.gz"
ancestrypath = "$(packageroot)/data/globalancestries.csv"

@testset "Unit" begin

    # Splitting of vectors into chunks
    @test LocalAncestry.vecsplit(1:7, 3) == [1:3, 4:5, 6:7]
    @test LocalAncestry.vecsplit(1:7, 1) == [1:7]
    @test LocalAncestry.vecsplit(1:7, 7) == [i:i for i in 1:7]
    @test reduce(vcat, LocalAncestry.vecsplit(collect(1:100), 3)) == collect(1:100)

    # Repeated elements (renamed from StatsBase.repeat)
    @test LocalAncestry.brepeat(["a", "b"], [2, 3]) == ["a", "a", "b", "b", "b"]
    @test LocalAncestry.brepeat([1], [0]) == Int[]
    @test_throws DimensionMismatch LocalAncestry.brepeat([1], [1, 2])

    # Reduction of vectors of block assignments to contiguous ranges
    @test LocalAncestry.vec_reduce([1, 1, 2, 2, 2, 3]) == [1:2, 3:5, 6:6]
    @test LocalAncestry.vec_reduce([5, 5]) == [1:2]

    # Pretty printing of assignments
    loci = DataFrame(position = [10, 20, 30, 40])
    assignments = DataFrame(individual = ["x"], haplotype = [1], block = [UnitRange(2:3)], ancestry = ["a"])
    LocalAncestry.pretty!(assignments, "chr1", loci)
    @test assignments.chromosome == ["chr1"]
    @test assignments.basepairs == [UnitRange(20, 30)]

    # Ancestry lookup at a position
    la = DataFrame(haplotype = [1, 1, 2],
        basepairs = [UnitRange(1, 100), UnitRange(101, 200), UnitRange(1, 200)],
        ancestry = ["A", "B", "A"])
    @test LocalAncestry.get_ancestry(la, 1, 150) == ["B"]
    @test LocalAncestry.get_ancestry(la, 2, 150) == ["A"]
    @test LocalAncestry.get_ancestry(la, 1, 300) == []

    # Hidden Markov model posterior decoding
    ancestry = zeros(Int8, 3)
    probabilities = [0.9 0.1 0.9; 0.1 0.9 0.1]
    LocalAncestry.hmm_all!(ancestry, probabilities, 0.9, 0.05)
    @test ancestry == Int8[1, 2, 1]

    # Block ranges and frequencies
    refdata = Int8[0 0 1 0 1 0;
                   0 0 0 0 0 1;
                   1 0 0 0 0 0;
                   0 0 0 0 0 0;
                   0 0 1 0 1 0;
                   0 0 0 0 0 1;
                   1 0 1 0 0 0;
                   1 0 0 1 1 0;
                   0 0 0 1 0 1;
                   0 0 1 1 1 0;
                   1 1 1 1 1 0;
                   1 1 1 1 0 0;
                   1 0 1 0 1 0;
                   1 0 1 1 1 1;
                   1 1 0 0 1 0;
                   1 1 1 1 1 1;
                   0 1 1 0 1 1;
                   1 1 1 1 1 0;
                   1 1 1 1 1 0;
                   1 1 1 1 1 1;
                   0 1 1 1 1 1;
                   1 1 1 0 1 0;
                   0 1 1 1 1 1;
                   1 1 1 1 1 0;
                   1 1 1 1 1 1;
                   1 0 1 1 1 1]
    popDict = Dict("a" => 1:10, "b" => 11:20, "c" => 21:26)
    refloci = DataFrame(chromosome = fill("chr1", 6), position = collect(1:6))
    library = LocalAncestry.get_haplotype_library(refdata, popDict, 0.2, refloci)

    ## Block ranges should cover all loci contiguously
    blockstarts = sort!(first.(collect(keys(library))))
    @test blockstarts[1] == 1
    @test sort!(vcat(collect.(keys(library))...)) == collect(1:6)

    ## Each haplotype in a block has a vector of ancestry frequencies summing to one
    for (block, blocklib) in library
        for frequencies in values(blocklib)
            @test length(frequencies) == 3
            @test sum(frequencies) ≈ 1.0
        end
    end
end

@testset "Accuracy" begin
    x = localancestry(referencepath, targetpath, ancestrypath)
    y = CSV.read("$(packageroot)/data/truelocalancestries.csv", DataFrame)
    y.basepairs = UnitRange.(parse.(Int, replace.(y.basepairs, r"\:.+" => "")), parse.(Int, replace.(y.basepairs, r".+\:" => "")))
    y.ancestry = String.(y.ancestry)
    a = LocalAncestry.evaluate(x, y)
    @test sum(a.pcorrect[a.ncorrect .> 0]/sum(a.ncorrect .> 0)) > 0.88
end

@testset "Parameters" begin
    # Omit both haplotypes of the first reference individual
    refsamples = LocalAncestry.QGIO.samples(referencepath)
    omits = DataFrame(individual = [refsamples[1], refsamples[1]], haplotype = [1, 2])
    omitpath = joinpath(mktempdir(), "omits.csv")
    CSV.write(omitpath, omits)

    x = localancestry(referencepath, targetpath, ancestrypath;
        omitpath = omitpath, chromosome = "chr22", maf = 0.01, totcm = 75.0,
        printlevel = "debug")
    @test sort(names(x)) == ["ancestry", "basepairs", "chromosome", "haplotype", "individual"]
    @test all(x.chromosome .== "chr22")
    @test nrow(x) > 0
end
