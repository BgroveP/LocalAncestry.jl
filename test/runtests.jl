using CSV
using LocalAncestry
using Test
using DataFrames
using Base.Threads

# Full
packageroot = dirname(dirname(pathof(LocalAncestry)))

@testset "Accuracy" begin
    x = localancestry("$(packageroot)/data/reference.vcf.gz", 
    "$(packageroot)/data/target.vcf.gz", 
    "$(packageroot)/data/globalancestries.csv")
    y = CSV.read("$(packageroot)/data/truelocalancestries.csv", DataFrame)
    y.basepairs = UnitRange.(parse.(Int, replace.(y.basepairs, r"\:.+" => "")), parse.(Int, replace.(y.basepairs, r".+\:" => "")))
    y.ancestry = String.(y.ancestry)
    a = LocalAncestry.evaluate(x, y)
    @test sum(a.pcorrect[a.ncorrect .> 0]/sum(a.ncorrect .> 0)) > 0.88
end


@testset "Time" begin
    dataset = "hapmap3"
    _, time, _, _, _ = @timed localancestry("data/current_$(dataset)_2000_1.vcf.gz", "data/target_$(dataset)_2000_1.vcf.gz", "data/ga_own_1.csv")
    @test time < 30.0
end


@testset "RAM" begin
    dataset = "hapmap3"
    _, _, bytes, _, _ = @timed localancestry("data/current_$(dataset)_2000_1.vcf.gz", "data/target_$(dataset)_2000_1.vcf.gz", "data/ga_own_1.csv")
    @test bytes < 3*10^9
end
