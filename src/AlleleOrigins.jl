module AlleleOrigins
using VariantCallFormat, BenchmarkTools, CSV, DataFrames, OrderedCollections

# Write your package code here.
include("assignMissing.jl")
include("haplotypeLibrary.jl")
include("haplotypeSearch.jl")
include("informativeness.jl")
include("logLikelihood.jl")
include("origins.jl")
include("misc.jl")
include("predict.jl")
include("priors.jl")
include("read.jl")
end
