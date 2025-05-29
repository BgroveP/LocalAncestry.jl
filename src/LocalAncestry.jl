module LocalAncestry
using VariantCallFormat, CSV, DataFrames, OrderedCollections, Tables


# Write your package code here.
include("assignMissing.jl")
include("checks.jl")
include("haplotypeLibrary.jl")
include("haplotypeSearch.jl")
include("informativeness.jl")
include("logLikelihood.jl")
include("origins.jl")
include("misc.jl")
include("predict.jl")
include("priors.jl")
include("read.jl")
include("evaluate.jl")
end
