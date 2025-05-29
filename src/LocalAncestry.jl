module LocalAncestry
using CSV 
using DataFrames 
using OrderedCollections 
using Tables
using VariantCallFormat

# Write your package code here.
include("assignMissing.jl")
include("checks.jl")
include("haplotypeLibrary.jl")
include("haplotypeSearch.jl")
include("informativeness.jl")
include("logLikelihood.jl")
include("interfaces.jl")
include("misc.jl")
include("predict.jl")
include("priors.jl")
include("read.jl")
include("evaluate.jl")

export getLocalAncestries

end
