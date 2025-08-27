module LocalAncestry
using CSV
using DataFrames
using GZip
using OrderedCollections
using Tables
using VariantCallFormat

# Constants


# Write your package code here.
include("assignMissing.jl")
include("checks.jl")
include("haplotypeLibrary.jl")
include("informativeness.jl")
include("interfaces.jl")
include("misc.jl")
include("predict.jl")
include("read.jl")

export get_local_ancestries
export get_local_ancestries2

end
