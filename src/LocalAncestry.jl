module LocalAncestry
using CSV
using DataFrames
using GZip
using OrderedCollections
using Tables
using VariantCallFormat

# Constants
READLINE_BUFFER_SIZE = 10000
PLOIDITY = 2

# Write your package code here.
include("assignMissing.jl")
include("checks.jl")
include("haplotypeLibrary.jl")
include("interfaces.jl")
include("misc.jl")
include("predict.jl")
include("read.jl")

export get_local_ancestries
export get_local_ancestries2

end
