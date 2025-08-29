module LocalAncestry
using CSV
using DataFrames
using GZip
using OrderedCollections
using Tables
using VariantCallFormat
using Base.Threads

# Constants
READLINE_BUFFER_SIZE = 10000
PLOIDITY = 2
NCHUNKS = nthreads()
NEARZERO_FLOAT::Float64 = 0.00000000000001
HMM_STATECHANGE_PROB::Float64 = NEARZERO_FLOAT

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
