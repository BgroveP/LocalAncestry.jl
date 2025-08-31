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
NCHUNKS = 9
NEARZERO_FLOAT::Float64 = 0.0000000000001
HMM_STATECHANGE_PROB::Float64 = 0.00001

# Init
function __init__()
    NCHUNKS = nthreads() 
end

# Write your package code here.
include("checks.jl")
include("library.jl")
include("interfaces.jl")
include("misc.jl")
include("predict.jl")
include("read.jl")

export get_local_ancestries

end
