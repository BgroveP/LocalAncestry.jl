module LocalAncestry

# Imports/uses
using Base.Threads
using CSV
using DataFrames
using Dates
using GZip
using OrderedCollections
using Tables
using VariantCallFormat

# Constants
READLINE_BUFFER_SIZE = 10000
PLOIDITY = 2
NEARZERO_FLOAT::Float64 = 0.0000000000001
HMM_STATECHANGE_PROB::Float64 = NEARZERO_FLOAT
MAX_BLOCK_FRACTION::Float64 = 0.1

# Write your package code here.
include("library.jl")
include("interfaces.jl")
include("misc.jl")
include("predict.jl")
include("QGIO.jl/QGIO.jl")
include("print.jl")
include("evaluate.jl")

export get_local_ancestries

end
