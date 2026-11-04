module LocalAncestry

# Imports/uses
using Base.Threads
using CSV
using DataFrames
using Dates
using GZip
using Sockets
using Statistics: mean

# Constants
PLOIDITY = 2
NEARZERO_FLOAT::Float64 = 0.0000000000001
HMM_STATECHANGE_PROB::Float64 = NEARZERO_FLOAT

# Include
include("library.jl")
include("interfaces.jl")
include("check.jl")
include("misc.jl")
include("predict.jl")
include("QGIO.jl/QGIO.jl")
include("print.jl")
include("evaluate.jl")
include("convert.jl")

export localancestry, read_rfmix, write_flare

end
