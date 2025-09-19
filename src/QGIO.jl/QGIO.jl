module QGIO

# Imports
using CSV
using DataFrames
using GZip
using PrettyTables

# Constants
READLINE_BUFFER_SIZE = 10000
VCF_HEADER = ["chromosome", "position", "identifier", "reference", "alternate", "quality", "filter", "format"]
NEARZERO_FLOAT = 0.000000000000001
PLOIDITY = 2
ALLELE_DICT_UINT8toINT8 = Dict{UInt8, Int8}(0x30 => 0, 0x31 => 1)

# Include
include("open.jl")
include("buffer.jl")
include("convert.jl")
include("header.jl")
include("loci.jl")
include("samples.jl")
include("read.jl")
include("readline.jl")
include("skipline.jl")
include("summary.jl")
include("print.jl")
include("haplotypes.jl")
end