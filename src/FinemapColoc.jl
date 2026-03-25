module FinemapColoc

using ArgParse
using CSV
using DataFrames
using DelimitedFiles
using Statistics
using Parquet2

include("cli.jl")
include("utils.jl")
include("finemap_gtex.jl")

end
