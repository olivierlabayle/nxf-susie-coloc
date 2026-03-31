module FinemapColoc

using ArgParse
using CSV
using DataFrames
using DelimitedFiles
using Statistics
using Parquet2

include("cli.jl")
include("utils.jl")
include("prepare_gwas_results.jl")
include("finemap_gwas.jl")
include("finemap_gtex.jl")
include("aggregate_coloc.jl")

export julia_main

end
