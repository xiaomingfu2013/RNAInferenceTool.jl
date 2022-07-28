module RNAInferenceTool

using Distances, Distributions, Random, BlackBoxOptim, LinearAlgebra
using SparseArrays, Statistics, StatsBase
using DSP, UnPack

include("utils.jl")
include("optim.jl")
include("models.jl")
include("dist.jl")
include("filter.jl")

export OptimStruct, optim_function
export compute_distribution, direct_err_function, signal_distribution
export sum_with_non, convert_histo
export G1, G2, Delay, Telegraph, Poisson, DelayComplete, DelaySync, TelegraphSync

end
