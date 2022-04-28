using RNAInferenceTool
using Test, SafeTestsets

@time begin
    @time @safetestset "compute_distribution" begin
        include("test_compute_dist.jl")
    end
    @time @safetestset "Inference" begin
        include("test_experiment.jl")
    end 
end
