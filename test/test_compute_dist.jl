using Test, RNAInferenceTool

τ = 0.5
params = [1.0526,8.2034,57.989,0,0,τ]
NT = 100
L1 = 862; L2 = 2200;

model_list = [Telegraph(),Delay()]
stage_list = [G1(), G2()]
filter_uniform = hcat(RNAInferenceTool.convolve_uniform([L1, L2], NT)...)

@testset for infer_counts in [true, false], model in model_list, stage in stage_list
    molecule_distribution = RNAInferenceTool.model_selection(params, NT, model) # select which model to use
    if infer_counts
        estimate_data = molecule_distribution[1:NT]
    else
        estimate_data = RNAInferenceTool.convolve_filter(filter_uniform,molecule_distribution) # pass the signal filter such that count -> signal 
    end
    dist0 = RNAInferenceTool.check_stage(stage, estimate_data, NT)
    @time dist1 = RNAInferenceTool.compute_distribution(params, NT, model, stage, infer_counts)
    @time dist2 = RNAInferenceTool.compute_distribution(params, NT, model, stage, infer_counts, filter = filter_uniform)

    @info infer_counts, model, stage
    @test sum(abs2, dist1 .- dist0) < 1e-10 
    @test sum(abs2, dist2 .- dist0) < 1e-10 
end

# RNAInferenceTool.compute_distribution(params, NT, model, stage, infer_counts, filter = filter_uniform)

# infer_counts = true
# model = Delay()
# stage = G2()

# molecule_distribution = RNAInferenceTool.model_selection(params, NT, model) # select which model to use
# if infer_counts
#     estimate_data = molecule_distribution[1:NT]
# else
#     estimate_data = RNAInferenceTool.convolve_filter(filter_uniform,molecule_distribution) # pass the signal filter such that count -> signal 
# end
# dist0 = RNAInferenceTool.check_stage!(stage, estimate_data, NT)
# dist1 = RNAInferenceTool.compute_distribution(params, NT, model, stage, infer_counts)
# # println(infer_counts, model, stage)
# using Plots
# plot(dist0)
# plot!(dist1)
# @test sum(abs2, dist1 .- dist0) < 1e-10 
