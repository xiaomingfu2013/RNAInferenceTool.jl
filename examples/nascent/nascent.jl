using DelimitedFiles, StatsBase, DataFrames
using RNAInferenceTool
# import nascent RNA data
include("import_nascentRNA.jl")

# set inference search range and initial values
L1, L2 = [862, 2200];
τ = (L1 + L2) / 65 / 60
SRange = [(0.0, 50.0), (0.0, 50.0), (0.0, 200.0), (0.0, 0.0), (0.0, 0.0), (τ, τ)] # For delay model the parameter order: σ_off, σ_on, ρ_on, ρ_off, d, τ



# construct data structure for inference
begin
    numofRun = 12
    data_sets = [1, 2, 3, 4]
    infer_struct = Any[]
    # 1:4 are G1 stage
    for i in data_sets
        push!(infer_struct, OptimStruct(histo_data_vec_normalized[i], G1(), Likelihood(), Delay(); infer_counts=false, L1=L1, L2=L2)) # G1 stage
    end
    # 5:8 are G2 stage
    for i in 5:8
        push!(infer_struct, OptimStruct(histo_data_vec_normalized[i], G2(), Likelihood(), Delay(); infer_counts=false, L1=L1, L2=L2)) # G2 stage
    end
    # last 4 sets are merged 
    for i in 1:4
        push!(infer_struct, OptimStruct(histo_data_vec_merged_normalized[i], G1(), Likelihood(), Delay(); infer_counts=false, L1=L1, L2=L2)) # G2 stage
    end
end

infer_struct

# Start inference
result_para = []
for iter = 1:numofRun
    push!(result_para, optim_function(SRange, infer_struct[iter], MaxFuncEvals=15000, TraceMode=:compact))
end


# Present results
parameter_set = hcat([result_para[j][1] for j in 1:numofRun]...) |> transpose
model1_result = [result_para[j][2] for j in 1:numofRun]
column1 = maximum([size(model1_result[j])[1] for j in 1:numofRun])
distribution_set = zeros(column1, 2 * numofRun)
for j in 1:numofRun
    distribution_set[1:size(result_para[j][2])[1], 2*(j-1)+1:2*j] = result_para[j][2]
end
parameter_set

function into_df(parameter_set)
    df = DataFrame()
    df.sigma_off = parameter_set[:, 1]
    df.sigma_on = parameter_set[:, 2]
    df.rho_on = parameter_set[:, 3]
    df.rho_off = parameter_set[:, 4]
    df.tau = parameter_set[:, 6]
    df.effective_ini = @. df.rho_on * df.sigma_on * df.tau / (df.sigma_on + df.sigma_off)
    df.burst_size = @. df.rho_on / df.sigma_off
    df.f_on = @. df.sigma_on / (df.sigma_on + df.sigma_off)
    df.loss = parameter_set[:, 7]
    return df
end


# order 1:4 G1 5:8 G2 9:12 merged
df1 = into_df(parameter_set)