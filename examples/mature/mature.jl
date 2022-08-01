using DelimitedFiles, StatsBase, DataFrames
using RNAInferenceTool
# import mature RNA data
include("import_matureRNA.jl")

# set inference search range
SRange = [(0.0, 50.0), (0.0, 50.0), (0.0, 200.0), (0.0, 0.0), (1.0, 1.0)] # For telegraph model the parameter order: σ_off, σ_on, ρ_on, ρ_off, d

# construct data structure for inference
begin
    numofRun = 12
    data_sets = [1, 2, 3, 4]
    infer_struct = Any[]
    # 1:4 are G1 stage
    for i in data_sets
        push!(infer_struct, OptimStruct(histo_data_vec[i], G1(), Likelihood(), Telegraph(); infer_counts=true)) # G1 stage
    end
    # 5:8 are G2 stage
    for i in 5:8
        push!(infer_struct, OptimStruct(histo_data_vec[i], G2(), Likelihood(), Telegraph(); infer_counts=true)) # G2 stage
    end
    # last 4 sets are merged 
    for i in 1:4
        push!(infer_struct, OptimStruct(histo_data_vec_merged[i], G1(), Likelihood_merge(), Telegraph(); infer_counts=true)) # G2 stage
    end
end

infer_struct
# Start inference
result_para = []
for iter = 1:numofRun
    push!(result_para, optim_function(SRange, infer_struct[iter], MaxFuncEvals=15000, TraceMode=:silent))
end

# Present results
parameter_set = hcat([result_para[j][1] for j in 1:numofRun]...) |> transpose
model1_result = [result_para[j][2] for j in 1:numofRun]
column1 = maximum([size(model1_result[j])[1] for j in 1:numofRun])
distribution_set = zeros(column1, 2 * numofRun)
result_para[1][2]
size(result_para[1][2])[1]
for j in 1:numofRun
    distribution_set[1:size(result_para[j][2])[1], 2*(j-1)+1:2*j] = result_para[j][2]
end

function into_df(parameter_set)
    df = DataFrame()
    df.sigma_off = parameter_set[:, 1]
    df.sigma_on = parameter_set[:, 2]
    df.rho_on = parameter_set[:, 3]
    df.rho_off = parameter_set[:, 4]
    df.d = parameter_set[:, 5]
    df.effective_ini = @. df.rho_on * df.sigma_on / (df.sigma_on + df.sigma_off)
    df.burst_size = @. df.rho_on / df.sigma_off
    df.f_on = @. df.sigma_on / (df.sigma_on + df.sigma_off)
    df.loss = parameter_set[:, 6]
    return df
end

df1 = into_df(parameter_set)