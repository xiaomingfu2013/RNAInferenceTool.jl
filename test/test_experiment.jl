using RNAInferenceTool
using DelimitedFiles
using Test
# histo_synthetic = readdlm("test/histo_synthetic.txt")|>vec
histo_synthetic = readdlm("histo_synthetic.txt")|>vec

# Inference
τ = 0.5
params = [1.0526,8.2034,57.989,0,0,τ]
#For delay model the parameter order: σ_off, σ_on, ρ_on, ρ_off, d, τ, fp
SRange = [(0.0,50.0),(0.0,50.0),(0.0,100.0),(0.0,0.0),(0.0,0.0),(τ,τ),]
infer_struct = OptimStruct(histo_synthetic, G1(), Likelihood(), Delay())
estimated_params, distributions = optim_function(SRange, infer_struct, MaxFuncEvals = 7000, TraceMode = :silent)

# estimated_params_ini, distributions_ini = optim_function(SRange, infer_struct, params, MaxFuncEvals = 7000, TraceMode = :compact, log_search_range = true)

estimated_params
# params
#Compare the paramters
reltol = 0.1
@test params[3] ≈ estimated_params[3] atol = reltol*params[3] 
@test params[2] ≈ estimated_params[2] atol = 2*reltol*params[2] 
@test params[1] ≈ estimated_params[1] atol = 3*reltol*params[1] 
