
function sample_moments(array, order::Int64)
    sum(array .^ order) / length(array)
end

function prob_moments(prob_array, order::Int64; start_num::Int64=0)
    N = length(prob_array)
    sum((start_num:start_num+N-1) .^ order .* prob_array)
end

function Err_func_moment_matching(histo_data, params, model; order_up::Int64=3)
    _, reference_data = Signal_data(histo_data)
    NT = length(reference_data)
    molecule_distribution = model_selection(params, NT, model)
    sample_moments_list = [sample_moments(histo_data, i) for i in 1:order_up]
    prob_moments_list = [prob_moments(molecule_distribution, i) for i in 1:order_up]
    error = sum(abs2, sample_moments_list .- prob_moments_list)
    return error, molecule_distribution
end

function optim_moment_matching(SRange, histo_data, model; order_up::Int=3, MaxFuncEvals::Int=20000, MaxTime=false, TraceMode=:compact)
    opts = bbsetup(params -> Err_func_moment_matching(histo_data, params, model; order_up=order_up)[1];
        Method=:adaptive_de_rand_1_bin_radiuslimited,
        SearchRange=SRange, TraceMode=TraceMode,
        MaxTime=MaxTime, MaxFuncEvals=MaxFuncEvals)
    res = bboptimize(opts)
    thetax = [best_candidate(res); best_fitness(res)]
    return thetax, Err_func_moment_matching(histo_data, best_candidate(res), model; order_up=order_up)[2]
end

