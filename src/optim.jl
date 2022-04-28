
""" 
Inference data struct
1. data: default type is histogram data; another supported type is signal distribution;
2. stage: G1 or G2; where G2 type data is inferred by setting G2 = G1*G1 (convolution);
3. dist: the distance function: likelihood, likelihood_fusion, likelihood_rejection and other distance functions in Distances.jl are supported;
4. model: telegraph model, delay telegraph model, bursty model and Poisson model are supported;
5. ishisto_data: if the data is of histogram type (true) esle dist type (false)
"""
struct OptimStruct{S,D,M}
    data::Vector
    stage::S
    dist::D
    model::M
    ishisto_data::Bool
end
function OptimStruct(data::Vector,stage::S,dist::D,model::M; ishisto_data::Bool=true) where {S,D,M}
    return OptimStruct(data, stage, dist, model, ishisto_data)
end

function inference_setup(data, model, dist, SRange; stage = G1(), ishisto_data::Bool = true,  kwargs...)
    optim_struct = OptimStruct(data, stage, dist, model, ishisto_data)
    optim_function(SRange, optim_struct; kwargs...)
end

struct OptimStructWrapper{S,D,M,F,FP,SR,EF}
    data::Vector
    stage::S
    dist::D
    model::M
    ishisto_data::Bool
	filter::F
	infer_counts::Bool
	fp::FP
	log_search_range::Bool
	SRange::SR
	err_func::EF
	reference_data::Vector
end

"""
	function optim_function(SRange, optim_struct; kwargs...)

SRange: search range
optim_struct: inferred data structure
important kwargs:
1. infer_counts: use only when the values are integer.
2. MaxFuncEvals: max runs
3. log_search_range: use when the search range is very large
4. TraceMode: :silent, :compact
"""
function optim_function(SRange, optim_struct::OptimStruct, args...; infer_counts::Bool=false, mom_match::Bool=false, log_search_range::Bool=false, L1 = 862, L2 = 2200, falsepositive = 0, NT = nothing, kwargs...)
	fp = falsepositive

	filter_uniform, err_func, reference_data, NT = ini_collection(optim_struct, falsepositive, infer_counts, mom_match, L1, L2, NT;kwargs...) 

	optim_struct_wrapper = OptimStructWrapper(optim_struct.data, optim_struct.stage, optim_struct.dist, optim_struct.model, optim_struct.ishisto_data, filter_uniform, infer_counts, fp, log_search_range, SRange, err_func, reference_data)

	thetax = start_BBO_optim(optim_struct_wrapper, args...; kwargs...)

	model, stage = optim_struct.model, optim_struct.stage
	params = thetax[1:end-1]
	estimate_data = compute_distribution(params, NT, model, stage, infer_counts, filter = filter_uniform)
    return thetax, hcat(estimate_data, reference_data)
end

function ini_collection(optim_struct::OptimStruct, falsepositive, infer_counts, mom_match, L1, L2, NT; kwargs...)
	fp = falsepositive
	if optim_struct.ishisto_data # if ssa_data or distribution data
		data = optim_struct.data
    	plotrange = Int(ceil(maximum(data)) + 10)
		err_func = err_func_histo 
		reference_data = signal_transform(data,fp)
	else
		@assert !(typeof(optim_struct.dist)<:AbstractHistogramDistance) "Do not support likelihood function as a distance function when the data is of distribution type!"
		@assert mom_match!=true "Do not support moment matching with distribution type!"
		plotrange = length(optim_struct.data) + 10
		if plotrange >= 200
			error("The plot range is too large, thus will greatly slow down the optimization!")
		end
		err_func = err_func_dist
		println("Variable [false positive] is not taken into account")
		reference_data = optim_struct.data
	end
	
	if NT === nothing 
		NT = copy(plotrange)
	else
		NT = max(plotrange, NT)
	end

	if infer_counts  # check if infer count data or signal data
    	filter_uniform = nothing;
	else
		filter_uniform = hcat(convolve_uniform([L1, L2], NT; kwargs...)...);
	end
	l_ref = length(reference_data)
	if l_ref < NT
		append!(reference_data, fill(0, NT-l_ref))
	end
	return filter_uniform, err_func, reference_data, NT
end

function start_BBO_optim(optim_struct_wrapper::OptimStructWrapper, args...; MaxFuncEvals::Int=20000, TraceMode=:compact, Method = :adaptive_de_rand_1_bin_radiuslimited, kwargs...)
	@unpack SRange, err_func, log_search_range = optim_struct_wrapper 
	if log_search_range # if use log search range, recommend when the search range is very large 
		SRange_log = [(log(SRange[i][1]+1),log(SRange[i][2]+1)) for i in eachindex(SRange)]
		opts = bbsetup(param->err_func(exp.(param) .- 1, optim_struct_wrapper; kwargs...);
        Method = Method,
        SearchRange = SRange_log, TraceMode = TraceMode, 
		MaxFuncEvals = MaxFuncEvals, kwargs...)
		if isempty(args)
			res = bboptimize(opts)
		else
			x0 = log.(args[1] .+ 1) 
			res = bboptimize(opts, x0)
		end
		thetax = [exp.(best_candidate(res)).-1; best_fitness(res)]
	else
		opts = bbsetup(param->err_func(param, optim_struct_wrapper; kwargs...);
        Method = Method,
        SearchRange = SRange, TraceMode = TraceMode, 
		MaxFuncEvals = MaxFuncEvals)
		res = bboptimize(opts, args...)
	    thetax = [best_candidate(res); best_fitness(res)]
	end
	return thetax
end

abstract type Stage end
struct G1<:Stage end 
struct G2<:Stage end 

function check_stage(::G1, estimate_data, N::Int64)
	estimate_data[1:N]
end
function check_stage(::G2, estimate_data, N::Int64)
	DSP.conv(estimate_data, estimate_data)[1:N]
end

function signal_transform(histo_normalized_data::Vector,falsepostive)
	temp = falsepostive > 0 ? filter(x->x.>=falsepostive,histo_normalized_data) : histo_normalized_data
	reference_data = convert_histo(temp)[2]
	return reference_data
end

function err_func_histo(params, optim_struct_wrapper::OptimStructWrapper; order_up::Int=4,  mom_match::Bool=false, NT = length(optim_struct_wrapper.reference_data), kwargs...)
	@unpack data, stage, dist, model, filter, infer_counts, fp, reference_data = optim_struct_wrapper
	histo_normalized_data = data
	filter_uniform = filter

	estimate_data = compute_distribution(params, NT, model, stage, infer_counts, filter = filter_uniform)
	if fp > 0 && !(typeof(dist)<:AbstractHistogramDistance)
		update_fp!(fp, estimate_data, reference_data)
	end
	error = prob_dist(estimate_data,reference_data,histo_normalized_data,dist; kwargs...)
	if mom_match # if add extra weigh with moment information
		sample_moments_list = [sample_moments(histo_normalized_data, i) for i in 1:order_up]
		prob_moments_list = [prob_moments(estimate_data, i) for i in 1:order_up]
		error_mom = sum(abs2, sample_moments_list.-prob_moments_list)
		error =  error + error_mom
	end    
	return error
end

function convolve_filter(filter::Matrix, count_data::Vector)
	NT = length(count_data)-1
	trim_filter =@view filter[1:NT,1:NT]
	signalv1 = zeros(NT)
	md_ =@view count_data[2:NT+1]
	mul!(signalv1,trim_filter,md_)
	[count_data[1]+signalv1[1]; signalv1[2:end]]
end


function signal_distribution(params, NT::Int, model::M; filter = nothing, L1=862, L2=2200, kwargs...) where {M}
	count_data = model_selection(params,NT,model)
	plotrange = length(count_data)+1
	if filter === nothing
		filter =  hcat(convolve_uniform([L1, L2], plotrange; kwargs...)...)
	end
	convolve_filter(filter, count_data)
end

function compute_distribution(params, NT::Int, model::M, stage::S, infer_counts::Bool; filter = nothing, L1 = 862, L2 = 2200, kwargs...) where {M, S}
    if !infer_counts
		est_data = signal_distribution(params, NT, model; filter = filter, L1 = L1, L2=L2, kwargs...)
    else
        est_data = model_selection(params, NT, model)
    end
    return check_stage(stage, est_data, NT)
end


function update_fp!(fp, estimate_data, reference_data)
	cut_off = Int(floor(fp+1))
	estimate_data =@view estimate_data[cut_off:end]
	reference_data =@view reference_data[cut_off:end]
	estimate_data = max.(estimate_data,0);
	reference_data = max.(reference_data,0);
end

function err_func_dist(params, optim_struct_wrapper::OptimStructWrapper; NT = length(reference_data), kwargs...)
	@unpack data, stage, dist, model, filter, infer_counts, fp = optim_struct_wrapper
	reference_data = data
	filter_uniform = filter
	@assert !(typeof(dist)<:AbstractHistogramDistance) 
	estimate_data = compute_distribution(params, NT, model, stage, infer_counts, filter = filter_uniform)
	P = max.(estimate_data,0); Q = max.(reference_data,0);
	if fp >0
		cut_off = Int(floor(fp+1))
		P =@view estimate_data[cut_off:end]
		Q =@view reference_data[cut_off:end]
	end
	error = prob_dist(P, Q, [], dist; kwargs...)
    return error
end

function direct_err_function(params, optim_struct::OptimStruct; infer_counts::Bool=false, L1 = 862, L2 = 2200, falsepositive = 0, log_search_range = false, SRange = nothing, NT = nothing, mom_match::Bool = false, kwargs...)
	
	filter_uniform, err_func, reference_data, NT = ini_collection(optim_struct, falsepositive, infer_counts, mom_match, L1, L2, NT; kwargs...) 
	
	optim_struct_wrapper = OptimStructWrapper(optim_struct.data, optim_struct.stage, optim_struct.dist, optim_struct.model, optim_struct.ishisto_data, filter_uniform, infer_counts, falsepositive, log_search_range, SRange, err_func, reference_data)

	return err_func(params, optim_struct_wrapper; kwargs...)
end