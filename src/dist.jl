export Likelihood, Likelihood_fusion, Likelihood_weight, Wasserstein

abstract type AbstractHistogramDistance end
abstract type AbstractDistributionDistance end

struct Likelihood <: AbstractHistogramDistance end
struct Likelihood_fusion <: AbstractHistogramDistance
    combine_bin::Int
end


struct Wasserstein <: AbstractDistributionDistance end

function prob_dist(estimate_signal::Vector, ref_signal::Vector, histo_normalized_data::Vector, dist::AbstractHistogramDistance; kwargs...)
    dist(estimate_signal, histo_normalized_data; kwargs...)
end

function prob_dist(estimate_signal::Vector, ref_signal::Vector, histo_normalized_data::Vector, dist::AbstractDistributionDistance; kwargs...)
    dist(estimate_signal, ref_signal; kwargs...)
end

function prob_dist(estimate_signal::Vector, ref_signal::Vector, histo_normalized_data::Vector, dist::Union{Distances.UnionMetric,Distances.UnionSemiMetric}; kwargs...)
    dist(estimate_signal, ref_signal; kwargs...)
end

function (f::Likelihood)(estimate_signal::Vector, histo_normalized_data::Vector; n_data=length(histo_normalized_data), kwargs...)
    estimate_signal_ = @. max(estimate_signal, 0)
    ind = @. Int(floor(histo_normalized_data) + 1) # because the index is from 0 
    -sum(log.(estimate_signal_[ind] .+ 0.1 / n_data)) # add a small quantity to avoid log(0)
end
function (f::Likelihood_fusion)(estimate_signal::Vector, histo_normalized_data::Vector; n_data=length(histo_normalized_data), kwargs...)
    estimate_signal_ = combine_dist(max.(estimate_signal, 0), f.combine_bin)
    ind = @. Int(floor(histo_normalized_data) + 1)
    ind_fusion = combine_ind(ind, f.combine_bin)
    -sum(log.(estimate_signal_[ind_fusion] .+ 0.1 / n_data))
end

function combine_dist(dist, combine::Int)
    dist_new_p1 = sum(dist[1:combine])
    dist_new_p2 = dist[combine+1:end]
    [dist_new_p1; dist_new_p2; zeros(combine - 1)]
end

function combine_ind(ind, combine::Int)
    ind_copy = copy(ind)
    ind_copy[ind_copy.<=combine] .= 1
    ind_copy[ind_copy.>combine] .= ind[ind_copy.>combine] .- combine .+ 1
    ind_copy
end


struct Likelihood_weight <: AbstractHistogramDistance
    ξ::Float64
end
Likelihood_weight() = Likelihood_weight(1e-3)
function (f::Likelihood_weight)(estimate_signal::Vector, histo_normalized_data::Vector; kwargs...)
    estimate_signal = max.(estimate_signal, 0)
    ind = Int.(floor.(histo_normalized_data) .+ 1)
    -sum(log.(estimate_signal[ind] .+ f.ξ)) # to avoid log(0)
end
# test=[1,1,2,3,4,5,4,3,2,1]
# combine_ind(test,3)
# combine_dist(test,3)

"
 Wasserstein distance for 1d case
"
function (f::Wasserstein)(P::Array, Q::Array)
    cumP = cumsum(P)
    cumQ = cumsum(Q)
    cityblock(cumP, cumQ)
end