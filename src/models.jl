
abstract type Model end
struct Delay<:Model end
struct DelaySync<:Model end
struct DelayComplete<:Model end
struct Telegraph<:Model end
struct TelegraphSync<:Model end
struct Poisson<:Model end
struct Bursty<:Model end

function model_selection(params,NT::Int64,::Delay)
	σ_off, σ_on, ρ_on, ρ_off, d, τ = params
	M_delay([σ_off, σ_on, ρ_on, ρ_off, d, τ], NT)
end
function model_selection(params,NT::Int64,::DelaySync)
	σ_off, σ_on, ρ_on, ρ_off, d, τ = params
	M_delay_sync([σ_off, σ_on, ρ_on, ρ_off, d, τ], NT)
end

function model_selection(params,NT::Int64,::DelayComplete)
	σ_off, σ_on, ρ_on, ρ_off, d, τ = params
	M_delay_complete([σ_off, σ_on, ρ_on, ρ_off, d, τ], NT)
end

function model_selection(params,NT::Int64,::Telegraph)
	σ_off, σ_on, ρ_on, ρ_off, d = params
	M_tele([σ_off, σ_on, ρ_on, ρ_off, d],NT)
end
function model_selection(params,NT::Int64,::TelegraphSync)
	σ_off, σ_on, ρ_on, ρ_off, d = params
	M_tele_sync([σ_off, σ_on, ρ_on, ρ_off, d],NT)
end

function model_selection(params,NT::Int64,::Poisson)
	ρ, d = params
	M_Poisson([ρ, d],NT)
end

function model_selection(params,NT::Int64,::Bursty)
	a, b, d = params
	M_bursty([a, b, d],NT)
end

# telegraph model
"""
    function M_tele(p, NT)

params: σ_off, σ_on, ρ_on, ρ_off, d
"""
function M_tele(p, NT::Int; kwargs...)
	σ_off, σ_on, ρ_on, ρ_off, d = p
    M = FSP_mat([σ_off, σ_on, ρ_on, ρ_off, d],NT; kwargs...)
    M[end,:]=ones(2*(NT+1))
    b = zeros(2*(NT+1))
    b[end] = 1
    sol = M\b
    sol[1:NT+1] + sol[NT+2:end]
end
function M_tele_sync(p, NT::Int; kwargs...)
	σ_off, σ_on, ρ_on, ρ_off, d = p
    M = FSP_mat([σ_off, σ_on, 2*ρ_on, ρ_off, d],NT; kwargs...)
    M[end,:]=ones(2*(NT+1))
    b = zeros(2*(NT+1))
    b[end] = 1
    sol = M\b
    sol[1:NT+1] + sol[NT+2:end]
end

# Poisson model
function M_Poisson(p,NT::Int)
	ρ, d = p
    λ = ρ/d
    sol = [exp(-λ)*λ^k/factorial(big(k)) for k in 0:NT]
    Float64.(sol)
end

#---
##======================================================
# telegraph model with delay
##======================================================
function FSP_mat(p, NT::Int; sparseM::Bool = true)
	# Extracting reaction rates
    σ_off, σ_on, ρ_on, ρ_off, d = p
    N = NT+1;

    IdM = Matrix(I,N,N)
    A11 = -σ_on*IdM
    A12 = σ_off*IdM
    if sparseM
        A = sparse([A11 A12; -A11 -A12])
        B0 = sparse([i==j ? -1 : i==j+1 ? 1 : 0 for i = 1:N, j=1:N])
        B11 = ρ_off*B0
        B22 = ρ_on*B0
        B = blockdiag(B11,B22)
    else
        A = [A11 A12; -A11 -A12]
        B = diagm(0=>[ρ_off*fill(-1,N);ρ_on*fill(-1,N)],-1=>[ρ_off*fill(1,N-1);0;ρ_on*fill(1,N-1)])
    end
    if d == 0
        M = A.+B
    else
        if sparseM
            C11 = d*sparse([i==j ? -i+1 : i==j-1 ? i : 0 for i = 1:N, j=1:N])
            C = blockdiag(C11, C11)
        else
            C11 = d*[i==j ? -i+1 : i==j-1 ? i : 0 for i = 1:N, j=1:N]
            C = [C11 0*I;0*I C11]
        end
        M = A.+ B.+ C    
    end
    return M
end


##======================================================
# Define FSP 

# Define FSP delay term
function delay_term(p)
    σ_off, σ_on, ρ_on, ρ_off, d, τ = p
    f0 = σ_off/(σ_off+σ_on)
    f1 = σ_on/(σ_off+σ_on)
    h0 = exp(-d*τ)*ρ_off*f0
    h1 = exp(-d*τ)*ρ_on*f1
    return [h0,h1]
end
"""
    function M_tele(p, NT)
        
params: σ_off, σ_on, ρ_on, ρ_off, d, τ
"""
function M_delay_complete(p, NT::Int; kwargs...)
	σ_off, σ_on, ρ_on, ρ_off, d, τ = p
    M = FSP_mat([σ_off, σ_on, ρ_on, ρ_off, d],NT; kwargs...)
    semigroup = exp(Matrix(M)*τ)

    sol1 = semigroup[:, 1]
    sol2 = semigroup[:, NT+2]

    sol1_dual = [sol1[1];diff(sol1)]
    sol2_dual = [sol2[1];diff(sol2)]

    a1, a2 = delay_term(p)

    h1 = a1*sol1_dual.+a2*sol2_dual
    sol = M\(-h1)
    return sol[1:NT+1]+sol[NT+2:end]
end

function M_delay_signal(p,NT::Int; L1=862,L2=2200, kwargs...)
    sol = M_delay(p, NT)
    filter_uniform = hcat(convolve_uniform([L1, L2], NT; kwargs...)...)
    convolve_filter(filter_uniform, sol)
end


function taylor_coefficients(NT::Int,at_x,gen::Function)
    Q = zeros(NT)
    taylor_gen = taylor_expand(u->gen(u),at_x,order=NT)
    for j in 1:NT
        Q[j] = taylor_gen[j-1]
    end
    return Q
end

function M_bursty(p, NT::Int)
	σ_on, b, d = p
    gen(u) = (1/(1-b*u))^(σ_on/d)
    return taylor_coefficients(NT+1,-1,gen)
end


# u0 = [(σ_off/(σ_off+σ_on));fill(0.,NT);(σ_on/(σ_off+σ_on));fill(0.,NT)]
# sol = semigroup*u0
function M_delay(p, NT::Int)
    σ_off, σ_on, _, _, _, τ = p
    M = FSP_mat(p, NT)
    semigroup = exp(Matrix(M)*τ)
    sol1 =@view semigroup[:, 1]
    sol2 =@view semigroup[:, NT+2]
    sol = (σ_off/(σ_off+σ_on))*sol1 + (σ_on/(σ_off+σ_on))*sol2
    sol[1:NT+1]+sol[NT+2:end]
end



function M_delay_sync(p, NT::Int)
    σ_off, σ_on, ρ_on, ρ_off, d, τ = p
    M = FSP_mat([σ_off, σ_on, 2*ρ_on, ρ_off, d, τ], NT)
    semigroup = exp(Matrix(M)*τ)
    sol1 =@view semigroup[:, 1]
    sol2 =@view semigroup[:, NT+2]
    sol = (σ_off/(σ_off+σ_on))*sol1 + (σ_on/(σ_off+σ_on))*sol2
    sol[1:NT+1]+sol[NT+2:end]
end