using DelaySSAToolkit # install package at https://github.com/palmtree2013/DelaySSAToolkit.jl
using JumpProcesses # add JumpProcesses
#σ_off: Gon -> Goff
#σ_on: Goff -> Gon
#ρ_on: Gon -> Gon + N, triggers N => τ 0
#ρ_off: Goff -> Goff + N, triggers N => τ 0
#d: P -> 0
# 1. Gon, 2. Goff, 3. P
function construct_prob_delaytelegraph(params)
    σ_off, σ_on, ρ_on, ρ_off, d, τ, tf = params
    rates = [σ_off, σ_on, ρ_on, ρ_off, d]
    react_stoich = [[1=>1],[2=>1],[1=>1],[2=>1],[3=>1]] 
    net_stoich = [[1=>-1,2=>1],[1=>1,2=>-1],[3=>1],[3=>1],[3=>-1]]
    mass_action_jump = MassActionJump(rates, react_stoich, net_stoich; scale_rates = false)
    jumpset = JumpSet((),(),nothing,mass_action_jump)
    delay_trigger = Dict(3=>[1=>τ],4=>[1=>τ])
    delay_complete = Dict(1=>[3=>-1])
    delayjumpset = DelayJumpSet(delay_trigger,delay_complete,Dict())
    u0 = [0,1,0] # initial condition
    de_chan0 = [[]] # initial delay channel
    tspan = (0.,tf)
    dprob = DiscreteProblem(u0, tspan)
    djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions = (false,false), save_delay_channel = true)
    return djprob
end

function delay_telegraph_SSA_signal(params; seed = nothing, L1 = 862, L = 3062)
    djprob = construct_prob_delaytelegraph(params)
    τ = params[end]
    sol = solve(djprob, SSAStepper(), seed = seed, save_delay_channel = true)    
    positions = sol.channel[end][1]
    filtered_pos = signal_function.(positions;τ = τ, L1 = L1, L = L)
    return filtered_pos
end

function signal_function(x; τ=τ , L1 = 862, L = 3062)
    return min(x/τ*L/L1,1)
end
# TEST

