function density_func_uniform(z, Δ::Float64, L1, L2)
    density1(z) = z>=0&&z<=1 ? 1 : 0
    density2(z) = z>1&&z<=1+Δ ? 1/Δ : 0
    L1/(L1+L2)*density1(z) .+ L2/(L1+L2)*density2(z)
end

function convolve_uniform(param, plotrange::Int; Δ=0.01, kwargs...)
    L1, L2  = param
    num_of_conv = plotrange
    Δx = Δ/10
    x = Δx:Δx:plotrange
    list_conv = density_func_uniform.(x, Δ, L1, L2)
    list_conv_save = Array{Any,1}(undef,num_of_conv)
    list_conv_save[1] = copy(list_conv)
    for i in 1:num_of_conv-1
        list_conv_save[i+1] = DSP.conv(list_conv_save[i],list_conv*Δx)[1:length(list_conv)]
    end
    base_save = Array{Any,1}(undef,plotrange)
    N = Int(1/Δx)
    for j in 1:num_of_conv
        base_save[j]=sum.([list_conv_save[j][1+i*N:(1+i)*N] for i in 0:plotrange-1])*Δx
    end
    base_save
end
