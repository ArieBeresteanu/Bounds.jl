module Bounds

using Statistics
using LinearAlgebra

### Silverman Rule of Thumb for bandwidth selection ##
function silverman(x::Vector{T}) where T<:Real
    # Silverman rule of thumb
    return 1.06*length(x)^(-0.2)*std(x)
end

### Epanechnikov Kernel function ###
function epa(x::T,x0::T=0.0,h::T=1.0) where T<:AbstractFloat 
    temp = min(abs((x-x0)/h),1.0)
    return 0.75*(1.0-temp^2)
end

function epa(x::Vector{T},x0::T=0.0,h::T=1.0) where T<:AbstractFloat 
    temp = min.(norm.((x.-x0)./h),1.0)
    return 0.75*(1.0.-temp.^2)
end


### Kernel density estimation ###
function kdens(X::Vector{T},x0::T,h::T) where T<:AbstractFloat 
    f = 0.0
    for x in x
        temp = (x-x0)/h
        f +=epa(temp)
    end
    return f/(length(X)*h)
end

end # module
