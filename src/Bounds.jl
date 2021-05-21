module Bounds

using Statistics
using LinearAlgebra

### Silverman Rule of Thumb for bandwidth selection ##
function silverman(x::Vector{T}) where T<:Real
    # Silverman rule of thumb
    return 1.06*length(x)^(-0.2)*std(x)
end

### Epanechnikov Kernel function ###
function epa(x::T,x0::T=0.0,h::T=1.0) where T<:Real 
    temp = min(abs((x-x0)/h),1.0)
    return 0.75*(1.0-temp^2)
end

function epa(x::Vector{T},x0::T=0.0,h::T=1.0) where T<:Real 
    temp = min.(norm.((x.-x0)./h),1.0)
    return 0.75*(1.0.-temp.^2)
end


### Kernel density estimation ###
function kdens(X::Vector{T},x0::T,h::T) where T<:Real 
    f = sum(epa.((X .- x0) ./ h))
    return f/(length(X)*h)
end

### Kernel regression estimation ###
function kreg(Y::Vector{T},X::Vector{T},x0::T,h::T) where T<:Real 
    f = 0.0
    m = 0.0
	nx=length(X)
	ny=length(Y)
	if nx != ny
		error("length of Y and X do not match")
	else
		df = epa.((X .- x0) ./ h)
		f = sum(df)
		m = sum(Y.*df)
	end
	if f >0 
        return m/f
    else
        return 0.0
    end
end


end # module
