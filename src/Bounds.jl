module Bounds

using Statistics
using LinearAlgebra

export SimpleBound,missingObs

### Silverman Rule of Thumb for bandwidth selection ##
function silverman(x::Vector{T}) where T<:Real
    # Silverman rule of thumb
    return 1.06*length(x)^(-0.2)*std(x)
end

### Epanechnikov Kernel function ###
#     Version 1: X,x0 and h are scalars
function epa(x::T,x0::T=0.0,h::T=1.0) where T<:Real 
    temp = min(abs((x-x0)/h),1.0)
    return 0.75*(1.0-temp^2)
end

#     Version 2: X is a vector, x0,h are scalars
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


function kreg(Y::Vector{T},X::Vector{T},x0::T,h::T,w::Vector{T}) where T<:AbstractFloat 
    # This is a version of kreg where we can introduce weights. 
    # It is done mostly so we can conditon on binary weights, i.e. exclude some observations

    f = 0.0
    m = 0.0
	nx=length(X)
	ny=length(Y)
    nw=length(w)
	
    if nx != ny || nx !=nw
		error("length of Y and X do not match")
	else
		df = w.*epa.((X .- x0) ./ h)
		f = sum(df)
		m = sum(w.*Y.*df)
	end
	if f > 0.000001 
        return m/f
    else
        return 0.0
    end
end

mutable struct assumptions
    tol   :: Float64
    Yₗ     :: Float64
    Yᵤ    :: Float64
    model :: String
end

### Worst Case Scenario Bounds 
mutable struct Results{T} 
    prob0 :: T
    prob1 :: T
    yhat0 :: T
    yhat1 :: T 
    bound0L :: T 
    bound0U :: T 
    bound1L :: T 
    bount1U :: T 
    treatL :: T 
    treatU :: T 
    model  :: String
end

function missingObs(y::Vector{T},z::Vector{T},x::Vector{T},x0::T,cont::Bool,h::T=1.0) where T<:Real
    # if z=0,  y is missing
    ny = length(y)
    nx = length(x)
    nz = length(z)
    if (ny !=nx) || (ny != nz) 
        error("vector length not matching")
    end
    zer0 ::T = 0
    res::Vector{T} = Results(zer0,0,0,0,0,0,0,0,0,0,"")
    if cont
        res.prob1 = kreg(z,x,x0,h)
        res.prob0 = 1- res.prob1
        res.yhat1 = kreg(y,x,x0,h,z)
        res.yhat0 = NaN 
    else
        nx0 = count(x .== x0)
        temp1 = z .* (x .== x0)
        nz1 = count(temp1)

        res.prob1 = nz1/nx0
        res.prob0 = 1 - res.prob1    
        res.yhat1 = count(y .* z .* (x .== x0))/nz1
        res.yhat0 = NaN 
    end
    res.bound0L = NaN
    res.bound0U = NaN
    res.bound1L = K0 * res.prob0 + res.yhat1 * res.prob1
    res.bound1U = K1 * res.prob0 + res.yhat1 * res.prob1
    res.treatL = NaN
    res.treatU = NaN
    res.model  = "Missing Y observations"
end

function treatmentEffect(y::Vector{T},z::Vector{T},x::Vector{T},x0::T,cont::Bool,h::T=1.0) where T<:Real
    ny = length(y)
    nx = length(x)
    nz = length(z)
    if (ny !=nx) || (ny != nz) 
        error("vector length not matching")
    end
    zer0 ::T = 0
    res::Vector{T} = Results(zer0,0,0,0,0,0,0,0,0,0,"")
    if cont
        res.prob1 = kreg(z,x,x0,h)
        res.prob0 = 1- res.prob1
        res.yhat1 = kreg(y,x,x0,h,z)
        res.yhat0 = kreg(y,x,x0,h,1 .-z) 
    else
        nx0 = count(x .== x0)
        temp1 = z .* (x .== x0)
        nz1 = count(temp1)

        res.prob1 = nz1/nx0
        res.prob0 = 1 - res.prob1    
        res.yhat1 = count(y .* z .* (x .== x0))/nz1
        res.yhat0 = count(y .* (1 .-z) .* (x .== x0))/nz1 
    end
    res.model  = "Treatment Effect, no assumptions"
end


struct SimpleBound
	LB :: Float64
	UB :: Float64
	method :: String
end

end # module
