module Bounds

using Statistics

function silverman(x::Vector{T}) where T<:Real
    return 1.06*length(x)^(-0.2)*std(x)
end

function epa(x::T) where T<:Real
    abs(x)<1 ? temp=x : temp=0.0
    return 0.75*(1-temp^2)
end

function kdens(X::Vector{T},x0::T,h::T) where T<:Real
    f = 0.0
    for x in x
        temp = (x-x0)/h
        f +=epa(temp)
    end
    return f/(length(X)*h)
end

end # module
