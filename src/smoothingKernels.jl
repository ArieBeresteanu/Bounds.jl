module smoothingKernels

# The unnormalized version of the kernel CAN be used in kernel regressions 
# because we take a ratio of two summations and the constant will cancel. 
# This may save some time and improve computational accuracy... maybe

# For kernel density we MUST use the notmalized (regular) version of the kernel
# functions because there is no ratio and the constant must be included!

uniform(u::T) where {T} = 1//2 * T(abs(u) <= 1.0)
uniform_unnormalized(u::T) where {T} = T(abs(u) <= 1.0)

triangular(u::T) where {T} = (1 - abs(u)) * T(abs(u) <= 1.0)
triangular_unnormalized(u::T) where {T} = (1 - abs(u)) * T(abs(u) <= 1.0)

epanechnikov(u::T) where {T} = 3//4 * (1 - u^2) * T(abs(u) <= 1.0)
epanechnikov_unnormalized(u::T) where {T} = (1 - u^2) * T(abs(u) <= 1.0)

biweight(u::T) where {T} = 15//16 * (1 - u^2)^2 * T(abs(u) <= 1.0)
biweight_unnormalized(u::T) where {T} = (1 - u^2)^2 * T(abs(u) <= 1.0)

triweight(u::T) where {T} = 35//32 * (1 - u^2)^3 * T(abs(u) <= 1.0)
triweight_unnormalized(u::T) where {T} = (1 - u^2)^3 * T(abs(u) <= 1.0)

tricube(u::T) where {T} = 70//81 * (1 - abs(u)^3)^3 * T(abs(u) <= 1.0)
tricube_unnormalized(u::T) where {T} = (1 - abs(u)^3)^3 * T(abs(u) <= 1.0)

gaussian(u::T) where {T} = (1 / sqrt(2 * pi)) *  exp(-1//2 * u^2)
gaussian_unnormalized(u::T) where {T} = exp(-1//2 * u^2)

cosine(u::T) where {T} = (pi / 4) * cos((pi / 2) * u) * T(abs(u) <= 1.0)
cosine_unnormalized(u::T) where {T} = cos((pi / 2) * u) * T(abs(u) <= 1.0)

logistic(u::T) where {T} = 1 / (exp(u) + 2 + exp(-u))
logistic_unnormalized(u::T) where {T} = 1 / (exp(u) + 2 + exp(-u))



kernels = Dict(
                :uniform => uniform,
                :triangular => triangular,
                :epanechnikov => epanechnikov,
                :gaussian => gaussian,
                :cosine => cosine,
                :logistic => logistic
            )  
unnormalized_kernels = Dict(
                :uniform => uniform_unnormalized,
                :triangular => triangular_unnormalized,
                :epanechnikov => epanechnikov_unnormalized,
                :gaussian => gaussian_unnormalized,
                :cosine => cosine_unnormalized,
                :logistic => logistic_unnormalized
            )


end