# The unnormalized version of the kernel CAN be used in kernel regressions 
# because we take a ratio of two summations and the constant will cancel. 
# This may save some time and improve computational accuracy... maybe

# For kernel density MUST use the notmalized (regular) version of the kernel
# functions because there is no ratio and the constant must be included!
import math 

uniform = lambda u: 1/2*(abs(u)<=1.0)
uniform_unnormalized = lambda u: 1*(abs(u)<=1.0) 

triangular = lambda u: (1-abs(u))*(abs(u)<=1.0)
triangular_unnormalized = lambda u: (1-abs(u))*(abs(u)<=1.0)

epanechnikov = lambda u: 3/4 * (1 - u**2) * (abs(u) <= 1.0)
epanechnikov_unnormalized = lambda u: (1 - u**2) * (abs(u) <= 1.0)

biweight = lambda u: 15/16 * (1 - u**2)**2 * (abs(u) <= 1.0)
biweight_unnormalized = lambda u: (1 - u**2)**2 * (abs(u) <= 1.0)

triweight= lambda u: 35/32 * (1 - u**2)**3 * (abs(u) <= 1.0)
triweight_unnormalized= lambda u: (1 - u**2)**3 * (abs(u) <= 1.0)

tricube= lambda u: 70/81 * (1 - abs(u)**3)**3 * (abs(u) <= 1.0)
tricube_unnormalized= lambda u: (1 - abs(u)**3)**3 * (abs(u) <= 1.0)

gaussian= lambda u: (1 / math.sqrt(2 * math.pi)) *  math.exp(-1/2 * u**2)
gaussian_unnormalized= lambda u: math.exp(-1/2 * u**2)

cosine= lambda u: (math.pi / 4) * math.cos((math.pi / 2) * u) * (abs(u) <= 1.0)
cosine_unnormalized= lambda u: math.cos((math.pi / 2) * u) * (abs(u) <= 1.0)

logistic= lambda u: 1 / (math.exp(u) + 2 + math.exp(-u))
logistic_unnormalized= lambda u: 1 / (math.exp(u) + 2 + math.exp(-u))

kernels = {
    "uniform": uniform,
    "triangular": triangular,
    "epanechnikov": epanechnikov,
    "gaussian": gaussian,
    "cosine": cosine,
    "logistic": logistic
}

unnormalized_kernels = {
    "uniform": uniform_unnormalized,
    "triangular": triangular_unnormalized,
    "epanechnikov": epanechnikov_unnormalized,
    "gaussian": gaussian_unnormalized,
    "cosine": cosine_unnormalized,
    "logistic": logistic_unnormalized
}