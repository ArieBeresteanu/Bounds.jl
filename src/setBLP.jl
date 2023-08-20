module setBLP

import LinearAlgebra ,Base
using Statistics, Random, Distributions
using DataFrames

####################
###   Includes:  ###
####################

include("vertices.jl") 
export Vertex, subVertex, addVertex, lambdaVertex, negVertex, xangle, fetchX, fetchY

include("segments.jl") 
export Segment, dotDist, xangle

include("polygons.jl")  
export Polygon, minkowskiSum, lambdaPolygon, dirHausdorff, hausdorff

###############################
###   Defined Structures:   ###
###############################  

mutable struct Options
	MC_iterations::Int64
	seed::Int64
	rng::AbstractRNG
	conf_level::Float64
end

function Base.show(o::Options; io::IO=stdout)
	println(io, "Options:")
	println(io, "  Number of MC iterations: ", o.MC_iterations)
	println(io, "  Seed: ", o.seed)
	println(io, "  Random Number Generator: ", o.rng)
	println(io, "  Confidence level: ", o.conf_level)
  end
  
mutable struct TestResults
	testStat :: Real
	criticalVal :: Real
	ConfidenceInterval :: Vector{Real}
end

mutable struct Results
	null  :: Vector{<:Real}
	bound :: Vector{<:Real}
	Htest :: TestResults
	dHtest :: TestResults
end

function Base.show(r::Results; io::IO=stdout, digits::Int=4)
    print(io, "Results: \n")
    print(io, "  Null: $(round.(r.null, digits=digits))\n") 
    print(io, "  Bound: $(round.(r.bound, digits=digits))\n") 
    print(io, "  Hausdorff based test: \n")
    print(io, "    Test Stat: $(round(r.Htest.testStat, digits=digits))\n")
    print(io, "    Critical Value: $(round(r.Htest.criticalVal, digits=digits))\n")
    print(io, "    Confidence Interval: $(round.(r.Htest.ConfidenceInterval, digits=digits))\n")
    print(io, "  directed Hausdorff test: \n")
    print(io, "    Test Stat: $(round(r.dHtest.testStat, digits=digits))\n")
    print(io, "    Critical Value: $(round(r.dHtest.criticalVal, digits=digits))\n")
    print(io, "    Confidence Interval: $(round.(r.dHtest.ConfidenceInterval, digits=digits))\n")
end

#####################
###   Constants   ###
#####################

const default_options = Options(2000,15217,MersenneTwister(),0.95)


#####################
###  Functions:   ###
#####################

plus(x::Real)=max(0.0,x)
minus(x::Real)=max(0.0,-x)

function HdistInterval(v1::Vector{<:Real},v2::Vector{<:Real})
    v = v1 - v2
    return maximum(abs.(v))
end
  
function dHdistInterval(v1::Vector{<:Real},v2::Vector{<:Real})
    v = v1 - v2
	return maximum([plus(v[1]),minus(v[2])])
end

## Plan: add DataFrame capabilities
function EY(yl::Vector{<:Real},yu::Vector{<:Real},H0::Vector{<:Real},options::Options=default_options,method="Asymptotic")
	#THis is the shell function that calls either the asymtotic distribution version or the bootstrap version of EY
	if method =="Asymptotic"
		EYasy(yl,yu,H0,options)
	else
		EYboot(yl,yu,H0,options)
	end
end


function EYboot(yl::Vector{<:Real},yu::Vector{<:Real},H0::Vector{<:Real},options::Options=default_options)
	#This function uses a bootstrap test. This option is not in BM(2008) for EY but it is proved for BLP in section 4
	LB = mean(yl)
	UB = mean(yu)
	bound = [LB,UB]

	# test Statistic
	n = length(yl)
	sqrt_n = sqrt(n)
	bound = vec(bound)
	testStat_H = sqrt_n*HdistInterval(bound,H0)
	testStat_dH = sqrt_n*dHdistInterval(bound,H0)

	B = options.MC_iterations #number of MC iterations to compute the critical value
	α = options.conf_level  #confidence level for the critical value1
	distribution = DiscreteUniform(1,n)

	r_H=zeros(n)
	r_dH = zeros(n)

	for i=1:B
		indx = rand(options.rng,distribution,n)
		yl_b = yl[indx]
		yu_b = yu[indx]
		bound_b = [mean(yl_b),mean(yu_b)]
		r_H[i] = sqrt_n * HdistInterval(bound_b,bound)
		r_dH[i] = sqrt_n * dHdistInterval(bound_b,bound)
	end
	sort!(r_H)
	c_H = r_H[floor(Int64,α*B)]
	CI_H = [LB-c_H/sqrt_n,UB+c_H/sqrt_n]
	Htest = TestResults(testStat_H,c_H,CI_H) 

	sort!(r_dH)
	c_dH = r_dH[floor(Int64,α*B)]
	CI_dH = [LB-c_dH/sqrt_n,UB+c_dH/sqrt_n]
	dHtest = TestResults(testStat_dH,c_dH,CI_dH)

	results = Results(bound,Htest,dHtest)

	return results
end

function EYasy(yl::Vector{<:Real},yu::Vector{<:Real},H0::Vector{<:Real},options::Options=default_options)
	#This function uses the test based on the asymptotic distributin as developed in BM(2008) pp. 778-779
    LB = mean(yl)
	UB = mean(yu)
	bound = [LB,UB]

	# test Statistic
	n = length(yl)
	sqrt_n = sqrt(n)
	testStat_H = sqrt_n*HdistInterval(bound,H0)
	testStat_dH = sqrt_n*dHdistInterval(bound,H0)

	#Simulating the asy. distribution using a MC method to establish a critical value (quantile):

	# Drawing pairs of bivariate normal r.v.'s 
	σ = cov(yl,yu)
	Pi = [var(yl) σ; σ var(yu)] #covariance matrix for yl,yu
	
	d = MvNormal([0, 0],Pi) #defining the joint normal distribution
	B = options.MC_iterations #number of MC iterations to compute the critical value
	α = options.conf_level  #confidence level for the critical value1

	## Following Algorithm on page 780 in BM2008:
	rr = (rand(d,B)); #drawing B pairs from a bivariate-normal distribution.
	
	## test based on Hausdorff distance:
	r_H = maximum(abs.(rr),dims=1);
	sort!(r_H,dims=2)
	c_H = r_H[floor(Int64,α*B)]
	CI_H = [LB-c_H/sqrt_n,UB+c_H/sqrt_n]
	Htest = TestResults(testStat_H,c_H,CI_H) 

	#test based on directed Hausdorff distance:
	r_dH = maximum([plus.(rr[1,:]) minus.(rr[2,:])],dims=2)
	sort!(r_dH,dims=1)
	c_dH = r_dH[floor(Int64,α*B)]
	CI_dH = [LB-c_dH/sqrt_n,UB+c_dH/sqrt_n]
	dHtest = TestResults(testStat_dH,c_dH,CI_dH)

	results = Results(bound,Htest,dHtest)

	return results
end

###########################
###  oneDproj functions ###
###########################

## Vector/Matrix versions ##

# 1. x is assumed to be one dimensional

function oneDproj(yl::Vector{<:Real},yu::Vector{<:Real},x::Vector{<:Real})
	x = x.-mean(x) #demean x
	M = [x.*yl x.*yu]
	s = sum(x.*x)
	lb = sum(minimum(M,dims=2)) / s
	ub = sum(maximum(M,dims=2)) / s 
	return [lb ub]
end

# 2. X is assumed to be a matrix of covariates and a single coordinate is specified

function oneDproj(yl::Vector{<:Real},yu::Vector{<:Real},x::Matrix{<:Real},cord::Int64)
    # The function assumes that the matrix x does not contain a 1 Vector
    our_x = x[:,cord] #taking out the coordinate of interest
    new_x = copy(x)
    new_x[:,cord] .= 1.0  #replacing the column with a vector of ones
    pred_x =new_x*(inv(new_x'*new_x)*new_x'*our_x)
	res_x = our_x - pred_x
    bound = oneDproj(yl,yu,res_x)
    return bound
end

# 3. X is assumed to be a matrix of covariates and a vector of coordinates is specified

function oneDproj(yl::Vector{<:Real},yu::Vector{<:Real},x::Matrix{<:Real},cords::Vector{Int64})
    # The function assumes that the matrix x does not contain a 1 Vector
	bounds = []
    for cord in cords
		our_x = x[:,cord]
		new_x = copy(x)
		new_x[:,cord] .= 1.0  #replacing the column with a vector of ones
		pred_x =new_x*(inv(new_x'*new_x)*new_x'*our_x)
		res_x = our_x - pred_x
    	bound = oneDproj(yl,yu,res_x)
		push!(bounds,bound)
	end
    return bound
end

# 4. X is assumed to be a matrix of covariates and a coordinate is NOT specified

function oneDproj(yl::Vector{<:Real},yu::Vector{<:Real},x::Matrix{<:Real})
    # The function assumes that the matrix x does not contain a 1 Vector
	ncols = size(x,2)
	bounds = []
    for cord in 1:ncols
		our_x = x[:,cord]
		new_x = copy(x)
		new_x[:,cord] .= 1.0  #replacing the column with a vector of ones
		pred_x =new_x*(inv(new_x'*new_x)*new_x'*our_x)
		res_x = our_x - pred_x
    	bound = oneDproj(yl,yu,res_x)
		push!(bounds,bound)
	end
    return bounds
end

## Data frame versions ##

function oneDproj(df::DataFrame, yl::Symbol,yu::Symbol,x::Symbol)
	y_l = copy(df[!,yl])
	y_u = copy(df[!,yu])
	new_x = Vector(df[!,x])
	bound = oneDproj(y_l,y_u,new_x)
	return bound
end

function oneDproj(df::DataFrame, yl::Symbol,yu::Symbol,x::Vector{Symbol})
	y_l = copy(df[!,yl])
	y_u = copy(df[!,yu])
	new_x = Matrix(df[!,x])
	bounds = oneDproj(y_l,y_u,new_x)
	return bounds
end

function oneDproj(df::DataFrame, yl::Symbol,yu::Symbol,x::Vector{Symbol},cord::Int64)
	y_l = copy(df[!,yl])
	y_u = copy(df[!,yu])
	new_x = Matrix(df[!,x])
	bounds = oneDproj(y_l,y_u,new_x,cord)
	return bounds
end

###################### End of oneDproj functions ######################################


function CI1d(yl::Vector{<:Real},yu::Vector{<:Real},x::Vector{<:Real},H0::Vector{<:Real},options::Options=default_options)
	## computes the 1D projection of the identification set on a specific dimesion of the explanatory variable

	#step 1: Compute the bounds on page 787 in BM2008
	
	
	bound = vec(oneDproj(yl,yu,x)) #vectorizing is necessary because the HdistInterval function wants two vectrs as input
	LB = bound[1]
	UB = bound[2]

	#step 2: Compute the test statisticss

	n = length(yl)
	sqrt_n = sqrt(n)
	testStat_H = sqrt_n*HdistInterval(bound,H0)
	testStat_dH = sqrt_n*dHdistInterval(bound,H0)

	#step 3: Bootstrap iterationss 

	B = options.MC_iterations #number of MC iterations to compute the critical value
	α = options.conf_level  #confidence level for the critical value1
	distribution = DiscreteUniform(1,n)

	r_H=zeros(B)
	r_dH = zeros(B)

	for i=1:B
		indx = rand(options.rng,distribution,n)
		yl_b = yl[indx]
		yu_b = yu[indx]
		x_b  = x[indx]
		bound_b = vec(oneDproj(yl_b,yu_b,x_b))
		r_H[i] = sqrt_n * HdistInterval(bound_b,bound)
		r_dH[i] = sqrt_n * dHdistInterval(bound_b,bound)
	end

	sort!(r_H)
	c_H = r_H[floor(Int64,α*B)]
	CI_H = [LB-c_H/sqrt_n,UB+c_H/sqrt_n]
	Htest = TestResults(testStat_H,c_H,CI_H) 

	sort!(r_dH)
	c_dH = r_dH[floor(Int64,α*B)]
	CI_dH = [LB-c_dH/sqrt_n,UB+c_dH/sqrt_n]
	dHtest = TestResults(testStat_dH,c_dH,CI_dH)

	results = Results(H0,bound,Htest,dHtest)

	return results #,r_H,r_dH
end
###########################
###  Export Statement:  ###
###########################

export Options,default_options, Results, TestResults, EY, CI1d, oneDproj

end #of module
