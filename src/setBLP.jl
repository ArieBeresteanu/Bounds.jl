module setBLP

import LinearAlgebra ,Base
using Statistics, Random

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

mutable struct testResults
	testStat :: Real
	criticalVal :: Real
	ConfidenceInterval :: Vector{<:Real}
end

mutable struct Results
	bound :: Vector{<:Real}
	Htest :: testResults
	dHtest :: testResults
end


#####################
###   Constants   ###
#####################

const default_options = Options(2000,15217,MersenneTwister(15217),0.95)


#####################
###  Functions:   ###
#####################

plus(::Real)=max(0.0,x)
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
	Htest = testResults(testStat_H,c_H,CI_H) 

	sort!(r_dH)
	c_dH = r_dH[floor(Int64,α*B)]
	CI_dH = [LB-c_dH/sqrt_n,UB+c_dH/sqrt_n]
	dHtest = testResults(testStat_dH,c_dH,CI_dH)

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
	sort!(r_H,dims=1)
	c_H = r_H[floor(Int64,α*length(rr))]
	CI_H = [LB-c_H/sqrt_n,UB+c_H/sqrt_n]
	Htest = testResults(testStat_H,c_H,CI_H) 

	#test based on directed Hausdorff distance:
	r_dH = maximum([plus.(rr[1,:]);minus.(rr[2,:])],dims=1)
	sort!(r_dH,dims=1)
	c_dH = r_dH[floor(Int64,α*length(rr))]
	CI_dH = [LB-c_dH/sqrt_n,UB+c_dH/sqrt_n]
	dHtest = testResults(testStat_dH,c_dH,CI_dH)

	results = Results(bound,Htest,dHtest)

	return results
end

###########################
###  Export Statement:  ###
###########################

export Options,default_options, Results, testResults, EY

end #of module
