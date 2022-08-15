module setBLP

import LinearAlgebra ,Base
using Statistics, Random

include("vertices.jl") 
export Vertex, subVertex, addVertex, lambdaVertex, negVertex, xangle, fetchX, fetchY

include("segments.jl") 
export Segment, dotDist, xangle

include("polygons.jl")  
export Polygon, minkowskiSum, lambdaPolygon, dirHausdorff, hausdorff

mutable struct Options
	MC_iterations::Int64
	seed::Int64
	rng::AbstractRNG
	conf_level::Float64
end

const default_options = Options(2000,15217,MersenneTwister(15217),0.95)


plus(::Real)=max(0.0,x)
minus(x::Real)=max(0.0,-x)


function HdistInterval(v1::Vector{Real},v2::Vector{Real})
    v = v1 - v2
    return maximum(abs.(v))
end
  
function dHdistInterval(v1::Vector{Real},v2::Vector{Real})
    v = v1 - v2
	return maximum([plus(v[1])+minus(v[2])])
end

function EY(yl::Vector{Float64},yu::Vector{Float64},H0::Vector{Float64},options::Options=default_options)
    LB = mean(yl)
	UB = mean(yu)
	bound = [LB,UB]

	# test Statistic
	n = length(yl)
	sqrt_n = sqrt(n)
	testStat_H = sqrt_n*distVertex(bound,H0)
	destStat_dH

	#critical value based on Hausdorff distance
	σ = cov(yl,yU)
	Pi = [var(yl) σ; σ var(yu)] #covariance matrix for yl,yu
	
	d = MvNormal([0, 0],Pi) #defining the joint normal distribution
	B= options.MC_iterations #number of MC iterations to compute the critical value
	α = options.conf_level  #confidence level for the critical value1

	## Following Algorithm on page 780 in BM2008:
	rr = (rand(d,B)); #drawing B pairs from a bivariate-normal distribution.
	
	## test based on Hausdorff distance:
	r_H = maximum(abs.(rr),dims=1);
	sort!(r_H,dims=1)
	c_H = r_H[floor(Int64,α*length(rr))]
	CI_H = [yl-c_H/sqrt_n,yu+c_H/sqrt_n]

	#test based on directed Hausdorff distance:
	r_dH = maximum([plus.(rr[1,:]);minus.(rr[2,:])],dims=1)
	sort!(r_dH,dims=1)
	c_dH = r_dH[floor(Int64,α*length(rr))]
	CI_dH = [yl-c_dH/sqrt_n,yu+c_dH/sqrt_n]


	return bound, testStat_H, c_H, CI_H, c_dH, CI_dH
end


export Options,EY

end #of module