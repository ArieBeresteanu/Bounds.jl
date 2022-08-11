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


function EY(yl::Vector{Float64},yu::Vector{Float64},H0::Vector{Float64},options::Options=default_options)
    LB = mean(yl)
	UB = mean(yu)
	bound = Vertex([LB,UB])

	# test Statistic
	n = length(yl)
	testStat = sqrt(n)*distVertex(bound,H0)

	#critical value based on Hausdorff distance
	σ = cov(yl,yU)
	Pi = [var(yl) σ; σ var(yu)] #covariance matrix for yl,yu
	
	d = MvNormal([0, 0],Pi) #defining the joint normal distribution
	B= options.MC_iterations #number of MC iterations to compute the critical value
	α = options.conf_level  #confidence level for the critical value1``

	rr = abs.(rand(d,B));
	r = maximum(rr,dims=1);
	sort!(r,dims=1)

	c_H = r[floor(Int64,α*length(r))]

	rr = abs.(rand(d,B));
	r = maximum(rr,dims=1);
	sort!(r,dims=1)
	

	c_H = r[floor(Int64,α*length(r))]

	
	return bound, testStat, c_H
end


export Options,EY

end #of module