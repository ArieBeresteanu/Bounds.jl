module setBLP

import LinearAlgebra ,Base
using Statistics

include("vertices.jl") 
export Vertex, subVertex, addVertex, lambdaVertex, negVertex, xangle, fetchX, fetchY

include("segments.jl") 
export Segment, dotDist, xangle

include("polygons.jl")  
export Polygon, minkowskiSum, lambdaPolygon, dirHausdorff, hausdorff


function EY(yl::Vector{Float64},yu::Vector{Float64},H0::Vector{Float64},B::Int64=1000,α::Float64=0.95;)
    LB = mean(yl)
	UB = mean(yu)
	bound = Vertex([LB,UB])

	# test Statistic
	n = length(yl)
	testStat = sqrt(n)*distVertex(bound,H0)

	#critical value based on 
	σ = cov(yl,yU)
	Pi = [var(yl) σ; σ var(yu)]
	
	d = MvNormal([0, 0],Pi)

	rr = abs.(rand(d,B));
	r = maximum(rr,dims=1);
	sort!(r,dims=1)

	c_H = r[floor(Int64,α*length(r))]

	rr = abs.(rand(d,10));
	r = maximum(rr,dims=1);
	sort!(r,dims=1)
	α = 0.95;

	c_H = r[floor(Int64,α*length(r))]

	
	return bound, testStat, c_H
end


export EY

end #of module