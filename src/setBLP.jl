module setBLP

import LinearAlgebra ,Base
using Statistics

include("vertices.jl") 
export Vertex, subVertex, addVertex, lambdaVertex, negVertex, xangle, fetchX, fetchY, distVertex

include("segments.jl") 
export Segment, dotDist, xangle

include("polygons.jl")  
export Polygon, minkowskiSum, lambdaPolygon, dirHausdorff, hausdorff

function EY(yl::Vector{Float64},yu::Vector{Float64},H0::Vector{Float64})
	LB = mean(yl)
	UB = mean(yu)
	bound = Vertex([LB,UB])
	return bound, distVertex(bound,H0)
end

export EY

end #of module