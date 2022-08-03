module setBLP

import LinearAlgebra ,Base
using Statistics

include("vertices.jl") 
export Vertex, subVertex, addVertex, lambdaVertex, negVertex, xangle, fetchX, fetchY

include("segments.jl") 
export Segment, dotDist, xangle

include("polygons.jl")  
export Polygon, minkowskiSum, lambdaPolygon, dirHausdorff, hausdorff

function EY(yl::Vector{Float64},yu::Vector{Float64})
	return mean(yl), mean(yu)
end

export EY

end #of module