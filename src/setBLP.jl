module setBLP

import LinearAlgebra

mutable struct Segment
    p1::Vector{Real}
    p2::Vector{Real}
    checkInput::Function
    length::Function
    dim::Function

    
    function Segment(p1,p2)
        this = new()

        this.p1 = p1
        this.p2 = p2

        this.checkInput = function()
            return size(this.p1) == size(this.p2)       
        end

        this.length = function()
            return norm(this.p1-this.p2)
        end

        this.dim = function()
            if this.checkInput()
                return length(this.p1)
            else
                return false
            end
        end

        return this        
    end
end


function dotDist(p::Vector{<:Real}, segment::Segment) 
    if segment.checkInput()
        if length(p) == segment.dim()
            p1_p2 = segment.p1 -segment.p2
            p_p2 = p -segment.p2

            λ = dot(p1_p2,p_p2)/dot(p1_p2,p1_p2)
            λ = max(min(λ,1),0)

            p0 = λ*segment.p1 + (1-λ)*segment.p2 

            return norm(p-p0)
        else
            return "dimention of p doesnt match dimention of segment"
        end
    else
        return "Segment has wrong dimentions"
    end
end

end #of module