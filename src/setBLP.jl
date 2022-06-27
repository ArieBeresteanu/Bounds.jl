module setBLP

import LinearAlgebra,Base

export Vertex,Segment,dotDist,Polygon, xangle

mutable struct Vertex
    v::Vector{Real}
end

Base.:(-)(v1::Vertex,v2::Vertex) = Vertex(v1.v-v2.v)

mutable struct Segment
    p1::Vertex
    p2::Vertex
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

function xangle(seg::Segment)
    Δ = seg.p2-seg.p1
    flag=false
    if Δ[2] < 0
        Δ[2] = -Δ[2]
        flag = true
    end
    xang =atan(Δ[2],abs(Δ[1]))
    if Δ[1]<0
        xang = pi-xang
    end
    if flag
        xang=2*pi-xang
    end
    return xang

end

function xangle(p1::Vertex,p2::Vertex)
    Δ = (p2-p1).v
    flag=false
    if Δ[2] < 0
        Δ[2] = -Δ[2]
        flag = true
    end
    xang =atan(Δ[2],abs(Δ[1]))
    if Δ[1]<0
        xang = pi-xang
    end
    if flag
        xang=2*pi-xang
    end
    return xang

end

function fetchY(ver::Vertex)
    return ver.v[2]
end


mutable struct Polygon
    vertices :: Vector{Vertex}
    sort :: Function

    function Polygon(vertices)
        this = new()

        this.vertices=vertices
        this.sort = function()
            n=length(this.vertices)
            #step 1: find the point with a minimal y coordinate and put it first.
            # comment: sorting is complexity nlog(n) but the following is just n
            #using sorting:
            #I = sortperm(fetchY.(this.vertices))
            #this.vertices = this.vertices[I]
            #going over the list
            m=fetchY(this.vertices[1])
            for i=2:n
                l=fetchY(this.vertices[i])
                if l<m #then swap
                    m=l
                    temp=this.vertices[i]
                    this.vertices[i]=this.vertices[1]
                    this.vertices[1]=temp
                end
            end
            #step 2: compute angles between the minimal vertex and all other vertices
            angs =zeros(n)
            angs[1]=-1
            v1 =this.vertices[1]
            for i=2:n
                angs[i] = xangle(v1,this.vertices[i])
            end
            #step 3: sort by angle
            I=sortperm(angs)
            this.vertices=this.vertices[I]   
        end
        return this
    end
end

end #of module
