# The polygon structure and its functions.
#module Polygons

#export Polygon, minkowskiSum, lambdaPolygon, dirHausdorff, hausdorff

mutable struct Polygon
    vertices :: Vector{Vertex}
    isSorted :: Bool 
    #sort :: Function
    #plot :: Function
    #angles :: Function
    #scatter :: Function

    function Polygon(vertices)
        this = new()

        this.vertices=vertices
        this.isSorted = false
        
        # ## Sorting function
        # this.sort = function()
        #     n=length(this.vertices)
        #     #step 1: find the point with a minimal y coordinate and put it first.
        #     # comment: sorting is complexity nlog(n) but the following is just n
        #     #using sorting:
        #     #I = sortperm(fetchY.(this.vertices))
        #     #this.vertices = this.vertices[I]
        #     #going over the list
        #     m=fetchY(this.vertices[1])
        #     for i=2:n
        #         l=fetchY(this.vertices[i])
        #         if l<m #then swap
        #             m=l
        #             temp=this.vertices[i]
        #             this.vertices[i]=this.vertices[1]
        #             this.vertices[1]=temp
        #         end
        #     end
        #     #step 2: compute angles between the minimal vertex and all other vertices
        #     angs =zeros(n) #first column for angles and second column for the x coordinate
        #     angs[1]=-1
        #     v1 =this.vertices[1]
        #     for i=2:n
        #         angs[i] = xangle(v1,this.vertices[i])
        #     end
        #     #step 3: sort by angle
        #     I=sortperm(angs)
        #     this.vertices=this.vertices[I] 
        #     this.isSorted = true
        # end
        
        # ## polygon angles function
        # this.angles = function()
        #     if this.isSorted == false
        #         this.sort()
        #     end
        #     #this function
        #     n=length(this.vertices)
        #     ang=zeros(n)
        #     for i=1:n
        #         i==n ? j=1 : j=i+1
        #         ang[i] =xangle(this.vertices[i],this.vertices[j])                
        #     end
        #     return ang
        # end
        
        
    #     ## Scatter plot function
    #     this.scatter = function()
    #         if this.isSorted == false
    #             this.sort()
    #         end
    #         n=length(this.vertices)
    #         x=zeros(n); y=zeros(n);
    #         for i=1:n
    #             x[i]=this.vertices[i].v[1]
    #             y[i]=this.vertices[i].v[2]
    #         end
    #         scatter(x,y,label="")
    #     end       
        return this
    end
end

## Scatter plot function
function scatterPolygon(p::Polygon)
    if p.isSorted == false
        p.sort()
    end
    n=length(p.vertices)
    x=zeros(n); y=zeros(n);
    for i=1:n
        x[i]=p.vertices[i].v[1]
        y[i]=p.vertices[i].v[2]
    end

    mx, Mx = minimum(x), maximum(x)
    my, My = minimum(y), maximum(y)

    mx -= 0.05*abs(mx)
    Mx += 0.05*abs(Mx)

    my -= 0.05*abs(my)
    My += 0.05*abs(My)

    scatter(x,y,label="")
end

## Line plot function
function plotPolygon(p::Polygon)
    if p.isSorted == false
        sortPolygon!(p)
    end
    n=length(p.vertices)
    x=zeros(n+1); y=zeros(n+1);
    for i=1:n
        x[i]=p.vertices[i].v[1]
        y[i]=p.vertices[i].v[2]
    end
    x[n+1]=p.vertices[1].v[1]
    y[n+1]=p.vertices[1].v[2]

    mx, Mx = minimum(x), maximum(x)
    my, My = minimum(y), maximum(y)

    mx -= 0.05*abs(mx)
    Mx += 0.05*abs(Mx)

    my -= 0.05*abs(my)
    My += 0.05*abs(My)

    plot(x,y,label="",fill=true,xlims = (mx,Mx), ylims=(my,My) )
end

## polygon angles function
function angles(p::Polygon)
    if p.isSorted == false
        sortPolygon!(p)
    end
    #this function
    n=length(p.vertices)
    ang=zeros(n)
    for i=1:n
        i==n ? j=1 : j=i+1
        ang[i] =xangle(p.vertices[i],p.vertices[j])                
    end
    return ang
end

# the y coordinate of a vertex
function fetchY(ver::Vertex)
    # taking the y-coordinate out of the vertex
    return ver.v[2]
end

function xangle(p1::Vertex,p2::Vertex)
    #computes the angle that the vector starting from vertex p1 and ending at vertex p2 makes with the x-axis
    Δ = (p2-p1).v
    flag=false
    if Δ[2] < 0
        Δ[2] = -Δ[2]
        flag = true
    end
    xang =atan(Δ[2],abs(Δ[1]))
    if Δ[1]<0
        xang = π-xang
    end
    if flag
        xang=2*π-xang
    end
    return xang

end

## Sorting a Polygon function
function sortPolygon!(P::Polygon)
    n=length(P.vertices)
    #step 1: find the point with a minimal y coordinate and put it first.
    # comment: sorting is complexity nlog(n) but the following is just n
    #using sorting:
    #I = sortperm(fetchY.(P.vertices))
    #P.vertices = P.vertices[I]
    #going over the list
    m=fetchY(P.vertices[1])
    for i=2:n
        l=fetchY(P.vertices[i])
        if l<m #then swap
            m=l
            temp=P.vertices[i]
            P.vertices[i]=P.vertices[1]
            P.vertices[1]=temp
        end
    end
    #step 2: compute angles between the minimal vertex and all other vertices
    angs =zeros(n) #first column for angles and second column for the x coordinate
    angs[1]=-1
    v1 =P.vertices[1]
    for i=2:n
        angs[i] = xangle(v1,P.vertices[i])
    end
    #step 3: sort by angle
    I=sortperm(angs)
    P.vertices=P.vertices[I] 
    P.isSorted = true
    return P
end

###################################################
############# Summation functions: ################
###################################################

function sumTwoSegments(s1::Segment,s2::Segment)
    vers = [s1.p1+s2.p1,
        s1.p1+s2.p2,
        s1.p2+s2.p1,
        s1.p2+s2.p2        
    ]
    poly = Polygon(vers)
    sortPolygon!(poly)
    return poly
end    

Base.:+(s1::Segment,s2::Segment) = sumTwoSegments(s1,s2)

function minkowskiSum(v::Vertex,P::Polygon)
    # this function adds v to every vertex of P
    n=length(P.vertices)
    poly=P #initial value
    for i=1:n
        poly.vertices[i] +=v
    end
    sortPolygon!(poly)
    return poly
end

function minkowskiSum(P::Polygon,v::Vertex)
    # this function adds v to every vertex of P
    n=length(P.vertices)
    poly=P #initial value
    for i=1:n
        poly.vertices[i] +=v
    end
    sortPolygon!(poly)
    return poly
end
    
Base.:(+)(v::Vertex,P::Polygon) = minkowskiSum(v,P)
Base.:(+)(P::Polygon,v::Vertex) = minkowskiSum(P,v)

function minkowskiSum(P::Polygon,Q::Polygon)
    # Computes the minkowski sum of two convex polygons: P and Q. The polygons
    # are represented by their vertices and are ordered counter clockwise such
    #* that the first vertex will be the one who has the smallest Y coordinate
    # (and smallest X coordinate in case of a tie).  This assumption is maintained
    # in twoDproj by conditions in BLPcalculator.
    
    m = length(P.vertices)
    n = length(Q.vertices)
    
# case 1: Both P and Q are length 1 (vertices)
    if m==1 && n==1 
        R = Polygon([P.vertices[1]+Q.vertices[1]])

# case 2: P is length 1 (a vertex) and Q is not    
    elseif m==1    
        R = P.vertices[1] + Q

# case 3: Q is length 1 (a vertex) and P is not
    elseif n==1
        R = Q.vertices[1] + P
    
# case 4: both Q and P have more than 1 vertex
    else
        angP=[angles(P); 100] # 100 is just a big number that we know is larger
        angQ=[angles(Q); 100] # than all the angles which are between 0 and 2π
    
        m = length(P.vertices)
        n = length(Q.vertices)
    
        PP = [P.vertices; P.vertices[1]]
        QQ = [Q.vertices; Q.vertices[1]]
    
        #println("m=",m," n=",n)
    
        #println("angP= ", angP)
        #println("angQ= ", angQ)
    
        i=1; j=1;
        #println("----- begin ----------")
    
        R =Polygon([PP[1]+QQ[1]]) # a polygon with the sum of the two lower points as the first vertex.
        #println("R vertices: ",R.vertices)
        while (i<m+1 || j<n+1)
            if j == n+1 #angP[i]<angQ[j] 
                #println("angP[i] is minimal")
                i +=1
            elseif i == m+1 #angQ[j]<angP[i]
                #println("angQ[j] is minimal")
                j +=1
            else
                dif = angP[i]-angQ[j]
                if dif ≤ 0
                    i +=1
                end
                if dif ≥ 0
                    j +=1
                end
            end
            R.vertices = [ R.vertices; PP[i]+QQ[j]]
            #println(i,j)
            #println("R vertices: ",R.vertices)
        end
    end
    sortPolygon!(R)
    return R
end

Base.:(+)(P::Polygon,Q::Polygon) = minkowskiSum(P,Q)

function minkowskiSum(P::Polygon,s::Segment)
    # first, convert the segment to a sorted polygon with two vertices
    Q = Polygon([s.p1, s.p2])
    sortPolygon!(Q)
    
    # second, sum the two polygons
    return P+Q
end

Base.:(+)(P::Polygon,s::Segment) = minkowskiSum(P,s)


function lambdaPolygon(P::Polygon,λ::Real)
    n = length(P.vertices)
    R=P
    for i=1:n
       R.vertices[i]=λ*P.vertices[i] 
    end
    return R
end

#question: can this function be written using a map() function?

function dirHausdorff(P::Polygon,Q::Polygon)
    n = length(P.vertices)
    m = length(Q.vertices)
    dist = - Inf
    for i=1:n
        d = Inf
        for j=1:m
            next_j = j+1>m ? 1 : j+1
            d = min(d, dotDist(P.vertices[i],Segment(Q.vertices[j],Q.vertices[next_j])))
        end
        dist = max(dist,d)
    end
    return dist
end

function hausdorff(P::Polygon,Q::Polygon)
    d1 = dirHausdorff(P,Q)
    d2 = dirHausdorff(Q,P)
    return max(d1,d2)
end

#end #of module Polygons 
