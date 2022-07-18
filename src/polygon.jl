# The polygon structure and its functions.

mutable struct Polygon
    vertices :: Vector{Vertex}
    isSorted :: Bool 
    sort :: Function
    plot :: Function
    angles :: Function
    scatter :: Function

    function Polygon(vertices)
        this = new()

        this.vertices=vertices
        this.isSorted = false
        
        ## Sorting function
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
            angs =zeros(n) #first column for angles and second column for the x coordinate
            angs[1]=-1
            v1 =this.vertices[1]
            for i=2:n
                angs[i] = xangle(v1,this.vertices[i])
            end
            #step 3: sort by angle
            I=sortperm(angs)
            this.vertices=this.vertices[I] 
            this.isSorted = true
        end
        
        ## polygon angles function
        this.angles = function()
            if this.isSorted == false
                this.sort()
            end
            #this function
            n=length(this.vertices)
            ang=zeros(n)
            for i=1:n
                i==n ? j=1 : j=i+1
                ang[i] =xangle(this.vertices[i],this.vertices[j])                
            end
            return ang
        end
        
        ## Line plot function
        this.plot = function()
            if this.isSorted == false
                this.sort()
            end
            n=length(this.vertices)
            x=zeros(n+1); y=zeros(n+1);
            for i=1:n
                x[i]=this.vertices[i].v[1]
                y[i]=this.vertices[i].v[2]
            end
            x[n+1]=this.vertices[1].v[1]
            y[n+1]=this.vertices[1].v[2]
            plot(x,y,label="",fill=true)
        end
        
        ## Scatter plot function
        this.scatter = function()
            if this.isSorted == false
                this.sort()
            end
            n=length(this.vertices)
            x=zeros(n); y=zeros(n);
            for i=1:n
                x[i]=this.vertices[i].v[1]
                y[i]=this.vertices[i].v[2]
            end
            scatter(x,y,label="")
        end
            
        return this
    end
end

function minkowskiSum(v::Vertex,P::Polygon)
    # this function adds v to every vertex of P
    n=length(P.vertices)
    R=P #initial value
    for i=1:n
        R.vertices[i] +=v
    end
    return R
end

function minkowskiSum(P::Polygon,v::Vertex)
    # this function adds v to every vertex of P
    n=length(P.vertices)
    R=P #initial value
    for i=1:n
        R.vertices[i] +=v
    end
    return R
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
        angP=[P.angles(); 100]
        angQ=[Q.angles(); 100]
    
    #m = length(angP)
    #n = length(angQ)
    
        PP = [P.vertices; P.vertices[1]]
        QQ = [Q.vertices; Q.vertices[1]]
    
        println("m=",m," n=",n)
    
        println("angP= ", angP)
        println("angQ= ", angQ)
    
        i=1; j=1;
        println("----- begin ----------")
    
        R =Polygon([PP[1]+QQ[1]]) # a polygon with the sum of the two lower points as the first vertex.
        println("R vertices: ",R.vertices)
        while (i<m+1 || j<n+1)
            if angP[i]<angQ[j] 
                println("angP[i] is minimal")
                i +=1
            elseif angQ[j]<angP[i]
                println("angQ[j] is minimal")
                j +=1
            else
                i +=1
                j +=1
            end
            R.vertices = [ R.vertices; PP[i]+QQ[j]]
            println(i,j)
            println("R vertices: ",R.vertices)
        end
    end
    return R
end

Base.:(+)(P::Polygon,Q::Polygon) = minkowskiSum(P,Q)

function lambdaPolygon(λ::Real,P::Polygon)
    n = length(P.vertices)
    R=P
    for i=1:n
       R.vertices[i]=λ*P.vertices[i] 
    end
end

#question: can this function be written using a map() function?

function dirHausdorff(P::Polygon,Q::Polygon)
    n = length(P.vertices)
    m = length(Q.vertices)
    dist_global = Inf
    for i=1:n
        dist = Inf
        for j=1:m
            next_j = j+1>m ? 1 : j+1
            d = dotDist(P.vertices[i],Segment(Q.vertices[j],Q.vertices[next_j]))
            dist = min(dist,d)
        end
        dist_global=max(dis_global,dist)
    end
    return dist
end

function hausdorff(P::Polygon,Q::Polygon)
    d1 = dirHausdorff(P,Q)
    d2 = dirHausdorff(Q,P)
    return max(d1,d2)
end
