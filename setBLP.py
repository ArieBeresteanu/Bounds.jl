#########################
### Imported packages ###
#########################
import pandas as pd
import math

#############
## dotDist ##
#############

def dotDist(p1, p2, p3):
    # Finds the minimal distance between a point p3 and
    # the line segment connecting points p1 and p2
    p31 = [a-b for (a,b) in zip(p3,p1)]
    p21 = [a-b for (a,b) in zip(p2,p1)]
    t = sum([a*b for (a,b) in zip(p31, p21)])/sum([a**2 for a in p21])
    t = min(1,max(t,0))
    p0 = [a+t*b for (a,b) in zip(p1,p21)]
    dist = math.dist(p3,p0)
    return dist

############
## xangle ##
############
def xangle(p1, p2):
    # Computes the angle between the vector from p1 to p2, and the X-axis
    theta = (p2[0]-p1[0])/math.dist(p1,p2)
    theta = math.acos(theta)
    theta = theta+2*(math.pi-theta)*(p2[1]<p1[1])
    return theta

# for 2-dimension only
class Polygons:
    # inputs are coordinates (x,y) of all vertices of a convex polygon in any order
    # x: x-coordinate of vertices ordered counter-clockwise, starting with the vertex with the smallest y-coordinate.
    # y: y-coordinate of vertices ordered counter-clockwise, starting with the vertex with the smallest y-coordinate.
    
    # theta: the angle between the line connecting the first vertex and other vertices in counter-clockwise order
    
    def __init__(self, *args):
        # find the starting point p1
        df = pd.DataFrame(args, columns=['x','y'])
        df = df.sort_values(by=['y','x'],ascending = [True, False])
        p1 = df.iloc[0,:]
        
        # calculate theta between p1 and other vertices
        theta = []
        for i in range(df.shape[0]):
            theta.append(xangle(p1, df.iloc[i,:]))
        df['theta'] = theta
        
        # sort by theta (counter-clockwise order)
        df.sort_values(by=['theta'], inplace = True)
        
        # get rid of non-vertex points:
        # caculate d, the distance between any vertex to the line segment connecting the previous and next vertices
        # if d is 0, drop the middle point
        d = [dotDist(df.iloc[-1,0:2], df.iloc[1,0:2], df.iloc[0,0:2])]
        for i in range(1,df.shape[0]-1):
            d.append(dotDist(df.iloc[i-1,0:2], df.iloc[i+1,0:2], df.iloc[i,0:2]))
        d.append(dotDist(df.iloc[-2,0:2], df.iloc[0,0:2], df.iloc[-1,0:2]))
        df['d'] = d
        df = df.loc[df['d']>0.0]
        
        self.x = list(df.loc[:,'x'])
        self.y = list(df.loc[:,'y'])
#         self.theta = list(df.loc[:,'theta'])  

#############################################
###      M A I N   F U N C T I O N S      ###
#############################################

#########################
# 1. Hausdorff Distance #
#########################

# directed Hausdorff distance
def dHausdorff(P1:Polygons, P2:Polygons):
    # directed Hausdorff distance between two polygons
    d_inf = []
    for i in range(len(P1.x)):
        P3 = [P1.x[i], P1.y[i]]
        d = dotDist([P2.x[0], P2.y[0]],[P2.x[-1], P2.y[-1]], P3)
        for j in range(1,len(P2.x)-1):
            d = min(d,dotDist([P2.x[j], P2.y[j]],[P2.x[j+1], P2.y[j+1]], P3))
        d_inf.append(d)
    return max(d_inf)

def Hausdorff(P1:Polygons, P2:Polygons):
    # Hausdorff distance between two polygons
    return max(dHausdorff(P1, P2), dHausdorff(P2, P1))

####################
# 2. Minkowski Sum #
####################

def minksum(P1:Polygons, P2:Polygons):
    P = ()
    for i in range(len(P1.x)):
        for j in range(len(P2.x)):
            P+=([P1.x[i]+P2.x[j],P1.y[i]+P2.y[j]],)
    return Polygons(*P)