import matplotlib.pyplot as plt 
import numpy as np

class Polygon:
	def __init__(self, *vertices):
		# inputs are coordinates (x,y) of all vertices of a convex polygon in any order
		self.vertices = list(vertices)
		self.isSorted = False

	## sorting
	def sort(self):
		n = len(self.vertices)
		# step 1: find the point with a minimal y coordinate 
		# (and smallest X coordinate in case of a tie) and put it first
		m = fetchY(self.vertices[0])
		for i in range(1,n):
			l = fetchY(self.vertices[i])
			if l < m or (l == m and fetchX(self.vertices[i])< fetchX(self.vertices[0])): # swap
				m = l
				temp = self.vertices[i]
				self.vertices[i] = self.vertices[0]
				self.vertices[0] = temp	

		#step 2: compute angles between the first vertex and all other vertices
		v1 = self.vertices[0]
		angs = [xangle(v1, ver) for ver in self.vertices]
		# sort by angs, then by X-coordinate (in case of collinearity on the first edge)
		I = np.lexsort((list(map(fetchX, self.vertices)),angs))
		self.vertices = [self.vertices[i] for i in I]
		self.isSorted = True

	## polygon angles
	def angles(self):
		if not self.isSorted:
			self.sort()
		n = len(self.vertices)
		ang = []
		for i in range(n-1):
			ang.append(xangle(self.vertices[i], self.vertices[i+1]))
		ang.append(xangle(self.vertices[-1], self.vertices[0]))
		return ang

	## line plot
	def plot(self):
		if not self.isSorted:
			self.sort()
		x = [ver.v[0] for ver in self.vertices]
		y = [ver.v[1] for ver in self.vertices]
		plt.fill(x,y)

	## scatter plot
	def scatter(self):
		if not self.isSorted:
			self.sort()
		x = [ver.v[0] for ver in self.vertices]
		y = [ver.v[1] for ver in self.vertices]
		plt.scatter(x,y)

def minkowskiSum(P:Polygon, Q:Polygon):
    # Computes the minkowski sum of two convex polygons: P and Q. The polygons
    # are represented by their vertices and are ordered counter clockwise such
    #* that the first vertex will be the one who has the smallest Y coordinate
    # (and smallest X coordinate in case of a tie).  This assumption is maintained
    # in twoDproj by conditions in BLPcalculator.
    m = len(P.vertices)
    n = len(Q.vertices)
    R = []
    if m ==1 or n == 1:
        for p in P.vertices:
            for q in Q.vertices:
                R+= [addVertex(p,q)]
    else:
        i = 0
        j = 0
        PP = P.vertices + [P.vertices[0]]
        QQ = Q.vertices + [Q.vertices[0]]
        angP = P.angles() + [P.angles()[0]]
        angQ = Q.angles() + [Q.angles()[0]]
        while i<m or j<n:
            R+= [addVertex(PP[i],QQ[j])]
            if i == m:
                j+=1
            elif j==n:
                i+=1
            else:        
                dif = angP[i]-angQ[j]
                if dif>=0:
                    j+=1
                if dif<=0:
                    i+=1
    return Polygon(*tuple(R))

def lambdaPolygon(c, P:Polygon):
	R = P
	R.vertices = [lmbdaVertex(c,ver) for ver in P.vertices]
	return R

# directed Hausdorff distance
def dHausdorff(P:Polygon, Q:Polygon):
	# directed Hausdorff distance between two polygons
	d_inf = []
	for i in range(len(P.vertices)):
		p = P.vertices[i]
		d = dotDist(p, Segment(Q.vertices[0],Q.vertices[-1]))
		for j in range(len(Q.vertices)-1):
			d = min(d,dotDist(p, Segment(Q.vertices[j],Q.vertices[j+1])))
		d_inf.append(d)
	return max(d_inf)

def hausdorff(P:Polygon, Q:Polygon):
    # Hausdorff distance between two polygons
    return max(dHausdorff(P,Q), dHausdorff(Q,P))