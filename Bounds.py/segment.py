class Segment:
	def __init__(self, p1:Vertex, p2:Vertex):
		self.p1 = p1
		self.p2 = p2

	def checkInput(self):
		return len(self.p1.v) == len(self.p2.v)

	def length(self):
		if self.checkInput():
			return norm(subVertex(self.p1,self.p2))
		else:
			return False 

	def dim(self):
		if self.checkInput():
			return len(self.p1.v)
		else:
			return False

def dotDist(p:Vertex, segment:Segment):
	# Finds the minimal distance between a point p and
    # the line segment connecting points p1 and p2
	if segment.checkInput():
		if len(p.v) == segment.dim():
			p31 = subVertex(p, segment.p1)
			p21 = subvertex(segment.p2, segment.p1)

			t = dot(p31, p21)/dot(p21, p21)
			t = min(1,max(t,0))

			p0 = addVertex(segment.p1, lambdaVertex(t, p21))

			return norm(subVertex(p, p0))
		else:
			raise ValueError("dimension of p doesnt match dimension of segment")
	else:
		raise ValueError("segment has wrong dimensions")