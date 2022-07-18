import numpy as np

class Vertex:
	'''This is a vertex class'''
	def __init__(self, lst:list):
		self.v = np.array(lst)

def subVertex(v1:Vertex, v2:Vertex):
	# substraction
	return Vertex(v1.v-v2.v) 

def addVertex(v1:Vertex, v2:Vertex):
	# summation
	return Vertex(v1.v+v2.v)

def lambdaVertex(c, ver:Vertex):
	# a constant times a vertex
	return Vertex(c*ver.v)

def negVertex(ver:Vertex):
	return Vertex(-ver.v)

def dot(v1:Vertex, v2:Vertex):
	return v1.v.dot(v2.v)

def norm(ver:Vertex):
	return np.sqrt(ver.v.dot(ver.v))

def xangle(p1:Vertex, p2:Vertex):
    # Computes the polar angle of the vector from p1 to p2
    d = p1.v-p2.v
    theta = (p2.v[0]-p1.v[0])/np.sqrt(d.dot(d))
    theta = np.arccos(theta)
    theta = theta+2*(np.pi-theta)*(p2.v[1]<p1.v[1])
    return theta

def fetchY(ver:Vertex):
    # taking the y-coordinate out of the vertex
    return ver.v[1]

def fetchX(ver:Vertex):
    # taking the x-coordinate out of the vertex
    return ver.v[0]