#########################
### Imported packages ###
#########################
import numpy as np
import pandas as pd

import vertex
# Class: Vertex
# Functions: subVertex, addVertex, lambdaVertex, negVertex, dot, norm, xangle, fetchY, fetchX, distVertex

import segment
# Class: Segment
##-methods: checkInput, length, dim
# Functions: dotDist

import polygons
# Class: Polygon
##-method: sort, angles, plot, scatter
# Functions: minkowskiSum, lambdaPolygon, dHausdorff, hausdorff

def EY(yl, yu, H0):
	LB = np.mean(yl)
	UB = np.mean(yu)
	bound = Vertex([LB, UB])
	return bound, distVertex(bound, Vertex(H0))

